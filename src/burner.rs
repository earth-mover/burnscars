//! Public entry point.

use std::ops::RangeInclusive;

use geo_traits::{GeometryTrait, GeometryType};

use crate::approx;
use crate::coverage;
use crate::error::BurnError;
use crate::grid::{GeomId, Grid, NorthUp, NorthUpGrid};
use crate::input::{Envelope, Loaded, for_each_xy};
use crate::output::{BurnResult, GeomStart, PolyRecs, Stretches};
use crate::scratch::Scratch;
use crate::{line, point};

/// What a cell of a polygon counts. Lines and points do not read it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rule {
    /// Exact coverage: an interior cell at 1, an edge cell at its fraction
    /// (controlledburn `BurnMode::Coverage`).
    Fraction,
    /// Every cell with coverage above the drop threshold, with no weight.
    Touched,
    /// Every cell whose centre is inside, with no weight
    /// (controlledburn `BurnMode::Approx`, fasterize semantics).
    Centre,
}

impl Rule {
    pub(crate) fn has_values(self) -> bool {
        self == Self::Fraction
    }
}

/// The origin of the row and column indices in [`Stretches`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Origin {
    /// Whole-grid indices.
    Grid,
    /// Indices from the window's first row and column.
    Window,
}

/// `Send` and `Sync`; use one per worker, because it owns reused scratch buffers.
#[derive(Debug)]
pub struct Burner {
    rule: Rule,
    scratch: Scratch,
    /// The polygon records of the current geometry.
    polys: PolyRecs,
}

impl Burner {
    #[must_use]
    pub fn new(rule: Rule) -> Self {
        Self { rule, scratch: Scratch::default(), polys: PolyRecs::default() }
    }

    #[must_use]
    pub fn rule(&self) -> Rule {
        self.rule
    }

    /// The capacity of the burner's own buffers now, in bytes:
    /// input copy, walk scratch, cell accumulator and the polygon records of `reserve_records`.
    #[must_use]
    pub fn scratch_bytes(&self) -> u64 {
        u64::try_from(self.scratch.bytes() + self.polys.bytes()).unwrap_or(u64::MAX)
    }

    /// Reserves room for `n` polygon records (16 B each), so the buffer does not grow during a burn.
    /// One geometry gives about one record per boundary cell and per interior run, per copy.
    pub fn reserve_records(&mut self, n: usize) {
        self.polys.reserve(n);
    }

    /// Appends the polygons' spans to `out`, with indices from `origin`.
    ///
    /// Refuses `out` of another rule (`values` `Some` against `None`).
    /// Refuses a geometry of a non-polygonal type, also an empty one, before it loads,
    /// so every window gives the same answer.
    /// Refuses an id that is not above the previous one in this call or above the last
    /// fiber already in `out`. Refuses more than `i32::MAX` values.
    /// Stops at the first error; `out` is then as it was before the failing geometry.
    pub fn burn_polygons<G, I>(
        &mut self,
        grid: &Grid,
        geoms: I,
        origin: Origin,
        out: &mut Stretches,
    ) -> Result<(), BurnError>
    where
        G: GeometryTrait<T = f64>,
        I: IntoIterator<Item = (GeomId, G)>,
    {
        out.check_rule(self.rule)?;
        let origin = match origin {
            Origin::Grid => (0, 0),
            Origin::Window => (grid.rows().start().0, grid.cols().start().0),
        };
        let mut previous = out.id.last().copied().map(GeomId);
        geoms.into_iter().try_for_each(|(id, g)| {
            if let Some(previous) = previous.filter(|&p| id <= p) {
                return Err(BurnError::IdOrder { id, previous });
            }
            previous = Some(id);
            self.burn_polygon(grid, id, &g, origin, out)
        })
    }

    /// Burns each `(id, line or point)` onto the window of `grid` and appends to `out`.
    /// A polygon is `BurnError::Polygonal`: use [`Burner::burn_polygons`].
    ///
    /// Stops at the first error. There is no rollback. Each geometry is validated before
    /// it writes, and its `GeomIndex` is pushed only when it is complete. After
    /// `WalkDiverged` or `Internal`, line and point records of the failing geometry can
    /// follow the last `GeomIndex`; truncate each table to the end of the last `GeomIndex`
    /// to drop them.
    pub fn burn_into<G, I>(
        &mut self,
        grid: &Grid,
        geoms: I,
        out: &mut BurnResult,
    ) -> Result<(), BurnError>
    where
        G: GeometryTrait<T = f64>,
        I: IntoIterator<Item = (GeomId, G)>,
    {
        geoms.into_iter().try_for_each(|(id, g)| self.burn_one(grid, id, &g, out))
    }

    /// [`Burner::burn_into`] into a new `BurnResult`.
    pub fn burn<G, I>(&mut self, grid: &Grid, geoms: I) -> Result<BurnResult, BurnError>
    where
        G: GeometryTrait<T = f64>,
        I: IntoIterator<Item = (GeomId, G)>,
    {
        let mut out = BurnResult::default();
        self.burn_into(grid, geoms, &mut out)?;
        Ok(out)
    }

    fn burn_polygon<G: GeometryTrait<T = f64>>(
        &mut self,
        grid: &Grid,
        id: GeomId,
        g: &G,
        origin: (u32, u32),
        out: &mut Stretches,
    ) -> Result<(), BurnError> {
        self.polys.clear();
        // The declared type, so an empty line or point is refused too. `load` refuses a collection.
        if is_polygonal(g) == Some(false) {
            return Err(BurnError::NotPolygonal { id });
        }
        match self.scratch.input.load(g, id)? {
            Loaded::Empty => return Ok(()),
            Loaded::Lines | Loaded::Points => return Err(BurnError::NotPolygonal { id }),
            Loaded::Polygons => {}
        }
        let n = grid.north_up();
        if n.is_identity() {
            self.burn_polys_on(&n.grid, id)?;
        } else {
            let (period, ks) = copies(&n, self.scratch.input.envelope(), Loaded::Polygons, id)?;
            for k in ks {
                self.burn_polys_on(&n.shifted(k, period), id)?;
            }
            self.polys.map_rows(n.rows);
        }
        self.polys.encode(id, origin, out)
    }

    /// The ported polygon algorithms, on a north-up grid, into the record buffer.
    fn burn_polys_on(&mut self, grid: &NorthUpGrid, id: GeomId) -> Result<(), BurnError> {
        let s = &mut self.scratch;
        match self.rule {
            Rule::Fraction | Rule::Touched => {
                coverage::burn(&s.input, grid, id, &mut s.walk, &mut s.cov, &mut self.polys)
            }
            Rule::Centre => approx::burn(&s.input, grid, id, &mut s.approx, &mut self.polys),
        }
    }

    fn burn_one<G: GeometryTrait<T = f64>>(
        &mut self,
        grid: &Grid,
        id: GeomId,
        g: &G,
        out: &mut BurnResult,
    ) -> Result<(), BurnError> {
        let start = out.start_geom();
        let loaded = self.scratch.input.load(g, id)?;
        if loaded == Loaded::Polygons {
            return Err(BurnError::Polygonal { id });
        }
        let n = grid.north_up();
        if n.is_identity() {
            self.burn_loaded(loaded, &n.grid, id, out)?;
        } else {
            self.burn_mapped(loaded, &n, id, start, out)?;
        }
        out.end_geom(id, start);
        Ok(())
    }

    /// The ported line and point algorithms, on a north-up grid.
    fn burn_loaded(
        &mut self,
        loaded: Loaded,
        grid: &NorthUpGrid,
        id: GeomId,
        out: &mut BurnResult,
    ) -> Result<(), BurnError> {
        let s = &mut self.scratch;
        match loaded {
            Loaded::Empty | Loaded::Polygons => {}
            Loaded::Lines => line::burn(&s.input, grid, id, &mut s.walk, out)?,
            Loaded::Points => point::burn(&s.input, grid, id, out),
        }
        Ok(())
    }

    /// A south-up grid or a period: burns on the north-up grid, one copy per period that
    /// meets the geometry, then maps the rows to the caller's indices.
    fn burn_mapped(
        &mut self,
        loaded: Loaded,
        n: &NorthUp,
        id: GeomId,
        start: GeomStart,
        out: &mut BurnResult,
    ) -> Result<(), BurnError> {
        let (period, ks) = copies(n, self.scratch.input.envelope(), loaded, id)?;
        let copies = ks.map(|k| n.shifted(k, period));
        if loaded == Loaded::Points {
            point::burn_shifted(&self.scratch.input, &copies, id, out);
        } else {
            let mut copies_hit = 0;
            for g in copies {
                let before = out.len();
                self.burn_loaded(loaded, &g, id, out)?;
                copies_hit += usize::from(out.len() > before);
            }
            if copies_hit > 1 {
                out.merge_tail(id, start);
            }
        }
        out.map_tail(start, n.rows);
        Ok(())
    }
}

/// The declared type is polygonal: `Some(true)` or `Some(false)`; `None` for a collection.
fn is_polygonal<G: GeometryTrait<T = f64>>(g: &G) -> Option<bool> {
    match g.as_type() {
        GeometryType::Polygon(_)
        | GeometryType::MultiPolygon(_)
        | GeometryType::Rect(_)
        | GeometryType::Triangle(_) => Some(true),
        GeometryType::GeometryCollection(_) => None,
        GeometryType::Point(_)
        | GeometryType::LineString(_)
        | GeometryType::MultiPoint(_)
        | GeometryType::MultiLineString(_)
        | GeometryType::Line(_) => Some(false),
    }
}

/// The period and the copies `k` (shifts by `k` periods) that meet the window.
/// No period: one copy, `k = 0`.
fn copies(
    n: &NorthUp,
    env: Option<Envelope>,
    loaded: Loaded,
    id: GeomId,
) -> Result<(f64, RangeInclusive<i64>), BurnError> {
    let (Some(period), Some(env)) = (n.period, env) else {
        return Ok((0.0, 0..=0));
    };
    let width = env.xmax - env.xmin;
    if width > period {
        return Err(BurnError::WiderThanPeriod { id, width, period });
    }
    let l = n.grid.lattice();
    // Points need every copy on the full grid: the seam rule compares the
    // hits of two copies, and a window can hold only one of them.
    let cols = if loaded == Loaded::Points { n.grid.full_cols() } else { n.grid.cols() };
    let west = l.x_edge(i64::from(cols.start().0));
    let east = l.x_edge(i64::from(cols.end().0));
    // Copy k covers [west + k p, east + k p]; keep each copy that meets the envelope.
    let lo = to_i64(((env.xmin - east) / period).ceil());
    let hi = to_i64(((env.xmax - west) / period).floor());
    Ok((period, lo..=hi))
}

/// Saturating: a copy that far away does not meet the window.
#[expect(clippy::cast_possible_truncation, reason = "saturating float to int cast")]
fn to_i64(v: f64) -> i64 {
    v as i64
}

/// Envelope of a geometry (finite coordinates only), as `[xmin, ymin, xmax, ymax]`.
/// Lets an orchestrator select the geometries for a chunk.
pub fn bbox<G: GeometryTrait<T = f64>>(geom: &G) -> Option<[f64; 4]> {
    let mut env: Option<[f64; 4]> = None;
    for_each_xy(geom, &mut |x, y| {
        if x.is_finite() && y.is_finite() {
            env = Some(match env {
                None => [x, y, x, y],
                Some([a, b, c, d]) => [a.min(x), b.min(y), c.max(x), d.max(y)],
            });
        }
    });
    env
}

#[cfg(test)]
mod tests {
    use std::ops::Range;

    use geo_types::{Geometry, LineString, Polygon, Rect};

    use super::*;
    use crate::coverage::fraction::to_f32;
    use crate::grid::{Col, Lattice, Row, Span, XAxis, YAxis, YOrigin};
    use crate::output::{Rec, StretchCounts};

    type TestResult = Result<(), Box<dyn std::error::Error>>;

    /// An edge record (test only).
    fn edge(row: u32, col: u32, w: f32) -> Rec {
        Rec::edge(Row(row), Col(col), w)
    }

    /// A run record (test only).
    fn run(row: u32, cols: Range<u32>) -> Rec {
        Rec::run(Row(row), Span::from_parts(Col(cols.start), Col(cols.end)))
    }

    /// Feeds one geometry's records, in the given order, to the merge and the encoder (test only).
    fn merge_one(rule: Rule, recs: &[Rec]) -> Result<Stretches, BurnError> {
        let mut p = PolyRecs::default();
        recs.iter().for_each(|&r| p.push(r));
        let mut out = Stretches::new(rule);
        p.encode(GeomId(0), (0, 0), &mut out)?;
        Ok(out)
    }

    fn fixture() -> [Rec; 5] {
        [edge(0, 9, 0.6), run(0, 4..9), edge(0, 2, 0.3), edge(0, 10, 0.1), edge(0, 3, 0.85)]
    }

    fn square(x0: f64, y0: f64) -> Geometry<f64> {
        Rect::new((x0, y0), (x0 + 2.0, y0 + 2.0)).into()
    }

    #[test]
    fn a_geometry_gives_edge_blocks_beside_the_interior_run() -> TestResult {
        // Row 0: edges at cols 2, 3 (.3, .85), a run 4..9, edges at 9, 10 (.6, .1), pushed out of order.
        let spans = merge_one(Rule::Fraction, &fixture())?;
        assert_eq!(spans.col_start, [2, 4, 9]);
        assert_eq!(spans.col_end, [4, 9, 11]);
        let values = spans.values.ok_or("no values")?;
        assert_eq!(values.value_ptr, [0, 2, 3, 5]);
        assert_eq!(values.weight, [0.3, 0.85, 1.0, 0.6, 0.1]);
        Ok(())
    }

    #[test]
    fn periodic_copies_merge_as_merge_tail_does() -> TestResult {
        // Copy A: run 2..5, edge 8 (.4). Copy B: run 4..7, edge 5 (.5) inside the joined run, edge 8 (.3).
        let recs = [run(0, 2..5), edge(0, 8, 0.4), run(0, 4..7), edge(0, 5, 0.5), edge(0, 8, 0.3)];
        let spans = merge_one(Rule::Fraction, &recs)?;
        assert_eq!(
            (spans.col_start.as_slice(), spans.col_end.as_slice()),
            (&[2, 8][..], &[7, 9][..])
        );
        let sum = to_f32(f64::from(0.4_f32) + f64::from(0.3_f32));
        assert_eq!(spans.values.as_ref().ok_or("no values")?.weight, [1.0, sum]);
        assert_eq!(spans.densify().len(), 6);
        Ok(())
    }

    #[test]
    fn a_summed_seam_edge_near_one_joins_the_run() -> TestResult {
        let spans =
            merge_one(Rule::Fraction, &[edge(0, 7, 0.5), run(0, 2..7), edge(0, 7, 0.499_999_5)])?;
        assert_eq!((spans.col_start, spans.col_end), (vec![2], vec![8]));
        assert_eq!(spans.values.ok_or("no values")?.weight, [1.0]);
        Ok(())
    }

    #[test]
    fn touched_merges_edges_into_runs() -> TestResult {
        let spans = merge_one(Rule::Touched, &fixture())?;
        assert_eq!((spans.col_start, spans.col_end), (vec![2], vec![11]));
        assert_eq!(spans.values, None);
        Ok(())
    }

    #[test]
    fn rows_of_one_geometry_are_fibers_in_row_order() -> TestResult {
        // A south-up grid reverses the rows.
        let spans = merge_one(Rule::Fraction, &[run(3, 0..2), run(1, 0..2), edge(2, 4, 0.5)])?;
        assert_eq!((spans.id, spans.row), (vec![0, 0, 0], vec![1, 2, 3]));
        assert_eq!(spans.fiber_ptr, [0, 1, 2, 3]);
        Ok(())
    }

    #[test]
    fn an_empty_burn_has_one_pointer() -> TestResult {
        let grid = Grid::from_extent(0.0, 0.0, 10.0, 10.0, 10, 10)?;
        let mut out = Stretches::new(Rule::Fraction);
        Burner::new(Rule::Fraction).burn_polygons(
            &grid,
            [(GeomId(0), square(20.0, 20.0))],
            Origin::Grid,
            &mut out,
        )?;
        assert_eq!(out.fiber_ptr, [0]);
        assert_eq!(out.values.ok_or("no values")?.value_ptr, [0]);
        Ok(())
    }

    #[test]
    fn burn_polygons_refuses_a_line_an_unsorted_id_and_a_repeated_id() -> TestResult {
        let grid = Grid::from_extent(0.0, 0.0, 10.0, 10.0, 10, 10)?;
        let mut b = Burner::new(Rule::Fraction);
        for x in [1.0, 50.0] {
            let line: Geometry<f64> = LineString::from(vec![(x, 1.0), (x + 2.0, 3.0)]).into();
            let r = b.burn_polygons(
                &grid,
                [(GeomId(0), line)],
                Origin::Grid,
                &mut Stretches::new(Rule::Fraction),
            );
            assert!(matches!(r, Err(BurnError::NotPolygonal { id: GeomId(0) })), "{r:?}");
        }
        for second in [1, 2] {
            let mut out = Stretches::new(Rule::Fraction);
            let items = [(GeomId(2), square(1.0, 1.0)), (GeomId(second), square(5.0, 5.0))];
            let r = b.burn_polygons(&grid, items, Origin::Grid, &mut out);
            assert!(matches!(r, Err(BurnError::IdOrder { .. })), "{r:?}");
            assert_eq!(out.id, [2, 2]);
            let again =
                b.burn_polygons(&grid, [(GeomId(2), square(5.0, 5.0))], Origin::Grid, &mut out);
            assert!(matches!(again, Err(BurnError::IdOrder { .. })), "{again:?}");
        }
        Ok(())
    }

    #[test]
    fn burn_into_refuses_a_polygon() -> TestResult {
        let grid = Grid::from_extent(0.0, 0.0, 10.0, 10.0, 10, 10)?;
        let p = Polygon::new(LineString::from(vec![(1.0, 1.0), (3.0, 1.0), (3.0, 3.0)]), vec![]);
        let r = Burner::new(Rule::Fraction).burn(&grid, [(GeomId(4), p)]);
        assert!(matches!(r, Err(BurnError::Polygonal { id: GeomId(4) })), "{r:?}");
        Ok(())
    }

    #[test]
    fn scratch_bytes_grows_with_a_burn() -> TestResult {
        let grid = Grid::from_extent(0.0, 0.0, 10.0, 10.0, 10, 10)?;
        let mut b = Burner::new(Rule::Fraction);
        assert_eq!(b.scratch_bytes(), 0);
        b.burn_polygons(
            &grid,
            [(GeomId(0), square(1.5, 1.5))],
            Origin::Grid,
            &mut Stretches::new(Rule::Fraction),
        )?;
        assert!(b.scratch_bytes() > 0);
        Ok(())
    }

    #[test]
    fn an_edge_at_a_run_start_goes_into_the_run() -> TestResult {
        let spans = merge_one(Rule::Fraction, &[edge(0, 4, 0.3), run(0, 4..9)])?;
        assert_eq!((spans.col_start, spans.col_end), (vec![4], vec![9]));
        assert_eq!(spans.values.ok_or("no values")?.weight, [1.0]);
        Ok(())
    }

    #[test]
    fn a_sum_below_the_full_threshold_stays_a_block_and_a_sum_above_one_clamps() -> TestResult {
        let below = [run(0, 2..7), edge(0, 7, 0.5), edge(0, 7, 0.499_998_5)];
        let spans = merge_one(Rule::Fraction, &below)?;
        assert_eq!((spans.col_start, spans.col_end), (vec![2, 7], vec![7, 8]));
        let sum = to_f32(f64::from(0.5_f32) + f64::from(0.499_998_5_f32));
        assert_eq!(spans.values.ok_or("no values")?.weight, [1.0, sum]);
        let over = merge_one(Rule::Fraction, &[edge(0, 7, 0.6), edge(0, 7, 0.6), run(0, 2..7)])?;
        assert_eq!((over.col_start, over.col_end), (vec![2], vec![8]));
        assert_eq!(over.values.ok_or("no values")?.weight, [1.0]);
        Ok(())
    }

    #[test]
    fn truncate_to_a_mark_restores_the_table() -> TestResult {
        let mut out = Stretches::new(Rule::Fraction);
        let mut p = PolyRecs::default();
        let mut encode = |id: u32, recs: &[Rec], out: &mut Stretches| {
            p.clear();
            recs.iter().for_each(|&r| p.push(r));
            p.encode(GeomId(id), (0, 0), out)
        };
        encode(0, &[run(0, 1..4), edge(1, 2, 0.5), edge(1, 3, 0.25)], &mut out)?;
        let first = out.clone();
        let mark = out.mark();
        encode(1, &[edge(0, 0, 0.5), run(0, 1..3), run(2, 0..9)], &mut out)?;
        assert_ne!(out, first);
        out.truncate(mark);
        assert_eq!(out, first);
        Ok(())
    }

    #[test]
    fn a_rule_must_match_the_table() -> TestResult {
        let grid = Grid::from_extent(0.0, 0.0, 10.0, 10.0, 10, 10)?;
        let n = StretchCounts { fibers: 4, stretches: 8, values: 16 };
        let mut touched = Stretches::with_capacity(Rule::Touched, n);
        let sq = [(GeomId(0), square(1.5, 1.5))];
        Burner::new(Rule::Touched).burn_polygons(&grid, sq.clone(), Origin::Grid, &mut touched)?;
        assert!(!touched.is_empty());
        assert_eq!(touched.values, None);

        for mut out in [Stretches::new(Rule::Fraction), touched] {
            let before = out.clone();
            let rule = if out.values.is_some() { Rule::Touched } else { Rule::Fraction };
            let r = Burner::new(rule).burn_polygons(&grid, sq.clone(), Origin::Grid, &mut out);
            assert!(matches!(r, Err(BurnError::RuleMismatch)), "{r:?}");
            assert_eq!(out, before);
        }
        Ok(())
    }

    /// One rect on a window at rows 3.., cols 5..: `Window` is `Grid` less `(3, 5)`.
    #[test]
    fn a_window_origin_subtracts_the_window_corner() -> TestResult {
        for origin in [YOrigin::NorthUp, YOrigin::SouthUp] {
            let y0 = if origin == YOrigin::NorthUp { 10.0 } else { 0.0 };
            let l = Lattice::new(XAxis::new(0.0, 1.0)?, YAxis::new(y0, 1.0, origin)?);
            let full = Grid::new(l, Span::new(Row(0), Row(10))?, Span::new(Col(0), Col(10))?)?;
            let w = full.window(Span::new(Row(3), Row(10))?, Span::new(Col(5), Col(10))?)?;
            let rect = [(GeomId(0), Geometry::from(Rect::new((4.5, 1.5), (8.5, 5.5))))];
            let mut b = Burner::new(Rule::Fraction);
            let mut grid = Stretches::new(Rule::Fraction);
            b.burn_polygons(&w, rect.clone(), Origin::Grid, &mut grid)?;
            let mut window = Stretches::new(Rule::Fraction);
            b.burn_polygons(&w, rect, Origin::Window, &mut window)?;
            assert!(!grid.is_empty());
            assert!(grid.row.iter().all(|&r| r >= 3) && grid.col_start.iter().all(|&c| c >= 5));
            let shift = |v: &[u32], d: u32| v.iter().map(|x| x - d).collect::<Vec<_>>();
            assert_eq!(window.row, shift(&grid.row, 3), "{origin:?}");
            assert_eq!(window.col_start, shift(&grid.col_start, 5), "{origin:?}");
            assert_eq!(window.col_end, shift(&grid.col_end, 5), "{origin:?}");
            assert_eq!((&window.fiber_ptr, &window.values), (&grid.fiber_ptr, &grid.values));
        }
        Ok(())
    }

    /// One copy on a south-up grid: rows strictly descend, so the encoder walks them in reverse
    /// and does not sort.
    #[test]
    fn descending_rows_are_walked_in_reverse_with_no_sort() -> TestResult {
        let recs = [run(3, 0..2), edge(3, 4, 0.5), edge(2, 1, 0.25), run(2, 2..5), run(1, 7..9)];
        let mut p = PolyRecs::default();
        recs.iter().for_each(|&r| p.push(r));
        let mut out = Stretches::new(Rule::Fraction);
        p.encode(GeomId(0), (0, 0), &mut out)?;
        assert_eq!(p.recs(), recs);
        assert_eq!(out.row, [1, 2, 3]);
        assert_eq!((out.col_start, out.col_end), (vec![7, 1, 2, 0, 4], vec![9, 2, 5, 2, 5]));
        assert_eq!(out.values.ok_or("no values")?.weight, [1.0, 0.25, 1.0, 1.0, 0.5]);
        Ok(())
    }

    #[test]
    fn reserve_records_reserves() {
        let mut b = Burner::new(Rule::Fraction);
        b.reserve_records(1000);
        assert!(b.scratch_bytes() >= 1000 * size_of::<Rec>() as u64);
    }
}
