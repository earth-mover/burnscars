//! The padded grid of one geometry's processing region.
//!
//! Layout of exactextract `Grid<infinite_extent>` and `grid_cell` (Apache-2.0, Copyright
//! 2018-2022 ISciences, LLC), commit 67969d99d1dfe1115be955087e2c48e1aab133e2,
//! `src/grid.h`, `src/grid.cpp`: local index 0 and n + 1 are padding, 1..=n are real cells.
//! Cell edges come from the lattice. The snap is controlledburn `snap_to_lattice`
//! (Apache-2.0, Copyright 2025 Michael Sumner), commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558,
//! `cpp/src/scanline.cpp`, without its range test.

use crate::geom::{box_::CellBox, coord::Coord, side::Side};
use crate::grid::{Col, Grid, Pos, Row, Span};
use crate::input::Envelope;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct LocalRow(pub(crate) usize);
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct LocalCol(pub(crate) usize);

#[derive(Clone, Copy, Debug)]
pub(crate) struct WalkerGrid {
    grid: Grid,
    rows: Span<Row>,
    cols: Span<Col>,
}

impl WalkerGrid {
    /// exactextract `processing_region` + `get_geometry_grid` (`src/raster_cell_intersection.cpp`)
    /// on the lattice, with a one-cell margin in place of `Grid::shrink_to_fit` fixes.
    pub(crate) fn for_envelope(grid: &Grid, env: &Envelope) -> Option<Self> {
        let l = grid.lattice();
        let (wr, wc) = (grid.rows(), grid.cols());
        let (c0, c1) = (i64::from(wc.start().0), i64::from(wc.end().0));
        let (r0, r1) = (i64::from(wr.start().0), i64::from(wr.end().0));
        // Closed intersection with the window: touching counts.
        if env.xmax < l.x_edge(c0)
            || env.xmin > l.x_edge(c1)
            || env.ymin > l.y_edge(r0)
            || env.ymax < l.y_edge(r1)
        {
            return None;
        }
        let (fc, fr) = (grid.full_cols(), grid.full_rows());
        let ci = |p: Pos<Col>| match p {
            Pos::Before => i64::from(fc.start().0) - 1,
            Pos::In(c) => i64::from(c.0),
            Pos::After => i64::from(fc.end().0),
        };
        let ri = |p: Pos<Row>| match p {
            Pos::Before => i64::from(fr.start().0) - 1,
            Pos::In(r) => i64::from(r.0),
            Pos::After => i64::from(fr.end().0),
        };
        let c_lo = (ci(grid.col_of(env.xmin)) - 1).max(c0);
        let c_hi = (ci(grid.col_of(env.xmax)) + 1).min(c1 - 1);
        let r_lo = (ri(grid.row_of(env.ymax)) - 1).max(r0); // ymax gives the low row index
        let r_hi = (ri(grid.row_of(env.ymin)) + 1).min(r1 - 1);
        if c_lo > c_hi || r_lo > r_hi {
            return None;
        }
        let u = |v: i64| u32::try_from(v).ok();
        Some(Self {
            grid: *grid,
            rows: Span::from_parts(Row(u(r_lo)?), Row(u(r_hi + 1)?)),
            cols: Span::from_parts(Col(u(c_lo)?), Col(u(c_hi + 1)?)),
        })
    }

    pub(crate) fn n_rows(&self) -> usize {
        self.rows.width() as usize
    }

    pub(crate) fn n_cols(&self) -> usize {
        self.cols.width() as usize
    }

    /// Spec walker rule 2: one formula, closed edges only at the full-grid edge.
    pub(crate) fn local_row_of(&self, y: f64) -> LocalRow {
        let n = self.n_rows();
        LocalRow(match self.grid.row_of(y) {
            Pos::Before => 0,
            Pos::After => n + 1,
            Pos::In(r) => local(i64::from(r.0) - i64::from(self.rows.start().0) + 1, n),
        })
    }

    pub(crate) fn local_col_of(&self, x: f64) -> LocalCol {
        let n = self.n_cols();
        LocalCol(match self.grid.col_of(x) {
            Pos::Before => 0,
            Pos::After => n + 1,
            Pos::In(c) => local(i64::from(c.0) - i64::from(self.cols.start().0) + 1, n),
        })
    }

    /// exactextract `grid_cell(const Grid<infinite_extent>&, ...)` (src/grid.cpp) with
    /// lattice edges; padding boxes reach `f64::MIN` / `f64::MAX`.
    pub(crate) fn cell_box(&self, r: LocalRow, c: LocalCol) -> CellBox {
        let l = self.grid.lattice();
        let (c0, r0) = (i64::from(self.cols.start().0), i64::from(self.rows.start().0));
        let (ci, ri) = (to_i64(c.0), to_i64(r.0));
        CellBox {
            xmin: if c.0 == 0 { f64::MIN } else { l.x_edge(c0 + ci - 1) },
            xmax: if c.0 > self.n_cols() { f64::MAX } else { l.x_edge(c0 + ci) },
            ymax: if r.0 == 0 { f64::MAX } else { l.y_edge(r0 + ri - 1) },
            ymin: if r.0 > self.n_rows() { f64::MIN } else { l.y_edge(r0 + ri) },
        }
    }

    pub(crate) fn row(&self, r: LocalRow) -> Option<Row> {
        (1..=self.n_rows())
            .contains(&r.0)
            .then(|| u32::try_from(r.0 - 1).ok())
            .flatten()
            .map(|k| Row(self.rows.start().0 + k))
    }

    pub(crate) fn col(&self, c: LocalCol) -> Option<Col> {
        (1..=self.n_cols())
            .contains(&c.0)
            .then(|| u32::try_from(c.0 - 1).ok())
            .flatten()
            .map(|k| Col(self.cols.start().0 + k))
    }

    /// controlledburn `walk_polyline` move; `None` if the walk leaves the padding ring.
    pub(crate) fn step(
        &self,
        r: LocalRow,
        c: LocalCol,
        side: Side,
    ) -> Option<(LocalRow, LocalCol)> {
        let (nr, nc) = (self.n_rows() + 1, self.n_cols() + 1);
        match side {
            Side::Top => Some((LocalRow(r.0.checked_sub(1)?), c)),
            Side::Bottom => (r.0 < nr).then_some((LocalRow(r.0 + 1), c)),
            Side::Left => Some((r, LocalCol(c.0.checked_sub(1)?))),
            Side::Right => (c.0 < nc).then_some((r, LocalCol(c.0 + 1))),
        }
    }

    /// controlledburn `snap_to_lattice` (cpp/src/scanline.cpp) on the global lattice.
    /// No range test: the subgrid range made tiles disagree (spec walker rule 1).
    pub(crate) fn snap(&self, c: Coord) -> Coord {
        let l = self.grid.lattice();
        let x = snap_axis(c.x, ((c.x - l.x0()) / l.dx()).round(), l.dx(), |k| l.x_edge(k));
        let y = snap_axis(c.y, ((l.y0() - c.y) / l.dy()).round(), l.dy(), |k| l.y_edge(k));
        Coord::new(x, y)
    }
}

fn snap_axis(v: f64, k: f64, step: f64, edge: impl Fn(i64) -> f64) -> f64 {
    if k.is_nan() || k.abs() >= 9.0e15 {
        return v; // far beyond any lattice index we can represent exactly
    }
    let e = edge(round_to_i64(k));
    if (v - e).abs() <= step * 1e-9 { e } else { v }
}

#[expect(clippy::cast_possible_truncation, reason = "|k| < 9e15 is checked by the caller")]
fn round_to_i64(k: f64) -> i64 {
    k as i64
}

fn to_i64(u: usize) -> i64 {
    i64::try_from(u).unwrap_or(i64::MAX)
}

fn local(k: i64, n: usize) -> usize {
    usize::try_from(k.clamp(0, to_i64(n) + 1)).unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g10() -> Grid {
        Grid::from_extent(0.0, 0.0, 10.0, 10.0, 10, 10).expect("grid")
    }
    fn env(xmin: f64, ymin: f64, xmax: f64, ymax: f64) -> Envelope {
        Envelope { xmin, ymin, xmax, ymax }
    }

    #[test]
    fn region_has_a_one_cell_margin_clamped_to_the_window() {
        let wg = WalkerGrid::for_envelope(&g10(), &env(2.5, 4.5, 6.5, 8.5)).expect("region");
        // cols 2..=6 plus margin -> 1..=7; rows 1..=5 plus margin -> 0..=6
        assert_eq!((wg.cols.range(), wg.rows.range()), (1..8, 0..7));
        let wg = WalkerGrid::for_envelope(&g10(), &env(-100.0, -100.0, 100.0, 100.0)).expect("all");
        assert_eq!((wg.cols.range(), wg.rows.range()), (0..10, 0..10));
    }

    #[test]
    fn degenerate_envelope_is_not_empty_and_outside_is_none() {
        assert!(WalkerGrid::for_envelope(&g10(), &env(1.0, 5.0, 9.0, 5.0)).is_some());
        assert!(WalkerGrid::for_envelope(&g10(), &env(11.0, 1.0, 12.0, 2.0)).is_none());
        assert!(WalkerGrid::for_envelope(&g10(), &env(1.0, -5.0, 2.0, -1.0)).is_none());
    }

    #[test]
    fn local_indices_and_padding() {
        let wg = WalkerGrid::for_envelope(&g10(), &env(2.5, 4.5, 6.5, 8.5)).expect("region");
        assert_eq!(wg.local_col_of(-5.0), LocalCol(0));
        assert_eq!(wg.local_col_of(0.5), LocalCol(0)); // left of the region
        assert_eq!(wg.local_col_of(1.5), LocalCol(1));
        assert_eq!(wg.local_col_of(7.5), LocalCol(7));
        assert_eq!(wg.local_col_of(8.5), LocalCol(8)); // right padding
        assert_eq!(wg.col(LocalCol(1)), Some(Col(1)));
        assert_eq!(wg.col(LocalCol(0)), None);
        assert_eq!(wg.col(LocalCol(8)), None);
        assert_eq!(wg.row(LocalRow(1)), Some(Row(0)));
    }

    #[test]
    fn cell_boxes_come_from_the_lattice() {
        let wg = WalkerGrid::for_envelope(&g10(), &env(2.5, 4.5, 6.5, 8.5)).expect("region");
        let b = wg.cell_box(LocalRow(1), LocalCol(1));
        assert_eq!((b.xmin, b.xmax, b.ymin, b.ymax), (1.0, 2.0, 9.0, 10.0));
        let p = wg.cell_box(LocalRow(1), LocalCol(0));
        assert_eq!((p.xmin, p.xmax), (f64::MIN, 1.0));
        let p = wg.cell_box(LocalRow(1), LocalCol(8));
        assert_eq!((p.xmin, p.xmax), (8.0, f64::MAX));
        let p = wg.cell_box(LocalRow(0), LocalCol(1));
        assert_eq!((p.ymin, p.ymax), (10.0, f64::MAX));
    }

    #[test]
    fn step_moves_and_stops_at_the_padding_ring() {
        let wg = WalkerGrid::for_envelope(&g10(), &env(2.5, 4.5, 6.5, 8.5)).expect("region");
        assert_eq!(wg.step(LocalRow(2), LocalCol(2), Side::Top), Some((LocalRow(1), LocalCol(2))));
        assert_eq!(
            wg.step(LocalRow(2), LocalCol(2), Side::Right),
            Some((LocalRow(2), LocalCol(3)))
        );
        assert_eq!(wg.step(LocalRow(0), LocalCol(2), Side::Top), None);
        assert_eq!(wg.step(LocalRow(2), LocalCol(8), Side::Right), None);
    }

    #[test]
    fn snap_moves_near_edge_vertices_onto_the_lattice() {
        let g = Grid::from_extent(0.0, 0.0, 10.0, 10.0, 100, 100).expect("grid");
        let wg = WalkerGrid::for_envelope(&g, &env(0.0, 0.0, 5.0, 5.0)).expect("region");
        let edge = g.lattice().x_edge(3); // 3 * 0.1 = 0.30000000000000004
        let s = wg.snap(Coord::new(0.3, 0.55));
        assert_eq!(s.x, edge);
        assert_eq!(s.y, 0.55);
        // Snap has no range test: an edge far outside the region still snaps.
        let far = g.lattice().x_edge(90);
        assert_eq!(wg.snap(Coord::new(far + 1e-12, 0.55)).x, far);
    }
}
