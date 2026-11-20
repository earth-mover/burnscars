//! Output: the polygon stretch table, and the line and point tables in `SoA` form.
//!
//! The line and point tables follow controlledburn (Apache-2.0, Copyright 2025 Michael Sumner),
//! commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558: `cpp/include/controlledburn/output.hpp`.
//! The polygon output is `Stretches`, which replaces its runs and edges tables.

use std::cmp::Reverse;
use std::ops::Range;

use crate::burner::Rule;
use crate::coverage::fraction::{CellClass, to_f32};
use crate::error::BurnError;
use crate::grid::{Col, Flip, GeomId, Row, Span};

/// Per fiber (one geometry, one row), the burnt cells sorted by column, as stretches of cells.
/// A stretch with one value is a run; a stretch with one value per cell is a block.
/// A one-cell stretch is both.
///
/// It is the sparse coverage matrix `W[geometry, cell]` in a CSR form whose row is the fiber.
/// Zonal statistics are `(W v) / (W 1)` for a raster `v`.
#[derive(Clone, Debug, PartialEq)]
pub struct Stretches {
    /// Fiber `f` is geometry `id[f]` in row `row[f]`. Fibers are strictly sorted by `(id, row)`.
    pub id: Vec<u32>,
    pub row: Vec<u32>,
    /// The stretches of fiber `f` are `fiber_ptr[f]..fiber_ptr[f + 1]`. `[0]` when empty.
    pub fiber_ptr: Vec<usize>,
    /// Stretch `s` covers columns `col_start[s]..col_end[s]`. Within a fiber, sorted and not overlapping.
    pub col_start: Vec<u32>,
    pub col_end: Vec<u32>,
    /// `None` under `Rule::Touched` and `Rule::Centre`.
    pub values: Option<StretchValues>,
}

/// The weights of the stretches. Stretch `s` holds `weight[value_ptr[s]..value_ptr[s + 1]]`:
/// one value (a run) or one value per cell (a block).
#[derive(Clone, Debug, PartialEq)]
pub struct StretchValues {
    /// One entry per stretch, plus one. `[0]` when empty.
    pub value_ptr: Vec<u32>,
    pub weight: Vec<f32>,
}

/// The capacity to reserve in a [`Stretches`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StretchCounts {
    pub fibers: usize,
    pub stretches: usize,
    /// Ignored for a rule with no values.
    pub values: usize,
}

impl Default for StretchValues {
    fn default() -> Self {
        Self { value_ptr: vec![0], weight: Vec::new() }
    }
}

/// Column lengths at the start of one geometry, to roll back a failed encode.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SpansMark {
    fibers: usize,
    stretches: usize,
    values: usize,
}

impl Stretches {
    /// An empty table for `rule`: `values` is `Some` exactly under `Rule::Fraction`.
    #[must_use]
    pub fn new(rule: Rule) -> Self {
        Self::with_capacity(rule, StretchCounts::default())
    }

    /// [`Stretches::new`] with every column reserved.
    #[must_use]
    pub fn with_capacity(rule: Rule, n: StretchCounts) -> Self {
        let mut fiber_ptr = Vec::with_capacity(n.fibers + 1);
        fiber_ptr.push(0);
        let values = rule.has_values().then(|| {
            let mut value_ptr = Vec::with_capacity(n.stretches + 1);
            value_ptr.push(0);
            StretchValues { value_ptr, weight: Vec::with_capacity(n.values) }
        });
        Self {
            id: Vec::with_capacity(n.fibers),
            row: Vec::with_capacity(n.fibers),
            fiber_ptr,
            col_start: Vec::with_capacity(n.stretches),
            col_end: Vec::with_capacity(n.stretches),
            values,
        }
    }

    /// Empties every column and keeps the capacity and the rule.
    pub fn clear(&mut self) {
        self.truncate(SpansMark { fibers: 0, stretches: 0, values: 0 });
    }

    /// The number of stretches, not of fibers.
    #[must_use]
    pub fn stretches(&self) -> usize {
        self.col_start.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.col_start.is_empty()
    }

    pub(crate) fn mark(&self) -> SpansMark {
        SpansMark {
            fibers: self.id.len(),
            stretches: self.col_start.len(),
            values: self.values.as_ref().map_or(0, |v| v.weight.len()),
        }
    }

    pub(crate) fn truncate(&mut self, m: SpansMark) {
        self.id.truncate(m.fibers);
        self.row.truncate(m.fibers);
        self.fiber_ptr.truncate(m.fibers + 1);
        self.fiber_ptr.resize(m.fibers + 1, 0);
        self.col_start.truncate(m.stretches);
        self.col_end.truncate(m.stretches);
        if let Some(v) = &mut self.values {
            v.value_ptr.truncate(m.stretches + 1);
            v.value_ptr.resize(m.stretches + 1, 0);
            v.weight.truncate(m.values);
        }
    }

    /// Refuses a rule whose value form is not the one of `self`. It never changes `self`.
    pub(crate) fn check_rule(&self, rule: Rule) -> Result<(), BurnError> {
        if self.values.is_some() == rule.has_values() {
            Ok(())
        } else {
            Err(BurnError::RuleMismatch)
        }
    }

    /// `(id, row, col) -> weight`, 1.0 with no values (test only).
    /// Panics on a cell that two stretches cover.
    #[cfg(test)]
    pub(crate) fn densify(&self) -> std::collections::BTreeMap<(u32, u32, u32), f32> {
        let mut m = std::collections::BTreeMap::new();
        for f in 0..self.id.len() {
            for s in self.fiber_ptr[f]..self.fiber_ptr[f + 1] {
                let cols = self.col_start[s]..self.col_end[s];
                let w = |i: usize| {
                    self.values.as_ref().map_or(1.0, |v| {
                        let vs = &v.weight[v.value_ptr[s] as usize..v.value_ptr[s + 1] as usize];
                        if vs.len() == 1 { vs[0] } else { vs[i] }
                    })
                };
                for (i, c) in cols.enumerate() {
                    if m.insert((self.id[f], self.row[f], c), w(i)).is_some() {
                        panic!("cell {:?} twice", (self.id[f], self.row[f], c));
                    }
                }
            }
        }
        m
    }
}

/// One polygon record of one geometry: columns `start..end` of `row`.
/// A run has weight 1.0; an edge has one column and its fraction. 16 B.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Rec {
    pub(crate) row: u32,
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) weight: f32,
}

impl Rec {
    pub(crate) fn run(row: Row, cols: Span<Col>) -> Self {
        Self { row: row.0, start: cols.start().0, end: cols.end().0, weight: 1.0 }
    }

    pub(crate) fn edge(row: Row, col: Col, fraction: f32) -> Self {
        // A full cell is a run; `is_edge` reads the weight.
        debug_assert!(fraction < 1.0, "edge weight {fraction}");
        Self { row: row.0, start: col.0, end: col.0 + 1, weight: fraction }
    }

    fn is_edge(&self) -> bool {
        self.weight < 1.0
    }
}

/// The polygon records of one geometry, before the encode. Reused across geometries.
#[derive(Debug, Default)]
pub(crate) struct PolyRecs {
    recs: Vec<Rec>,
}

impl PolyRecs {
    pub(crate) fn clear(&mut self) {
        self.recs.clear();
    }

    pub(crate) fn push(&mut self, r: Rec) {
        self.recs.push(r);
    }

    pub(crate) fn bytes(&self) -> usize {
        self.recs.capacity() * size_of::<Rec>()
    }

    pub(crate) fn reserve(&mut self, n: usize) {
        self.recs.reserve(n);
    }

    /// Maps the rows to the caller's rows, in place.
    pub(crate) fn map_rows(&mut self, rows: Flip) {
        if rows.is_flip() {
            self.recs.iter_mut().for_each(|r| r.row = rows.index(r.row));
        }
    }

    /// Merges the records into one per cell and appends the fibers of `id` to `out`,
    /// less `origin` (`(row, col)`). On error, `out` is as it was before.
    ///
    /// The merge, in sorted order: overlapping or touching runs join; an edge inside a run
    /// goes; the edges of one cell add in `f64`, at most 1, and a full sum joins the runs.
    pub(crate) fn encode(
        &mut self,
        id: GeomId,
        origin: (u32, u32),
        out: &mut Stretches,
    ) -> Result<(), BurnError> {
        // Runs before edges at one start: an edge then meets every run that can hold it.
        let key = |r: &Rec| (r.row, r.start, r.is_edge());
        let same_row = |a: &Rec, b: &Rec| a.row == b.row;
        let mark = out.mark();
        let mut fiber = |row: &[Rec]| Fiber::new(id, origin, out).encode(row);
        let r = if self.recs.is_sorted_by_key(key) {
            self.recs.chunk_by(same_row).try_for_each(&mut fiber)
        } else if self.recs.is_sorted_by_key(|r| (Reverse(r.row), r.start, r.is_edge())) {
            // One copy on a south-up grid: the sweep order with the rows flipped.
            self.recs.chunk_by(same_row).rev().try_for_each(&mut fiber)
        } else {
            // Two or more copies. The weight bits make the order total, so the sums do not
            // depend on the burn order of the copies.
            self.recs.sort_unstable_by_key(|r| (r.row, r.start, r.is_edge(), r.weight.to_bits()));
            self.recs.chunk_by(same_row).try_for_each(&mut fiber)
        };
        if r.is_err() {
            out.truncate(mark);
        }
        r
    }
}

/// The encoder of one fiber. At most one of `run` and `block` is open.
struct Fiber<'o> {
    id: GeomId,
    origin: (u32, u32),
    out: &'o mut Stretches,
    run: Option<(u32, u32)>,
    block: Option<(u32, u32)>,
}

impl<'o> Fiber<'o> {
    fn new(id: GeomId, origin: (u32, u32), out: &'o mut Stretches) -> Self {
        Self { id, origin, out, run: None, block: None }
    }

    /// `recs` is one row, sorted by `(start, is_edge)`.
    fn encode(mut self, recs: &[Rec]) -> Result<(), BurnError> {
        let Some(row) = recs.first().map(|r| r.row) else { return Ok(()) };
        let values = self.out.values.is_some();
        let first = self.out.col_start.len();
        for g in recs.chunk_by(|a, b| a.is_edge() && b.is_edge() && a.start == b.start) {
            let Some(&r) = g.first() else { continue };
            if !r.is_edge() {
                self.add_run(r.start, r.end)?;
                continue;
            }
            if self.run.is_some_and(|(_, end)| r.start < end) {
                continue;
            }
            let sum = g.iter().map(|e| f64::from(e.weight)).sum::<f64>().min(1.0);
            match CellClass::of(to_f32(sum)) {
                CellClass::Empty => {}
                CellClass::Edge(w) if values => self.add_edge(r.start, w)?,
                CellClass::Edge(_) | CellClass::Full => self.add_run(r.start, r.start + 1)?,
            }
        }
        self.close_run()?;
        self.close_block()?;
        if self.out.col_start.len() > first {
            let row = self.local(row, self.origin.0)?;
            self.out.id.push(self.id.0);
            self.out.row.push(row);
            self.out.fiber_ptr.push(self.out.col_start.len());
        }
        Ok(())
    }

    /// `v - origin`. A record lies in the window, so this does not underflow.
    fn local(&self, v: u32, origin: u32) -> Result<u32, BurnError> {
        v.checked_sub(origin)
            .ok_or(BurnError::Internal { id: self.id, detail: "record before the window origin" })
    }

    fn add_run(&mut self, a: u32, b: u32) -> Result<(), BurnError> {
        self.close_block()?;
        match &mut self.run {
            Some((_, end)) if a <= *end => *end = (*end).max(b),
            _ => {
                self.close_run()?;
                self.run = Some((a, b));
            }
        }
        Ok(())
    }

    fn add_edge(&mut self, c: u32, w: f32) -> Result<(), BurnError> {
        self.close_run()?;
        match &mut self.block {
            Some((_, end)) if *end == c => *end += 1,
            _ => {
                self.close_block()?;
                self.block = Some((c, c + 1));
            }
        }
        if let Some(v) = &mut self.out.values {
            v.weight.push(w);
        }
        Ok(())
    }

    fn close_run(&mut self) -> Result<(), BurnError> {
        let Some((a, b)) = self.run.take() else { return Ok(()) };
        if let Some(v) = &mut self.out.values {
            v.weight.push(1.0);
        }
        self.push_stretch(a, b)
    }

    /// The block's weights are already in `weight`.
    fn close_block(&mut self) -> Result<(), BurnError> {
        let Some((a, b)) = self.block.take() else { return Ok(()) };
        self.push_stretch(a, b)
    }

    fn push_stretch(&mut self, a: u32, b: u32) -> Result<(), BurnError> {
        let (a, b) = (self.local(a, self.origin.1)?, self.local(b, self.origin.1)?);
        self.out.col_start.push(a);
        self.out.col_end.push(b);
        if let Some(v) = &mut self.out.values {
            let n = u32::try_from(v.weight.len())
                .map_err(|_| BurnError::TooManyValues { id: self.id })?;
            v.value_ptr.push(n);
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Lines {
    pub row: Vec<u32>,
    pub col: Vec<u32>,
    pub length: Vec<f32>,
    pub id: Vec<u32>,
}
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Points {
    pub row: Vec<u32>,
    pub col: Vec<u32>,
    pub id: Vec<u32>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LineRef {
    pub row: Row,
    pub col: Col,
    pub length: f32,
    pub id: GeomId,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointRef {
    pub row: Row,
    pub col: Col,
    pub id: GeomId,
}

impl Lines {
    #[must_use]
    pub fn len(&self) -> usize {
        self.row.len()
    }
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.row.is_empty()
    }
    pub fn iter(&self) -> impl ExactSizeIterator<Item = LineRef> + '_ {
        self.row.iter().zip(&self.col).zip(&self.length).zip(&self.id).map(
            |(((&r, &c), &l), &id)| LineRef { row: Row(r), col: Col(c), length: l, id: GeomId(id) },
        )
    }
}
impl Points {
    #[must_use]
    pub fn len(&self) -> usize {
        self.row.len()
    }
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.row.is_empty()
    }
    pub fn iter(&self) -> impl ExactSizeIterator<Item = PointRef> + '_ {
        self.row.iter().zip(&self.col).zip(&self.id).map(|((&r, &c), &id)| PointRef {
            row: Row(r),
            col: Col(c),
            id: GeomId(id),
        })
    }
}

/// The table rows of one input item (the CSC index). Ranges can be empty.
#[derive(Clone, Debug, PartialEq)]
pub struct GeomIndex {
    pub id: GeomId,
    pub lines: Range<usize>,
    pub points: Range<usize>,
}

/// Table lengths at the start of one geometry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct GeomStart {
    lines: usize,
    points: usize,
}

/// One table record. The only way to write the tables.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Emit {
    Line { row: Row, col: Col, length: f32 },
    Point { row: Row, col: Col },
}

/// The result of a line and point burn. Polygons go to [`Stretches`].
///
/// - `lines`: per cell, the length of the line in the cell.
/// - `points`: per point, its cell.
/// - The burn order is geometry-major, so `geom_ptr` is a CSC index with the geometry
///   as the column.
/// - `to_row_major()` builds the CSR index.
///
/// One geometry has at most one line record per cell.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BurnResult {
    pub lines: Lines,
    pub points: Points,
    pub geom_ptr: Vec<GeomIndex>,
}

/// CSR view. For row r, `lines[line_ptr[r - r0]..line_ptr[r - r0 + 1]]` are indices into
/// `BurnResult::lines`; the same for points. Within a row, the entries
/// keep burn order; they are not sorted by column.
#[derive(Clone, Debug, PartialEq)]
pub struct RowMajor {
    pub rows: Span<Row>,
    pub line_ptr: Vec<usize>,
    pub lines: Vec<usize>,
    pub point_ptr: Vec<usize>,
    pub points: Vec<usize>,
}

#[derive(Debug, thiserror::Error)]
#[error("table row {row:?} is outside rows {rows:?}")]
pub struct RowMajorError {
    pub row: Row,
    pub rows: Span<Row>,
}

impl BurnResult {
    /// Clears every table and keeps the capacity.
    pub fn clear(&mut self) {
        let Self { lines, points, geom_ptr } = self;
        for v in [&mut lines.row, &mut lines.col, &mut lines.id] {
            v.clear();
        }
        lines.length.clear();
        for v in [&mut points.row, &mut points.col, &mut points.id] {
            v.clear();
        }
        geom_ptr.clear();
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.lines.len() + self.points.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub(crate) fn push(&mut self, id: GeomId, e: Emit) {
        match e {
            Emit::Line { row, col, length } => {
                self.lines.row.push(row.0);
                self.lines.col.push(col.0);
                self.lines.length.push(length);
                self.lines.id.push(id.0);
            }
            Emit::Point { row, col } => {
                self.points.row.push(row.0);
                self.points.col.push(col.0);
                self.points.id.push(id.0);
            }
        }
    }

    pub(crate) fn start_geom(&self) -> GeomStart {
        GeomStart { lines: self.lines.len(), points: self.points.len() }
    }

    /// Pushed only after the geometry is complete (spec: state of `out` after an error).
    pub(crate) fn end_geom(&mut self, id: GeomId, s: GeomStart) {
        self.geom_ptr.push(GeomIndex {
            id,
            lines: s.lines..self.lines.len(),
            points: s.points..self.points.len(),
        });
    }

    /// Maps the rows of the records after `s` to the caller's rows, in place.
    pub(crate) fn map_tail(&mut self, s: GeomStart, rows: Flip) {
        if !rows.is_flip() {
            return;
        }
        for (row, from) in [(&mut self.lines.row, s.lines), (&mut self.points.row, s.points)] {
            for r in row.get_mut(from..).into_iter().flatten() {
                *r = rows.index(*r);
            }
        }
    }

    /// Makes the line records after `s` one record per cell, when more than one burn of one
    /// geometry wrote them: the lengths add. The polygon merge is in `PolyRecs::encode`.
    pub(crate) fn merge_tail(&mut self, id: GeomId, s: GeomStart) {
        let l = &self.lines;
        let mut lines: Vec<(u32, u32, f32)> =
            (s.lines..l.len()).map(|i| (l.row[i], l.col[i], l.length[i])).collect();
        let Self { lines: t, .. } = self;
        for v in [&mut t.row, &mut t.col, &mut t.id] {
            v.truncate(s.lines);
        }
        t.length.truncate(s.lines);
        // Stable sort: the sums keep the burn order of the copies.
        lines.sort_by_key(|&(r, c, _)| (r, c));
        for g in lines.chunk_by(|a, b| (a.0, a.1) == (b.0, b.1)) {
            let Some(&(row, col, _)) = g.first() else { continue };
            let length = to_f32(g.iter().map(|x| f64::from(x.2)).sum::<f64>());
            self.push(id, Emit::Line { row: Row(row), col: Col(col), length });
        }
    }

    /// Row-major (CSR) index over `rows`.
    #[must_use = "the index is the only result"]
    pub fn to_row_major(&self, rows: Span<Row>) -> Result<RowMajor, RowMajorError> {
        let (line_ptr, lines) = csr(&self.lines.row, rows)?;
        let (point_ptr, points) = csr(&self.points.row, rows)?;
        Ok(RowMajor { rows, line_ptr, lines, point_ptr, points })
    }
}

/// Stable counting sort of table indices by row: O(n + rows).
fn csr(row: &[u32], rows: Span<Row>) -> Result<(Vec<usize>, Vec<usize>), RowMajorError> {
    let r0 = rows.start().0;
    let mut ptr = vec![0_usize; rows.width() as usize + 1];
    for &r in row {
        if !rows.contains(Row(r)) {
            return Err(RowMajorError { row: Row(r), rows });
        }
        ptr[(r - r0) as usize + 1] += 1;
    }
    for i in 1..ptr.len() {
        ptr[i] += ptr[i - 1];
    }
    let mut next = ptr.clone();
    let mut idx = vec![0_usize; row.len()];
    for (i, &r) in row.iter().enumerate() {
        let k = (r - r0) as usize;
        idx[next[k]] = i;
        next[k] += 1;
    }
    Ok((ptr, idx))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> BurnResult {
        let mut r = BurnResult::default();
        let s = r.start_geom();
        r.push(GeomId(7), Emit::Line { row: Row(2), col: Col(1), length: 0.5 });
        r.push(GeomId(7), Emit::Line { row: Row(1), col: Col(3), length: 0.25 });
        r.end_geom(GeomId(7), s);
        let s = r.start_geom();
        r.push(GeomId(2), Emit::Line { row: Row(2), col: Col(0), length: 1.5 });
        r.push(GeomId(2), Emit::Point { row: Row(1), col: Col(9) });
        r.end_geom(GeomId(2), s);
        r
    }

    #[test]
    fn push_iter_round_trip() {
        let r = sample();
        assert_eq!(r.len(), 4);
        let lines: Vec<LineRef> = r.lines.iter().collect();
        assert_eq!(lines[0], LineRef { row: Row(2), col: Col(1), length: 0.5, id: GeomId(7) });
        assert_eq!(r.points.iter().next().map(|p| p.col), Some(Col(9)));
    }

    #[test]
    fn geom_ptr_is_the_csc_index() {
        let r = sample();
        assert_eq!(r.geom_ptr.len(), 2);
        assert_eq!(r.geom_ptr[0].id, GeomId(7));
        assert_eq!(r.geom_ptr[0].lines, 0..2);
        assert_eq!(r.geom_ptr[1].lines, 2..3);
        assert_eq!(r.geom_ptr[1].points, 0..1);
    }

    #[test]
    fn row_major_is_a_stable_counting_sort() {
        let r = sample();
        let rows = Span::new(Row(0), Row(3)).expect("rows");
        let m = r.to_row_major(rows).expect("csr");
        assert_eq!(m.line_ptr, vec![0, 0, 1, 3]);
        assert_eq!(m.lines, vec![1, 0, 2]); // row 1 holds line 1, row 2 holds lines 0 and 2
        assert_eq!(m.point_ptr, vec![0, 0, 1, 1]);
    }

    #[test]
    fn row_major_rejects_rows_outside() {
        let r = sample();
        let rows = Span::new(Row(0), Row(2)).expect("rows");
        let e = r.to_row_major(rows).expect_err("row 2 is outside");
        assert_eq!(e.row, Row(2));
    }

    #[test]
    fn clear_keeps_capacity() {
        let mut r = sample();
        let cap = r.lines.row.capacity();
        r.clear();
        assert!(r.is_empty() && r.geom_ptr.is_empty());
        assert_eq!(r.lines.row.capacity(), cap);
    }

    #[test]
    fn spans_clear_resets_to_default_and_keeps_capacity() {
        let n = StretchCounts { fibers: 4, stretches: 8, values: 16 };
        let mut s = Stretches::with_capacity(Rule::Fraction, n);
        s.id.push(1);
        s.row.push(0);
        s.col_start.push(2);
        s.col_end.push(3);
        s.fiber_ptr.push(1);
        let cap = s.col_start.capacity();
        s.clear();
        assert_eq!(s, Stretches::new(Rule::Fraction));
        assert_eq!(s.col_start.capacity(), cap);
    }
}
