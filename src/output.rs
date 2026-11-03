//! Sparse output tables: the controlledburn four-table contract in `SoA` form.
//!
//! Contract from controlledburn (Apache-2.0, Copyright 2025 Michael Sumner),
//! commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558:
//! `cpp/include/controlledburn/output.hpp`. See the crate docs for the CSR/CSC relation.

use std::ops::Range;

use crate::coverage::fraction::to_f32;
use crate::grid::{Col, Flip, GeomId, Row, Span};

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Runs {
    pub row: Vec<u32>,
    pub col_start: Vec<u32>,
    pub col_end: Vec<u32>,
    pub id: Vec<u32>,
}
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Edges {
    pub row: Vec<u32>,
    pub col: Vec<u32>,
    pub fraction: Vec<f32>,
    pub id: Vec<u32>,
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
pub struct RunRef {
    pub row: Row,
    pub cols: Span<Col>,
    pub id: GeomId,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EdgeRef {
    pub row: Row,
    pub col: Col,
    pub fraction: f32,
    pub id: GeomId,
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

impl Runs {
    #[must_use]
    pub fn len(&self) -> usize {
        self.row.len()
    }
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.row.is_empty()
    }
    pub fn iter(&self) -> impl ExactSizeIterator<Item = RunRef> + '_ {
        self.row.iter().zip(&self.col_start).zip(&self.col_end).zip(&self.id).map(
            |(((&r, &s), &e), &id)| RunRef {
                row: Row(r),
                cols: Span::from_parts(Col(s), Col(e)),
                id: GeomId(id),
            },
        )
    }
}
impl Edges {
    #[must_use]
    pub fn len(&self) -> usize {
        self.row.len()
    }
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.row.is_empty()
    }
    pub fn iter(&self) -> impl ExactSizeIterator<Item = EdgeRef> + '_ {
        self.row.iter().zip(&self.col).zip(&self.fraction).zip(&self.id).map(
            |(((&r, &c), &f), &id)| EdgeRef {
                row: Row(r),
                col: Col(c),
                fraction: f,
                id: GeomId(id),
            },
        )
    }
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
    pub runs: Range<usize>,
    pub edges: Range<usize>,
    pub lines: Range<usize>,
    pub points: Range<usize>,
}

/// Table lengths at the start of one geometry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct GeomStart {
    runs: usize,
    edges: usize,
    lines: usize,
    points: usize,
}

/// One table record. The only way to write the tables.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Emit {
    Run { row: Row, cols: Span<Col> },
    Edge { row: Row, col: Col, fraction: f32 },
    Line { row: Row, col: Col, length: f32 },
    Point { row: Row, col: Col },
}

/// The result of a burn: the sparse coverage matrix `A[cell, geometry]`.
///
/// The value is the fraction of the cell that the geometry covers (1 for interior cells).
/// Zonal statistics are `(Aᵀv) / (Aᵀ1)` for a raster `v`.
///
/// - `runs` is a run-length layer: one interval of columns in one row, all with value 1
///   ("RLE-CSR").
/// - `edges` is COO: explicit `(row, col, value)` triples.
/// - The burn order is geometry-major, so `geom_ptr` is a CSC index with the geometry
///   as the column.
/// - `to_row_major()` builds the CSR index.
///
/// Two geometries can cover the same cell, so the full object is a sparse `(row, col, id)`
/// tensor. One geometry has at most one record per cell.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BurnResult {
    pub runs: Runs,
    pub edges: Edges,
    pub lines: Lines,
    pub points: Points,
    pub geom_ptr: Vec<GeomIndex>,
}

/// CSR view. For row r, `runs[run_ptr[r - r0]..run_ptr[r - r0 + 1]]` are indices into
/// `BurnResult::runs`; the same for edges, lines and points. Within a row, the entries
/// keep burn order; they are not sorted by column.
#[derive(Clone, Debug, PartialEq)]
pub struct RowMajor {
    pub rows: Span<Row>,
    pub run_ptr: Vec<usize>,
    pub runs: Vec<usize>,
    pub edge_ptr: Vec<usize>,
    pub edges: Vec<usize>,
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
        let Self { runs, edges, lines, points, geom_ptr } = self;
        for v in [&mut runs.row, &mut runs.col_start, &mut runs.col_end, &mut runs.id] {
            v.clear();
        }
        for v in [&mut edges.row, &mut edges.col, &mut edges.id] {
            v.clear();
        }
        edges.fraction.clear();
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
        self.runs.len() + self.edges.len() + self.lines.len() + self.points.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub(crate) fn push(&mut self, id: GeomId, e: Emit) {
        match e {
            Emit::Run { row, cols } => {
                self.runs.row.push(row.0);
                self.runs.col_start.push(cols.start().0);
                self.runs.col_end.push(cols.end().0);
                self.runs.id.push(id.0);
            }
            Emit::Edge { row, col, fraction } => {
                self.edges.row.push(row.0);
                self.edges.col.push(col.0);
                self.edges.fraction.push(fraction);
                self.edges.id.push(id.0);
            }
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
        GeomStart {
            runs: self.runs.len(),
            edges: self.edges.len(),
            lines: self.lines.len(),
            points: self.points.len(),
        }
    }

    /// Pushed only after the geometry is complete (spec: state of `out` after an error).
    pub(crate) fn end_geom(&mut self, id: GeomId, s: GeomStart) {
        self.geom_ptr.push(GeomIndex {
            id,
            runs: s.runs..self.runs.len(),
            edges: s.edges..self.edges.len(),
            lines: s.lines..self.lines.len(),
            points: s.points..self.points.len(),
        });
    }

    /// Maps the records after `s` to other indices, in place.
    pub(crate) fn map_tail(&mut self, s: GeomStart, rows: Flip, cols: Flip) {
        let r = &mut self.runs;
        for row in r.row.get_mut(s.runs..).into_iter().flatten() {
            *row = rows.index(*row);
        }
        let ends = r.col_start.iter_mut().zip(r.col_end.iter_mut()).skip(s.runs);
        for (a, b) in ends {
            (*a, *b) = cols.range(*a, *b);
        }
        let cells = [
            (&mut self.edges.row, &mut self.edges.col, s.edges),
            (&mut self.lines.row, &mut self.lines.col, s.lines),
            (&mut self.points.row, &mut self.points.col, s.points),
        ];
        for (row, col, from) in cells {
            for (r, c) in row.iter_mut().zip(col.iter_mut()).skip(from) {
                (*r, *c) = (rows.index(*r), cols.index(*c));
            }
        }
    }

    /// Makes the records after `s` one record per cell, when more than one burn of one
    /// geometry wrote them: runs join, edge fractions and line lengths add.
    /// An edge in a run cell goes: the run already holds the full cell.
    pub(crate) fn merge_tail(&mut self, id: GeomId, s: GeomStart) {
        let r = &self.runs;
        let mut runs: Vec<(u32, u32, u32)> =
            (s.runs..r.len()).map(|i| (r.row[i], r.col_start[i], r.col_end[i])).collect();
        let e = &self.edges;
        let mut edges: Vec<(u32, u32, f32)> =
            (s.edges..e.len()).map(|i| (e.row[i], e.col[i], e.fraction[i])).collect();
        let l = &self.lines;
        let mut lines: Vec<(u32, u32, f32)> =
            (s.lines..l.len()).map(|i| (l.row[i], l.col[i], l.length[i])).collect();
        self.truncate(s);

        runs.sort_unstable();
        let mut joined: Vec<(u32, u32, u32)> = Vec::with_capacity(runs.len());
        for (row, a, b) in runs {
            match joined.last_mut() {
                Some(last) if last.0 == row && a <= last.2 => last.2 = last.2.max(b),
                _ => joined.push((row, a, b)),
            }
        }
        for &(row, a, b) in &joined {
            self.push(id, Emit::Run { row: Row(row), cols: Span::from_parts(Col(a), Col(b)) });
        }
        let in_run = |row: u32, col: u32| {
            let i = joined.partition_point(|&(r, a, _)| (r, a) <= (row, col));
            i > 0 && joined.get(i - 1).is_some_and(|&(r, _, b)| r == row && col < b)
        };
        // Stable sorts: the sums keep the burn order of the copies.
        edges.sort_by_key(|&(r, c, _)| (r, c));
        for g in edges.chunk_by(|a, b| (a.0, a.1) == (b.0, b.1)) {
            let Some(&(row, col, _)) = g.first() else { continue };
            if !in_run(row, col) {
                let f = g.iter().map(|x| f64::from(x.2)).sum::<f64>().min(1.0);
                self.push(id, Emit::Edge { row: Row(row), col: Col(col), fraction: to_f32(f) });
            }
        }
        lines.sort_by_key(|&(r, c, _)| (r, c));
        for g in lines.chunk_by(|a, b| (a.0, a.1) == (b.0, b.1)) {
            let Some(&(row, col, _)) = g.first() else { continue };
            let length = to_f32(g.iter().map(|x| f64::from(x.2)).sum::<f64>());
            self.push(id, Emit::Line { row: Row(row), col: Col(col), length });
        }
    }

    /// Drops the run, edge and line records after `s`.
    fn truncate(&mut self, s: GeomStart) {
        let Self { runs, edges, lines, .. } = self;
        for v in [&mut runs.row, &mut runs.col_start, &mut runs.col_end, &mut runs.id] {
            v.truncate(s.runs);
        }
        for v in [&mut edges.row, &mut edges.col, &mut edges.id] {
            v.truncate(s.edges);
        }
        edges.fraction.truncate(s.edges);
        for v in [&mut lines.row, &mut lines.col, &mut lines.id] {
            v.truncate(s.lines);
        }
        lines.length.truncate(s.lines);
    }

    /// Row-major (CSR) index over `rows`.
    #[must_use = "the index is the only result"]
    pub fn to_row_major(&self, rows: Span<Row>) -> Result<RowMajor, RowMajorError> {
        let (run_ptr, runs) = csr(&self.runs.row, rows)?;
        let (edge_ptr, edges) = csr(&self.edges.row, rows)?;
        let (line_ptr, lines) = csr(&self.lines.row, rows)?;
        let (point_ptr, points) = csr(&self.points.row, rows)?;
        Ok(RowMajor { rows, run_ptr, runs, edge_ptr, edges, line_ptr, lines, point_ptr, points })
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

    fn span(a: u32, b: u32) -> Span<Col> {
        Span::new(Col(a), Col(b)).expect("span")
    }

    fn sample() -> BurnResult {
        let mut r = BurnResult::default();
        let s = r.start_geom();
        r.push(GeomId(7), Emit::Run { row: Row(2), cols: span(1, 4) });
        r.push(GeomId(7), Emit::Edge { row: Row(1), col: Col(3), fraction: 0.5 });
        r.push(GeomId(7), Emit::Run { row: Row(1), cols: span(4, 6) });
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
        assert_eq!(r.len(), 5);
        let runs: Vec<RunRef> = r.runs.iter().collect();
        assert_eq!(runs[0], RunRef { row: Row(2), cols: span(1, 4), id: GeomId(7) });
        assert_eq!(r.edges.iter().next().map(|e| e.fraction), Some(0.5));
        assert_eq!(r.lines.iter().next().map(|l| l.length), Some(1.5));
        assert_eq!(r.points.iter().next().map(|p| p.col), Some(Col(9)));
    }

    #[test]
    fn geom_ptr_is_the_csc_index() {
        let r = sample();
        assert_eq!(r.geom_ptr.len(), 2);
        assert_eq!(r.geom_ptr[0].id, GeomId(7));
        assert_eq!(r.geom_ptr[0].runs, 0..2);
        assert_eq!(r.geom_ptr[0].edges, 0..1);
        assert_eq!(r.geom_ptr[1].runs, 2..2);
        assert_eq!(r.geom_ptr[1].lines, 0..1);
    }

    #[test]
    fn row_major_is_a_stable_counting_sort() {
        let r = sample();
        let rows = Span::new(Row(0), Row(3)).expect("rows");
        let m = r.to_row_major(rows).expect("csr");
        assert_eq!(m.run_ptr, vec![0, 0, 1, 2]);
        assert_eq!(m.runs, vec![1, 0]); // row 1 holds run 1, row 2 holds run 0
        assert_eq!(m.edge_ptr, vec![0, 0, 1, 1]);
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
        let cap = r.runs.row.capacity();
        r.clear();
        assert!(r.is_empty() && r.geom_ptr.is_empty());
        assert_eq!(r.runs.row.capacity(), cap);
    }
}
