//! Lattice, grid window and index types.
//!
//! Index rules follow exactextract (Apache-2.0, Copyright 2018-2022 ISciences, LLC),
//! commit 67969d99d1dfe1115be955087e2c48e1aab133e2: `src/grid.h` (`Grid::get_column`,
//! `Grid::get_row`). Cell edges come from the lattice origin, not from a window corner.

use std::fmt;
use std::ops::Range;

/// Global lattice row. Row 0 is at the top.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Row(pub u32);

/// Global lattice column.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Col(pub u32);

/// Caller-supplied geometry id. Ids do not have to be unique.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GeomId(pub u32);

impl From<u32> for GeomId {
    fn from(v: u32) -> Self {
        Self(v)
    }
}
impl From<Row> for u32 {
    fn from(v: Row) -> Self {
        v.0
    }
}
impl From<Col> for u32 {
    fn from(v: Col) -> Self {
        v.0
    }
}
impl fmt::Display for GeomId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Non-empty half-open index range. `start < end` by construction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Span<T> {
    start: T,
    end: T,
}

impl<T: Copy + Ord + Into<u32>> Span<T> {
    /// Errors if `start >= end`.
    pub fn new(start: T, end: T) -> Result<Self, GridError> {
        if start < end {
            Ok(Self { start, end })
        } else {
            Err(GridError::EmptySpan { start: start.into(), end: end.into() })
        }
    }
    /// For code that already holds the invariant.
    #[expect(dead_code, reason = "first used by a later task")]
    pub(crate) const fn from_parts(start: T, end: T) -> Self {
        Self { start, end }
    }
    #[must_use]
    pub fn start(&self) -> T {
        self.start
    }
    #[must_use]
    pub fn end(&self) -> T {
        self.end
    }
    #[must_use]
    pub fn width(&self) -> u32 {
        self.end.into() - self.start.into()
    }
    #[must_use]
    pub fn contains(&self, t: T) -> bool {
        self.start <= t && t < self.end
    }
    #[must_use]
    pub fn contains_span(&self, o: Span<T>) -> bool {
        self.start <= o.start && o.end <= self.end
    }
    #[must_use]
    pub fn range(&self) -> Range<u32> {
        self.start.into()..self.end.into()
    }
}

/// Position of a coordinate against the full grid, on one axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pos<T> {
    Before,
    In(T),
    After,
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum GridError {
    #[error("non-finite lattice value: {0}")]
    NonFinite(&'static str),
    #[error("resolution must be > 0: dx = {dx}, dy = {dy}")]
    Resolution { dx: f64, dy: f64 },
    #[error("empty span: start {start}, end {end}")]
    EmptySpan { start: u32, end: u32 },
    #[error("window rows {rows:?}, cols {cols:?} is outside the full grid")]
    Outside { rows: Span<Row>, cols: Span<Col> },
    #[error("extent is invalid: {0}")]
    Extent(&'static str),
}

/// Infinite regular lattice. Column edge k is `x0 + k*dx`; row edge k is `y0 - k*dy`.
/// Row 0 is at the top. Cell (r, c) is `[x_edge(c), x_edge(c+1)] x [y_edge(r+1), y_edge(r)]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lattice {
    x0: f64,
    y0: f64,
    dx: f64,
    dy: f64,
}

#[expect(clippy::cast_precision_loss, reason = "lattice indices stay below 2^53")]
pub(crate) fn to_f64(i: i64) -> f64 {
    i as f64
}

impl Lattice {
    /// Errors: non-finite value, `dx <= 0`, `dy <= 0`.
    pub fn new(x0: f64, y0: f64, dx: f64, dy: f64) -> Result<Self, GridError> {
        for (name, v) in [("x0", x0), ("y0", y0), ("dx", dx), ("dy", dy)] {
            if !v.is_finite() {
                return Err(GridError::NonFinite(name));
            }
        }
        if !(dx > 0.0 && dy > 0.0) {
            return Err(GridError::Resolution { dx, dy });
        }
        Ok(Self { x0, y0, dx, dy })
    }
    #[must_use]
    pub fn x0(&self) -> f64 {
        self.x0
    }
    #[must_use]
    pub fn y0(&self) -> f64 {
        self.y0
    }
    #[must_use]
    pub fn dx(&self) -> f64 {
        self.dx
    }
    #[must_use]
    pub fn dy(&self) -> f64 {
        self.dy
    }
    #[must_use]
    pub fn x_edge(&self, col: i64) -> f64 {
        self.x0 + to_f64(col) * self.dx
    }
    #[must_use]
    pub fn y_edge(&self, row: i64) -> f64 {
        self.y0 - to_f64(row) * self.dy
    }
    /// Column centre, controlledburn formula `xmin + (col + 0.5) * dx`.
    #[expect(dead_code, reason = "first used by a later task")]
    pub(crate) fn x_mid(&self, col: i64) -> f64 {
        self.x0 + (to_f64(col) + 0.5) * self.dx
    }
    /// Row centre, controlledburn formula `ymax - (r + 0.5) * dy`.
    #[expect(dead_code, reason = "first used by a later task")]
    pub(crate) fn y_mid(&self, row: i64) -> f64 {
        self.y0 - (to_f64(row) + 0.5) * self.dy
    }
}

/// A window of a full extent on a lattice. All spans are global lattice indices.
/// `window()` keeps the full extent, so the edge rules of the full grid hold in every chunk.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid {
    lattice: Lattice,
    full_rows: Span<Row>,
    full_cols: Span<Col>,
    rows: Span<Row>,
    cols: Span<Col>,
}

impl Grid {
    /// Full grid = window = `rows` x `cols`.
    /// Do not use this for a chunk: a chunk made here closes its own right and bottom
    /// edges, so a point on a seam goes into two chunks. Use [`Grid::window`] for chunks.
    pub fn new(lattice: Lattice, rows: Span<Row>, cols: Span<Col>) -> Result<Self, GridError> {
        let right = lattice.x_edge(i64::from(cols.end().0));
        let bottom = lattice.y_edge(i64::from(rows.end().0));
        if !(right.is_finite() && bottom.is_finite()) {
            return Err(GridError::Extent("the far grid edge is not finite"));
        }
        Ok(Self { lattice, full_rows: rows, full_cols: cols, rows, cols })
    }
    /// Same meaning as controlledburn `GridSpec`: x0 = xmin, y0 = ymax,
    /// dx = (xmax - xmin) / ncol, dy = (ymax - ymin) / nrow.
    #[expect(clippy::neg_cmp_op_on_partial_ord, reason = "the negated form also rejects NaN")]
    pub fn from_extent(
        xmin: f64,
        ymin: f64,
        xmax: f64,
        ymax: f64,
        ncol: u32,
        nrow: u32,
    ) -> Result<Self, GridError> {
        for (name, v) in [("xmin", xmin), ("ymin", ymin), ("xmax", xmax), ("ymax", ymax)] {
            if !v.is_finite() {
                return Err(GridError::NonFinite(name));
            }
        }
        if !(xmax > xmin) {
            return Err(GridError::Extent("xmax must be > xmin"));
        }
        if !(ymax > ymin) {
            return Err(GridError::Extent("ymax must be > ymin"));
        }
        let rows = Span::new(Row(0), Row(nrow))?;
        let cols = Span::new(Col(0), Col(ncol))?;
        let lattice = Lattice::new(
            xmin,
            ymax,
            (xmax - xmin) / f64::from(ncol),
            (ymax - ymin) / f64::from(nrow),
        )?;
        Self::new(lattice, rows, cols)
    }
    /// A sub-window with the same lattice and the same full grid.
    pub fn window(&self, rows: Span<Row>, cols: Span<Col>) -> Result<Grid, GridError> {
        if !(self.full_rows.contains_span(rows) && self.full_cols.contains_span(cols)) {
            return Err(GridError::Outside { rows, cols });
        }
        Ok(Self { rows, cols, ..*self })
    }
    #[must_use]
    pub fn lattice(&self) -> &Lattice {
        &self.lattice
    }
    #[must_use]
    pub fn rows(&self) -> Span<Row> {
        self.rows
    }
    #[must_use]
    pub fn cols(&self) -> Span<Col> {
        self.cols
    }
    #[must_use]
    pub fn full_rows(&self) -> Span<Row> {
        self.full_rows
    }
    #[must_use]
    pub fn full_cols(&self) -> Span<Col> {
        self.full_cols
    }
    /// exactextract `Grid<infinite_extent>::get_column` on the lattice; padding is
    /// `Before` / `After`. The cell edges are the authority.
    #[must_use]
    #[expect(clippy::neg_cmp_op_on_partial_ord, reason = "the negated form also rejects NaN")]
    pub fn col_of(&self, x: f64) -> Pos<Col> {
        let l = &self.lattice;
        let (s, e) = (i64::from(self.full_cols.start.0), i64::from(self.full_cols.end.0));
        if !(x >= l.x_edge(s)) {
            return Pos::Before; // also NaN
        }
        let right = l.x_edge(e);
        if x > right {
            return Pos::After;
        }
        if x == right {
            return Pos::In(Col(to_u32(e - 1))); // closed far edge of the full grid
        }
        let mut c = floor_clamped((x - l.x0) / l.dx, s, e - 1);
        while c > s && x < l.x_edge(c) {
            c -= 1;
        }
        while c + 1 < e && x >= l.x_edge(c + 1) {
            c += 1;
        }
        Pos::In(Col(to_u32(c)))
    }

    /// exactextract `Grid<infinite_extent>::get_row` on the lattice. Cell r is
    /// `y_edge(r + 1) < y <= y_edge(r)`.
    #[must_use]
    #[expect(clippy::neg_cmp_op_on_partial_ord, reason = "the negated form also rejects NaN")]
    pub fn row_of(&self, y: f64) -> Pos<Row> {
        let l = &self.lattice;
        let (s, e) = (i64::from(self.full_rows.start.0), i64::from(self.full_rows.end.0));
        if !(y <= l.y_edge(s)) {
            return Pos::Before; // above the grid, or NaN
        }
        let bottom = l.y_edge(e);
        if y < bottom {
            return Pos::After;
        }
        if y == bottom {
            return Pos::In(Row(to_u32(e - 1))); // closed far edge of the full grid
        }
        let mut r = floor_clamped((l.y0 - y) / l.dy, s, e - 1);
        while r > s && y > l.y_edge(r) {
            r -= 1;
        }
        while r + 1 < e && y <= l.y_edge(r + 1) {
            r += 1;
        }
        Pos::In(Row(to_u32(r)))
    }
}

#[expect(clippy::cast_possible_truncation, reason = "clamped to [lo, hi] before the cast")]
fn floor_clamped(v: f64, lo: i64, hi: i64) -> i64 {
    v.floor().clamp(to_f64(lo), to_f64(hi)) as i64
}

/// Callers pass values inside a `Span<u32>`, so the fallback never applies.
fn to_u32(i: i64) -> u32 {
    u32::try_from(i).unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g10() -> Grid {
        Grid::from_extent(0.0, 0.0, 10.0, 10.0, 10, 10).expect("grid")
    }

    #[test]
    fn lattice_rejects_bad_values() {
        assert!(matches!(Lattice::new(f64::NAN, 0.0, 1.0, 1.0), Err(GridError::NonFinite("x0"))));
        assert!(matches!(Lattice::new(0.0, 0.0, 0.0, 1.0), Err(GridError::Resolution { .. })));
        assert!(matches!(Lattice::new(0.0, 0.0, 1.0, -1.0), Err(GridError::Resolution { .. })));
        assert!(matches!(
            Grid::from_extent(0.0, 0.0, -1.0, 10.0, 10, 10),
            Err(GridError::Extent(_))
        ));
        assert!(matches!(
            Grid::from_extent(0.0, 0.0, 10.0, 10.0, 0, 10),
            Err(GridError::EmptySpan { .. })
        ));
    }

    #[test]
    fn from_extent_matches_gridspec() {
        let g = g10();
        let l = g.lattice();
        assert_eq!((l.x_edge(0), l.x_edge(10)), (0.0, 10.0));
        assert_eq!((l.y_edge(0), l.y_edge(10)), (10.0, 0.0));
        assert_eq!(g.rows().range(), 0..10);
    }

    #[test]
    fn col_of_maps_each_edge_to_its_column() {
        let g = Grid::new(
            Lattice::new(0.0, 0.0, 0.1, 0.1).expect("lattice"),
            Span::new(Row(0), Row(1000)).expect("rows"),
            Span::new(Col(0), Col(1000)).expect("cols"),
        )
        .expect("grid");
        let l = *g.lattice();
        for k in 0..1000_i64 {
            let c = u32::try_from(k).expect("k");
            assert_eq!(g.col_of(l.x_edge(k)), Pos::In(Col(c)), "edge {k}");
            assert_eq!(g.col_of(l.x_edge(k) + 0.05), Pos::In(Col(c)), "inside {k}");
            if k > 0 {
                assert_eq!(g.col_of(l.x_edge(k).next_down()), Pos::In(Col(c - 1)), "below {k}");
            }
            assert_eq!(g.row_of(l.y_edge(k)), Pos::In(Row(c)), "row edge {k}");
            if k > 0 {
                assert_eq!(g.row_of(l.y_edge(k).next_up()), Pos::In(Row(c - 1)), "row above {k}");
            }
        }
    }

    #[test]
    fn closed_far_edges_and_outside() {
        let g = g10();
        assert_eq!(g.col_of(10.0), Pos::In(Col(9)));
        assert_eq!(g.row_of(0.0), Pos::In(Row(9)));
        assert_eq!(g.col_of(0.0), Pos::In(Col(0)));
        assert_eq!(g.row_of(10.0), Pos::In(Row(0)));
        assert_eq!(g.col_of(-0.1), Pos::Before);
        assert_eq!(g.col_of(10.1), Pos::After);
        assert_eq!(g.row_of(10.1), Pos::Before);
        assert_eq!(g.row_of(-0.1), Pos::After);
        assert_eq!(g.col_of(f64::NAN), Pos::Before);
    }

    #[test]
    fn large_finite_coordinates_do_not_overflow() {
        let g = g10();
        assert_eq!(g.col_of(1e300), Pos::After);
        assert_eq!(g.col_of(-1e300), Pos::Before);
        assert_eq!(g.row_of(1e300), Pos::Before);
        assert_eq!(g.row_of(-1e300), Pos::After);
    }

    #[test]
    fn window_keeps_full_grid_edges() {
        let g = g10();
        let w = g
            .window(Span::new(Row(0), Row(5)).expect("r"), Span::new(Col(0), Col(5)).expect("c"))
            .expect("window");
        assert_eq!(w.cols().range(), 0..5);
        assert_eq!(w.full_cols().range(), 0..10);
        // x = 5 is a seam: the column to its right, not the closed edge of the window.
        assert_eq!(w.col_of(5.0), Pos::In(Col(5)));
        assert!(matches!(
            g.window(Span::new(Row(0), Row(11)).expect("r"), Span::new(Col(0), Col(5)).expect("c")),
            Err(GridError::Outside { .. })
        ));
    }
}
