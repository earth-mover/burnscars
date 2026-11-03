//! Lattice, grid window and index types.
//!
//! Index rules follow exactextract (Apache-2.0, Copyright 2018-2022 ISciences, LLC),
//! commit 67969d99d1dfe1115be955087e2c48e1aab133e2: `src/grid.h` (`Grid::get_column`,
//! `Grid::get_row`). Cell edges come from the lattice origin, not from a window corner.

use std::fmt;
use std::ops::Range;

/// Global lattice row. Row 0 is next to the `y0` edge of the lattice.
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
        self.end.into().saturating_sub(self.start.into())
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
    #[error("resolution must be > 0: {axis} = {step}")]
    Resolution { axis: &'static str, step: f64 },
    #[error("geotransform is rotated: rx = {rx}, ry = {ry}")]
    Rotated { rx: f64, ry: f64 },
    #[error("period must be finite and > 0: {0}")]
    Period(f64),
    #[error("empty span: start {start}, end {end}")]
    EmptySpan { start: u32, end: u32 },
    #[error("window rows {rows:?}, cols {cols:?} is outside the full grid")]
    Outside { rows: Span<Row>, cols: Span<Col> },
    #[error("extent is invalid: {0}")]
    Extent(&'static str),
}

/// Which edge `y0` is, and the direction in which the rows run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum YOrigin {
    /// `y0` is the north edge. Row 0 is the north row and rows run south (GDAL `dy < 0`).
    NorthUp,
    /// `y0` is the south edge. Row 0 is the south row and rows run north (GDAL `dy > 0`),
    /// e.g. ascending latitude.
    SouthUp,
}

/// The x axis of a [`Lattice`]. Columns always run east: `x0` is the west edge of column 0.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct XAxis {
    x0: f64,
    dx: f64,
    period: Option<f64>,
}

impl XAxis {
    /// Errors: non-finite value, `dx <= 0`.
    pub fn new(x0: f64, dx: f64) -> Result<Self, GridError> {
        finite(&[("x0", x0), ("dx", dx)])?;
        positive("dx", dx)?;
        Ok(Self { x0, dx, period: None })
    }
    /// Makes x periodic with `period` (e.g. 360 for a longitude). Geometry x and cell x
    /// then compare modulo the period: a lattice that runs over 0..360 burns geometries
    /// in [-180, 180), and a window can cross the seam.
    ///
    /// A geometry must not be wider in x than one period ([`BurnError::WiderThanPeriod`]).
    /// When the full grid is wider than one period, a place on the earth falls in more
    /// than one column, and each of them gets its record.
    ///
    /// Errors: `period` is not finite or not > 0.
    ///
    /// [`BurnError::WiderThanPeriod`]: crate::BurnError::WiderThanPeriod
    pub fn with_period(self, period: f64) -> Result<Self, GridError> {
        if !(period.is_finite() && period > 0.0) {
            return Err(GridError::Period(period));
        }
        Ok(Self { period: Some(period), ..self })
    }
    #[must_use]
    pub fn x0(&self) -> f64 {
        self.x0
    }
    #[must_use]
    pub fn dx(&self) -> f64 {
        self.dx
    }
    #[must_use]
    pub fn period(&self) -> Option<f64> {
        self.period
    }
}

/// The y axis of a [`Lattice`]: the edge of row 0, a positive step, and the row direction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YAxis {
    y0: f64,
    dy: f64,
    origin: YOrigin,
}

impl YAxis {
    /// Errors: non-finite value, `dy <= 0`.
    pub fn new(y0: f64, dy: f64, origin: YOrigin) -> Result<Self, GridError> {
        finite(&[("y0", y0), ("dy", dy)])?;
        positive("dy", dy)?;
        Ok(Self { y0, dy, origin })
    }
    #[must_use]
    pub fn y0(&self) -> f64 {
        self.y0
    }
    #[must_use]
    pub fn dy(&self) -> f64 {
        self.dy
    }
    #[must_use]
    pub fn origin(&self) -> YOrigin {
        self.origin
    }
}

fn finite(vals: &[(&'static str, f64)]) -> Result<(), GridError> {
    match vals.iter().find(|(_, v)| !v.is_finite()) {
        Some(&(name, _)) => Err(GridError::NonFinite(name)),
        None => Ok(()),
    }
}

#[expect(clippy::neg_cmp_op_on_partial_ord, reason = "the negated form also rejects NaN")]
fn positive(axis: &'static str, step: f64) -> Result<(), GridError> {
    if !(step > 0.0) {
        return Err(GridError::Resolution { axis, step });
    }
    Ok(())
}

/// Infinite regular lattice. Column edge k is `x0 + k*dx`. Row edge k is `y0 - k*dy` for
/// [`YOrigin::NorthUp`] and `y0 + k*dy` for [`YOrigin::SouthUp`]. Both steps are > 0.
///
/// Cell (r, c) is the box between `x_edge(c)`, `x_edge(c+1)`, `y_edge(r)` and
/// `y_edge(r+1)`. Every output record is in these indices.
///
/// A GDAL geotransform `[x0, dx, 0, y0, 0, dy]` maps to `XAxis::new(x0, dx)` and
/// `YAxis::new(y0, |dy|, origin)`, with `NorthUp` for `dy < 0` and `SouthUp` for
/// `dy > 0`; see [`Lattice::from_geotransform`].
///
/// A burn on a `SouthUp` lattice gives the records of the `NorthUp` lattice over the same
/// cells, with the rows flipped. The edge rules are geographic on both (a point on a shared
/// edge goes into the cell to its east and to its south).
///
/// [`XAxis::with_period`] makes x periodic (a longitude). Upstream controlledburn keeps
/// row order and the antimeridian in a front end; burnscars does both in the engine, so a
/// cell that the seam cuts gets exact coverage.
///
/// # Example
///
/// A brick of 3 rows by 4 columns: latitude centres 10, 20, 30 (ascending, so row 0 is
/// south), longitude centres 170, 180, 190, 200 (a 0..360 grid). The zone is in
/// [-180, 180).
///
/// ```
/// use burnscars::{Burner, Col, GeomId, Grid, Lattice, Mode, Row, Span, XAxis, YAxis, YOrigin};
/// use geo_types::Rect;
///
/// let x = XAxis::new(165.0, 10.0)?.with_period(360.0)?;
/// let y = YAxis::new(5.0, 10.0, YOrigin::SouthUp)?; // y0 = 5 is the south edge
/// let grid = Grid::new(Lattice::new(x, y), Span::new(Row(0), Row(3))?, Span::new(Col(0), Col(4))?)?;
/// let zone = Rect::new((-175.0, 15.0), (-165.0, 35.0));
/// let out = Burner::new(Mode::Approx).burn(&grid, [(GeomId(0), zone)])?;
/// let runs: Vec<(u32, std::ops::Range<u32>)> =
///     out.runs.iter().map(|r| (r.row.0, r.cols.range())).collect();
/// // The cells at longitude 190, latitude 30 and 20.
/// assert_eq!(runs, [(2, 2..3), (1, 2..3)]);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lattice {
    x0: f64,
    y0: f64,
    dx: f64,
    dy: f64,
    /// 1 / dx and 1 / dy: the first guess of `col_of` / `row_of`, which the edge loops
    /// then correct, so a rounding difference from a division cannot change the result.
    inv_dx: f64,
    inv_dy: f64,
    period: Option<f64>,
    origin: YOrigin,
    /// Row edge k is `y0 + (y_base - k) * dy` (`NorthUp`) or `y0 + (y_base + k) * dy`
    /// (`SouthUp`). 0 on a caller's lattice. On the north-up form of a `SouthUp` grid it is
    /// the flip sum, so the burn uses the caller's edges bit for bit.
    y_base: i64,
}

#[expect(clippy::cast_precision_loss, reason = "lattice indices stay below 2^53")]
pub(crate) fn to_f64(i: i64) -> f64 {
    i as f64
}

impl Lattice {
    #[must_use]
    pub fn new(x: XAxis, y: YAxis) -> Self {
        Self {
            x0: x.x0,
            y0: y.y0,
            dx: x.dx,
            dy: y.dy,
            inv_dx: 1.0 / x.dx,
            inv_dy: 1.0 / y.dy,
            period: x.period,
            origin: y.origin,
            y_base: 0,
        }
    }
    /// From a GDAL geotransform `[x0, dx, rx, y0, ry, dy]`: `x = x0 + col*dx + row*rx`,
    /// `y = y0 + col*ry + row*dy`. `dy < 0` is [`YOrigin::NorthUp`], `dy > 0` is
    /// [`YOrigin::SouthUp`].
    ///
    /// Errors: non-finite value, `rx != 0` or `ry != 0`, `dx <= 0` (columns run east),
    /// `dy == 0`.
    pub fn from_geotransform(gt: [f64; 6]) -> Result<Self, GridError> {
        let [x0, dx, rx, y0, ry, dy] = gt;
        finite(&[("x0", x0), ("dx", dx), ("rx", rx), ("y0", y0), ("ry", ry), ("dy", dy)])?;
        if rx != 0.0 || ry != 0.0 {
            return Err(GridError::Rotated { rx, ry });
        }
        let origin = if dy < 0.0 { YOrigin::NorthUp } else { YOrigin::SouthUp };
        Ok(Self::new(XAxis::new(x0, dx)?, YAxis::new(y0, dy.abs(), origin)?))
    }
    /// The GDAL geotransform `[x0, dx, 0, y0, 0, dy]`. It does not hold the period.
    #[must_use]
    pub fn to_geotransform(&self) -> [f64; 6] {
        let y = self.y();
        let dy = match y.origin {
            YOrigin::NorthUp => -y.dy,
            YOrigin::SouthUp => y.dy,
        };
        [self.x0, self.dx, 0.0, y.y0, 0.0, dy]
    }
    #[must_use]
    pub fn x(&self) -> XAxis {
        XAxis { x0: self.x0, dx: self.dx, period: self.period }
    }
    #[must_use]
    pub fn y(&self) -> YAxis {
        YAxis { y0: self.y0, dy: self.dy, origin: self.origin }
    }
    #[must_use]
    pub fn x_edge(&self, col: i64) -> f64 {
        self.x0 + to_f64(col) * self.dx
    }
    #[must_use]
    pub fn y_edge(&self, row: i64) -> f64 {
        // NorthUp with y_base 0: y0 + (-k)*dy, which is y0 - k*dy bit for bit.
        self.y0 + to_f64(self.y_base + self.y_sign() * row) * self.dy
    }
    /// -1 for `NorthUp`, 1 for `SouthUp`.
    fn y_sign(&self) -> i64 {
        match self.origin {
            YOrigin::NorthUp => -1,
            YOrigin::SouthUp => 1,
        }
    }
    pub(crate) fn x0(&self) -> f64 {
        self.x0
    }
    pub(crate) fn dx(&self) -> f64 {
        self.dx
    }
    pub(crate) fn dy(&self) -> f64 {
        self.dy
    }
    /// Fractional row index of `y` (north-up only): `y_edge(y_index(y)) == y` up to rounding.
    /// controlledburn `(ymax - y) / dy`.
    pub(crate) fn y_index(&self, y: f64) -> f64 {
        debug_assert!(self.origin == YOrigin::NorthUp, "y_index on a south-up lattice");
        to_f64(self.y_base) + (self.y0 - y) / self.dy
    }
    /// Column centre, controlledburn formula `xmin + (col + 0.5) * dx`.
    pub(crate) fn x_mid(&self, col: i64) -> f64 {
        self.x0 + (to_f64(col) + 0.5) * self.dx
    }
    /// Row centre (north-up only), controlledburn formula `ymax - (r + 0.5) * dy`.
    pub(crate) fn y_mid(&self, row: i64) -> f64 {
        debug_assert!(self.origin == YOrigin::NorthUp, "y_mid on a south-up lattice");
        // With y_base 0: y0 + (-r - 0.5)*dy, which is y0 - (r + 0.5)*dy bit for bit.
        self.y0 + (to_f64(self.y_base - row) - 0.5) * self.dy
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
    /// Same meaning as controlledburn `GridSpec` ([cpp/include/controlledburn/burn.hpp#L40](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/include/controlledburn/burn.hpp#L40)): x0 = xmin, y0 = ymax,
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
            XAxis::new(xmin, (xmax - xmin) / f64::from(ncol))?,
            YAxis::new(ymax, (ymax - ymin) / f64::from(nrow), YOrigin::NorthUp)?,
        );
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
    /// With a period, `x` first moves by whole periods into the period that starts at the
    /// west edge of the full grid.
    ///
    /// exactextract `Grid<infinite_extent>::get_column`
    /// ([src/grid.h#L60](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/grid.h#L60)) on the lattice;
    /// padding is `Before` / `After`. The cell edges are the authority.
    #[must_use]
    #[expect(clippy::neg_cmp_op_on_partial_ord, reason = "the negated form also rejects NaN")]
    pub fn col_of(&self, x: f64) -> Pos<Col> {
        let l = &self.lattice;
        let (s, e) = (i64::from(self.full_cols.start.0), i64::from(self.full_cols.end.0));
        let x = match l.period {
            Some(p) => {
                let west = l.x_edge(s);
                west + (x - west).rem_euclid(p)
            }
            None => x,
        };
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
        let mut c = floor_clamped((x - l.x0) * l.inv_dx, s, e - 1);
        while c > s && x < l.x_edge(c) {
            c -= 1;
        }
        while c + 1 < e && x >= l.x_edge(c + 1) {
            c += 1;
        }
        Pos::In(Col(to_u32(c)))
    }

    /// exactextract `Grid<infinite_extent>::get_row`
    /// ([src/grid.h#L85](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/grid.h#L85)) on the lattice.
    /// Cell r is `y_edge(r + 1) < y <= y_edge(r)`.
    #[must_use]
    #[expect(clippy::neg_cmp_op_on_partial_ord, reason = "the negated form also rejects NaN")]
    pub fn row_of(&self, y: f64) -> Pos<Row> {
        if self.lattice.origin == YOrigin::SouthUp {
            let n = self.north_up();
            return n.rows.pos(n.grid.row_of(y));
        }
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
        let mut r = floor_clamped(to_f64(l.y_base) + (l.y0 - y) * l.inv_dy, s, e - 1);
        while r > s && y > l.y_edge(r) {
            r -= 1;
        }
        while r + 1 < e && y <= l.y_edge(r + 1) {
            r += 1;
        }
        Pos::In(Row(to_u32(r)))
    }
}

impl Grid {
    /// This grid in the layout of the ported algorithms (north-up, no period), and the
    /// row map back to the caller's indices.
    pub(crate) fn north_up(&self) -> NorthUp {
        let l = self.lattice;
        let rows = Flip::new(l.origin == YOrigin::SouthUp, self.full_rows);
        let lattice = Lattice {
            origin: YOrigin::NorthUp,
            y_base: rows.sum.unwrap_or(l.y_base),
            period: None,
            ..l
        };
        let grid = NorthUpGrid(Grid { lattice, rows: rows.span(self.rows), ..*self });
        NorthUp { grid, period: l.period, rows }
    }
}

/// A grid in the only layout that the ported algorithms take: north-up, no period.
/// [`Grid::north_up`] is the only constructor.
#[derive(Clone, Copy, Debug)]
pub(crate) struct NorthUpGrid(Grid);

impl std::ops::Deref for NorthUpGrid {
    type Target = Grid;
    fn deref(&self) -> &Grid {
        &self.0
    }
}

/// A grid as a [`NorthUpGrid`], and the map back to the caller's indices.
#[derive(Clone, Copy, Debug)]
pub(crate) struct NorthUp {
    pub(crate) grid: NorthUpGrid,
    pub(crate) period: Option<f64>,
    pub(crate) rows: Flip,
}

impl NorthUp {
    /// No flip and no period: the records need no map.
    pub(crate) fn is_identity(&self) -> bool {
        !self.rows.is_flip() && self.period.is_none()
    }
    /// The north-up grid moved east by `k` periods.
    pub(crate) fn shifted(&self, k: i64, period: f64) -> NorthUpGrid {
        let mut g = self.grid;
        if k != 0 {
            g.0.lattice.x0 += to_f64(k) * period;
        }
        g
    }
}

/// A row flip over a full span `[s, e)`: `k` maps to `s + e - 1 - k`. It is its own inverse.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Flip {
    /// `s + e`, or `None` for no flip.
    sum: Option<i64>,
}

impl Flip {
    fn new(flip: bool, full: Span<Row>) -> Self {
        Self { sum: flip.then(|| i64::from(full.start.0) + i64::from(full.end.0)) }
    }
    pub(crate) fn is_flip(self) -> bool {
        self.sum.is_some()
    }
    pub(crate) fn index(self, k: u32) -> u32 {
        self.sum.map_or(k, |sum| to_u32(sum - 1 - i64::from(k)))
    }
    fn span(self, s: Span<Row>) -> Span<Row> {
        match self.sum {
            None => s,
            Some(sum) => Span::from_parts(
                Row(to_u32(sum - i64::from(s.end.0))),
                Row(to_u32(sum - i64::from(s.start.0))),
            ),
        }
    }
    fn pos(self, p: Pos<Row>) -> Pos<Row> {
        match (self.is_flip(), p) {
            (_, Pos::In(r)) => Pos::In(Row(self.index(r.0))),
            (false, p) => p,
            (true, Pos::Before) => Pos::After,
            (true, Pos::After) => Pos::Before,
        }
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

    fn north_up(x0: f64, y0: f64, dx: f64, dy: f64) -> Lattice {
        let x = XAxis::new(x0, dx).expect("x");
        Lattice::new(x, YAxis::new(y0, dy, YOrigin::NorthUp).expect("y"))
    }

    #[test]
    fn geotransform_round_trips_and_rejects_rotation() {
        for gt in [[10.0, 0.5, 0.0, 60.0, 0.0, -0.25], [-180.0, 1.0, 0.0, -90.0, 0.0, 2.0]] {
            let l = Lattice::from_geotransform(gt).expect("gt");
            assert_eq!(l.to_geotransform(), gt);
            assert_eq!((l.x_edge(3), l.y_edge(3)), (gt[0] + 3.0 * gt[1], gt[3] + 3.0 * gt[5]));
        }
        let south = Lattice::from_geotransform([0.0, 1.0, 0.0, 5.0, 0.0, 10.0]).expect("gt");
        assert_eq!(south.y().origin(), YOrigin::SouthUp);
        let bad = |gt| Lattice::from_geotransform(gt).expect_err("bad");
        assert!(matches!(bad([0.0, 1.0, 0.1, 0.0, 0.0, -1.0]), GridError::Rotated { .. }));
        assert!(matches!(bad([0.0, 1.0, 0.0, 0.0, 0.1, -1.0]), GridError::Rotated { .. }));
        assert!(matches!(
            bad([0.0, -1.0, 0.0, 0.0, 0.0, -1.0]),
            GridError::Resolution { axis: "dx", .. }
        ));
        assert!(matches!(
            bad([0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
            GridError::Resolution { axis: "dy", .. }
        ));
    }

    /// The north-up form of a south-up grid has the caller's row edges and centres, bit for bit.
    #[test]
    fn south_up_burns_on_the_callers_edges() {
        let x = XAxis::new(0.0, 0.1).expect("x");
        let l = Lattice::new(x, YAxis::new(-0.3, 0.1, YOrigin::SouthUp).expect("y"));
        let g = Grid::new(
            l,
            Span::new(Row(0), Row(700)).expect("r"),
            Span::new(Col(0), Col(3)).expect("c"),
        )
        .and_then(|g| {
            g.window(Span::new(Row(40), Row(90)).expect("r"), Span::new(Col(0), Col(3)).expect("c"))
        })
        .expect("grid");
        let n = g.north_up();
        assert_eq!(n.grid.rows().range(), 610..660);
        let nl = *n.grid.lattice();
        for j in 0..700_i64 {
            let k = 699 - j; // north-up row of caller row j
            assert_eq!(nl.y_edge(k + 1), l.y_edge(j), "south edge of row {j}");
            assert_eq!(nl.y_edge(k), l.y_edge(j + 1), "north edge of row {j}");
            assert_eq!(nl.y_mid(k), -0.3 + (to_f64(j) + 0.5) * 0.1, "centre of row {j}");
            let r = Row(u32::try_from(j).expect("j"));
            let south = Row(r.0.saturating_sub(1)); // a shared edge goes south
            assert_eq!(g.row_of(l.y_edge(j)), Pos::In(south), "south edge {j}");
            assert_eq!(g.row_of(l.y_edge(j).next_up()), Pos::In(r), "above south edge {j}");
        }
    }

    #[test]
    fn lattice_rejects_bad_values() {
        assert!(matches!(XAxis::new(f64::NAN, 1.0), Err(GridError::NonFinite("x0"))));
        assert!(matches!(
            YAxis::new(0.0, f64::INFINITY, YOrigin::NorthUp),
            Err(GridError::NonFinite("dy"))
        ));
        for step in [0.0, -1.0] {
            assert!(matches!(XAxis::new(0.0, step), Err(GridError::Resolution { axis: "dx", .. })));
            assert!(matches!(
                YAxis::new(0.0, step, YOrigin::SouthUp),
                Err(GridError::Resolution { axis: "dy", .. })
            ));
        }
        for p in [0.0, -360.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                XAxis::new(0.0, 1.0).and_then(|x| x.with_period(p)),
                Err(GridError::Period(_))
            ));
        }
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
            north_up(0.0, 0.0, 0.1, 0.1),
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
