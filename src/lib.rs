//! Sparse rasterization of polygons, lines and points onto a regular grid.
//!
//! A port of the controlledburn core (Apache-2.0, Copyright 2025 Michael Sumner),
//! commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558, with the exact coverage arithmetic
//! of exactextract (Apache-2.0, Copyright 2018-2022 ISciences, LLC),
//! commit 67969d99d1dfe1115be955087e2c48e1aab133e2.
//!
//! # Output as a sparse matrix
//!
//! The polygon output is the sparse coverage matrix `W[geometry, cell]`, as a [`Stretches`]
//! table. The value is the fraction of the cell that the geometry covers, and it is 1 for
//! interior cells. Zonal statistics are `(W v) / (W 1)` for a raster `v`. This is what
//! exactextract computes.
//!
//! - A fiber is one geometry in one row. Fibers are sorted by `(id, row)`, so `fiber_ptr`
//!   is a CSR index with the fiber as the row.
//! - A fiber holds stretches of columns, sorted and not overlapping. A run is a stretch
//!   with one value (the interior, at 1); a block is a stretch with one value per cell
//!   (the boundary). The runs keep the storage at O(perimeter), not O(area).
//! - Under [`Rule::Touched`] and [`Rule::Centre`] the cells have no weight, and every
//!   stretch is a run.
//! - [`Origin::Grid`] gives whole-grid row and column indices;
//!   [`Origin::Window`] gives them from the window's first row and column.
//!
//! Lines and points go to [`BurnResult`] through [`Burner::burn_into`].
//!
//! Input must be valid OGC geometry. The fill rule is nonzero winding.
//!
//! # Example
//!
//! ```
//! use burnscars::{Burner, GeomId, Grid, Origin, Rule, Stretches};
//! use geo_types::polygon;
//!
//! // A 4 x 4 grid over [0, 4] x [0, 4], and a square that covers 2 x 2 cells.
//! let grid = Grid::from_extent(0.0, 0.0, 4.0, 4.0, 4, 4)?;
//! let square = polygon![(x: 1.0, y: 1.0), (x: 3.0, y: 1.0), (x: 3.0, y: 3.0), (x: 1.0, y: 3.0)];
//! let mut out = Stretches::new(Rule::Fraction);
//! Burner::new(Rule::Fraction).burn_polygons(&grid, [(GeomId(0), &square)], Origin::Grid, &mut out)?;
//! // Two fibers (rows 1 and 2), each one run over columns 1..3 at weight 1.
//! assert_eq!(out.row, [1, 2]);
//! assert_eq!((out.col_start.as_slice(), out.col_end.as_slice()), (&[1, 1][..], &[3, 3][..]));
//! assert_eq!(out.values.map(|v| v.weight), Some(vec![1.0, 1.0]));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod approx;
mod burner;
mod coverage;
mod error;
mod geom;
mod grid;
mod input;
mod line;
mod output;
mod point;
mod scratch;
mod sort;

pub use burner::{Burner, Origin, Rule, bbox};
pub use error::{At, BurnError, GeomKind};
pub use grid::{Col, GeomId, Grid, GridError, Lattice, Pos, Row, Span, XAxis, YAxis, YOrigin};
pub use output::{
    BurnResult, GeomIndex, LineRef, Lines, PointRef, Points, RowMajor, RowMajorError,
    StretchValues, Stretches,
};
