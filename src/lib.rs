//! Sparse rasterization of polygons, lines and points onto a regular grid.
//!
//! A port of the controlledburn core (Apache-2.0, Copyright 2025 Michael Sumner),
//! commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558, with the exact coverage arithmetic
//! of exactextract (Apache-2.0, Copyright 2018-2022 ISciences, LLC),
//! commit 67969d99d1dfe1115be955087e2c48e1aab133e2.
//!
//! # Output as a sparse matrix
//!
//! The polygon output is the sparse coverage matrix `A[cell, geometry]`. The value is
//! the fraction of the cell that the geometry covers, and it is 1 for interior cells.
//! Zonal statistics are `(Aᵀv) / (Aᵀ1)` for a raster `v`. This is what exactextract
//! computes.
//!
//! - `runs` is a run-length layer: each entry is one interval of columns in one row,
//!   all with the value 1 ("RLE-CSR"). Plain CSR stores one entry per cell, which is
//!   O(area). The runs keep the storage at O(perimeter).
//! - `edges` is COO: explicit `(row, col, value)` triples.
//! - The burn order is geometry-major, so `geom_ptr` is a CSC index with the geometry
//!   as the column.
//! - `to_row_major()` builds the CSR index (`row_ptr`).
//! - Two geometries can cover the same cell, so the full object is a sparse
//!   `(row, col, id)` tensor. One geometry has at most one record per cell.
//!
//! Input must be valid OGC geometry. The fill rule is nonzero winding.
//!
//! # Example
//!
//! ```
//! use burnscars::{Burner, GeomId, Grid, Mode};
//! use geo_types::polygon;
//!
//! // A 4 x 4 grid over [0, 4] x [0, 4], and a square that covers 2 x 2 cells.
//! let grid = Grid::from_extent(0.0, 0.0, 4.0, 4.0, 4, 4)?;
//! let square = polygon![(x: 1.0, y: 1.0), (x: 3.0, y: 1.0), (x: 3.0, y: 3.0), (x: 1.0, y: 3.0)];
//! let out = Burner::new(Mode::Coverage).burn(&grid, [(GeomId(0), &square)])?;
//! // Full cells are in runs. Partial cells are in edges. Here all 4 cells are full.
//! let edge_cells = out.edges.iter().count();
//! let run_cells: u32 = out.runs.iter().map(|r| r.cols.width()).sum();
//! assert_eq!(edge_cells + run_cells as usize, 4);
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

pub use burner::{Burner, Mode, bbox};
pub use error::{At, BurnError, GeomKind};
pub use grid::{Col, GeomId, Grid, GridError, Lattice, Pos, Row, Span, XAxis, YAxis, YOrigin};
pub use output::{
    BurnResult, EdgeRef, Edges, GeomIndex, LineRef, Lines, PointRef, Points, RowMajor,
    RowMajorError, RunRef, Runs,
};
