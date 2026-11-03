//! Points: cell membership.
//!
//! From controlledburn (Apache-2.0, Copyright 2025 Michael Sumner), commit
//! 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558, `cpp/src/scanline.cpp` (`process_point`),
//! with exactextract `Grid::get_row` / `get_column` edge rules at the full-grid edges.
//!
//! - controlledburn `process_point`: <https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp#L884>
//! - exactextract `get_column`: <https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/grid.h#L60>
//! - exactextract `get_row`: <https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/grid.h#L85>

use crate::grid::Pos;
use crate::grid::{Col, GeomId, NorthUpGrid, Row};
use crate::input::InputBuf;
use crate::output::BurnResult;
use crate::output::Emit;

pub(crate) fn burn(input: &InputBuf, grid: &NorthUpGrid, id: GeomId, out: &mut BurnResult) {
    for p in input.points() {
        // Closed edges apply only at the full grid (spec: Closed edges of the full grid).
        if let (Pos::In(row), Pos::In(col)) = (grid.row_of(p.y), grid.col_of(p.x))
            && grid.rows().contains(row)
            && grid.cols().contains(col)
        {
            out.push(id, Emit::Point { row, col });
        }
    }
}

/// Points on copies of one window, one period apart: each point once per cell.
/// A hit on the closed far x edge goes when another copy holds the point:
/// on a full-period grid that edge is the seam, and the other copy's near edge wins.
/// The hits are on the full grid and the window filter comes last, so every tile
/// agrees on which copy wins.
pub(crate) fn burn_shifted<I>(input: &InputBuf, grids: &I, id: GeomId, out: &mut BurnResult)
where
    I: Iterator<Item = NorthUpGrid> + Clone,
{
    let mut hits: Vec<(Row, Col, bool, bool)> = Vec::new();
    for p in input.points() {
        hits.clear();
        for grid in grids.clone() {
            if let (Pos::In(row), Pos::In(col)) = (grid.row_of(p.y), grid.col_of(p.x)) {
                let far = grid.lattice().x_edge(i64::from(grid.full_cols().end().0));
                let inside = grid.rows().contains(row) && grid.cols().contains(col);
                hits.push((row, col, p.x == far, inside));
            }
        }
        let near = hits.iter().any(|h| !h.2);
        hits.sort_unstable_by_key(|h| (h.0, h.1, h.2)); // a near hit first in each cell
        hits.dedup_by_key(|h| (h.0, h.1));
        for &(row, col, far, inside) in &hits {
            if inside && !(far && near) {
                out.push(id, Emit::Point { row, col });
            }
        }
    }
}
