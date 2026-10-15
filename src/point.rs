//! Points: cell membership.
//!
//! From controlledburn (Apache-2.0, Copyright 2025 Michael Sumner), commit
//! 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558, `cpp/src/scanline.cpp` (`process_point`),
//! with exactextract `Grid::get_row` / `get_column` edge rules at the full-grid edges.

use crate::grid::Pos;
use crate::grid::{GeomId, Grid};
use crate::input::InputBuf;
use crate::output::BurnResult;
use crate::output::Emit;

pub(crate) fn burn(input: &InputBuf, grid: &Grid, id: GeomId, out: &mut BurnResult) {
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
