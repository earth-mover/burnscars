//! Lines: length in cell.
//!
//! Ported from controlledburn (Apache-2.0, Copyright 2025 Michael Sumner),
//! commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558: `cpp/src/scanline.cpp` (`process_line`).
//! One record per cell for all parts follows exactextract `add_ring_results`
//! (Apache-2.0, Copyright 2018-2022 ISciences, LLC), commit
//! 67969d99d1dfe1115be955087e2c48e1aab133e2, `src/raster_cell_intersection.cpp`.
//! `add_ring_results`: <https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/raster_cell_intersection.cpp#L479>

use crate::coverage::{fraction::to_f32, walk::walk};
use crate::error::At;
use crate::error::BurnError;
use crate::geom::walker_grid::WalkerGrid;
use crate::grid::{GeomId, Grid};
use crate::input::InputBuf;
use crate::output::BurnResult;
use crate::output::Emit;
use crate::scratch::WalkScratch;

/// From controlledburn `process_line`
/// ([cpp/src/scanline.cpp](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp)):
/// open walks, all traversals of a cell summed in visit order.
pub(crate) fn burn(
    input: &InputBuf,
    grid: &Grid,
    id: GeomId,
    s: &mut WalkScratch,
    out: &mut BurnResult,
) -> Result<(), BurnError> {
    let Some(env) = input.envelope() else {
        return Ok(());
    };
    // controlledburn walks the full grid; burnscars walks the window region.
    let Some(wg) = WalkerGrid::for_envelope(grid, &env) else {
        return Ok(());
    };
    s.arena.clear();
    for (part, line) in input.lines() {
        s.ring.clear();
        s.ring.extend(line.iter().map(|c| wg.snap(*c)));
        let at = At { id, part, ring: 0, vertex: 0 };
        walk(&mut s.ring, &wg, false, &mut s.arena).map_err(|e| e.at(at))?;
    }
    let arena = &mut s.arena;
    arena.travs.sort_unstable_by_key(|t| (t.row, t.col, t.seq));
    for group in arena.travs.chunk_by(|a, b| (a.row, a.col) == (b.row, b.col)) {
        let Some(first) = group.first() else { continue };
        let (Some(row), Some(col)) = (wg.row(first.row), wg.col(first.col)) else { continue };
        // Visit order, one running sum: the same order as the C++ loop.
        let len = group
            .iter()
            .filter_map(|t| arena.coords.get(t.range()))
            .flat_map(|cs| cs.windows(2))
            .fold(0.0, |sum, w| sum + seg_len(w));
        if len > 0.0 {
            out.push(id, Emit::Line { row, col, length: to_f32(len) });
        }
    }
    Ok(())
}

fn seg_len(w: &[crate::geom::coord::Coord]) -> f64 {
    match w {
        [a, b] => {
            let (sx, sy) = (b.x - a.x, b.y - a.y);
            (sx * sx + sy * sy).sqrt() // controlledburn uses sqrt, not hypot
        }
        _ => 0.0,
    }
}
