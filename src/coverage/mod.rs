//! Polygon coverage path.
//!
//! Ported from controlledburn (Apache-2.0, Copyright 2025 Michael Sumner),
//! commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558: `cpp/src/scanline.cpp`
//! (`walk_ring`, `process_polygon`). One accumulator and one sweep per geometry follows
//! exactextract `RasterCellIntersection::process` (Apache-2.0, Copyright 2018-2022
//! ISciences, LLC), commit 67969d99d1dfe1115be955087e2c48e1aab133e2,
//! `src/raster_cell_intersection.cpp`.

pub(crate) mod fraction;
pub(crate) mod sweep;
pub(crate) mod walk;

use std::ops::Range;

use crate::error::At;
use crate::error::BurnError;
use crate::geom::{box_::CellBox, walker_grid::WalkerGrid};
use crate::geom::{
    coord::Coord,
    walker_grid::{LocalCol, LocalRow},
};
use crate::grid::{GeomId, Grid};
use crate::input::InputBuf;
use crate::input::RingRole;
use crate::output::BurnResult;
use crate::scratch::WalkScratch;
use fraction::Chain;
use fraction::{analytical_fraction, closed_ring_fraction, left_hand_area, to_f32};
use walk::Valid;
use walk::{Trav, walk};

/// One ring's contribution to one cell (spec: Storage, O2).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct CellAcc {
    pub(crate) row: LocalRow,
    pub(crate) col: LocalCol,
    pub(crate) ring: u32,
    pub(crate) cov: f32,
    pub(crate) wind: i32,
}

#[derive(Debug, Default)]
pub(crate) struct CoverageScratch {
    valid: Vec<(Range<usize>, Valid)>,
    ranges: Vec<Range<usize>>,
    chains: Vec<Chain>,
    lha: Vec<Coord>,
    acc: Vec<CellAcc>,
}

/// controlledburn `process_polygon` ([cpp/src/scanline.cpp](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp)) for all rings of one
/// geometry at once, as exactextract `RasterCellIntersection::process`
/// ([src/raster_cell_intersection.cpp](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/raster_cell_intersection.cpp)); one sweep per geometry.
pub(crate) fn burn(
    input: &InputBuf,
    grid: &Grid,
    id: GeomId,
    walk_s: &mut WalkScratch,
    s: &mut CoverageScratch,
    out: &mut BurnResult,
) -> Result<(), BurnError> {
    let Some(env) = input.envelope() else {
        return Ok(());
    };
    let Some(wg) = WalkerGrid::for_envelope(grid, &env) else {
        return Ok(());
    };
    walk_s.arena.clear();
    s.acc.clear();
    let internal = |detail| BurnError::Internal { id, detail };
    for (ring_idx, ring) in input.rings().enumerate() {
        let at = At { id, part: ring.part, ring: ring.ring, vertex: 0 };
        walk_s.ring.clear();
        walk_s.ring.extend(ring.coords.iter().map(|c| wg.snap(*c)));
        let mark = walk_s.arena.travs.len();
        walk(&mut walk_s.ring, &wg, true, &mut walk_s.arena).map_err(|e| e.at(at))?;
        let ring_no = u32::try_from(ring_idx).map_err(|_| internal("more than 2^32 rings"))?;
        let arena = &mut walk_s.arena;
        let travs = arena.travs.get_mut(mark..).ok_or_else(|| internal("arena mark"))?;
        travs.sort_unstable_by_key(|t| (t.row, t.col, t.seq));
        for group in travs.chunk_by(|a, b| (a.row, a.col) == (b.row, b.col)) {
            if let Some(acc) =
                cell_record(group, &arena.coords, &wg, ring.role, ring_no, s).map_err(internal)?
            {
                s.acc.push(acc);
            }
        }
    }
    // The ring number makes the key unique, so the float sums keep one order everywhere.
    s.acc.sort_unstable_by_key(|a| (a.row, a.col, a.ring));
    sweep::emit(&s.acc, &wg, id, out);
    Ok(())
}

/// One ring's coverage and winding for one cell (controlledburn `walk_ring`,
/// [cpp/src/scanline.cpp](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp)).
fn cell_record(
    group: &[Trav],
    coords: &[Coord],
    wg: &WalkerGrid,
    role: RingRole,
    ring: u32,
    s: &mut CoverageScratch,
) -> Result<Option<CellAcc>, &'static str> {
    let Some(first) = group.first() else {
        return Ok(None);
    };
    let (r, c) = (first.row, first.col);
    if wg.row(r).is_none() {
        return Ok(None); // padding rows carry no winding for any grid row
    }
    s.valid.clear();
    s.valid.extend(group.iter().filter_map(|t| t.classify(coords).map(|v| (t.range(), v))));
    if s.valid.is_empty() {
        return Ok(None);
    }
    let bx = wg.cell_box(r, c);
    let frac = if wg.col(c).is_some() { to_f32(fraction_of(&bx, coords, s)?) } else { 0.0 };

    let y_mid = (bx.ymin + bx.ymax) / 2.0;
    let mut wind = 0;
    let mut crossed = false;
    for (range, v) in &s.valid {
        let Valid::Through { .. } = v else { continue }; // closed rings give no winding
        let Some(cs) = coords.get(range.clone()) else { continue };
        let (Some(a), Some(b)) = (cs.first(), cs.last()) else { continue };
        if cs.len() < 2 {
            continue;
        }
        // Half-open at y_mid: controlledburn fix 0efe32f.
        if (a.y >= y_mid && b.y < y_mid) || (a.y < y_mid && b.y >= y_mid) {
            crossed = true;
            wind += if a.y >= y_mid { -1 } else { 1 };
        }
    }
    if frac == 0.0 && !crossed {
        return Ok(None);
    }
    Ok(Some(CellAcc {
        row: r,
        col: c,
        ring,
        cov: role.sign_f32() * frac,
        wind: wind * role.sign_i32(),
    }))
}

fn fraction_of(
    bx: &CellBox,
    coords: &[Coord],
    s: &mut CoverageScratch,
) -> Result<f64, &'static str> {
    let err = "traversal end is not on the cell boundary";
    match s.valid.as_slice() {
        [(r, Valid::ClosedInCell)] => {
            Ok(closed_ring_fraction(bx, coords.get(r.clone()).ok_or(err)?))
        }
        [(r, Valid::Through { .. })] => {
            analytical_fraction(bx, coords.get(r.clone()).ok_or(err)?).ok_or(err)
        }
        many => {
            s.ranges.clear();
            s.ranges.extend(many.iter().map(|(r, _)| r.clone()));
            let area = bx.area();
            if area > 0.0 {
                Ok(left_hand_area(bx, coords, &s.ranges, &mut s.chains, &mut s.lha).ok_or(err)?
                    / area)
            } else {
                Ok(0.0)
            }
        }
    }
}
