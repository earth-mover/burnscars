//! Polyline walker: one traversal record per cell visit, in a flat arena.
//!
//! Ported from controlledburn (Apache-2.0, Copyright 2025 Michael Sumner),
//! commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558: `cpp/src/scanline.cpp`
//! (`walk_polyline`, `LightTraversal`). The arena replaces `CellMap`.

use std::ops::Range;

use crate::error::{At, BurnError};
use crate::geom::{
    box_::{CellBox, Location, clamp},
    coord::Coord,
    crossing::Crossing,
    side::Side,
    walker_grid::{LocalCol, LocalRow, WalkerGrid},
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Valid {
    Through { entry: Side, exit: Side },
    ClosedInCell,
}

/// One traversal of one cell. From controlledburn `LightTraversal`
/// ([cpp/src/scanline.cpp#L68](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp#L68)).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Trav {
    pub(crate) row: LocalRow,
    pub(crate) col: LocalCol,
    pub(crate) seq: u32,
    start: usize,
    len: usize,
    pub(crate) entry: Option<Side>,
    pub(crate) exit: Option<Side>,
}

impl Trav {
    pub(crate) fn range(&self) -> Range<usize> {
        self.start..self.start + self.len
    }
    /// controlledburn `walk_ring` filter ([cpp/src/scanline.cpp#L338](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp#L338)): traversed with > 1 unique coordinate, or a
    /// closed ring with no entry.
    pub(crate) fn classify(&self, arena: &[Coord]) -> Option<Valid> {
        let cs = arena.get(self.range())?;
        let first = *cs.first()?;
        match (self.entry, self.exit) {
            (Some(entry), Some(exit)) if cs.iter().any(|c| *c != first) => {
                Some(Valid::Through { entry, exit })
            }
            (None, _) if cs.len() >= 3 && cs.last() == Some(&first) => Some(Valid::ClosedInCell),
            _ => None,
        }
    }
}

#[derive(Debug, Default)]
pub(crate) struct Arena {
    pub(crate) coords: Vec<Coord>,
    pub(crate) travs: Vec<Trav>,
    next_seq: u32,
}

impl Arena {
    pub(crate) fn clear(&mut self) {
        self.coords.clear();
        self.travs.clear();
        self.next_seq = 0;
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum WalkError {
    Diverged { coords: usize },
    Internal(&'static str),
}

impl WalkError {
    pub(crate) fn at(self, at: At) -> BurnError {
        match self {
            Self::Diverged { coords } => BurnError::WalkDiverged { at, coords },
            Self::Internal(detail) => BurnError::Internal { id: at.id, detail },
        }
    }
}

/// controlledburn rounding guard in `walk_polyline` ([cpp/src/scanline.cpp#L240](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp#L240)): a crossing
/// further than `next` is replaced by `next` clamped onto the box.
pub(crate) fn rounding_guard(bx: &CellBox, from: Coord, next: Coord, x: Crossing) -> Crossing {
    let (dxn, dyn_) = (next.x - from.x, next.y - from.y);
    let (dxx, dyx) = (x.coord.x - from.x, x.coord.y - from.y);
    let d_next = dxn * dxn + dyn_ * dyn_;
    let d_exit = dxx * dxx + dyx * dyx;
    if x.coord.x.is_finite() && x.coord.y.is_finite() && d_exit <= d_next * (1.0 + 1e-9) {
        return x;
    }
    let coord = Coord::new(clamp(next.x, bx.xmin, bx.xmax), clamp(next.y, bx.ymin, bx.ymax));
    let side = if next.x < bx.xmin {
        Side::Left
    } else if next.x > bx.xmax {
        Side::Right
    } else if next.y < bx.ymin {
        Side::Bottom
    } else {
        Side::Top
    };
    Crossing { side, coord }
}

/// Walker rule 9 (burnscars only): a crossing within `tol` of a box corner moves onto
/// the corner. The crossing is computed from the segment and the corner from the lattice,
/// so they can differ by an ulp. Unsnapped, the walk enters the corner-adjacent cell for a
/// one-ulp traversal that `left_hand_area` mis-chains. Snapped, that traversal has one
/// unique coordinate and `classify` drops it, whichever side the crossing reports.
pub(crate) fn corner_snap(bx: &CellBox, x: Crossing, tol: (f64, f64)) -> Crossing {
    let near = |v: f64, lo: f64, hi: f64, t: f64| {
        if (v - lo).abs() <= t {
            lo
        } else if (v - hi).abs() <= t {
            hi
        } else {
            v
        }
    };
    let c = x.coord;
    let coord = match x.side {
        Side::Left | Side::Right => Coord::new(c.x, near(c.y, bx.ymin, bx.ymax, tol.1)),
        Side::Top | Side::Bottom => Coord::new(near(c.x, bx.xmin, bx.xmax, tol.0), c.y),
    };
    Crossing { side: x.side, coord }
}

/// From controlledburn `walk_polyline`
/// ([cpp/src/scanline.cpp](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp));
/// arena storage replaces `CellMap`. `coords` must already be snapped (rule 1); the
/// reprocess branch (rule 7) appends to it.
pub(crate) fn walk(
    coords: &mut Vec<Coord>,
    wg: &WalkerGrid,
    closed: bool,
    arena: &mut Arena,
) -> Result<(), WalkError> {
    let Some(&first) = coords.first() else {
        return Ok(());
    };
    let n_input = coords.len();
    // Deliberate difference (spec table): a bound on cell visits. C++ would loop forever;
    // here the walk returns `Diverged`.
    let max_visits =
        n_input.saturating_add(64).saturating_mul(wg.n_rows() + wg.n_cols() + 4).saturating_mul(2);
    let mut visits = 0_usize;
    let mut pos = 0_usize;
    let mut row = wg.local_row_of(first.y);
    let mut col = wg.local_col_of(first.x);
    let mut last_exit: Option<Coord> = None;

    while pos < coords.len() {
        visits += 1;
        if visits > max_visits {
            return Err(WalkError::Diverged { coords: coords.len() });
        }
        let bx = wg.cell_box(row, col);
        let start = arena.coords.len();
        let mut entry: Option<Side> = None;
        let mut exit: Option<Side> = None;

        while pos < coords.len() {
            let from_exit = last_exit.is_some();
            let next = last_exit.unwrap_or(coords[pos]);
            let prev_original = pos.checked_sub(1).and_then(|i| coords.get(i)).copied();

            if arena.coords.len() == start {
                // First coordinate of this visit: enter the cell.
                entry = bx.side(next);
                arena.coords.push(next);
                if from_exit {
                    last_exit = None;
                } else {
                    pos += 1;
                }
                continue;
            }
            if bx.location(next) != Location::Outside {
                arena.coords.push(next);
                if from_exit {
                    last_exit = None;
                } else {
                    pos += 1;
                }
            } else {
                let from = prev_original
                    .or_else(|| arena.coords.last().copied())
                    .ok_or(WalkError::Internal("walker has no previous coordinate"))?;
                let x = bx
                    .crossing(from, next)
                    .ok_or(WalkError::Internal("segment does not leave the cell"))?;
                let x = rounding_guard(&bx, from, next, x);
                // Deliberate difference (spec table): walker rule 9.
                let x = corner_snap(&bx, x, wg.snap_tol());
                arena.coords.push(x.coord);
                exit = Some(x.side);
                if x.coord != next {
                    last_exit = Some(x.coord);
                }
                break;
            }
        }

        // Rule 6: force an exit if the visit ends on the boundary.
        if exit.is_none()
            && let Some(&last) = arena.coords.get(start..).and_then(<[Coord]>::last)
            && bx.location(last) == Location::Boundary
        {
            exit = bx.side(last);
        }

        // Rule 7: reprocess the start cell of a closed ring.
        if closed && exit.is_some() && entry.is_none() {
            if coords.len() > 4 * n_input + 64 {
                return Err(WalkError::Diverged { coords: coords.len() });
            }
            let seg = arena.coords.get(start..).unwrap_or(&[]);
            coords.extend_from_slice(seg);
        }

        arena.travs.push(Trav {
            row,
            col,
            seq: arena.next_seq,
            start,
            len: arena.coords.len() - start,
            entry,
            exit,
        });
        arena.next_seq = arena.next_seq.wrapping_add(1);

        if let Some(side) = exit {
            match wg.step(row, col, side) {
                Some(rc) => (row, col) = rc,
                // Deliberate difference (spec table): C++ increments and stops; here the walk
                // ends without error after the last input. Unspent input is an error.
                None if pos >= coords.len() && last_exit.is_none() => return Ok(()),
                None => return Err(WalkError::Internal("walk left the walker grid")),
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::Grid;
    use crate::input::Envelope;

    fn wg(grid: &Grid, e: (f64, f64, f64, f64)) -> WalkerGrid {
        WalkerGrid::for_envelope(
            &grid.north_up().grid,
            &Envelope { xmin: e.0, ymin: e.1, xmax: e.2, ymax: e.3 },
        )
        .expect("region")
    }
    fn ring(pts: &[(f64, f64)]) -> Vec<Coord> {
        pts.iter().map(|&(x, y)| Coord::new(x, y)).collect()
    }

    #[test]
    fn offset_square_visits_its_boundary_cells() {
        let g = Grid::from_extent(0.0, 0.0, 10.0, 10.0, 10, 10).expect("grid");
        let w = wg(&g, (2.5, 4.5, 6.5, 8.5));
        let mut r = ring(&[(2.5, 4.5), (6.5, 4.5), (6.5, 8.5), (2.5, 8.5), (2.5, 4.5)]);
        let mut a = Arena::default();
        walk(&mut r, &w, true, &mut a).expect("walk");
        let mut cells: Vec<(usize, usize)> = a
            .travs
            .iter()
            .filter(|t| t.classify(&a.coords).is_some())
            .map(|t| (t.row.0, t.col.0))
            .collect();
        cells.sort_unstable();
        cells.dedup();
        assert_eq!(cells.len(), 16); // 5x5 ring of boundary cells
    }

    // Rule 7: a ring that starts strictly inside a cell gets a proper entry on the second pass.
    #[test]
    fn reprocess_branch_closes_the_start_cell() {
        let g = Grid::from_extent(0.0, 0.0, 10.0, 10.0, 10, 10).expect("grid");
        let w = wg(&g, (2.5, 4.5, 6.5, 8.5));
        let mut r = ring(&[(2.5, 4.5), (6.5, 4.5), (6.5, 8.5), (2.5, 8.5), (2.5, 4.5)]);
        let n0 = r.len();
        let mut a = Arena::default();
        walk(&mut r, &w, true, &mut a).expect("walk");
        assert!(r.len() > n0, "the first traversal was appended");
        let start = (w.local_row_of(4.5), w.local_col_of(2.5));
        let through = a.travs.iter().any(|t| {
            (t.row, t.col) == start && matches!(t.classify(&a.coords), Some(Valid::Through { .. }))
        });
        assert!(through);
    }

    // Rule 6: a line that ends on a cell wall gets a forced exit side.
    #[test]
    fn forced_exit_on_a_wall() {
        let g = Grid::from_extent(0.0, 0.0, 3.0, 3.0, 3, 3).expect("grid");
        let w = wg(&g, (0.5, 0.5, 1.0, 0.5));
        let mut l = ring(&[(0.5, 0.5), (1.0, 0.5)]);
        let mut a = Arena::default();
        walk(&mut l, &w, false, &mut a).expect("walk");
        assert_eq!(a.travs.first().and_then(|t| t.exit), Some(Side::Right));
    }

    // Rule 5: an extrapolated crossing beyond the segment end is clamped onto the box.
    #[test]
    fn rounding_guard_clamps_extrapolated_crossings() {
        let bx = CellBox { xmin: 0.0, ymin: 0.0, xmax: 0.3, ymax: 1.0 };
        let from = Coord::new(0.1, 0.5);
        let next = Coord::new(0.30000000000000004, 0.5);
        let bad = Crossing { side: Side::Top, coord: Coord::new(5.0, 5.0) };
        let g = rounding_guard(&bx, from, next, bad);
        assert_eq!(g.side, Side::Right);
        assert_eq!(g.coord, Coord::new(0.3, 0.5));
        let ok = Crossing { side: Side::Right, coord: Coord::new(0.3, 0.5) };
        assert_eq!(rounding_guard(&bx, from, next, ok), ok);
    }

    // Rule 9: a crossing one ulp from a corner moves onto the corner; others stay.
    #[test]
    fn corner_snap_moves_near_corner_crossings() {
        let bx = CellBox { xmin: 0.7, ymin: 0.4, xmax: 0.75, ymax: 0.45000000000000007 };
        let tol = (0.05e-9, 0.05e-9);
        let near = Crossing { side: Side::Right, coord: Coord::new(0.75, 0.45) };
        assert_eq!(corner_snap(&bx, near, tol).coord, Coord::new(0.75, 0.45000000000000007));
        let mid = Crossing { side: Side::Right, coord: Coord::new(0.75, 0.43) };
        assert_eq!(corner_snap(&bx, mid, tol), mid);
        let bottom = Crossing { side: Side::Bottom, coord: Coord::new(0.7 + 1e-12, 0.4) };
        assert_eq!(corner_snap(&bx, bottom, tol).coord, Coord::new(0.7, 0.4));
    }

    // a line that ends on a padding wall at f64::MAX is not an error.
    #[test]
    fn line_ending_on_the_outer_wall_walks() {
        let g = Grid::from_extent(0.0, 0.0, 3.0, 3.0, 3, 3).expect("grid");
        let w = wg(&g, (0.5, 0.5, 1.5, 0.5));
        let mut l = ring(&[(0.5, 0.5), (f64::MAX, 0.5)]);
        let mut a = Arena::default();
        assert_eq!(walk(&mut l, &w, false, &mut a), Ok(()));
    }

    #[test]
    fn classify_rules() {
        let pts = ring(&[(0.0, 0.5), (1.0, 0.5), (0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.2)]);
        let mk = |start, len, entry, exit| Trav {
            row: LocalRow(1),
            col: LocalCol(1),
            seq: 0,
            start,
            len,
            entry,
            exit,
        };
        assert_eq!(
            mk(0, 2, Some(Side::Left), Some(Side::Right)).classify(&pts),
            Some(Valid::Through { entry: Side::Left, exit: Side::Right })
        );
        assert_eq!(mk(2, 4, None, None).classify(&pts), Some(Valid::ClosedInCell));
        assert_eq!(mk(0, 2, None, Some(Side::Right)).classify(&pts), None);
    }
}
