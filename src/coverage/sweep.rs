//! Winding sweep: runs and edges for one geometry.
//!
//! Ported from controlledburn (Apache-2.0, Copyright 2025 Michael Sumner),
//! commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558: `cpp/src/scanline.cpp`
//! (`process_polygon`, the sweep loop). Local padded indices follow exactextract
//! `infinite_extent`; `Option<LocalCol>` replaces the `prev_col = -2` sentinel.

use crate::coverage::CellAcc;
use crate::coverage::fraction::CellClass;
use crate::geom::walker_grid::{LocalCol, WalkerGrid};
use crate::grid::GeomId;
use crate::grid::{Col, Row, Span};
use crate::output::BurnResult;
use crate::output::Emit;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct MergedCell {
    pub(crate) col: LocalCol,
    pub(crate) cov: f32,
    pub(crate) wind: i32,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct SweepState {
    winding: i32,
    prev: Option<LocalCol>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum LocalEmit {
    Run(LocalCol, LocalCol),
    Edge(LocalCol, f32),
}

/// One merged cell of the controlledburn `process_polygon` sweep loop
/// ([cpp/src/scanline.cpp](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp)): the run between cells, then the cell.
pub(crate) fn step(st: &mut SweepState, cell: MergedCell, n_cols: usize) -> [Option<LocalEmit>; 2] {
    let between = match st.prev {
        Some(p) if st.winding != 0 && cell.col.0 > p.0 + 1 => {
            Some(LocalEmit::Run(LocalCol(p.0 + 1), cell.col))
        }
        _ => None,
    };
    let own = if (1..=n_cols).contains(&cell.col.0) {
        match CellClass::of(cell.cov) {
            CellClass::Edge(w) => Some(LocalEmit::Edge(cell.col, w)),
            CellClass::Full => Some(LocalEmit::Run(cell.col, LocalCol(cell.col.0 + 1))),
            CellClass::Empty => None,
        }
    } else {
        None
    };
    st.winding += cell.wind;
    st.prev = Some(cell.col);
    [between, own]
}

/// Sweeps each row of `acc` (sorted by row, col, ring) and writes global records.
pub(crate) fn emit(acc: &[CellAcc], wg: &WalkerGrid, id: GeomId, out: &mut BurnResult) {
    for row_cells in acc.chunk_by(|a, b| a.row == b.row) {
        let Some(row) = row_cells.first().and_then(|c| wg.row(c.row)) else { continue };
        row_cells
            .chunk_by(|a, b| a.col == b.col)
            .filter_map(|cs| {
                Some(MergedCell {
                    col: cs.first()?.col,
                    cov: cs.iter().fold(0.0, |s, c| s + c.cov),
                    wind: cs.iter().map(|c| c.wind).sum(),
                })
            })
            .scan(SweepState::default(), |st, c| Some(step(st, c, wg.n_cols())))
            .flatten()
            .flatten()
            .filter_map(|e| to_global(e, row, wg))
            .for_each(|e| out.push(id, e));
    }
}

fn to_global(e: LocalEmit, row: Row, wg: &WalkerGrid) -> Option<Emit> {
    match e {
        LocalEmit::Run(a, b) => {
            let (a, b) = (a.0.max(1), b.0.min(wg.n_cols() + 1));
            if a >= b {
                return None;
            }
            let start = wg.col(LocalCol(a))?;
            let last = wg.col(LocalCol(b - 1))?;
            Some(Emit::Run { row, cols: Span::from_parts(start, Col(last.0 + 1)) })
        }
        LocalEmit::Edge(c, w) => Some(Emit::Edge { row, col: wg.col(c)?, fraction: w }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn step_emits_the_gap_run_and_the_cell() {
        let mut st = SweepState::default();
        let left_pad = MergedCell { col: LocalCol(0), cov: 0.0, wind: 1 };
        assert_eq!(step(&mut st, left_pad, 8), [None, None]);
        let edge = MergedCell { col: LocalCol(4), cov: 0.5, wind: -1 };
        assert_eq!(
            step(&mut st, edge, 8),
            [
                Some(LocalEmit::Run(LocalCol(1), LocalCol(4))),
                Some(LocalEmit::Edge(LocalCol(4), 0.5))
            ]
        );
        let full = MergedCell { col: LocalCol(6), cov: 1.0, wind: 0 };
        assert_eq!(step(&mut st, full, 8), [None, Some(LocalEmit::Run(LocalCol(6), LocalCol(7)))]);
    }
}
