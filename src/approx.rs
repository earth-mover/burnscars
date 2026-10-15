//! Approx path: the cell-centre rule by an intercept sweep.
//!
//! Ported from controlledburn (Apache-2.0, Copyright 2025 Michael Sumner),
//! commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558: `cpp/src/scanline.cpp`
//! (`process_polygon_approx`). The rule matches GDAL `GDALdllImageFilledPolygon`
//! (`alg/llrasterize.cpp`): centre sampling, half-open edge test. One sweep per geometry
//! (all components) replaces controlledburn's sweep per component.

use crate::error::BurnError;
use crate::grid::{Col, GeomId, Grid, Pos, Span, to_f64};
use crate::input::InputBuf;
use crate::output::{BurnResult, Emit};

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Icpt {
    pub(crate) row: u32,
    pub(crate) x: f64,
    pub(crate) delta: i32,
    pub(crate) seq: u32,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct ApproxState {
    winding: i32,
    run_start: Option<u32>,
}

/// One intercept of the controlledburn approx sweep. Left of centre is `x <= x_mid`.
pub(crate) fn step(st: &mut ApproxState, ix: &Icpt, grid: &Grid) -> Option<Span<Col>> {
    let wc = grid.cols();
    let (c0, c1) = (wc.start().0, wc.end().0);
    let col = match grid.col_of(ix.x) {
        Pos::Before => c0,
        Pos::After => c1 - 1,
        Pos::In(c) => c.0.clamp(c0, c1 - 1),
    };
    let x_mid = grid.lattice().x_mid(i64::from(col));
    let new = st.winding + ix.delta;
    let mut run = None;
    if ix.x <= x_mid {
        if st.winding != 0 && new == 0 {
            // Leave before the centre: close the run before this cell.
            if let Some(s) = st.run_start.take() {
                run = (s < col).then(|| Span::from_parts(Col(s), Col(col)));
            }
        } else if st.winding == 0 && new != 0 {
            st.run_start = Some(col); // enter before the centre: this cell is in
        }
    } else if st.winding != 0 && new == 0 {
        // Leave after the centre: include this cell.
        let s = st.run_start.take().unwrap_or(col);
        let e = (col + 1).min(c1);
        run = (s < e).then(|| Span::from_parts(Col(s), Col(e)));
    } else if st.winding == 0 && new != 0 {
        st.run_start = Some(col + 1); // enter after the centre: start after this cell
    }
    st.winding = new;
    run
}

/// From controlledburn `process_polygon_approx`
/// ([cpp/src/scanline.cpp](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp)),
/// on the lattice, all components in one sweep.
pub(crate) fn burn(
    input: &InputBuf,
    grid: &Grid,
    id: GeomId,
    icpt: &mut Vec<Icpt>,
    out: &mut BurnResult,
) -> Result<(), BurnError> {
    let Some(env) = input.envelope() else {
        return Ok(());
    };
    let l = *grid.lattice();
    let (wr, wc) = (grid.rows(), grid.cols());
    // Deliberate difference (spec table): x pre-test against the window columns (R8).
    if env.xmax < l.x_edge(i64::from(wc.start().0)) || env.xmin > l.x_edge(i64::from(wc.end().0)) {
        return Ok(());
    }
    let row_lo = match grid.row_of(env.ymax) {
        Pos::Before => i64::from(wr.start().0),
        Pos::In(r) => i64::from(r.0).max(i64::from(wr.start().0)),
        Pos::After => return Ok(()),
    };
    let row_hi = match grid.row_of(env.ymin) {
        Pos::Before => return Ok(()),
        Pos::In(r) => i64::from(r.0).min(i64::from(wr.end().0) - 1),
        Pos::After => i64::from(wr.end().0) - 1,
    };
    if row_lo > row_hi {
        return Ok(());
    }

    icpt.clear();
    let mut seq = 0_u32;
    let (y0, dy) = (l.y0(), l.dy());
    for ring in input.rings() {
        let factor = ring.role.sign_i32(); // rings are CCW after normalization
        for w in ring.coords.windows(2) {
            let (p, q) = (w[0], w[1]);
            if p.y == q.y {
                continue; // horizontal edge
            }
            let (ya, yb) = if p.y < q.y { (p.y, q.y) } else { (q.y, p.y) };
            // Rows with ya < y_mid <= yb (top-inclusive, bottom-exclusive).
            let e_lo = ceil_clamped((y0 - yb) / dy - 0.5 - 1e-10, row_lo, row_hi + 1).max(row_lo);
            let e_hi = (ceil_clamped((y0 - ya) / dy - 0.5 - 1e-10, row_lo - 1, row_hi + 1) - 1)
                .min(row_hi);
            if e_lo > e_hi {
                continue;
            }
            let inv = 1.0 / (q.y - p.y);
            // Sign from the edge endpoints: controlledburn fix 65af8d7.
            let delta = if q.y < p.y { -1 } else { 1 } * factor;
            for r in e_lo..=e_hi {
                let t = (l.y_mid(r) - p.y) * inv;
                let row = u32::try_from(r)
                    .map_err(|_| BurnError::Internal { id, detail: "approx row" })?;
                icpt.push(Icpt { row, x: p.x + t * (q.x - p.x), delta, seq });
                seq = seq.wrapping_add(1);
            }
        }
    }
    icpt.sort_unstable_by(|a, b| {
        a.row.cmp(&b.row).then(a.x.total_cmp(&b.x)).then(a.seq.cmp(&b.seq))
    });

    let c1 = wc.end().0;
    for row_icpts in icpt.chunk_by(|a, b| a.row == b.row) {
        let Some(row) = row_icpts.first().map(|i| crate::grid::Row(i.row)) else { continue };
        let mut st = ApproxState::default();
        row_icpts
            .iter()
            .filter_map(|ix| step(&mut st, ix, grid))
            .for_each(|cols| out.push(id, Emit::Run { row, cols }));
        // The polygon goes past the right window edge.
        if st.winding != 0
            && let Some(s) = st.run_start.filter(|&s| s < c1)
        {
            out.push(id, Emit::Run { row, cols: Span::from_parts(Col(s), Col(c1)) });
        }
    }
    Ok(())
}

#[expect(clippy::cast_possible_truncation, reason = "clamped to [lo, hi] before the cast")]
fn ceil_clamped(v: f64, lo: i64, hi: i64) -> i64 {
    v.ceil().clamp(to_f64(lo), to_f64(hi)) as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn step_follows_the_centre_rule() {
        let g = Grid::from_extent(0.0, 0.0, 10.0, 10.0, 10, 10).expect("grid");
        let mut st = ApproxState::default();
        let enter = Icpt { row: 0, x: 2.5, delta: 1, seq: 0 };
        assert_eq!(step(&mut st, &enter, &g), None);
        let leave = Icpt { row: 0, x: 6.5, delta: -1, seq: 1 };
        assert_eq!(step(&mut st, &leave, &g), Some(Span::from_parts(Col(2), Col(6))));
        let enter_right = Icpt { row: 0, x: 2.6, delta: 1, seq: 2 };
        let leave_right = Icpt { row: 0, x: 6.6, delta: -1, seq: 3 };
        assert_eq!(step(&mut st, &enter_right, &g), None);
        assert_eq!(step(&mut st, &leave_right, &g), Some(Span::from_parts(Col(3), Col(7))));
    }
}
