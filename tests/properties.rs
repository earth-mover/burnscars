//! Property tests on lattice-biased random rings on non-dyadic grids.
//!
//! Ported from controlledburn-rs `tests/properties.rs` (Apache-2.0, Copyright 2026
//! Michael Sumner), commit 6198e44f. The generator puts vertices on cell corners, walls,
//! row and column centre lines, the grid edge and just outside it: the places where a
//! typed decimal and a computed cell edge differ by an ulp.
//!
//! Properties:
//! 1. Coverage conserves area: the burned area equals the ring area clipped to the grid.
//! 2. The tables are well formed: fractions in (0, 1], runs inside the grid, each cell
//!    at most once.
//! 3. The per-cell cover does not change with ring orientation or start vertex.
//! 4. The per-cell cover equals an exact per-cell clip (Sutherland-Hodgman, shoelace)
//!    on grids up to 16 x 16.
//! 5. A full coverage cell is in approx runs; an untouched coverage cell is not.
//! 6. A burn at 2x resolution, averaged over 2 x 2 blocks, equals the 1x burn.
//!
//! Long run: `PROPTEST_CASES=100000 cargo test --release --test properties`.

#![expect(
    clippy::cast_precision_loss,
    clippy::expect_used,
    clippy::panic,
    reason = "test code: small indices, a failed burn fails the test"
)]

use burnscars::{BurnResult, Burner, GeomId, Grid, Mode};
use geo_types::{LineString, Polygon, coord};
use proptest::prelude::*;

type P = (f64, f64);

#[derive(Clone, Copy, Debug)]
struct Ext {
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64,
}

/// A test grid: the extent as typed, and the cell counts.
#[derive(Clone, Copy, Debug)]
struct G {
    ext: Ext,
    res: f64,
    ncol: u32,
    nrow: u32,
}

impl G {
    fn grid(&self) -> Grid {
        let e = self.ext;
        Grid::from_extent(e.xmin, e.ymin, e.xmax, e.ymax, self.ncol, self.nrow).expect("grid")
    }
    fn doubled(&self) -> G {
        G { res: self.res / 2.0, ncol: 2 * self.ncol, nrow: 2 * self.nrow, ..*self }
    }
    fn nc(&self) -> usize {
        self.ncol as usize
    }
    fn nr(&self) -> usize {
        self.nrow as usize
    }
}

fn mk_grid(res: f64, x0: f64, y0: f64, ncol: u32, nrow: u32) -> G {
    let ext = Ext {
        xmin: x0,
        ymin: y0,
        xmax: x0 + f64::from(ncol) * res,
        ymax: y0 + f64::from(nrow) * res,
    };
    G { ext, res, ncol, nrow }
}

/// Non-dyadic steps on several origins, so `origin + k * res` is seldom the typed decimal.
fn grid_up_to(max_dim: u32) -> impl Strategy<Value = G> {
    let res = prop::sample::select(vec![0.1, 0.05, 0.3, 1.0 / 3.0, 1.0 / 7.0, 0.7, 1.0, 2.5]);
    let origin = prop::sample::select(vec![
        (0.0, 0.0),
        (-180.0, -90.0),
        (0.3, 0.3),
        (-1.7, 2.9),
        (1000.1, -2000.7),
    ]);
    (res, origin, 4_u32..max_dim, 4_u32..max_dim)
        .prop_map(|(res, (x0, y0), nc, nr)| mk_grid(res, x0, y0, nc, nr))
}

#[derive(Clone, Copy, Debug)]
enum Snap {
    None,
    Corner,
    XWall,
    YWall,
    RowCentre,
    ColCentre,
    Centre,
    GridEdge,
    Outside,
}

fn snap_kind() -> impl Strategy<Value = Snap> {
    prop_oneof![
        3 => Just(Snap::None),
        3 => Just(Snap::Corner),
        2 => Just(Snap::XWall),
        2 => Just(Snap::YWall),
        2 => Just(Snap::RowCentre),
        2 => Just(Snap::ColCentre),
        1 => Just(Snap::Centre),
        1 => Just(Snap::GridEdge),
        1 => Just(Snap::Outside),
    ]
}

/// Lattice positions as a user types them: round the multiplier, recompute from the
/// typed origin. This keeps the "0.3 vs 3 * 0.1" difference.
fn apply_snap(g: &G, x: f64, y: f64, s: Snap) -> P {
    let (e, r) = (g.ext, g.res);
    let kx = ((x - e.xmin) / r).round();
    let ky = ((e.ymax - y) / r).round();
    let (xl, yl) = (e.xmin + kx * r, e.ymax - ky * r);
    let (xc, yc) = (e.xmin + (kx + 0.5) * r, e.ymax - (ky + 0.5) * r);
    match s {
        Snap::None => (x, y),
        Snap::Corner => (xl, yl),
        Snap::XWall => (xl, y),
        Snap::YWall => (x, yl),
        Snap::RowCentre => (x, yc),
        Snap::ColCentre => (xc, y),
        Snap::Centre => (xc, yc),
        Snap::GridEdge => {
            let (dl, dr) = ((x - e.xmin).abs(), (x - e.xmax).abs());
            let (db, dt) = ((y - e.ymin).abs(), (y - e.ymax).abs());
            let m = dl.min(dr).min(db).min(dt);
            if m == dl {
                (e.xmin, y)
            } else if m == dr {
                (e.xmax, y)
            } else if m == db {
                (x, e.ymin)
            } else {
                (x, e.ymax)
            }
        }
        Snap::Outside => {
            let (w, h) = (e.xmax - e.xmin, e.ymax - e.ymin);
            if (x - e.xmin) < w / 2.0 {
                (e.xmin - (x - e.xmin).abs() - 0.37 * r, y)
            } else {
                (x, e.ymax + (y - e.ymin).abs() * 0.5 + 0.37 * h / f64::from(g.nrow))
            }
        }
    }
}

/// A closed simple ring: vertices sorted by angle about a centre, then snapped. Snaps
/// that make the ring touch itself are filtered (invalid input).
fn ring_on(g: G) -> impl Strategy<Value = (G, Vec<P>)> {
    let e = g.ext;
    let (w, h) = (e.xmax - e.xmin, e.ymax - e.ymin);
    (3_usize..12)
        .prop_flat_map(|n| {
            (
                (0.05_f64..0.95, 0.05_f64..0.95),
                prop::collection::vec((0.0_f64..1.0, 0.05_f64..1.0, snap_kind()), n),
            )
        })
        .prop_map(move |((cx, cy), vs)| {
            let (cx, cy) = (e.xmin + cx * w, e.ymin + cy * h);
            let mut pts: Vec<(f64, f64, Snap)> = vs
                .iter()
                .enumerate()
                .map(|(i, &(jit, rad, s))| {
                    let ang = (i as f64 + jit * 0.9) / vs.len() as f64 * std::f64::consts::TAU;
                    (cx + ang.cos() * rad * w * 0.6, cy + ang.sin() * rad * h * 0.6, s)
                })
                .collect();
            pts.sort_by(|a, b| (a.1 - cy).atan2(a.0 - cx).total_cmp(&(b.1 - cy).atan2(b.0 - cx)));
            let mut ring: Vec<P> = pts.iter().map(|&(x, y, s)| apply_snap(&g, x, y, s)).collect();
            ring.dedup();
            if ring.first() != ring.last() {
                ring.push(ring[0]);
            }
            (g, ring)
        })
        .prop_filter("simple ring with area", |(_, r)| {
            r.len() >= 4 && shoelace(r).abs() > 1e-12 && is_simple(r)
        })
}

fn ring() -> impl Strategy<Value = (G, Vec<P>)> {
    grid_up_to(40).prop_flat_map(ring_on)
}

/// Crossing-number test; `p` is not on the ring.
fn inside(r: &[P], p: P) -> bool {
    r.windows(2).fold(false, |acc, w| {
        let (a, b) = (w[0], w[1]);
        let crosses =
            (a.1 > p.1) != (b.1 > p.1) && p.0 < a.0 + (p.1 - a.1) / (b.1 - a.1) * (b.0 - a.0);
        acc != crosses
    })
}

/// A shell from `ring_on` and one snapped star-shaped hole strictly inside it, so the
/// hole's boundary cells lie inside the shell and off its edges.
fn holed_on(g: G) -> impl Strategy<Value = (G, Vec<P>, Vec<P>)> {
    let star =
        (0.02_f64..0.25, prop::collection::vec((0.0_f64..1.0, 0.3_f64..1.0, snap_kind()), 3..8));
    (ring_on(g), star)
        .prop_map(|((g, r), (scale, vs))| {
            let n = (r.len() - 1) as f64;
            let cx = r[..r.len() - 1].iter().map(|p| p.0).sum::<f64>() / n;
            let cy = r[..r.len() - 1].iter().map(|p| p.1).sum::<f64>() / n;
            let rad = scale * (g.ext.xmax - g.ext.xmin).min(g.ext.ymax - g.ext.ymin);
            let mut h: Vec<P> = vs
                .iter()
                .enumerate()
                .map(|(i, &(jit, k, s))| {
                    let a = (i as f64 + jit * 0.9) / vs.len() as f64 * std::f64::consts::TAU;
                    apply_snap(&g, cx + a.cos() * k * rad, cy + a.sin() * k * rad, s)
                })
                .collect();
            h.dedup();
            if h.first() != h.last() {
                h.push(h[0]);
            }
            (g, r, h)
        })
        .prop_filter("simple hole strictly inside the shell", |(_, r, h)| {
            h.len() >= 4
                && shoelace(h).abs() > 1e-12
                && is_simple(h)
                && h.iter().all(|&p| inside(r, p))
                && r.windows(2).all(|e| h.windows(2).all(|f| !seg_touch(e[0], e[1], f[0], f[1])))
        })
}

/// Signed area of a closed ring.
fn shoelace(r: &[P]) -> f64 {
    r.windows(2).map(|w| w[0].0 * w[1].1 - w[1].0 * w[0].1).sum::<f64>() / 2.0
}

fn orient(p: P, q: P, r: P) -> f64 {
    (q.0 - p.0) * (r.1 - p.1) - (q.1 - p.1) * (r.0 - p.0)
}

fn on_segment(p: P, q: P, r: P) -> bool {
    r.0 >= p.0.min(q.0) && r.0 <= p.0.max(q.0) && r.1 >= p.1.min(q.1) && r.1 <= p.1.max(q.1)
}

/// Closed test: segments ab and cd cross or touch.
fn seg_touch(a: P, b: P, c: P, d: P) -> bool {
    let (o1, o2, o3, o4) = (orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b));
    let nonzero = o1 != 0.0 && o2 != 0.0 && o3 != 0.0 && o4 != 0.0;
    if nonzero && (o1 > 0.0) != (o2 > 0.0) && (o3 > 0.0) != (o4 > 0.0) {
        return true;
    }
    (o1 == 0.0 && on_segment(a, b, c))
        || (o2 == 0.0 && on_segment(a, b, d))
        || (o3 == 0.0 && on_segment(c, d, a))
        || (o4 == 0.0 && on_segment(c, d, b))
}

/// Strictly simple: no repeated vertex, no touching non-adjacent edges, no fold-back.
fn is_simple(r: &[P]) -> bool {
    let n = r.len() - 1;
    for i in 0..n {
        for j in i + 1..n {
            if r[i] == r[j] {
                return false;
            }
        }
    }
    for i in 0..n {
        for j in i + 1..n {
            let (a, b, c, d) = (r[i], r[i + 1], r[j], r[j + 1]);
            let adjacent = j == i + 1 || (i == 0 && j == n - 1);
            if !adjacent {
                if seg_touch(a, b, c, d) {
                    return false;
                }
            } else {
                // Fold-back: adjacent edges leave the shared vertex at an angle < 1e-9 rad.
                // A spike a few ulps wide is valid for GEOS, but exactextract burns its
                // cell as 0 (burnscars matches); upstream filters only the exact case.
                let shared = if j == i + 1 { b } else { a };
                let o_ab = if shared == b { a } else { b };
                let o_cd = if shared == c { d } else { c };
                let (u, v) = (
                    (o_ab.0 - shared.0, o_ab.1 - shared.1),
                    (o_cd.0 - shared.0, o_cd.1 - shared.1),
                );
                let dot = u.0 * v.0 + u.1 * v.1;
                let cross = u.0 * v.1 - u.1 * v.0;
                if dot > 0.0 && cross.abs() <= 1e-9 * u.0.hypot(u.1) * v.0.hypot(v.1) {
                    return false;
                }
            }
        }
    }
    true
}

/// Sutherland-Hodgman clip against one half-plane. `axis` 0 is x, 1 is y; `keep_ge`
/// keeps the side `>= v`.
fn clip_plane(input: &[P], axis: usize, v: f64, keep_ge: bool) -> Vec<P> {
    let coord = |p: P| if axis == 0 { p.0 } else { p.1 };
    let inside = |p: P| if keep_ge { coord(p) >= v } else { coord(p) <= v };
    let isect = |a: P, b: P| {
        let t = (v - coord(a)) / (coord(b) - coord(a));
        if axis == 0 { (v, a.1 + t * (b.1 - a.1)) } else { (a.0 + t * (b.0 - a.0), v) }
    };
    let mut out = Vec::with_capacity(input.len() + 4);
    for i in 0..input.len() {
        let (prev, cur) = (input[(i + input.len() - 1) % input.len()], input[i]);
        match (inside(prev), inside(cur)) {
            (true, true) => out.push(cur),
            (false, true) => {
                out.push(isect(prev, cur));
                out.push(cur);
            }
            (true, false) => out.push(isect(prev, cur)),
            (false, false) => {}
        }
    }
    out
}

/// Area of a simple closed ring inside a rectangle.
fn clipped_area(r: &[P], e: &Ext) -> f64 {
    let mut c: Vec<P> = r[..r.len() - 1].to_vec();
    for (axis, v, ge) in
        [(0, e.xmin, true), (0, e.xmax, false), (1, e.ymin, true), (1, e.ymax, false)]
    {
        if c.is_empty() {
            return 0.0;
        }
        c = clip_plane(&c, axis, v, ge);
    }
    if c.len() < 3 {
        return 0.0;
    }
    c.push(c[0]);
    shoelace(&c).abs()
}

fn burn1(r: &[P], g: &G, mode: Mode) -> BurnResult {
    burn_holed(r, &[], g, mode)
}

fn burn_holed(r: &[P], holes: &[Vec<P>], g: &G, mode: Mode) -> BurnResult {
    let ls =
        |r: &[P]| -> LineString<f64> { r.iter().map(|&(x, y)| coord! { x: x, y: y }).collect() };
    let p = Polygon::new(ls(r), holes.iter().map(|h| ls(h)).collect());
    Burner::new(mode)
        .burn(&g.grid(), [(GeomId(0), &p)])
        .unwrap_or_else(|e| panic!("burn failed: {e}\n{}\n{}", fmt_grid(g), fmt_ring(r)))
}

/// Dense row-major cover from the tables.
fn cover(res: &BurnResult, g: &G) -> Vec<f64> {
    let nc = g.nc();
    let mut m = vec![0.0; nc * g.nr()];
    for x in res.runs.iter() {
        for c in x.cols.range() {
            m[x.row.0 as usize * nc + c as usize] += 1.0;
        }
    }
    for e in res.edges.iter() {
        m[e.row.0 as usize * nc + e.col.0 as usize] += f64::from(e.fraction);
    }
    m
}

fn area(res: &BurnResult, g: &G) -> f64 {
    let l = *g.grid().lattice();
    let runs: f64 = res.runs.iter().map(|x| f64::from(x.cols.width())).sum();
    let edges: f64 = res.edges.iter().map(|e| f64::from(e.fraction)).sum();
    l.dx() * l.dy() * (runs + edges)
}

fn fmt_ring(r: &[P]) -> String {
    let body: Vec<String> = r.iter().map(|p| format!("{:?} {:?}", p.0, p.1)).collect();
    format!("POLYGON (({}))", body.join(", "))
}

fn fmt_grid(g: &G) -> String {
    let e = g.ext;
    format!(
        "from_extent({:?}, {:?}, {:?}, {:?}, {}, {})",
        e.xmin, e.ymin, e.xmax, e.ymax, g.ncol, g.nrow
    )
}

fn max_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

/// Property 4 on one case; also the body of the regression tests.
fn check_per_cell(g: &G, r: &[P]) -> Result<(), TestCaseError> {
    check_per_cell_holed(g, r, &[])
}

/// Property 4 for a shell with holes strictly inside it.
fn check_per_cell_holed(g: &G, r: &[P], holes: &[Vec<P>]) -> Result<(), TestCaseError> {
    let m = cover(&burn_holed(r, holes, g, Mode::Coverage), g);
    let (e, res) = (g.ext, g.res);
    for row in 0..g.nr() {
        for col in 0..g.nc() {
            let cell = Ext {
                xmin: e.xmin + col as f64 * res,
                xmax: e.xmin + (col + 1) as f64 * res,
                ymax: e.ymax - row as f64 * res,
                ymin: e.ymax - (row + 1) as f64 * res,
            };
            let hole_area: f64 = holes.iter().map(|h| clipped_area(h, &cell)).sum();
            let truth = (clipped_area(r, &cell) - hole_area) / (res * res);
            let got = m[row * g.nc() + col];
            prop_assert!(
                (got - truth).abs() <= 2e-6,
                "row {row} col {col}: burned {got} vs clipped {truth}\n{}\n{}\nholes {:?}",
                fmt_grid(g),
                fmt_ring(r),
                holes.iter().map(|h| fmt_ring(h)).collect::<Vec<_>>()
            );
        }
    }
    Ok(())
}

fn cases() -> u32 {
    std::env::var("PROPTEST_CASES").ok().and_then(|s| s.parse().ok()).unwrap_or(1000)
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: cases(),
        // The hole generator rejects about 2 of 3 draws.
        max_local_rejects: cases().saturating_mul(4).max(65_536),
        max_shrink_iters: 4000,
        ..ProptestConfig::default()
    })]

    /// A ring R and its complement (a box around R and the grid, with R as a hole) add up
    /// to 1 in every cell, so to the grid area in total. Each of the three rings (R alone,
    /// the box, R as the hole) gets a random orientation: the fill must not depend on it.
    #[test]
    fn ring_and_complement_sum_to_one((g, r) in ring(), flip in any::<[bool; 3]>()) {
        let e = g.ext;
        let (w, h) = (e.xmax - e.xmin, e.ymax - e.ymin);
        let lo_x = r.iter().map(|p| p.0).fold(e.xmin, f64::min) - w;
        let hi_x = r.iter().map(|p| p.0).fold(e.xmax, f64::max) + w;
        let lo_y = r.iter().map(|p| p.1).fold(e.ymin, f64::min) - h;
        let hi_y = r.iter().map(|p| p.1).fold(e.ymax, f64::max) + h;
        let bx = vec![(lo_x, lo_y), (hi_x, lo_y), (hi_x, hi_y), (lo_x, hi_y), (lo_x, lo_y)];
        let orient = |v: &[P], rev: bool| -> Vec<P> {
            if rev { v.iter().rev().copied().collect() } else { v.to_vec() }
        };
        let (ring, shell, hole) = (orient(&r, flip[0]), orient(&bx, flip[1]), orient(&r, flip[2]));
        let a = cover(&burn1(&ring, &g, Mode::Coverage), &g);
        let b = cover(&burn_holed(&shell, std::slice::from_ref(&hole), &g, Mode::Coverage), &g);
        for (i, (x, y)) in a.iter().zip(&b).enumerate() {
            prop_assert!((x + y - 1.0).abs() <= 2e-6,
                "cell {i}: ring {x} + complement {y} != 1 (flip {flip:?})\n{}\n{}", fmt_grid(&g), fmt_ring(&r));
        }
    }

    #[test]
    fn hole_conserves_area((g, r, h) in grid_up_to(40).prop_flat_map(holed_on)) {
        let res = burn_holed(&r, std::slice::from_ref(&h), &g, Mode::Coverage);
        let (got, truth) = (area(&res, &g), clipped_area(&r, &g.ext) - clipped_area(&h, &g.ext));
        let cell = g.res * g.res;
        let tol = 1e-6 * truth.max(cell) + res.edges.len() as f64 * cell * 2e-7;
        prop_assert!((got - truth).abs() <= tol,
            "area {got} vs {truth}\n{}\n{}\nhole {}", fmt_grid(&g), fmt_ring(&r), fmt_ring(&h));
    }

    #[test]
    fn hole_matches_exact_per_cell_oracle((g, r, h) in grid_up_to(16).prop_flat_map(holed_on)) {
        check_per_cell_holed(&g, &r, std::slice::from_ref(&h))?;
    }

    #[test]
    fn coverage_conserves_area((g, r) in ring()) {
        let res = burn1(&r, &g, Mode::Coverage);
        let (got, truth) = (area(&res, &g), clipped_area(&r, &g.ext));
        let cell = g.res * g.res;
        // f32 fractions: a few ulp of slack per edge cell.
        let tol = 1e-6 * truth.max(cell) + res.edges.len() as f64 * cell * 2e-7;
        prop_assert!((got - truth).abs() <= tol,
            "area {got} vs {truth} (diff {:+e}, tol {tol:e})\n{}\n{}", got - truth, fmt_grid(&g), fmt_ring(&r));
    }

    #[test]
    fn tables_are_well_formed((g, r) in ring()) {
        let res = burn1(&r, &g, Mode::Coverage);
        let mut seen = std::collections::HashSet::new();
        for x in res.runs.iter() {
            prop_assert!(x.row.0 < g.nrow && x.cols.end().0 <= g.ncol, "bad run {x:?}");
            for c in x.cols.range() {
                prop_assert!(seen.insert((x.row.0, c)), "cell twice (run) {x:?}\n{}\n{}", fmt_grid(&g), fmt_ring(&r));
            }
        }
        for e in res.edges.iter() {
            prop_assert!(e.row.0 < g.nrow && e.col.0 < g.ncol, "bad edge {e:?}");
            prop_assert!(e.fraction > 0.0 && e.fraction <= 1.0 + 1e-6,
                "fraction out of range {e:?}\n{}\n{}", fmt_grid(&g), fmt_ring(&r));
            prop_assert!(seen.insert((e.row.0, e.col.0)), "cell twice (edge) {e:?}\n{}\n{}", fmt_grid(&g), fmt_ring(&r));
        }
    }

    #[test]
    fn invariant_under_orientation_and_start((g, r) in ring(), k in 0_usize..8) {
        let base = cover(&burn1(&r, &g, Mode::Coverage), &g);
        let rev: Vec<P> = r.iter().rev().copied().collect();
        let worst = max_diff(&base, &cover(&burn1(&rev, &g, Mode::Coverage), &g));
        prop_assert!(worst <= 1e-5, "reversed ring differs by {worst}\n{}\n{}", fmt_grid(&g), fmt_ring(&r));
        let n = r.len() - 1;
        let s = k % n;
        let mut rot: Vec<P> = r[s..n].to_vec();
        rot.extend_from_slice(&r[..s]);
        rot.push(rot[0]);
        let worst = max_diff(&base, &cover(&burn1(&rot, &g, Mode::Coverage), &g));
        prop_assert!(worst <= 1e-5, "start-rotated ring differs by {worst}\n{}\n{}", fmt_grid(&g), fmt_ring(&r));
    }

    #[test]
    fn matches_exact_per_cell_oracle((g, r) in grid_up_to(16).prop_flat_map(ring_on)) {
        check_per_cell(&g, &r)?;
    }

    #[test]
    fn approx_agrees_with_coverage_on_full_and_empty_cells((g, r) in ring()) {
        let cov = cover(&burn1(&r, &g, Mode::Coverage), &g);
        let apx = cover(&burn1(&r, &g, Mode::Approx), &g);
        for (i, (c, a)) in cov.iter().zip(&apx).enumerate() {
            if *c >= 1.0 - 1e-6 {
                prop_assert!(*a > 0.0, "cell {i} full in coverage, absent in approx\n{}\n{}", fmt_grid(&g), fmt_ring(&r));
            }
            if *c == 0.0 {
                prop_assert!(*a == 0.0, "cell {i} untouched in coverage, in approx\n{}\n{}", fmt_grid(&g), fmt_ring(&r));
            }
        }
    }

    #[test]
    fn doubling_resolution_aggregates_exactly((g, r) in ring()) {
        let m1 = cover(&burn1(&r, &g, Mode::Coverage), &g);
        let g2 = g.doubled();
        let m2 = cover(&burn1(&r, &g2, Mode::Coverage), &g2);
        let (nc, nc2) = (g.nc(), g2.nc());
        for row in 0..g.nr() {
            for col in 0..nc {
                let at = |dr: usize, dc: usize| m2[(2 * row + dr) * nc2 + 2 * col + dc];
                let agg = (at(0, 0) + at(0, 1) + at(1, 0) + at(1, 1)) / 4.0;
                let got = m1[row * nc + col];
                prop_assert!((got - agg).abs() <= 2e-6,
                    "row {row} col {col}: 1x {got} vs 2x aggregate {agg}\n{}\n{}", fmt_grid(&g), fmt_ring(&r));
            }
        }
    }
}

/// The shrunken failures in controlledburn-rs `tests/properties.proptest-regressions`,
/// and the controlledburn#16 reproducer. Each must match the per-cell oracle. Upstream seeds 3
/// and 4 are invalid input (a repeated vertex; a vertex on another edge) and are left out.
#[test]
fn upstream_regressions() {
    let cases: [(G, Vec<P>); 3] = [
        (
            G {
                ext: Ext { xmin: -180.0, ymin: -90.0, xmax: -175.0, ymax: -80.666_666_666_666_67 },
                res: 1.0 / 3.0,
                ncol: 15,
                nrow: 28,
            },
            vec![
                (-179.896_547_563_642_14, -88.833_333_333_333_34),
                (-180.940_765_599_607_9, -92.368_089_775_040_1),
                (-179.6, -88.518_152_273_978),
                (-179.896_547_563_642_14, -88.833_333_333_333_34),
            ],
        ),
        (
            mk_grid(0.1, 0.0, 0.0, 4, 4),
            vec![
                (0.25, -0.045_985_306_711_012_02),
                (0.255_693_165_682_278_6, -0.089_198_004_038_124_35),
                (0.403_982_260_541_639, 0.098_403_822_466_420_17),
                (0.25, -0.045_985_306_711_012_02),
            ],
        ),
        (mk_grid(0.1, 0.0, 0.0, 100, 100), vec![(0.3, 0.3), (4.7, 0.3), (2.5, 6.2), (0.3, 0.3)]),
    ];
    for (g, r) in &cases {
        assert!(is_simple(r), "regression ring is not simple: {}", fmt_ring(r));
        if let Err(e) = check_per_cell(g, r) {
            panic!("{e}");
        }
    }
}

/// Found by this suite at 100k cases. The edge (0.7, 0.5) -> (0.775, 0.425) passes one
/// ulp from the lattice corner (0.75, 0.45000000000000007). Without walker rule 9 the
/// walk enters cell (8, 8) for a one-ulp traversal and `left_hand_area` gives 1.165.
#[test]
fn crossing_one_ulp_from_a_corner() {
    let g = mk_grid(0.05, 0.3, 0.3, 14, 11);
    let r = vec![
        (0.716_341_696_372_258_5, 0.438_511_550_039_540_6),
        (0.771_628_165_166_780_6, 0.3),
        (0.784_606_878_930_528_4, 0.458_084_146_815_100_5),
        (0.795_106_878_930_528_3, 0.464_078_122_671_144_73),
        (0.799_117_522_048_654_5, 0.473_776_579_333_970_55),
        (0.795_106_878_930_528_3, 0.483_475_035_996_796_37),
        (0.784_606_878_930_528_4, 0.489_469_011_852_840_6),
        (0.677_661_037_238_211_7, 0.716_698_433_943_938_6),
        (0.7, 0.5),
        (0.775, 0.425_000_000_000_000_04),
        (0.716_341_696_372_258_5, 0.438_511_550_039_540_6),
    ];
    assert!(is_simple(&r));
    if let Err(e) = check_per_cell(&g, &r) {
        panic!("{e}");
    }
}
