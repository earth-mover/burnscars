//! Partition equality (spec: Output order and determinism): feature-split and tile-split
//! burns give the same expanded form as the full burn, bit for bit.

mod common;

use burnscars::{BurnResult, Burner, Col, GeomId, Grid, Lattice, Mode, Row, Span};
use common::{expanded, ls};
use geo_types::{Geometry, LineString, MultiPoint, Point, Polygon};
use proptest::prelude::*;

/// Snap to the lattice (x0 = 0, y0 = 32): mode bit 0 snaps x, bit 1 snaps y.
fn snap_xy(dx: f64, dy: f64, x: f64, y: f64, mode: u8) -> (f64, f64) {
    let sx = if mode & 1 == 0 { x } else { (x / dx).round() * dx };
    let sy = if mode & 2 == 0 { y } else { 32.0 - ((32.0 - y) / dy).round() * dy };
    (sx, sy)
}

/// A random coordinate with a random snap mode (none, x seam, y seam, corner).
fn snapped_pt(dx: f64, dy: f64) -> impl Strategy<Value = (f64, f64)> {
    (0.0..32.0_f64, 0.0..32.0_f64, 0_u8..4).prop_map(move |(x, y, m)| snap_xy(dx, dy, x, y, m))
}

/// A star polygon around (cx, cy); some vertices snap to the lattice (x0 = 0, y0 = 32).
/// The filter keeps it simple: angles about the centre strictly increase and every
/// angular gap is < pi, so the polygon is valid (spec: valid input only).
#[expect(clippy::cast_precision_loss, reason = "test indices")]
fn star_polygon(dx: f64, dy: f64) -> impl Strategy<Value = Polygon<f64>> {
    (4.0..28.0_f64, 4.0..28.0_f64, prop::collection::vec((1.0..6.0_f64, any::<bool>()), 3..24))
        .prop_map(move |(cx, cy, radii)| {
            let n = radii.len();
            let pts: Vec<(f64, f64)> = radii
                .iter()
                .enumerate()
                .map(|(i, &(r, snap))| {
                    let a = std::f64::consts::TAU * (i as f64) / (n as f64);
                    let (x, y) = (cx + r * a.cos(), cy + r * a.sin());
                    if snap { snap_xy(dx, dy, x, y, 3) } else { (x, y) }
                })
                .collect();
            (cx, cy, pts)
        })
        .prop_filter("simple star", |(cx, cy, pts)| {
            let ang: Vec<f64> = pts
                .iter()
                .map(|&(x, y)| (y - cy).atan2(x - cx).rem_euclid(std::f64::consts::TAU))
                .collect();
            let increasing = ang.windows(2).all(|w| w[1] > w[0]);
            let gaps = ang
                .windows(2)
                .map(|w| w[1] - w[0])
                .chain(std::iter::once(ang[0] + std::f64::consts::TAU - ang[ang.len() - 1]));
            increasing && gaps.into_iter().all(|g| g < std::f64::consts::PI - 1e-9)
        })
        .prop_map(|(_, _, pts)| Polygon::new(LineString::from(pts), vec![]))
}

/// A polyline with some vertices on seams. An optional first segment lies along a
/// lattice wall (both ends on one x edge, or one y edge). Repeated vertices are removed.
fn polyline(dx: f64, dy: f64) -> impl Strategy<Value = Option<LineString<f64>>> {
    let wall = prop::option::of((any::<bool>(), 0.0..32.0_f64, 0.0..32.0_f64, 0.0..32.0_f64));
    (prop::collection::vec(snapped_pt(dx, dy), 1..8), wall).prop_map(move |(rest, wall)| {
        let mut v: Vec<(f64, f64)> = match wall {
            Some((true, x, a, b)) => vec![snap_xy(dx, dy, x, a, 1), snap_xy(dx, dy, x, b, 1)],
            Some((false, y, a, b)) => vec![snap_xy(dx, dy, a, y, 2), snap_xy(dx, dy, b, y, 2)],
            None => vec![],
        };
        v.extend(rest);
        v.dedup();
        (v.len() >= 2).then(|| ls(&v))
    })
}

fn chunks(n: u32, size: u32) -> Vec<(u32, u32)> {
    (0..n).step_by(size as usize).map(|s| (s, (s + size).min(n))).collect()
}

#[expect(clippy::cast_possible_truncation, clippy::expect_used, reason = "test helper")]
fn burn_all(mode: Mode, g: &Grid, geoms: &[Geometry<f64>]) -> BurnResult {
    let items = geoms.iter().enumerate().map(|(i, x)| (GeomId(i as u32), x));
    Burner::new(mode).burn(g, items).expect("burn")
}

/// One drawn case: lattice step, grid size, polygons, lines and points.
type Case = (f64, f64, u32, u32, Vec<Polygon<f64>>, Vec<LineString<f64>>, Vec<(f64, f64)>);

/// Lattice step first, so the geometry can snap to it and the grid (at most
/// ceil(32 / step) cells per axis, full size in 3 of 4 draws) usually covers the geometry.
#[expect(clippy::cast_possible_truncation, clippy::cast_sign_loss, reason = "small positive")]
fn case() -> impl Strategy<Value = Case> {
    (0.3..2.0_f64, 0.3..2.0_f64).prop_flat_map(|(dx, dy)| {
        let (nc, nr) = ((32.0 / dx).ceil() as u32, (32.0 / dy).ceil() as u32);
        (
            Just(dx),
            Just(dy),
            prop_oneof![3 => Just(nc), 1 => 1..=nc],
            prop_oneof![3 => Just(nr), 1 => 1..=nr],
            prop::collection::vec(star_polygon(dx, dy), 1..4),
            prop::collection::vec(polyline(dx, dy), 0..3)
                .prop_map(|v| v.into_iter().flatten().collect::<Vec<_>>()),
            prop::collection::vec(snapped_pt(dx, dy), 0..6),
        )
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]
    #[test]
    #[expect(clippy::cast_possible_truncation, reason = "test indices")]
    fn tiles_and_feature_splits_equal_the_full_burn(
        (dx, dy, ncol, nrow, polys, lines, pts) in case(),
        rchunk in 1_u32..12, cchunk in 1_u32..12, split in 0_usize..8,
    ) {
        let lat = Lattice::new(0.0, 32.0, dx, dy).expect("lattice");
        let g = Grid::new(lat, Span::new(Row(0), Row(nrow)).expect("r"), Span::new(Col(0), Col(ncol)).expect("c")).expect("grid");
        let mut geoms: Vec<Geometry<f64>> = polys.into_iter().map(Geometry::Polygon).collect();
        geoms.extend(lines.into_iter().map(Geometry::LineString));
        geoms.push(Geometry::MultiPoint(MultiPoint::new(pts.iter().map(|&(x, y)| Point::new(x, y)).collect())));

        for mode in [Mode::Coverage, Mode::Approx] {
            let full = expanded(&burn_all(mode, &g, &geoms));

            // Feature split: two workers, ids kept, results concatenated.
            let k = split.min(geoms.len());
            let mut b = Burner::new(mode);
            let mut out = BurnResult::default();
            let items = geoms.iter().enumerate().map(|(i, x)| (GeomId(i as u32), x));
            b.burn_into(&g, items.clone().take(k), &mut out).expect("a");
            Burner::new(mode).burn_into(&g, items.skip(k), &mut out).expect("b");
            prop_assert_eq!(&expanded(&out), &full);

            // Tile split: every chunk made with window(), outputs merged.
            let mut tiled = std::collections::BTreeMap::new();
            for (r0, r1) in chunks(nrow, rchunk) {
                for (c0, c1) in chunks(ncol, cchunk) {
                    let w = g.window(Span::new(Row(r0), Row(r1)).expect("r"), Span::new(Col(c0), Col(c1)).expect("c")).expect("w");
                    for (key, mut v) in expanded(&burn_all(mode, &w, &geoms)) {
                        tiled.entry(key).or_insert_with(Vec::new).append(&mut v);
                    }
                }
            }
            for v in tiled.values_mut() { v.sort_unstable(); }
            prop_assert_eq!(&tiled, &full);
        }
    }
}
