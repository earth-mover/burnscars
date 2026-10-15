//! Partition equality (spec: Output order and determinism): feature-split and tile-split
//! burns give the same expanded form as the full burn, bit for bit.

mod common;

use burnscars::{BurnResult, Burner, Col, GeomId, Grid, Lattice, Mode, Row, Span};
use common::{expanded, ls};
use geo_types::{Geometry, LineString, MultiPoint, Point, Polygon};
use proptest::prelude::*;

/// Lattice origin and step: x0 in a small range, y0 = 32.
#[derive(Clone, Copy, Debug)]
struct Lat {
    x0: f64,
    dx: f64,
    dy: f64,
}

/// Snap to the lattice (y0 = 32): mode bit 0 snaps x, bit 1 snaps y.
fn snap_xy(l: Lat, x: f64, y: f64, mode: u8) -> (f64, f64) {
    let Lat { x0, dx, dy } = l;
    let sx = if mode & 1 == 0 { x } else { x0 + ((x - x0) / dx).round() * dx };
    let sy = if mode & 2 == 0 { y } else { 32.0 - ((32.0 - y) / dy).round() * dy };
    (sx, sy)
}

/// A random coordinate with a random snap mode (none, x seam, y seam, corner).
fn snapped_pt(l: Lat) -> impl Strategy<Value = (f64, f64)> {
    (0.0..32.0_f64, 0.0..32.0_f64, 0_u8..4).prop_map(move |(x, y, m)| snap_xy(l, x, y, m))
}

/// True if the angles about (cx, cy) strictly increase and every gap is < pi: the ring
/// is then a simple star (spec: valid input only).
fn is_star(cx: f64, cy: f64, pts: &[(f64, f64)]) -> bool {
    let tau = std::f64::consts::TAU;
    let ang: Vec<f64> = pts.iter().map(|&(x, y)| (y - cy).atan2(x - cx).rem_euclid(tau)).collect();
    let increasing = ang.windows(2).all(|w| w[1] > w[0]);
    let wrap = ang[0] + tau - ang[ang.len() - 1];
    let gaps = ang.windows(2).map(|w| w[1] - w[0]).chain(std::iter::once(wrap));
    increasing && gaps.into_iter().all(|g| g < std::f64::consts::PI - 1e-9)
}

/// Distance from (cx, cy) to the nearest edge of the closed ring `pts` (first == last
/// not required).
fn centre_gap(cx: f64, cy: f64, pts: &[(f64, f64)]) -> f64 {
    (0..pts.len())
        .map(|i| {
            let ((ax, ay), (bx, by)) = (pts[i], pts[(i + 1) % pts.len()]);
            let (ex, ey) = (bx - ax, by - ay);
            let t = (((cx - ax) * ex + (cy - ay) * ey) / (ex * ex + ey * ey)).clamp(0.0, 1.0);
            (cx - (ax + t * ex)).hypot(cy - (ay + t * ey))
        })
        .fold(f64::INFINITY, f64::min)
}

/// A star polygon around (cx, cy), some vertices on the lattice. Half of the draws get
/// one hole: a smaller star about the same centre, kept inside by the filter.
#[expect(clippy::cast_precision_loss, reason = "test indices")]
fn star_polygon(l: Lat) -> impl Strategy<Value = Polygon<f64>> {
    let ring = move |cx: f64, cy: f64, radii: &[(f64, bool)]| -> Vec<(f64, f64)> {
        let n = radii.len();
        radii
            .iter()
            .enumerate()
            .map(|(i, &(r, snap))| {
                let a = std::f64::consts::TAU * (i as f64) / (n as f64);
                let (x, y) = (cx + r * a.cos(), cy + r * a.sin());
                if snap { snap_xy(l, x, y, 3) } else { (x, y) }
            })
            .collect()
    };
    (
        4.0..28.0_f64,
        4.0..28.0_f64,
        prop::collection::vec((1.0..6.0_f64, any::<bool>()), 3..24),
        prop::option::of(prop::collection::vec((0.2..0.9_f64, any::<bool>()), 3..12)),
    )
        .prop_filter_map("simple star with an inner hole", move |(cx, cy, radii, hole)| {
            let pts = ring(cx, cy, &radii);
            if !is_star(cx, cy, &pts) {
                return None;
            }
            let Some(h) = hole else {
                return Some(Polygon::new(LineString::from(pts), vec![]));
            };
            // Hole radii are fractions of the largest disc about the centre inside the star.
            let d = centre_gap(cx, cy, &pts);
            let hole_radii: Vec<(f64, bool)> = h.iter().map(|&(f, s)| (f * d * 0.9, s)).collect();
            let hole_pts = ring(cx, cy, &hole_radii);
            let inside = hole_pts.iter().all(|&(x, y)| (x - cx).hypot(y - cy) < d * 0.95);
            (d > 0.3 && inside && is_star(cx, cy, &hole_pts))
                .then(|| Polygon::new(LineString::from(pts), vec![LineString::from(hole_pts)]))
        })
}

/// A polyline with some vertices on seams. An optional first segment lies along a
/// lattice wall (both ends on one x edge, or one y edge). Repeated vertices are removed.
fn polyline(l: Lat) -> impl Strategy<Value = Option<LineString<f64>>> {
    let wall = prop::option::of((any::<bool>(), 0.0..32.0_f64, 0.0..32.0_f64, 0.0..32.0_f64));
    (prop::collection::vec(snapped_pt(l), 1..8), wall).prop_map(move |(rest, wall)| {
        let mut v: Vec<(f64, f64)> = match wall {
            Some((true, x, a, b)) => vec![snap_xy(l, x, a, 1), snap_xy(l, x, b, 1)],
            Some((false, y, a, b)) => vec![snap_xy(l, a, y, 2), snap_xy(l, b, y, 2)],
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
type Case = (f64, f64, f64, u32, u32, Vec<Polygon<f64>>, Vec<LineString<f64>>, Vec<(f64, f64)>);

/// Lattice origin x0 and step first, so the geometry can snap to it and the grid (at most
/// ceil(32 / step) cells per axis, full size in 3 of 4 draws) usually covers the geometry.
#[expect(clippy::cast_possible_truncation, clippy::cast_sign_loss, reason = "small positive")]
fn case() -> impl Strategy<Value = Case> {
    (0.0..1.5_f64, 0.3..2.0_f64, 0.3..2.0_f64).prop_flat_map(|(x0, dx, dy)| {
        let l = Lat { x0, dx, dy };
        let (nc, nr) = ((32.0 / dx).ceil() as u32, (32.0 / dy).ceil() as u32);
        (
            Just(x0),
            Just(dx),
            Just(dy),
            prop_oneof![3 => Just(nc), 1 => 1..=nc],
            prop_oneof![3 => Just(nr), 1 => 1..=nr],
            prop::collection::vec(star_polygon(l), 1..4),
            prop::collection::vec(polyline(l), 0..3)
                .prop_map(|v| v.into_iter().flatten().collect::<Vec<_>>()),
            prop::collection::vec(snapped_pt(l), 0..6),
        )
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]
    #[test]
    #[expect(clippy::cast_possible_truncation, reason = "test indices")]
    fn tiles_and_feature_splits_equal_the_full_burn(
        (x0, dx, dy, ncol, nrow, polys, lines, pts) in case(),
        rchunk in 1_u32..12, cchunk in 1_u32..12, split in 0_usize..8,
    ) {
        let lat = Lattice::new(x0, 32.0, dx, dy).expect("lattice");
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
