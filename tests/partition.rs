//! Partition equality (spec: Output order and determinism): feature-split and tile-split
//! burns give the same expanded form as the full burn, bit for bit.

mod common;

use burnscars::{BurnResult, Burner, Col, GeomId, Grid, Lattice, Mode, Row, Span};
use common::{expanded, ls};
use geo_types::{Geometry, LineString, MultiPoint, Point, Polygon};
use proptest::prelude::*;

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
                    if snap {
                        ((x / dx).round() * dx, 32.0 - ((32.0 - y) / dy).round() * dy)
                    } else {
                        (x, y)
                    }
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

fn polyline() -> impl Strategy<Value = LineString<f64>> {
    prop::collection::vec((0.0..32.0_f64, 0.0..32.0_f64), 2..8).prop_map(|v| ls(&v))
}

fn chunks(n: u32, size: u32) -> Vec<(u32, u32)> {
    (0..n).step_by(size as usize).map(|s| (s, (s + size).min(n))).collect()
}

#[expect(clippy::cast_possible_truncation, clippy::expect_used, reason = "test helper")]
fn burn_all(mode: Mode, g: &Grid, geoms: &[Geometry<f64>]) -> BurnResult {
    let items = geoms.iter().enumerate().map(|(i, x)| (GeomId(i as u32), x));
    Burner::new(mode).burn(g, items).expect("burn")
}

/// Lattice resolution first, so the polygons can snap to that lattice.
fn case() -> impl Strategy<Value = (f64, f64, Vec<Polygon<f64>>)> {
    (0.3..2.0_f64, 0.3..2.0_f64).prop_flat_map(|(dx, dy)| {
        (Just(dx), Just(dy), prop::collection::vec(star_polygon(dx, dy), 1..4))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]
    #[test]
    #[expect(clippy::cast_possible_truncation, reason = "test indices")]
    fn tiles_and_feature_splits_equal_the_full_burn(
        ncol in 1_u32..40, nrow in 1_u32..40, (dx, dy, polys) in case(),
        lines in prop::collection::vec(polyline(), 0..3),
        pts in prop::collection::vec((0.0..32.0_f64, 0.0..32.0_f64), 0..6),
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
