//! Benchmarks (spec: Benchmarks). Synthetic cases, both modes.

#![expect(clippy::expect_used, reason = "bench code; a failed burn must stop the run")]

use burnscars::{Burner, GeomId, Grid, Mode};
use criterion::{Criterion, criterion_group, criterion_main};
use geo_types::{LineString, Polygon, coord};
use std::hint::black_box;

#[expect(clippy::cast_precision_loss, reason = "bench indices")]
fn circle(n: usize, cx: f64, cy: f64, r: f64) -> Polygon<f64> {
    let pts: LineString<f64> = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64) / (n as f64);
            coord! { x: cx + r * a.cos(), y: cy + r * a.sin() }
        })
        .collect();
    Polygon::new(pts, vec![])
}

#[expect(clippy::cast_possible_truncation, reason = "bench indices")]
fn bench(c: &mut Criterion) {
    let big = circle(1_000_000, 50.0, 50.0, 40.0);
    let small: Vec<Polygon<f64>> = (0..10_000)
        .map(|i| circle(16, f64::from(i % 100) + 0.5, f64::from(i / 100) + 0.5, 0.4))
        .collect();
    for side in [100_u32, 1_000, 10_000, 31_623] {
        let g = Grid::from_extent(0.0, 0.0, 100.0, 100.0, side, side).expect("grid");
        for mode in [Mode::Coverage, Mode::Approx] {
            c.bench_function(&format!("circle_1e6/{side}/{mode:?}"), |b| {
                let mut burner = Burner::new(mode);
                b.iter(|| black_box(burner.burn(&g, [(GeomId(0), &big)]).expect("burn")));
            });
            c.bench_function(&format!("small_1e4/{side}/{mode:?}"), |b| {
                let mut burner = Burner::new(mode);
                b.iter(|| {
                    let items = small.iter().enumerate().map(|(i, p)| (GeomId(i as u32), p));
                    black_box(burner.burn(&g, items).expect("burn"))
                });
            });
        }
    }
}

criterion_group!(benches, bench);
criterion_main!(benches);
