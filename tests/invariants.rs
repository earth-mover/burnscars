//! Ported from controlledburn `cpp/tests/test_burn.cpp`
//! (commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558, Apache-2.0).

mod common;

use burnscars::{Burner, GeomId, Mode};
use common::{covered_area, grid, poly};
use geo_types::{MultiPolygon, Polygon};

#[expect(clippy::expect_used, reason = "test helper")]
fn burn1(p: &Polygon<f64>, g: &burnscars::Grid) -> burnscars::BurnResult {
    Burner::new(Mode::Coverage).burn(g, [(GeomId(0), p)]).expect("burn")
}

fn shoelace(pts: &[(f64, f64)]) -> f64 {
    pts.windows(2).map(|w| w[0].0 * w[1].1 - w[1].0 * w[0].1).sum::<f64>().abs() / 2.0
}

#[test]
fn aligned_rectangle() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let r = burn1(&poly(&[(2.0, 4.0), (6.0, 4.0), (6.0, 8.0), (2.0, 8.0), (2.0, 4.0)], &[]), &g);
    assert!(r.edges.is_empty());
    assert!(!r.runs.is_empty());
    assert!((covered_area(&r, &g) - 16.0).abs() < 1e-9);
    assert!(r.runs.iter().all(|x| x.id == GeomId(0)));
}

#[test]
fn offset_rectangle_quarters_and_halves() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let r = burn1(&poly(&[(2.5, 4.5), (6.5, 4.5), (6.5, 8.5), (2.5, 8.5), (2.5, 4.5)], &[]), &g);
    assert!((covered_area(&r, &g) - 16.0).abs() < 1e-6);
    let quarters = r.edges.iter().filter(|e| (e.fraction - 0.25).abs() < 1e-6).count();
    let halves = r.edges.iter().filter(|e| (e.fraction - 0.5).abs() < 1e-6).count();
    assert_eq!((quarters, halves), (4, 12));
}

#[test]
fn triangle_on_an_awkward_grid() {
    let g = grid(0.0, 0.0, 100.0, 100.0, 37, 41);
    let pts = [(13.3, 17.7), (88.1, 22.4), (41.9, 79.2), (13.3, 17.7)];
    let area = shoelace(&pts);
    let r = burn1(&poly(&pts, &[]), &g);
    assert!((covered_area(&r, &g) - area).abs() < area * 1e-4);
}

#[test]
fn hole_in_both_orientations() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 20, 20);
    let outer = [(1.0, 1.0), (9.0, 1.0), (9.0, 9.0), (1.0, 9.0), (1.0, 1.0)];
    let ccw = [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)];
    let cw = [(3.0, 3.0), (3.0, 7.0), (7.0, 7.0), (7.0, 3.0), (3.0, 3.0)];
    for hole in [&ccw[..], &cw[..]] {
        let r = burn1(&poly(&outer, &[hole]), &g);
        assert!((covered_area(&r, &g) - 48.0).abs() < 1e-6);
    }
}

#[test]
fn disjoint_multipolygon() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let mp = MultiPolygon::new(vec![
        poly(&[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)], &[]),
        poly(&[(6.0, 6.0), (9.0, 6.0), (9.0, 9.0), (6.0, 9.0), (6.0, 6.0)], &[]),
    ]);
    let r = Burner::new(Mode::Coverage).burn(&g, [(GeomId(0), &mp)]).expect("burn");
    assert!((covered_area(&r, &g) - 13.0).abs() < 1e-6);
}

#[test]
fn beyond_extent_fills_the_grid() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 5, 5);
    let r = burn1(
        &poly(
            &[(-100.0, -100.0), (100.0, -100.0), (100.0, 100.0), (-100.0, 100.0), (-100.0, -100.0)],
            &[],
        ),
        &g,
    );
    assert!(r.edges.is_empty());
    assert!((covered_area(&r, &g) - 100.0).abs() < 1e-9);
}

#[test]
fn straddle_is_clipped() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let r = burn1(&poly(&[(-5.0, 2.5), (4.5, 2.5), (4.5, 6.5), (-5.0, 6.5), (-5.0, 2.5)], &[]), &g);
    assert!((covered_area(&r, &g) - 18.0).abs() < 1e-6);
}

#[test]
fn polygon_inside_one_cell() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let r = burn1(&poly(&[(3.2, 3.2), (3.4, 3.2), (3.4, 3.4), (3.2, 3.4), (3.2, 3.2)], &[]), &g);
    assert!(r.runs.is_empty());
    assert_eq!(r.edges.len(), 1);
    assert!((covered_area(&r, &g) - 0.04).abs() < 1e-6);
}

// Review focus 1.
#[test]
fn outside_window_gives_no_records() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    for (dx, dy) in [(-20.0, 0.0), (20.0, 0.0), (0.0, 20.0), (0.0, -20.0)] {
        let p = poly(
            &[
                (2.0 + dx, 2.0 + dy),
                (4.0 + dx, 2.0 + dy),
                (4.0 + dx, 4.0 + dy),
                (2.0 + dx, 4.0 + dy),
                (2.0 + dx, 2.0 + dy),
            ],
            &[],
        );
        let r = burn1(&p, &g);
        assert!(r.is_empty());
        assert_eq!(r.geom_ptr.len(), 1);
        assert!(r.geom_ptr[0].runs.is_empty());
    }
}

// Review focus 5.
#[test]
fn ids_follow_caller_order() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let a = poly(&[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)], &[]);
    let b = poly(&[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 6.0), (5.0, 5.0)], &[]);
    let r = Burner::new(Mode::Coverage)
        .burn(&g, [(GeomId(9), &a), (GeomId(2), &b), (GeomId(9), &b)])
        .expect("burn");
    let ids: Vec<u32> = r.geom_ptr.iter().map(|gi| gi.id.0).collect();
    assert_eq!(ids, vec![9, 2, 9]);
    assert_eq!(r.runs.id, vec![9, 2, 9]);
}
