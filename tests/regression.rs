//! Regression cases from controlledburn `tests/testthat/test-lattice-snap.R`,
//! `test-edge-zoo.R`, `test-shared-boundary.R` (commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558).

mod common;

use burnscars::{Burner, GeomId, Mode};
use common::{covered_area, grid, poly};
use geo_types::{MultiPolygon, Polygon};

#[expect(clippy::expect_used, reason = "test helper")]
fn burn(p: &Polygon<f64>, g: &burnscars::Grid) -> burnscars::BurnResult {
    Burner::new(Mode::Coverage).burn(g, [(GeomId(0), p)]).expect("burn")
}

fn tri(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs() / 2.0
}

// controlledburn fix 1f9a53c: vertices a few ulp off a cell edge.
#[test]
fn lattice_snap_keeps_the_walk_on_the_grid() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 100, 100);
    let cases = [
        ([(0.3, 0.3), (4.7, 0.3), (2.5, 6.2), (0.3, 0.3)], tri((0.3, 0.3), (4.7, 0.3), (2.5, 6.2))),
        ([(4.7, 0.3), (2.5, 6.2), (0.3, 0.3), (4.7, 0.3)], tri((0.3, 0.3), (4.7, 0.3), (2.5, 6.2))),
        (
            [(0.3, 0.3), (4.7, 0.3), (2.5, 6.25), (0.3, 0.3)],
            tri((0.3, 0.3), (4.7, 0.3), (2.5, 6.25)),
        ),
    ];
    for (pts, area) in cases {
        let r = burn(&poly(&pts, &[]), &g);
        assert!(!r.runs.is_empty());
        assert!((covered_area(&r, &g) - area).abs() < 1e-6 * area, "{pts:?}");
    }
}

// controlledburn fix 0efe32f: a row-centre crossing on a column wall.
#[test]
fn row_centre_crossing_on_a_wall_keeps_its_winding() {
    let pts = [
        (0.8, 1.15),
        (2.0, 1.75),
        (2.85, 1.5),
        (3.45, 1.9),
        (4.25, 1.45),
        (3.95, 1.1),
        (2.6, 0.8),
        (2.15, 1.15),
        (1.3, 0.7),
        (0.8, 1.15),
    ];
    for n in [20, 100, 200, 1000] {
        let g = grid(0.0, 0.0, 10.0, 10.0, n, n);
        let r = burn(&poly(&pts, &[]), &g);
        assert!((covered_area(&r, &g) - 2.12875).abs() < 1e-6 * 2.12875, "ncol {n}");
    }
}

#[test]
fn edge_on_a_row_boundary() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let r = burn(&poly(&[(2.0, 2.0), (8.0, 2.0), (8.0, 5.0), (2.0, 5.0), (2.0, 2.0)], &[]), &g);
    assert!(r.edges.is_empty());
    assert!((covered_area(&r, &g) - 18.0).abs() < 1e-9);
}

#[test]
fn vertex_on_a_grid_node_and_on_an_edge_midpoint() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let r = burn(&poly(&[(2.0, 2.0), (8.0, 2.0), (5.0, 5.0), (2.0, 2.0)], &[]), &g);
    assert!((covered_area(&r, &g) - 9.0).abs() < 1e-9);
    let r = burn(&poly(&[(2.0, 2.0), (8.0, 2.0), (5.0, 3.5), (2.0, 2.0)], &[]), &g);
    assert!((covered_area(&r, &g) - 4.5).abs() < 1e-9);
}

#[test]
fn corner_touch_multipolygon_is_not_over_counted() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let mp = MultiPolygon::new(vec![
        poly(&[(2.0, 2.0), (5.0, 2.0), (5.0, 5.0), (2.0, 5.0), (2.0, 2.0)], &[]),
        poly(&[(5.0, 5.0), (8.0, 5.0), (8.0, 8.0), (5.0, 8.0), (5.0, 5.0)], &[]),
    ]);
    let r = Burner::new(Mode::Coverage).burn(&g, [(GeomId(0), &mp)]).expect("burn");
    assert!((covered_area(&r, &g) - 18.0).abs() < 1e-9);
}

#[test]
fn collinear_vertices_do_not_change_the_output() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let a = burn(
        &poly(&[(2.0, 2.0), (5.0, 2.0), (8.0, 2.0), (8.0, 8.0), (2.0, 8.0), (2.0, 2.0)], &[]),
        &g,
    );
    let b = burn(&poly(&[(2.0, 2.0), (8.0, 2.0), (8.0, 8.0), (2.0, 8.0), (2.0, 2.0)], &[]), &g);
    assert_eq!(common::expanded(&a), common::expanded(&b));
}

// Review focus 3.
#[test]
fn repeated_vertices_do_not_change_coverage() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let a = burn(
        &poly(
            &[(2.5, 2.5), (6.5, 2.5), (6.5, 2.5), (6.5, 6.5), (2.5, 6.5), (2.5, 6.5), (2.5, 2.5)],
            &[],
        ),
        &g,
    );
    let b = burn(&poly(&[(2.5, 2.5), (6.5, 2.5), (6.5, 6.5), (2.5, 6.5), (2.5, 2.5)], &[]), &g);
    assert!((covered_area(&a, &g) - covered_area(&b, &g)).abs() < 1e-12);
}

// Review focus 4.
#[test]
fn sliver_fractions_stay_in_range() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let pts = [(1.1, 1.1), (8.9, 1.1 + 1e-9), (8.9, 1.1 + 2e-9), (1.1, 1.1)];
    let r = burn(&poly(&pts, &[]), &g);
    assert!(r.edges.iter().all(|e| e.fraction > 0.0 && e.fraction < 1.0));
    let sub = poly(&[(1.1, 1.1), (1.2, 1.1), (1.15, 1.2), (1.1, 1.1)], &[]);
    let r = burn(&sub, &g);
    assert_eq!(r.edges.len(), 1);
    assert!(r.edges.iter().all(|e| e.fraction > 0.0 && e.fraction < 1.0));
}

#[test]
fn huge_finite_coordinates_do_not_panic() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let p = poly(
        &[(-1e300, -1e300), (1e300, -1e300), (1e300, 1e300), (-1e300, 1e300), (-1e300, -1e300)],
        &[],
    );
    let r = Burner::new(Mode::Coverage).burn(&g, [(GeomId(0), &p)]);
    assert!(r.is_ok(), "{r:?}");
}
