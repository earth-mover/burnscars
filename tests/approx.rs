//! Approx cases from controlledburn `cpp/tests/test_burn.cpp` and
//! `tests/testthat/test-lattice-snap.R` (commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558).

mod common;

use burnscars::{Burner, GeomId, Mode};
use common::{covered_area, grid, poly};
use geo_types::Polygon;

#[expect(clippy::expect_used, reason = "test helper")]
fn burn(p: &Polygon<f64>, g: &burnscars::Grid) -> burnscars::BurnResult {
    Burner::new(Mode::Approx).burn(g, [(GeomId(0), p)]).expect("burn")
}

#[test]
fn aligned_rectangle() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let r = burn(&poly(&[(2.0, 4.0), (6.0, 4.0), (6.0, 8.0), (2.0, 8.0), (2.0, 4.0)], &[]), &g);
    assert!(r.edges.is_empty());
    assert!((covered_area(&r, &g) - 16.0).abs() < 1e-9);
}

#[test]
fn offset_rectangle_is_left_and_top_inclusive() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let r = burn(&poly(&[(2.5, 4.5), (6.5, 4.5), (6.5, 8.5), (2.5, 8.5), (2.5, 4.5)], &[]), &g);
    assert!(r.edges.is_empty());
    // Cols 2..6 (centre 2.5 on the left edge is in, 6.5 on the right edge is out);
    // rows 1..5 (centre 8.5 on the top edge is in, 4.5 on the bottom edge is out).
    let cells: Vec<(u32, u32)> =
        r.runs.iter().flat_map(|x| x.cols.range().map(move |c| (x.row.0, c))).collect();
    assert_eq!(cells.len(), 16);
    assert!(cells.iter().all(|&(row, col)| (1..5).contains(&row) && (2..6).contains(&col)));
}

#[test]
fn beyond_extent_fills_the_grid() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 5, 5);
    let r = burn(
        &poly(
            &[(-100.0, -100.0), (100.0, -100.0), (100.0, 100.0), (-100.0, 100.0), (-100.0, -100.0)],
            &[],
        ),
        &g,
    );
    assert!((covered_area(&r, &g) - 100.0).abs() < 1e-9);
}

#[test]
fn hole_is_respected() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 20, 20);
    let r = burn(
        &poly(
            &[(1.0, 1.0), (9.0, 1.0), (9.0, 9.0), (1.0, 9.0), (1.0, 1.0)],
            &[&[(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0), (3.0, 3.0)]],
        ),
        &g,
    );
    let a = covered_area(&r, &g);
    assert!(a > 40.0 && a < 56.0, "{a}");
    assert!((a / 0.25 - (a / 0.25).round()).abs() < 1e-9);
}

// controlledburn test-lattice-snap.R: a vertex on a row centre line.
#[test]
fn vertex_on_a_row_centre_does_not_fill_the_row() {
    let g = grid(-180.0, -90.0, -175.0, -80.66666666666667, 15, 28);
    let p = poly(
        &[
            (-179.89654756364214, -88.83333333333334),
            (-180.94076559960791, -92.3680897750401),
            (-179.6, -88.518152273978),
            (-179.89654756364214, -88.83333333333334),
        ],
        &[],
    );
    let r = burn(&p, &g);
    assert!(r.runs.col_end.iter().all(|&e| e <= 2), "{:?}", r.runs);
    assert_eq!(r.runs.len(), 2);
}

#[test]
fn outside_window_gives_no_records() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let r =
        burn(&poly(&[(20.0, 2.0), (24.0, 2.0), (24.0, 4.0), (20.0, 4.0), (20.0, 2.0)], &[]), &g);
    assert!(r.is_empty());
    assert_eq!(r.geom_ptr.len(), 1);
}
