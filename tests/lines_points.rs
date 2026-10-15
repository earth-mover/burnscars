//! Line and point cases from controlledburn `cpp/tests/test_burn.cpp` and
//! `tests/testthat/test-edge-zoo.R` (commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558).

mod common;

use burnscars::{BurnError, Burner, Col, GeomId, Mode, Row, Span};
use common::{grid, ls};
use geo_types::{MultiLineString, MultiPoint, Point};

fn total_length(r: &burnscars::BurnResult) -> f64 {
    r.lines.iter().map(|l| f64::from(l.length)).sum()
}

#[expect(clippy::expect_used, reason = "test helper: a failed burn is a test bug")]
fn burn_line(pts: &[(f64, f64)], g: &burnscars::Grid, mode: Mode) -> burnscars::BurnResult {
    Burner::new(mode).burn(g, [(GeomId(0), &ls(pts))]).expect("burn")
}

#[test]
fn diagonal_line_length_is_conserved_in_both_modes() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    for mode in [Mode::Coverage, Mode::Approx] {
        let r = burn_line(&[(0.5, 0.5), (9.5, 7.5)], &g, mode);
        assert!(r.runs.is_empty() && r.edges.is_empty());
        assert!((total_length(&r) - 9.0_f64.hypot(7.0)).abs() < 1e-4);
    }
}

#[test]
fn lines_on_cell_walls_and_through_corners() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    assert!(
        (total_length(&burn_line(&[(1.0, 5.0), (9.0, 5.0)], &g, Mode::Coverage)) - 8.0).abs()
            < 1e-6
    );
    assert!(
        (total_length(&burn_line(&[(5.0, 1.0), (5.0, 9.0)], &g, Mode::Coverage)) - 8.0).abs()
            < 1e-6
    );
    let g3 = grid(0.0, 0.0, 3.0, 3.0, 3, 3);
    assert!(
        (total_length(&burn_line(&[(0.5, 0.5), (2.5, 2.5)], &g3, Mode::Coverage)) - 8.0_f64.sqrt())
            .abs()
            < 1e-6
    );
}

#[test]
fn short_lines() {
    let g3 = grid(0.0, 0.0, 3.0, 3.0, 3, 3);
    let r = burn_line(&[(0.2, 0.2), (0.8, 0.8)], &g3, Mode::Coverage);
    assert_eq!(r.lines.len(), 1);
    assert!((total_length(&r) - 0.72_f64.sqrt()).abs() < 1e-6);
    let r = burn_line(&[(0.5, 0.5), (0.501, 0.501)], &g3, Mode::Coverage);
    assert!((total_length(&r) - 2e-6_f64.sqrt()).abs() < 1e-7);
    let r = burn_line(&[(0.5, 0.5), (0.5, 0.5)], &g3, Mode::Coverage);
    assert!(r.lines.is_empty());
}

#[test]
fn long_lines() {
    let g = grid(0.0, 0.0, 1000.0, 1000.0, 1000, 1000);
    assert!(
        (total_length(&burn_line(&[(0.5, 500.5), (999.5, 500.5)], &g, Mode::Coverage)) - 999.0)
            .abs()
            < 1e-3
    );
    let d = burn_line(&[(0.5, 0.5), (999.5, 999.5)], &g, Mode::Coverage);
    assert!((total_length(&d) - 999.0 * 2.0_f64.sqrt()).abs() < 1e-2);
}

#[test]
fn line_parts_merge_into_one_record_per_cell() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let m =
        MultiLineString::new(vec![ls(&[(0.2, 0.5), (0.8, 0.5)]), ls(&[(0.5, 0.2), (0.5, 0.8)])]);
    let r = Burner::new(Mode::Coverage).burn(&g, [(GeomId(0), &m)]).expect("burn");
    assert_eq!(r.lines.len(), 1);
    assert!((total_length(&r) - 1.2).abs() < 1e-6);
}

#[test]
fn short_line_is_an_error() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let r = Burner::new(Mode::Coverage).burn(&g, [(GeomId(4), &ls(&[(1.0, 1.0)]))]);
    assert!(matches!(r, Err(BurnError::ShortLine { .. })));
}

#[test]
fn points_bin_into_cells() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let mp =
        MultiPoint::new(vec![Point::new(0.5, 9.5), Point::new(9.5, 0.5), Point::new(15.0, 5.0)]);
    let r = Burner::new(Mode::Coverage).burn(&g, [(GeomId(0), &mp)]).expect("burn");
    let cells: Vec<(u32, u32)> = r.points.iter().map(|p| (p.row.0, p.col.0)).collect();
    assert_eq!(cells, vec![(0, 0), (9, 9)]);
}

#[test]
fn point_on_the_far_edge_follows_exactextract_and_chunks_agree() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let p = Point::new(10.0, 0.0); // right and bottom edge of the full grid
    let r = Burner::new(Mode::Coverage).burn(&g, [(GeomId(0), &p)]).expect("burn");
    assert_eq!(r.points.iter().map(|q| (q.row.0, q.col.0)).collect::<Vec<_>>(), vec![(9, 9)]);
    // A point on an interior seam goes into exactly one chunk.
    let seam = Point::new(5.0, 5.0);
    let mut n = 0;
    for (rs, cs) in [((0, 5), (0, 5)), ((0, 5), (5, 10)), ((5, 10), (0, 5)), ((5, 10), (5, 10))] {
        let w = g
            .window(
                Span::new(Row(rs.0), Row(rs.1)).expect("r"),
                Span::new(Col(cs.0), Col(cs.1)).expect("c"),
            )
            .expect("window");
        n += Burner::new(Mode::Coverage).burn(&w, [(GeomId(0), &seam)]).expect("burn").points.len();
    }
    assert_eq!(n, 1);
}
