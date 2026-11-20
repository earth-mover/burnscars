//! Multipolygon decisions (spec: Multipolygons): one record per (row, col, id).

mod common;

use burnscars::{GeomId, Rule, Spans};
use common::{burn_spans, covered_area, densify, grid, poly};
use geo_types::MultiPolygon;

/// `densify` panics on a cell that two stretches cover.
fn one_record_per_cell(r: &Spans) -> bool {
    densify(r).len() == r.col_start.iter().zip(&r.col_end).map(|(a, b)| (b - a) as usize).sum()
}

#[test]
fn components_sharing_a_cell_give_one_record() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    // Two triangles split cell (5, 5) along its diagonal.
    let mp = MultiPolygon::new(vec![
        poly(&[(5.0, 4.0), (6.0, 4.0), (6.0, 5.0), (5.0, 4.0)], &[]),
        poly(&[(5.0, 4.0), (6.0, 5.0), (5.0, 5.0), (5.0, 4.0)], &[]),
    ]);
    let r = burn_spans(Rule::Fraction, &g, [(GeomId(0), &mp)]).expect("burn");
    assert!(one_record_per_cell(&r));
    assert!((covered_area(&r, &g) - 1.0).abs() < 1e-9);
}

#[test]
fn components_on_both_sides_of_a_window() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let w = g
        .window(
            burnscars::Span::new(burnscars::Row(0), burnscars::Row(10)).expect("r"),
            burnscars::Span::new(burnscars::Col(3), burnscars::Col(7)).expect("c"),
        )
        .expect("window");
    let mp = MultiPolygon::new(vec![
        poly(&[(0.5, 1.5), (3.5, 1.5), (3.5, 3.5), (0.5, 3.5), (0.5, 1.5)], &[]),
        poly(&[(6.5, 5.5), (9.5, 5.5), (9.5, 7.5), (6.5, 7.5), (6.5, 5.5)], &[]),
    ]);
    let r = burn_spans(Rule::Fraction, &w, [(GeomId(0), &mp)]).expect("burn");
    assert!(one_record_per_cell(&r));
    assert!(r.col_start.iter().zip(&r.col_end).all(|(&a, &b)| a >= 3 && b <= 7));
    // In the window: 0.5 col-width of the left part and 0.5 of the right part, 2 rows each.
    assert!((covered_area(&r, &g) - 2.0).abs() < 1e-9);
}

#[test]
fn hole_spanning_two_components_cells() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 20, 20);
    let mp = MultiPolygon::new(vec![
        poly(
            &[(1.0, 1.0), (5.0, 1.0), (5.0, 5.0), (1.0, 5.0), (1.0, 1.0)],
            &[&[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0), (2.0, 2.0)]],
        ),
        poly(&[(5.0, 1.0), (9.0, 1.0), (9.0, 5.0), (5.0, 5.0), (5.0, 1.0)], &[]),
    ]);
    let r = burn_spans(Rule::Fraction, &g, [(GeomId(0), &mp)]).expect("burn");
    assert!(one_record_per_cell(&r));
    assert!((covered_area(&r, &g) - (16.0 - 4.0 + 16.0)).abs() < 1e-9);
}

#[test]
fn approx_adjacent_components_give_one_record_per_cell() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let mp = MultiPolygon::new(vec![
        poly(&[(1.0, 1.0), (5.0, 1.0), (5.0, 5.0), (1.0, 5.0), (1.0, 1.0)], &[]),
        poly(&[(5.0, 1.0), (9.0, 1.0), (9.0, 5.0), (5.0, 5.0), (5.0, 1.0)], &[]),
    ]);
    let r = burn_spans(Rule::Centre, &g, [(GeomId(0), &mp)]).expect("burn");
    assert!(one_record_per_cell(&r));
    assert!((covered_area(&r, &g) - 32.0).abs() < 1e-9);
}
