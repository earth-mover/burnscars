//! Multipolygon decisions (spec: Multipolygons): one record per (row, col, id).

mod common;

use std::collections::BTreeSet;

use burnscars::{BurnResult, Burner, GeomId, Mode};
use common::{covered_area, grid, poly};
use geo_types::MultiPolygon;

fn one_record_per_cell(r: &BurnResult) -> bool {
    let mut seen = BTreeSet::new();
    let runs = r.runs.iter().flat_map(|x| x.cols.range().map(move |c| (x.id.0, x.row.0, c)));
    let edges = r.edges.iter().map(|e| (e.id.0, e.row.0, e.col.0));
    runs.chain(edges).all(|k| seen.insert(k))
}

#[test]
fn components_sharing_a_cell_give_one_record() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    // Two triangles split cell (5, 5) along its diagonal.
    let mp = MultiPolygon::new(vec![
        poly(&[(5.0, 4.0), (6.0, 4.0), (6.0, 5.0), (5.0, 4.0)], &[]),
        poly(&[(5.0, 4.0), (6.0, 5.0), (5.0, 5.0), (5.0, 4.0)], &[]),
    ]);
    let r = Burner::new(Mode::Coverage).burn(&g, [(GeomId(0), &mp)]).expect("burn");
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
    let r = Burner::new(Mode::Coverage).burn(&w, [(GeomId(0), &mp)]).expect("burn");
    assert!(one_record_per_cell(&r));
    assert!(r.runs.iter().all(|x| x.cols.start().0 >= 3 && x.cols.end().0 <= 7));
    assert!(r.edges.iter().all(|e| (3..7).contains(&e.col.0)));
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
    let r = Burner::new(Mode::Coverage).burn(&g, [(GeomId(0), &mp)]).expect("burn");
    assert!(one_record_per_cell(&r));
    assert!((covered_area(&r, &g) - (16.0 - 4.0 + 16.0)).abs() < 1e-9);
}
