//! Error contract: no records of a failed geometry, error kinds and `Display` text.

mod common;

use burnscars::{At, BurnError, BurnResult, Burner, GeomId, GeomKind, Mode};
use common::{grid, ls, poly};
use geo_types::{Geometry, GeometryCollection};

#[test]
fn failed_geometry_leaves_no_records() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let sq = poly(&[(1.0, 1.0), (4.0, 1.0), (4.0, 4.0), (1.0, 4.0), (1.0, 1.0)], &[]);
    let bad = poly(&[(1.0, 1.0), (f64::NAN, 1.0), (4.0, 4.0), (1.0, 4.0), (1.0, 1.0)], &[]);
    let mut out = BurnResult::default();
    let r =
        Burner::new(Mode::Coverage).burn_into(&g, [(GeomId(0), &sq), (GeomId(1), &bad)], &mut out);
    assert!(matches!(r, Err(BurnError::NonFinite { .. })));
    assert_eq!(out.geom_ptr.len(), 1);
    let first = &out.geom_ptr[0];
    assert_eq!(first.id, GeomId(0));
    assert_eq!(out.runs.len(), first.runs.end);
    assert_eq!(out.edges.len(), first.edges.end);
    assert_eq!(out.lines.len(), first.lines.end);
    assert_eq!(out.points.len(), first.points.end);
}

#[test]
fn short_ring_short_line_and_unsupported() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 10, 10);
    let mut b = Burner::new(Mode::Coverage);
    let tri = poly(&[(1.0, 1.0), (4.0, 1.0), (1.0, 1.0)], &[]);
    assert!(matches!(b.burn(&g, [(GeomId(0), &tri)]), Err(BurnError::ShortRing { .. })));
    let one = ls(&[(1.0, 1.0)]);
    assert!(matches!(b.burn(&g, [(GeomId(0), &one)]), Err(BurnError::ShortLine { .. })));
    let gc: Geometry<f64> = Geometry::GeometryCollection(GeometryCollection::new_from(vec![]));
    let r = b.burn(&g, [(GeomId(0), &gc)]);
    assert!(matches!(r, Err(BurnError::Unsupported { kind: GeomKind::GeometryCollection, .. })));
}

#[test]
fn display_text() {
    let at = At { id: GeomId(7), part: 1, ring: 2, vertex: 0 };
    let s = BurnError::WalkDiverged { at, coords: 99 }.to_string();
    assert!(s.contains("geometry 7, part 1, ring 2"), "{s}");
    assert!(s.contains("walk did not converge after 99"), "{s}");
    let s = BurnError::Internal { id: GeomId(7), detail: "x" }.to_string();
    assert!(s.contains("geometry 7: internal error: x"), "{s}");
}

/// Compile-time check of the thread-safety contract (test only).
#[test]
fn thread_safety_contract() {
    fn copy_send_sync<T: Copy + Send + Sync>() {}
    fn send<T: Send>() {}
    copy_send_sync::<burnscars::Grid>();
    send::<Burner>();
    send::<BurnResult>();
}
