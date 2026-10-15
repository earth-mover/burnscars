//! Steady state makes no allocation (spec: Storage, O6). One test per binary, so no
//! other test changes the counter.
#![expect(unsafe_code, reason = "test-only counting allocator")]

mod common;

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

use burnscars::{BurnResult, Burner, GeomId, Mode};
use common::{grid, ls, poly};
use geo_types::{Geometry, MultiPoint, Point};

struct Counting;
static ALLOCS: AtomicUsize = AtomicUsize::new(0);

// SAFETY: forwards every call to `System` unchanged; only counts.
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        ALLOCS.fetch_add(1, Ordering::Relaxed);
        // SAFETY: same contract as the caller's.
        unsafe { System.alloc(l) }
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        // SAFETY: same contract as the caller's.
        unsafe { System.dealloc(p, l) }
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, n: usize) -> *mut u8 {
        ALLOCS.fetch_add(1, Ordering::Relaxed);
        // SAFETY: same contract as the caller's.
        unsafe { System.realloc(p, l, n) }
    }
}

#[global_allocator]
static A: Counting = Counting;

#[test]
fn second_burn_makes_no_allocation() {
    let g = grid(0.0, 0.0, 10.0, 10.0, 50, 50);
    let geoms: Vec<Geometry<f64>> = vec![
        Geometry::Polygon(poly(
            &[(1.3, 1.3), (8.7, 2.1), (7.9, 8.4), (2.2, 7.7), (1.3, 1.3)],
            &[&[(3.0, 3.0), (5.0, 3.1), (4.2, 5.0), (3.0, 3.0)]],
        )),
        Geometry::LineString(ls(&[(0.5, 0.5), (9.5, 7.5), (2.0, 9.0)])),
        Geometry::MultiPoint(MultiPoint::new(vec![Point::new(1.0, 1.0), Point::new(9.0, 9.0)])),
    ];
    for mode in [Mode::Coverage, Mode::Approx] {
        let mut b = Burner::new(mode);
        let mut out = BurnResult::default();
        let items =
            || geoms.iter().enumerate().map(|(i, x)| (GeomId(u32::try_from(i).unwrap_or(0)), x));
        b.burn_into(&g, items(), &mut out).expect("warm-up");
        out.clear();
        let before = ALLOCS.load(Ordering::Relaxed);
        b.burn_into(&g, items(), &mut out).expect("burn");
        let after = ALLOCS.load(Ordering::Relaxed);
        assert_eq!(after - before, 0, "{mode:?}");
    }
}
