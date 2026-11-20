//! The stretch table on identity, south-up and periodic (0..360) grids, for random star polygons.
//!
//! - The layout holds: column lengths, pointers, sorted fibers, one value or one per cell.
//! - Compression: runs and blocks are maximal, so no two runs touch and no two blocks touch.
//! - `Touched` burns the cells of `Fraction`.
//! - A full periodic burn conserves the area of each star, and of a multipolygon with a hole.
//! - Oracles: a south-up burn is the north-up burn with the rows flipped, bit for bit;
//!   a periodic burn (`Fraction` and `Centre`) is a burn on a plain lattice over several
//!   periods, folded by the period.
//!   `properties.rs` and `oracle.rs` check `Fraction` on an identity grid against an exact clip
//!   and exactextract; `oracle.rs` checks `Centre` against the GDAL centre rule.

#![expect(clippy::expect_used, reason = "test code: a failed burn fails the test")]

mod common;

use std::collections::{BTreeMap, BTreeSet};

use burnscars::{Col, GeomId, Grid, Lattice, Row, Rule, Span, Stretches, XAxis, YAxis, YOrigin};
use common::{burn_spans, densify, offset};
use geo_traits::GeometryTrait;
use geo_types::{LineString, MultiPolygon, Polygon};
use proptest::prelude::*;

const NCOL: u32 = 36;
const NROW: u32 = 18;
const RULES: [Rule; 3] = [Rule::Fraction, Rule::Touched, Rule::Centre];

/// A star about `(cx, cy)`: angles strictly increase, every gap is below pi, radius >= 12.
/// A vertex within 1 of a lattice line (a multiple of 10) moves onto it; at radius >= 12 the
/// move keeps the angle order, so the ring stays simple.
fn star() -> impl Strategy<Value = Polygon<f64>> {
    (-200.0..560.0_f64, -25.0..25.0_f64, star_ring(12.0..60.0, true))
        .prop_map(|(cx, cy, ring)| Polygon::new(ring.at(cx, cy), vec![]))
}

/// A star's ring, placed at a centre by `at`.
#[derive(Clone, Debug)]
struct Ring {
    /// Per vertex: radius, angle jitter, snap.
    vs: Vec<(f64, f64, bool)>,
    snaps: bool,
}

impl Ring {
    fn at(&self, cx: f64, cy: f64) -> LineString<f64> {
        let n = self.vs.len();
        let snap = |v: f64| {
            let s = (v / 10.0).round() * 10.0;
            if (s - v).abs() < 1.0 { s } else { v }
        };
        #[expect(clippy::cast_precision_loss, reason = "small n")]
        let pts: Vec<(f64, f64)> = self
            .vs
            .iter()
            .enumerate()
            .map(|(i, &(r, jitter, s))| {
                let a = std::f64::consts::TAU * (i as f64 + jitter) / n as f64;
                let (x, y) = (cx + r * a.cos(), cy + r * a.sin());
                if s && self.snaps { (snap(x), snap(y)) } else { (x, y) }
            })
            .collect();
        LineString::from(pts)
    }
}

/// A star's ring with radii in `radii`; `snaps` lets vertices move onto lattice lines.
fn star_ring(radii: std::ops::Range<f64>, snaps: bool) -> impl Strategy<Value = Ring> {
    prop::collection::vec((radii, 0.0..0.5_f64, any::<bool>()), 5..14)
        .prop_map(move |vs| Ring { vs, snaps })
}

/// A star with a hole, and a second star 180 east: a valid multipolygon, narrower than a period.
/// The shell's edges stay over 12 cos(0.3 pi) > 7 from the centre, so a hole of radius < 6 is inside.
fn holed_pair() -> impl Strategy<Value = MultiPolygon<f64>> {
    (
        -200.0..380.0_f64,
        -25.0..25.0_f64,
        star_ring(12.0..60.0, true),
        star_ring(2.0..6.0, false),
        star_ring(12.0..60.0, true),
    )
        .prop_map(|(cx, cy, shell, hole, other)| {
            let mut hole = hole.at(cx, cy);
            hole.0.reverse();
            MultiPolygon::new(vec![
                Polygon::new(shell.at(cx, cy), vec![hole]),
                Polygon::new(other.at(cx + 180.0, cy), vec![]),
            ])
        })
}

fn span<T: Copy + Ord + Into<u32>>(a: T, b: T) -> Span<T> {
    Span::new(a, b).expect("span")
}

/// `NROW` x `ncol` cells of 10 x 10, from x = `x0` and between y = -90 and 90,
/// cut to `rows` x `cols`. Every star lies inside the rows.
fn grid(x: XAxis, origin: YOrigin, ncol: u32, rows: (u32, u32), cols: (u32, u32)) -> Grid {
    let y0 = if origin == YOrigin::NorthUp { 90.0 } else { -90.0 };
    let y = YAxis::new(y0, 10.0, origin).expect("y");
    Grid::new(Lattice::new(x, y), span(Row(0), Row(NROW)), span(Col(0), Col(ncol)))
        .and_then(|g| g.window(span(Row(rows.0), Row(rows.1)), span(Col(cols.0), Col(cols.1))))
        .expect("grid")
}

fn periodic(cols: (u32, u32)) -> Grid {
    let x = XAxis::new(0.0, 10.0).and_then(|x| x.with_period(360.0)).expect("x");
    grid(x, YOrigin::NorthUp, NCOL, (0, NROW), cols)
}

/// Five periods from x = -720, with no period.
fn plain() -> Grid {
    grid(XAxis::new(-720.0, 10.0).expect("x"), YOrigin::NorthUp, 5 * NCOL, (0, NROW), (0, 5 * NCOL))
}

/// A sorted window `a..=b` from two draws.
fn window((a, b): (u32, u32)) -> (u32, u32) {
    (a.min(b), a.max(b) + 1)
}

#[expect(clippy::cast_possible_truncation, reason = "test indices")]
fn burn<P>(rule: Rule, g: &Grid, polys: &[P]) -> Stretches
where
    for<'p> &'p P: GeometryTrait<T = f64>,
{
    burn_spans(rule, g, polys.iter().enumerate().map(|(i, p)| (GeomId(i as u32), p))).expect("burn")
}

/// The plain burn folded by the period, cut to `cols`: `(id, row, col) -> f64` sum.
fn folded<P>(rule: Rule, polys: &[P], cols: (u32, u32)) -> BTreeMap<(u32, u32, u32), f64>
where
    for<'p> &'p P: GeometryTrait<T = f64>,
{
    let mut want: BTreeMap<(u32, u32, u32), f64> = BTreeMap::new();
    for ((id, r, c), w) in densify(&burn(rule, &plain(), polys)) {
        *want.entry((id, r, c % NCOL)).or_default() += f64::from(w);
    }
    want.retain(|&(_, _, c), _| (cols.0..cols.1).contains(&c));
    want
}

/// The layout and the compression rules of `Stretches`.
fn check_layout(s: &Stretches) -> Result<(), TestCaseError> {
    let (nf, ns) = (s.id.len(), s.col_start.len());
    prop_assert_eq!(s.row.len(), nf);
    prop_assert_eq!(s.col_end.len(), ns);
    prop_assert_eq!(s.fiber_ptr.len(), nf + 1);
    prop_assert_eq!((s.fiber_ptr[0], s.fiber_ptr[nf]), (0, ns));
    prop_assert!(s.fiber_ptr.windows(2).all(|w| w[0] < w[1]), "a fiber with no stretch");
    let keys: Vec<(u32, u32)> = s.id.iter().copied().zip(s.row.iter().copied()).collect();
    prop_assert!(keys.windows(2).all(|w| w[0] < w[1]), "fibers not strictly sorted");
    prop_assert!((0..ns).all(|k| s.col_start[k] < s.col_end[k]), "an empty stretch");
    // A run has weight 1 (or no values); a block holds one weight below 1 per cell.
    let weights = |k: usize| -> Option<&[f32]> {
        s.values.as_ref().map(|v| &v.weight[offset(v.value_ptr[k])..offset(v.value_ptr[k + 1])])
    };
    let is_run = |k: usize| weights(k).is_none_or(|w| w == [1.0]);
    if let Some(v) = &s.values {
        prop_assert_eq!(v.value_ptr.len(), ns + 1);
        prop_assert_eq!(v.value_ptr[0], 0);
        prop_assert_eq!(offset(v.value_ptr[ns]), v.weight.len());
        prop_assert!(v.value_ptr.windows(2).all(|w| w[0] <= w[1]), "value_ptr not monotonic");
        for k in 0..ns {
            let n = weights(k).map_or(0, <[f32]>::len);
            let width = (s.col_end[k] - s.col_start[k]) as usize;
            prop_assert!(n == 1 || n == width, "stretch {k}: {n} values for {width} cells");
            prop_assert!(
                is_run(k) || weights(k).is_some_and(|w| w.iter().all(|&x| x > 0.0 && x < 1.0)),
                "stretch {k}: block weights {:?}",
                weights(k)
            );
        }
    }
    for f in 0..nf {
        for k in s.fiber_ptr[f]..s.fiber_ptr[f + 1] - 1 {
            prop_assert!(s.col_end[k] <= s.col_start[k + 1], "stretches overlap or unsorted");
            if s.col_end[k] == s.col_start[k + 1] {
                prop_assert!(is_run(k) != is_run(k + 1), "two runs or two blocks touch at {k}");
            }
        }
    }
    Ok(())
}

/// Shoelace area of the exterior ring.
fn area(p: &Polygon<f64>) -> f64 {
    ring_area(p.exterior())
}

fn ring_area(r: &LineString<f64>) -> f64 {
    let c: Vec<_> = r.coords().collect();
    c.windows(2).map(|w| w[0].x * w[1].y - w[1].x * w[0].y).sum::<f64>().abs() / 2.0
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn spans_hold_their_layout_on_every_grid(
        polys in prop::collection::vec(star(), 1..4),
        cols in (0..NCOL, 0..NCOL),
    ) {
        let cols = window(cols);
        let x = XAxis::new(0.0, 10.0).expect("x");
        let grids = [
            grid(x, YOrigin::NorthUp, NCOL, (0, NROW), cols),
            grid(x, YOrigin::SouthUp, NCOL, (0, NROW), cols),
            periodic(cols),
        ];
        for g in &grids {
            let fraction = burn(Rule::Fraction, g, &polys);
            for rule in RULES {
                let s = burn(rule, g, &polys);
                check_layout(&s)?;
                prop_assert_eq!(s.values.is_some(), rule == Rule::Fraction);
            }
            let touched: BTreeSet<_> = densify(&burn(Rule::Touched, g, &polys)).into_keys().collect();
            let want: BTreeSet<_> = densify(&fraction).into_keys().collect();
            prop_assert_eq!(touched, want);
        }
    }

    /// The full periodic grid holds every copy of a star, so its cells add up to its area,
    /// and every rule gives each star a fiber.
    #[test]
    fn a_full_periodic_burn_conserves_area(polys in prop::collection::vec(star(), 1..4)) {
        let g = periodic((0, NCOL));
        let cells = densify(&burn(Rule::Fraction, &g, &polys));
        for (i, p) in (0_u32..).zip(&polys) {
            let got = 100.0
                * cells.iter().filter(|(k, _)| k.0 == i).map(|(_, &w)| f64::from(w)).sum::<f64>();
            let want = area(p);
            prop_assert!((got - want).abs() <= 1e-5 * want, "geometry {i}: {got} vs {want}");
        }
        for rule in RULES {
            let ids: BTreeSet<u32> = burn(rule, &g, &polys).id.into_iter().collect();
            prop_assert_eq!(ids.len(), polys.len(), "{:?}", rule);
        }
    }

    /// A multipolygon with a hole: the full periodic burn conserves its area, and its spans hold their layout.
    #[test]
    fn a_holed_multipolygon_conserves_area(mp in holed_pair()) {
        let g = periodic((0, NCOL));
        let s = burn(Rule::Fraction, &g, std::slice::from_ref(&mp));
        check_layout(&s)?;
        let got = 100.0 * densify(&s).values().map(|&w| f64::from(w)).sum::<f64>();
        let want: f64 = mp
            .0
            .iter()
            .map(|p| area(p) - p.interiors().iter().map(ring_area).sum::<f64>())
            .sum();
        prop_assert!((got - want).abs() <= 1e-5 * want, "{got} vs {want}");
        for rule in [Rule::Touched, Rule::Centre] {
            check_layout(&burn(rule, &g, std::slice::from_ref(&mp)))?;
        }
    }

    /// Rows of 10 at integer edges: both row orders burn the same cells, bit for bit,
    /// on any window.
    #[test]
    fn south_up_is_north_up_with_the_rows_flipped(
        polys in prop::collection::vec(star(), 1..4),
        rows in (0..NROW, 0..NROW),
        cols in (0..NCOL, 0..NCOL),
    ) {
        let x = XAxis::new(0.0, 10.0).expect("x");
        let ((r0, r1), cols) = (window(rows), window(cols));
        let gn = grid(x, YOrigin::NorthUp, NCOL, (r0, r1), cols);
        let gs = grid(x, YOrigin::SouthUp, NCOL, (NROW - r1, NROW - r0), cols);
        for rule in RULES {
            let north = densify(&burn(rule, &gn, &polys));
            let south: BTreeMap<_, _> = densify(&burn(rule, &gs, &polys))
                .into_iter()
                .map(|((id, r, c), w)| ((id, NROW - 1 - r, c), w))
                .collect();
            prop_assert_eq!(&south, &north, "{:?}", rule);
        }
    }

    /// `Centre` on a periodic grid burns the cells of the folded plain burn.
    #[test]
    fn a_periodic_centre_burn_is_the_folded_plain_burn(
        polys in prop::collection::vec(star(), 1..4),
        cols in (0..NCOL, 0..NCOL),
    ) {
        let cols = window(cols);
        let got: BTreeSet<_> = densify(&burn(Rule::Centre, &periodic(cols), &polys)).into_keys().collect();
        let want: BTreeSet<_> = folded(Rule::Centre, &polys, cols).into_keys().collect();
        prop_assert_eq!(got, want);
    }

    /// The copies of a geometry on a periodic grid sum per cell, as the plain burn folded by the period.
    #[test]
    fn a_periodic_burn_is_the_folded_plain_burn(
        polys in prop::collection::vec(star(), 1..4),
        cols in (0..NCOL, 0..NCOL),
    ) {
        let cols = window(cols);
        let want = folded(Rule::Fraction, &polys, cols);
        let got = densify(&burn(Rule::Fraction, &periodic(cols), &polys));
        prop_assert_eq!(got.keys().collect::<Vec<_>>(), want.keys().collect::<Vec<_>>());
        for (k, w) in &got {
            let v = want[k].min(1.0);
            prop_assert!((f64::from(*w) - v).abs() <= 2e-6, "cell {k:?}: {w} vs {v}");
        }
    }
}
