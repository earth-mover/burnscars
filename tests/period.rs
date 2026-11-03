//! South-up rows and a longitude period: lon/lat bricks in the caller's own indices.

#![expect(clippy::expect_used, reason = "test code: a failed burn fails the test")]

mod common;

use std::collections::BTreeMap;

use burnscars::{
    BurnError, BurnResult, Burner, Col, GeomId, Grid, Lattice, Mode, Pos, Row, Span, XAxis, YAxis,
    YOrigin,
};
use common::{expanded, ls, poly};
use geo_types::{Geometry, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Rect};

type Expanded = BTreeMap<(u8, u32, u32, u32), Vec<u32>>;

fn span<T: Copy + Ord + Into<u32>>(a: T, b: T) -> Span<T> {
    Span::new(a, b).expect("span")
}

/// A lattice on a full grid of `nrow` x `ncol`, cut to the window `rows` x `cols`.
fn window(l: Lattice, nrow: u32, ncol: u32, rows: (u32, u32), cols: (u32, u32)) -> Grid {
    Grid::new(l, span(Row(0), Row(nrow)), span(Col(0), Col(ncol)))
        .and_then(|g| g.window(span(Row(rows.0), Row(rows.1)), span(Col(cols.0), Col(cols.1))))
        .expect("grid")
}

#[expect(clippy::cast_possible_truncation, reason = "test indices")]
fn burn(mode: Mode, g: &Grid, geoms: &[Geometry<f64>]) -> BurnResult {
    let items = geoms.iter().enumerate().map(|(i, x)| (GeomId(i as u32), x));
    Burner::new(mode).burn(g, items).expect("burn")
}

fn runs(r: &BurnResult) -> Vec<(u32, std::ops::Range<u32>)> {
    r.runs.iter().map(|x| (x.row.0, x.cols.range())).collect()
}

/// One cell per key, with the value as an f32.
fn cells(m: &Expanded) -> Vec<((u8, u32, u32), f32)> {
    m.iter().map(|(&(t, _, r, c), v)| ((t, r, c), f32::from_bits(v[0]))).collect()
}

/// Longitude step 10 from the west edge `x0`, period 360; latitude centres 10, 20, 30
/// (row 0 south): y0 = 5 is the south edge.
fn lonlat(x0: f64) -> Lattice {
    let x = XAxis::new(x0, 10.0).and_then(|x| x.with_period(360.0)).expect("x");
    Lattice::new(x, YAxis::new(5.0, 10.0, YOrigin::SouthUp).expect("y"))
}

fn zone() -> Geometry<f64> {
    Rect::new((-175.0, 15.0), (-165.0, 35.0)).into()
}

/// The brick: lon centres 170..200 on a 0..360 grid; the zone is in [-180, 180).
#[test]
fn the_brick_burns_in_its_own_indices() {
    let l = lonlat(165.0);
    let g = window(l, 3, 4, (0, 3), (0, 4));
    assert_eq!(runs(&burn(Mode::Approx, &g, &[zone()])), [(2, 2..3), (1, 2..3)]);
    let cov = cells(&expanded(&burn(Mode::Coverage, &g, &[zone()])));
    assert_eq!(cov, [((0, 1, 2), 1.0), ((0, 2, 2), 1.0)]);
}

/// The same brick as a window of a global 0..360 lattice (centres 0..350), across 180 or not.
#[test]
fn a_window_of_the_global_lattice_gives_global_indices() {
    let l = lonlat(-5.0);
    let across = window(l, 3, 36, (0, 3), (17, 21));
    assert_eq!(runs(&burn(Mode::Approx, &across, &[zone()])), [(2, 19..20), (1, 19..20)]);
    let east = window(l, 3, 36, (0, 3), (19, 23));
    assert_eq!(runs(&burn(Mode::Approx, &east, &[zone()])), [(2, 19..20), (1, 19..20)]);
    let west = window(l, 3, 36, (0, 3), (10, 19));
    assert!(burn(Mode::Approx, &west, &[zone()]).is_empty());
}

/// A zone split at 180 (two parts) fills the cell that the seam cuts: 0.5 from each part.
#[test]
fn a_zone_split_at_the_seam_fills_the_cell_across_it() {
    let l = lonlat(-5.0);
    let g = window(l, 3, 36, (0, 3), (17, 21));
    let split = Geometry::MultiPolygon(MultiPolygon::new(vec![
        poly(&[(175.0, 15.0), (180.0, 15.0), (180.0, 35.0), (175.0, 35.0), (175.0, 15.0)], &[]),
        poly(
            &[(-180.0, 15.0), (-175.0, 15.0), (-175.0, 35.0), (-180.0, 35.0), (-180.0, 15.0)],
            &[],
        ),
    ]));
    let cov = cells(&expanded(&burn(Mode::Coverage, &g, std::slice::from_ref(&split))));
    assert_eq!(cov, [((0, 1, 18), 1.0), ((0, 2, 18), 1.0)]);
    assert_eq!(runs(&burn(Mode::Approx, &g, &[split])), [(2, 18..19), (1, 18..19)]);

    let line = Geometry::MultiLineString(MultiLineString::new(vec![
        ls(&[(176.0, 20.0), (180.0, 20.0)]),
        ls(&[(-180.0, 20.0), (-176.0, 20.0)]),
    ]));
    assert_eq!(cells(&expanded(&burn(Mode::Coverage, &g, &[line]))), [((1, 1, 18), 8.0)]);
}

/// A zone in 0..360 on a [-180, 180) lattice: its copy one period west fills column 0.
#[test]
fn a_zone_past_180_burns_on_a_lattice_from_minus_180() {
    let l = lonlat(-185.0);
    let g = window(l, 3, 36, (0, 3), (0, 36));
    let z: Geometry<f64> = Rect::new((175.0, 15.0), (185.0, 35.0)).into();
    let cov = cells(&expanded(&burn(Mode::Coverage, &g, &[z])));
    assert_eq!(cov, [((0, 1, 0), 1.0), ((0, 2, 0), 1.0)]);
    let line = Geometry::LineString(ls(&[(176.0, 20.0), (184.0, 20.0)]));
    assert_eq!(cells(&expanded(&burn(Mode::Coverage, &g, &[line]))), [((1, 1, 0), 8.0)]);
}

/// On a lattice of exactly one period, the far x edge is the seam: a point on it goes
/// into column 0 only, as a point on the near edge does.
#[test]
fn a_point_on_the_seam_is_in_one_cell() {
    let l = lonlat(0.0);
    let g = window(l, 3, 36, (0, 3), (0, 36));
    let pts = |x: f64| Geometry::MultiPoint(MultiPoint::new(vec![Point::new(x, 20.0)]));
    for (x, col) in [(0.0, 0), (360.0, 0), (-360.0, 0), (-180.0, 18), (185.0, 18), (-5.0, 35)] {
        let got = cells(&expanded(&burn(Mode::Coverage, &g, &[pts(x)])));
        assert_eq!(got, [((2, 1, col), 0.0)], "x = {x}");
    }
}

/// A window east of 180 on a 0..360 lattice equals the window 18 columns to the west on a
/// [-180, 180) lattice with no period, bit for bit.
#[test]
fn a_window_east_of_180_matches_a_shifted_lattice() {
    let y = YAxis::new(40.0, 10.0, YOrigin::NorthUp).expect("y");
    let x = XAxis::new(0.0, 10.0).and_then(|x| x.with_period(360.0)).expect("x");
    let wrapped = Lattice::new(x, y);
    let plain = Lattice::new(XAxis::new(-180.0, 10.0).expect("x"), y);
    let a = window(wrapped, 4, 36, (0, 4), (19, 25));
    let b = window(plain, 4, 36, (0, 4), (1, 7));
    let geoms: Vec<Geometry<f64>> = vec![
        poly(
            &[(-163.3, 3.1), (-131.7, 9.9), (-140.2, 37.4), (-155.0, 24.25), (-163.3, 3.1)],
            &[&[(-150.0, 15.0), (-142.5, 17.5), (-147.0, 22.2), (-150.0, 15.0)]],
        )
        .into(),
        Geometry::LineString(LineString::from(vec![(-171.1, 2.2), (-150.0, 30.0), (-128.3, 11.7)])),
        Geometry::MultiPoint(MultiPoint::new(vec![
            Point::new(-160.0, 20.0),
            Point::new(-133.3, 1.0),
        ])),
    ];
    for mode in [Mode::Coverage, Mode::Approx] {
        let shifted: Expanded = expanded(&burn(mode, &a, &geoms))
            .into_iter()
            .map(|((t, id, r, c), v)| ((t, id, r, c - 18), v))
            .collect();
        let want = expanded(&burn(mode, &b, &geoms));
        assert!(!want.is_empty());
        assert_eq!(shifted, want, "{mode:?}");
    }
}

#[test]
fn a_geometry_wider_than_the_period_is_refused() {
    let l = lonlat(0.0);
    let g = window(l, 3, 36, (0, 3), (0, 36));
    let wide: Geometry<f64> = Rect::new((-180.0, 15.0), (181.0, 35.0)).into();
    let e = Burner::new(Mode::Coverage).burn(&g, [(GeomId(4), &wide)]).expect_err("too wide");
    assert!(matches!(e, BurnError::WiderThanPeriod { id: GeomId(4), .. }), "{e}");
    let world: Geometry<f64> = Rect::new((-180.0, 15.0), (180.0, 35.0)).into();
    let out = burn(Mode::Approx, &g, &[world]);
    assert_eq!(runs(&out), [(2, 0..36), (1, 0..36)]);
}

/// `col_of` and `row_of` answer in the caller's indices.
#[test]
fn col_of_and_row_of_use_the_callers_indices() {
    let l = lonlat(0.0);
    let g = window(l, 3, 36, (0, 3), (0, 36));
    assert_eq!(g.col_of(-175.0), Pos::In(Col(18)));
    assert_eq!(g.col_of(5.0), Pos::In(Col(0)));
    assert_eq!(g.row_of(10.0), Pos::In(Row(0)));
    assert_eq!(g.row_of(30.0), Pos::In(Row(2)));
    assert_eq!(g.row_of(40.0), Pos::After);
    assert_eq!(g.row_of(0.0), Pos::Before);
}

/// A point on the seam of a full-period grid goes into one tile only: column 0, as in the
/// full burn, and not also the last column of the tile that holds the far edge.
#[test]
fn a_point_on_the_seam_is_in_one_tile() {
    let pts = Geometry::MultiPoint(MultiPoint::new(vec![Point::new(360.0, 20.0)]));
    let full = cells(&expanded(&burn(
        Mode::Coverage,
        &window(lonlat(0.0), 3, 36, (0, 3), (0, 36)),
        std::slice::from_ref(&pts),
    )));
    assert_eq!(full, [((2, 1, 0), 0.0)]);
    let west = burn(
        Mode::Coverage,
        &window(lonlat(0.0), 3, 36, (0, 3), (0, 1)),
        std::slice::from_ref(&pts),
    );
    let east = burn(Mode::Coverage, &window(lonlat(0.0), 3, 36, (0, 3), (35, 36)), &[pts]);
    assert_eq!(cells(&expanded(&west)), full);
    assert!(east.is_empty(), "{east:?}");
}
