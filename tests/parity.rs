//! controlledburn parity fixtures (`fixtures/`, commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558).
//! Port of `cpp/tests/test_parity_fixtures.cpp`: `edges_empty` is checked only when 1.

mod common;

use burnscars::{Burner, GeomId, Grid, Mode};
use common::{covered_area, hex, known_diffs, read_csv};

const GEOMS: &str = include_str!("../fixtures/geometries.csv");
const EXPECTED: &str = include_str!("../fixtures/expected.csv");

fn num(s: &str) -> Option<f64> {
    (s != "NA").then(|| s.parse().ok()).flatten()
}

#[test]
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    reason = "fixture values"
)]
fn parity_with_controlledburn_fixtures() {
    let skip = known_diffs("parity");
    let geoms = read_csv(GEOMS);
    let expected = read_csv(EXPECTED);
    let col = |hdr: &[String], name: &str| hdr.iter().position(|h| h == name).expect(name);
    let (gh, eh) = (&geoms[0], &expected[0]);
    let mut checked = 0;
    for e in &expected[1..] {
        let case = &e[col(eh, "case")];
        if skip.iter().any(|(c, _)| c == case) {
            continue;
        }
        let g = geoms.iter().skip(1).find(|r| &r[col(gh, "case")] == case).expect("case row");
        let f = |n: &str| g[col(gh, n)].parse::<f64>().expect(n);
        let grid = Grid::from_extent(
            f("xmin"),
            f("ymin"),
            f("xmax"),
            f("ymax"),
            f("ncol") as u32,
            f("nrow") as u32,
        )
        .expect("grid");
        let bytes = hex(&g[col(gh, "wkb_hex")]);
        let geom = wkb::reader::read_wkb(&bytes).expect("wkb");
        let r = Burner::new(Mode::Coverage).burn(&grid, [(GeomId(0), &geom)]).expect("burn");
        let tol = num(&e[col(eh, "tol_rel")]).unwrap_or(1e-9);
        if let Some(a) = num(&e[col(eh, "covered_area")]) {
            assert!((covered_area(&r, &grid) - a).abs() <= tol * a.max(1.0), "{case}: area");
        }
        if num(&e[col(eh, "edges_empty")]) == Some(1.0) {
            assert!(r.edges.is_empty(), "{case}: edges");
        }
        if let Some(l) = num(&e[col(eh, "line_length")]) {
            let got: f64 = r.lines.iter().map(|x| f64::from(x.length)).sum();
            assert!((got - l).abs() <= tol * l, "{case}: length");
        }
        if let Some(n) = num(&e[col(eh, "n_points")]) {
            assert_eq!(r.points.len() as f64, n, "{case}: points");
        }
        checked += 1;
    }
    assert!(checked >= 10 - skip.len());
}
