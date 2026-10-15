//! Oracle comparison (spec: Oracle): exactextract coverage and GDAL centre rule.

mod common;

use std::collections::{BTreeMap, BTreeSet};

use burnscars::{BurnResult, Burner, GeomId, Grid, Mode};
use common::{hex, known_diff_cells, read_csv};

const CASES: &str = include_str!("../fixtures/oracle/cases.csv");
const COVERAGE: &str = include_str!("../fixtures/oracle/coverage.csv");
const APPROX: &str = include_str!("../fixtures/oracle/approx.csv");
const TOL: f64 = 1e-6;

#[expect(clippy::expect_used, reason = "test helper: a bad fixture row is a test bug")]
fn grid_of(r: &[String]) -> Grid {
    let f = |i: usize| r[i].parse::<f64>().expect("grid value");
    Grid::from_extent(
        f(2),
        f(3),
        f(4),
        f(5),
        r[6].parse().expect("ncol"),
        r[7].parse().expect("nrow"),
    )
    .expect("grid")
}

#[expect(clippy::expect_used, reason = "test helper: a bad fixture row is a test bug")]
fn burn(case: &[String], mode: Mode) -> BurnResult {
    let bytes = hex(&case[1]);
    let g = wkb::reader::read_wkb(&bytes).expect("wkb");
    Burner::new(mode).burn(&grid_of(case), [(GeomId(0), &g)]).expect("burn")
}

/// Per-cell fraction, runs as 1, cells below TOL removed.
fn cells(r: &BurnResult) -> BTreeMap<(u32, u32), f64> {
    let mut m = BTreeMap::new();
    for x in r.runs.iter() {
        for c in x.cols.range() {
            m.insert((x.row.0, c), 1.0);
        }
    }
    for e in r.edges.iter() {
        m.insert((e.row.0, e.col.0), f64::from(e.fraction));
    }
    m.retain(|_, f| *f >= TOL);
    m
}

#[expect(clippy::expect_used, reason = "test helper: a bad fixture row is a test bug")]
fn oracle_cells(rows: &[Vec<String>], case: &str) -> BTreeMap<(u32, u32), f64> {
    rows.iter()
        .filter(|r| r[0] == case)
        .map(|r| {
            (
                (r[1].parse().expect("row"), r[2].parse().expect("col")),
                r.get(3).map_or(1.0, |f| f.parse().expect("f")),
            )
        })
        .filter(|(_, f)| *f >= TOL)
        .collect()
}

/// Cells where burnscars and the oracle disagree: a key on one side only, or |Δf| > TOL.
fn diff_cells(
    got: &BTreeMap<(u32, u32), f64>,
    want: &BTreeMap<(u32, u32), f64>,
) -> BTreeSet<(u32, u32)> {
    let keys: BTreeSet<(u32, u32)> = got.keys().chain(want.keys()).copied().collect();
    keys.into_iter()
        .filter(|k| match (got.get(k), want.get(k)) {
            (Some(a), Some(b)) => (a - b).abs() > TOL,
            _ => true,
        })
        .collect()
}

/// Case names whose difference set is not the listed one (empty if not listed).
fn mismatches(table: &str, mode: Mode, cases: &str, oracle: &str) -> Vec<String> {
    let known = known_diff_cells(table);
    let rows = read_csv(oracle);
    let mut failures = Vec::new();
    for case in &read_csv(cases)[1..] {
        let got = cells(&burn(case, mode));
        let want = oracle_cells(&rows[1..], &case[0]);
        let diff = diff_cells(&got, &want);
        if known.get(&case[0]).map_or(!diff.is_empty(), |listed| *listed != diff) {
            // Print in the known_diffs.toml form, for review.
            let list: Vec<String> = diff.iter().map(|(r, c)| format!("[{r}, {c}]")).collect();
            eprintln!("{table} {}: cells = [{}]", case[0], list.join(", "));
            failures.push(case[0].clone());
        }
    }
    failures
}

fn compare(table: &str, mode: Mode, oracle: &str) {
    let failures = mismatches(table, mode, CASES, oracle);
    assert!(failures.is_empty(), "{table} differs from the oracle for: {failures:?}");
}

#[test]
fn coverage_matches_exactextract() {
    compare("coverage", Mode::Coverage, COVERAGE);
}

#[test]
fn approx_matches_gdal_centre_rule() {
    compare("approx", Mode::Approx, APPROX);
}

/// Differential fuzz against the live oracle: `just fuzz-oracle` (needs uv).
#[test]
#[ignore = "needs uv and the Python oracle"]
fn fuzz_against_live_oracle() {
    let dir = std::env::temp_dir().join(format!("burnscars-fuzz-oracle-{}", std::process::id()));
    let status = std::process::Command::new("uv")
        .args(["run", "--locked", "--project", "oracle", "python", "oracle/generate.py", "--out"])
        .arg(&dir)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .status()
        .expect("run uv");
    assert!(status.success());
    let cases = std::fs::read_to_string(dir.join("cases.csv")).expect("cases");
    let cov = std::fs::read_to_string(dir.join("coverage.csv")).expect("coverage");
    let failures = mismatches("coverage", Mode::Coverage, &cases, &cov);
    assert!(failures.is_empty(), "coverage differs from the live oracle for: {failures:?}");
    let approx = std::fs::read_to_string(dir.join("approx.csv")).expect("approx");
    let failures = mismatches("approx", Mode::Approx, &cases, &approx);
    assert!(failures.is_empty(), "approx differs from the live oracle for: {failures:?}");
}
