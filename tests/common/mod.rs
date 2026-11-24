//! Test-only helpers shared by the integration tests.
#![allow(dead_code)]

use std::collections::{BTreeMap, BTreeSet};

use burnscars::{BurnError, BurnResult, Burner, GeomId, Grid, Origin, Rule, Stretches};
use geo_traits::GeometryTrait;
use geo_types::{Geometry, LineString, Polygon, coord};

#[expect(clippy::panic, reason = "test helper: a bad test grid is a test bug")]
pub(crate) fn grid(xmin: f64, ymin: f64, xmax: f64, ymax: f64, ncol: u32, nrow: u32) -> Grid {
    match Grid::from_extent(xmin, ymin, xmax, ymax, ncol, nrow) {
        Ok(g) => g,
        Err(e) => panic!("bad test grid: {e}"),
    }
}

pub(crate) fn ls(pts: &[(f64, f64)]) -> LineString<f64> {
    pts.iter().map(|&(x, y)| coord! { x: x, y: y }).collect()
}

pub(crate) fn poly(ext: &[(f64, f64)], holes: &[&[(f64, f64)]]) -> Polygon<f64> {
    Polygon::new(ls(ext), holes.iter().map(|h| ls(h)).collect())
}

/// Burns polygons into a new `Stretches`.
pub(crate) fn burn_spans<G, I>(rule: Rule, g: &Grid, items: I) -> Result<Stretches, BurnError>
where
    G: GeometryTrait<T = f64>,
    I: IntoIterator<Item = (GeomId, G)>,
{
    let mut out = Stretches::new(rule);
    Burner::new(rule).burn_polygons(g, items, Origin::Grid, &mut out)?;
    Ok(out)
}

/// `(id, row, col) -> weight`, 1.0 with no values (test only).
/// Panics on a cell that two stretches cover.
#[expect(clippy::panic, reason = "test helper: a cell twice is a test failure")]
pub(crate) fn densify(s: &Stretches) -> BTreeMap<(u32, u32, u32), f32> {
    let mut m = BTreeMap::new();
    for st in s.iter() {
        for (i, c) in st.cols.range().enumerate() {
            let w = st.weight.map_or(1.0, |ws| if ws.len() == 1 { ws[0] } else { ws[i] });
            if m.insert((st.id.0, st.row.0, c), w).is_some() {
                panic!("cell {:?} twice", (st.id.0, st.row.0, c));
            }
        }
    }
    m
}

/// A value offset as an index.
#[expect(clippy::expect_used, reason = "test helper: a negative offset is a test failure")]
pub(crate) fn offset(o: i32) -> usize {
    usize::try_from(o).expect("value offset >= 0")
}

/// Sum of the cell weights times the cell area.
pub(crate) fn covered_area(s: &Stretches, g: &Grid) -> f64 {
    let l = g.lattice();
    let cells: f64 = densify(s).values().map(|&w| f64::from(w)).sum();
    l.x().dx() * l.y().dy() * cells
}

/// The cells with a weight below 1: the boundary cells under `Rule::Fraction`.
pub(crate) fn edge_cells(s: &Stretches) -> Vec<((u32, u32, u32), f32)> {
    densify(s).into_iter().filter(|&(_, w)| w < 1.0).collect()
}

/// Expanded form: (table, id, row, col) -> value bits.
pub(crate) type Expanded = BTreeMap<(u8, u32, u32, u32), Vec<u32>>;

/// Expanded form of the polygon spans, as table 0.
pub(crate) fn expanded_spans(s: &Stretches) -> Expanded {
    densify(s).into_iter().map(|((id, r, c), w)| ((0, id, r, c), vec![w.to_bits()])).collect()
}

/// Expanded form of the lines (table 1) and the points (table 2).
pub(crate) fn expanded(r: &BurnResult) -> Expanded {
    let mut m = Expanded::new();
    for l in r.lines.iter() {
        m.entry((1, l.id.0, l.row.0, l.col.0)).or_default().push(l.length.to_bits());
    }
    for p in r.points.iter() {
        m.entry((2, p.id.0, p.row.0, p.col.0)).or_default().push(0);
    }
    for v in m.values_mut() {
        v.sort_unstable();
    }
    m
}

/// True for a polygon, a multipolygon, a rect or a triangle.
pub(crate) fn is_polygonal(g: &Geometry<f64>) -> bool {
    matches!(
        g,
        Geometry::Polygon(_)
            | Geometry::MultiPolygon(_)
            | Geometry::Rect(_)
            | Geometry::Triangle(_)
    )
}

/// Burns a mixed list, id = index: polygons into `Stretches`, lines and points into `BurnResult`.
/// The expanded forms of both, in one map.
#[expect(clippy::cast_possible_truncation, reason = "test indices")]
pub(crate) fn burn_mixed(
    rule: Rule,
    g: &Grid,
    geoms: &[Geometry<f64>],
) -> Result<Expanded, BurnError> {
    let items = || geoms.iter().enumerate().map(|(i, x)| (GeomId(i as u32), x));
    let polys = burn_spans(rule, g, items().filter(|(_, x)| is_polygonal(x)))?;
    let rest = Burner::new(rule).burn(g, items().filter(|(_, x)| !is_polygonal(x)))?;
    let mut m = expanded(&rest);
    m.extend(expanded_spans(&polys));
    Ok(m)
}

/// Minimal CSV reader: comma separator, double-quoted fields, no escaped quotes.
pub(crate) fn read_csv(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            let mut fields = Vec::new();
            let mut cur = String::new();
            let mut quoted = false;
            for ch in line.chars() {
                match (ch, quoted) {
                    ('"', _) => quoted = !quoted,
                    (',', false) => fields.push(std::mem::take(&mut cur)),
                    _ => cur.push(ch),
                }
            }
            fields.push(cur);
            fields
        })
        .collect()
}

pub(crate) fn hex(s: &str) -> Vec<u8> {
    (0..s.len() / 2).filter_map(|i| u8::from_str_radix(s.get(2 * i..2 * i + 2)?, 16).ok()).collect()
}

const KNOWN_DIFFS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/oracle/known_diffs.toml");

/// All `[[diff]]` entries of `oracle/known_diffs.toml`. A missing or bad file fails the test.
#[expect(clippy::panic, reason = "test helper: a bad known_diffs.toml is a test bug")]
fn diff_entries() -> Vec<toml::Value> {
    let path = KNOWN_DIFFS;
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    let doc: toml::Table = text.parse().unwrap_or_else(|e| panic!("cannot parse {path}: {e}"));
    match doc.get("diff") {
        None => Vec::new(),
        Some(v) => v
            .as_array()
            .unwrap_or_else(|| panic!("{path}: `diff` is not an array of tables"))
            .clone(),
    }
}

/// A non-blank string field of a `[[diff]]` entry.
#[expect(clippy::panic, reason = "test helper: a bad known_diffs.toml is a test bug")]
fn diff_field(d: &toml::Value, k: &str) -> String {
    match d.get(k).and_then(toml::Value::as_str) {
        Some(s) if !s.trim().is_empty() => s.to_owned(),
        _ => panic!("{KNOWN_DIFFS}: a [[diff]] entry has no `{k}`: {d}"),
    }
}

/// Whole-case `[[diff]]` entries for `table`, as (case, reason). Each entry needs a
/// `table`, a `case` and a `reason`.
pub(crate) fn known_diffs(table: &str) -> Vec<(String, String)> {
    diff_entries()
        .iter()
        .filter_map(|d| {
            let (t, case, reason) =
                (diff_field(d, "table"), diff_field(d, "case"), diff_field(d, "reason"));
            (t == table).then_some((case, reason))
        })
        .collect()
}

/// Per-cell `[[diff]]` entries for `table`: case -> the exact set of differing (row, col).
/// An entry for `table` with no `cells` list, or a case listed two times, fails the test.
#[expect(clippy::panic, reason = "test helper: a bad known_diffs.toml is a test bug")]
pub(crate) fn known_diff_cells(table: &str) -> BTreeMap<String, BTreeSet<(u32, u32)>> {
    // Checks `table`, `case` and `reason` of all entries.
    let _ = known_diffs(table);
    let mut out = BTreeMap::new();
    for d in diff_entries().iter().filter(|d| diff_field(d, "table") == table) {
        let case = diff_field(d, "case");
        let cells = d
            .get("cells")
            .and_then(toml::Value::as_array)
            .unwrap_or_else(|| panic!("{KNOWN_DIFFS}: {table}/{case} has no `cells` list"));
        let set: BTreeSet<(u32, u32)> = cells
            .iter()
            .map(|rc| {
                let n = |i: usize| {
                    rc.get(i)
                        .and_then(toml::Value::as_integer)
                        .and_then(|v| u32::try_from(v).ok())
                        .unwrap_or_else(|| panic!("{KNOWN_DIFFS}: {table}/{case}: bad cell {rc}"))
                };
                (n(0), n(1))
            })
            .collect();
        assert!(!set.is_empty(), "{KNOWN_DIFFS}: {table}/{case}: `cells` is empty");
        assert!(
            out.insert(case.clone(), set).is_none(),
            "{KNOWN_DIFFS}: {table}/{case} is listed two times"
        );
    }
    out
}
