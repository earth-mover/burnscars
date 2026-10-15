//! Test-only helpers shared by the integration tests.
#![allow(dead_code)]

use std::collections::BTreeMap;

use burnscars::{BurnResult, Grid};
use geo_types::{LineString, Polygon, coord};

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

/// Sum of runs and edge fractions times the cell area.
pub(crate) fn covered_area(r: &BurnResult, g: &Grid) -> f64 {
    let l = g.lattice();
    let runs: f64 = r.runs.iter().map(|x| f64::from(x.cols.width())).sum();
    let edges: f64 = r.edges.iter().map(|e| f64::from(e.fraction)).sum();
    l.dx() * l.dy() * (runs + edges)
}

/// Expanded form: (table, id, row, col) -> value bits. Runs give 1.0 per cell.
pub(crate) fn expanded(r: &BurnResult) -> BTreeMap<(u8, u32, u32, u32), Vec<u32>> {
    let mut m: BTreeMap<(u8, u32, u32, u32), Vec<u32>> = BTreeMap::new();
    for x in r.runs.iter() {
        for c in x.cols.range() {
            m.entry((0, x.id.0, x.row.0, c)).or_default().push(1.0_f32.to_bits());
        }
    }
    for e in r.edges.iter() {
        m.entry((0, e.id.0, e.row.0, e.col.0)).or_default().push(e.fraction.to_bits());
    }
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
