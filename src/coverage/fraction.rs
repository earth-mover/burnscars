//! Coverage fraction of one cell.
//!
//! `analytical_fraction` and `closed_ring_fraction`: controlledburn
//! `analytical_covered_fraction` and `closed_ring_covered_fraction` (Apache-2.0, Copyright
//! 2025 Michael Sumner), commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558,
//! `cpp/src/ee/analytical_coverage.h`. `left_hand_area`: exactextract (Apache-2.0,
//! Copyright 2018 ISciences, LLC) `src/traversal_areas.cpp`, the GEOS-free version that
//! controlledburn vendors as `cpp/src/ee/traversal_areas.cpp`.

use std::ops::Range;

use crate::geom::{
    box_::CellBox,
    coord::Coord,
    measures::area_signed_ee,
    perimeter::{perimeter_distance, perimeter_distance_ccw},
};

/// Classification of an accumulated cell fraction, with the controlledburn tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum CellClass {
    Empty,
    Edge(f32),
    Full,
}

#[expect(clippy::cast_possible_truncation, reason = "f32 output, same as controlledburn")]
pub(crate) fn to_f32(v: f64) -> f32 {
    v as f32
}

#[derive(Clone, Debug)]
pub(crate) struct Chain {
    start: f64,
    stop: f64,
    coords: ChainCoords,
    visited: bool,
}

#[derive(Clone, Debug)]
enum ChainCoords {
    Range(Range<usize>),
    Corner(Coord),
}

impl CellClass {
    pub(crate) fn of(w: f32) -> Self {
        let tol = 1e-6_f32;
        if w > tol && w < 1.0 - tol {
            Self::Edge(w)
        } else if w >= 1.0 - tol {
            Self::Full
        } else {
            Self::Empty
        }
    }
}

/// Port of controlledburn `analytical_covered_fraction`
/// ([cpp/src/ee/analytical_coverage.h](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/ee/analytical_coverage.h)).
/// O3: the polygon is an iterator chain, not a Vec; the bits equal the C++ result.
pub(crate) fn analytical_fraction(b: &CellBox, coords: &[Coord]) -> Option<f64> {
    let cell_area = b.area();
    if cell_area <= 0.0 || coords.len() < 2 {
        return Some(0.0);
    }
    let (first, last) = (*coords.first()?, *coords.last()?);
    let perim = b.perimeter();
    let exit_pd = perimeter_distance(b, last)?;
    let entry_pd = perimeter_distance(b, first)?;
    let arc = if exit_pd > entry_pd + 1e-12 {
        exit_pd - entry_pd
    } else if entry_pd > exit_pd + 1e-12 {
        perim - entry_pd + exit_pd
    } else {
        // Entry and exit coincide: the area of the closed path.
        let close = (first != last).then_some(first);
        return Some(area_signed_ee(coords.iter().copied().chain(close)).abs() / cell_area);
    };
    let (h, w) = (b.height(), b.width());
    let corners = [
        (Coord::new(b.xmin, b.ymin), 0.0),
        (Coord::new(b.xmin, b.ymax), h),
        (Coord::new(b.xmax, b.ymax), h + w),
        (Coord::new(b.xmax, b.ymin), 2.0 * h + w),
    ];
    let mut in_arc = [(Coord::new(0.0, 0.0), 0.0_f64); 4];
    let mut n = 0;
    for (cc, pd) in corners {
        let mut d = exit_pd - pd;
        if d < 0.0 {
            d += perim;
        }
        if d > 1e-12 && d < arc - 1e-12 {
            in_arc[n] = (cc, d);
            n += 1;
        }
    }
    // Corner distances are distinct, so any sort gives the C++ insertion-sort order.
    in_arc[..n].sort_unstable_by(|a, b| a.1.total_cmp(&b.1));
    let poly =
        coords.iter().copied().chain(in_arc[..n].iter().map(|t| t.0)).chain(std::iter::once(first));
    Some(area_signed_ee(poly).abs() / cell_area)
}

/// Port of controlledburn `closed_ring_covered_fraction`
/// ([cpp/src/ee/analytical_coverage.h#L164](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/ee/analytical_coverage.h#L164)).
pub(crate) fn closed_ring_fraction(b: &CellBox, coords: &[Coord]) -> f64 {
    let cell_area = b.area();
    if cell_area <= 0.0 {
        return 0.0;
    }
    area_signed_ee(coords.iter().copied()).abs() / cell_area
}

impl Chain {
    fn is_point(&self) -> bool {
        match &self.coords {
            ChainCoords::Range(r) => r.len() == 1,
            ChainCoords::Corner(_) => true,
        }
    }
}

/// Port of exactextract `left_hand_area`
/// ([src/traversal_areas.cpp](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/traversal_areas.cpp));
/// the GEOS-free version vendored by controlledburn. `chains` and `buf` are scratch.
pub(crate) fn left_hand_area(
    b: &CellBox,
    coords: &[Coord],
    lists: &[Range<usize>],
    chains: &mut Vec<Chain>,
    buf: &mut Vec<Coord>,
) -> Option<f64> {
    chains.clear();
    for r in lists {
        let s = coords.get(r.clone())?;
        chains.push(Chain {
            start: perimeter_distance(b, *s.first()?)?,
            stop: perimeter_distance(b, *s.last()?)?,
            coords: ChainCoords::Range(r.clone()),
            visited: false,
        });
    }
    let (h, w, perim) = (b.height(), b.width(), b.perimeter());
    for (cc, pd) in [
        (Coord::new(b.xmin, b.ymin), 0.0),
        (Coord::new(b.xmin, b.ymax), h),
        (Coord::new(b.xmax, b.ymax), h + w),
        (Coord::new(b.xmax, b.ymin), 2.0 * h + w),
    ] {
        chains.push(Chain { start: pd, stop: pd, coords: ChainCoords::Corner(cc), visited: false });
    }

    let mut sum = 0.0;
    for i in 0..chains.len() {
        if chains[i].visited || chains[i].is_point() {
            continue;
        }
        buf.clear();
        let mut cur = i;
        loop {
            chains[cur].visited = true;
            match &chains[cur].coords {
                ChainCoords::Range(r) => buf.extend_from_slice(coords.get(r.clone())?),
                ChainCoords::Corner(cc) => buf.push(*cc),
            }
            cur = next_chain(chains, cur, i, perim)?;
            if cur == i {
                break;
            }
        }
        buf.push(*buf.first()?);
        sum += area_signed_ee(buf.iter().copied()).abs();
    }
    Some(sum)
}

/// Port of exactextract `next_chain` ([src/traversal_areas.cpp#L48](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/traversal_areas.cpp#L48)).
fn next_chain(chains: &[Chain], cur: usize, kill: usize, perim: f64) -> Option<usize> {
    let stop = chains.get(cur)?.stop;
    let mut best = None;
    let mut min_d = f64::MAX;
    for (j, c) in chains.iter().enumerate() {
        if c.visited && j != kill {
            continue;
        }
        let d = perimeter_distance_ccw(stop, c.start, perim);
        if d < min_d {
            min_d = d;
            best = Some(j);
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::perimeter::perimeter_distance;
    use proptest::prelude::*;

    fn c(x: f64, y: f64) -> Coord {
        Coord::new(x, y)
    }
    fn b10() -> CellBox {
        CellBox { xmin: 0.0, ymin: 0.0, xmax: 10.0, ymax: 10.0 }
    }
    fn lha(lists: &[&[Coord]]) -> Option<f64> {
        let mut coords = Vec::new();
        let mut ranges = Vec::new();
        for l in lists {
            let s = coords.len();
            coords.extend_from_slice(l);
            ranges.push(s..coords.len());
        }
        left_hand_area(&b10(), &coords, &ranges, &mut Vec::new(), &mut Vec::new())
    }

    /// Test only: controlledburn `polygon_signed_area` ([cpp/src/ee/analytical_coverage.h#L31](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/ee/analytical_coverage.h#L31)) as an index loop.
    fn area_signed_reference(ring: &[Coord]) -> f64 {
        if ring.len() < 3 {
            return 0.0;
        }
        let x0 = ring[0].x;
        let mut sum = 0.0;
        for i in 1..ring.len() - 1 {
            sum += (ring[i].x - x0) * (ring[i - 1].y - ring[i + 1].y);
        }
        sum / 2.0
    }

    /// Test only: controlledburn `analytical_covered_fraction` ([cpp/src/ee/analytical_coverage.h#L58](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/ee/analytical_coverage.h#L58)) with a Vec.
    fn analytical_fraction_reference(b: &CellBox, coords: &[Coord]) -> Option<f64> {
        let cell_area = b.area();
        if cell_area <= 0.0 || coords.len() < 2 {
            return Some(0.0);
        }
        let perim = b.perimeter();
        let exit_pd = perimeter_distance(b, *coords.last()?)?;
        let entry_pd = perimeter_distance(b, coords[0])?;
        let arc = if exit_pd > entry_pd + 1e-12 {
            exit_pd - entry_pd
        } else if entry_pd > exit_pd + 1e-12 {
            perim - entry_pd + exit_pd
        } else {
            let mut poly = coords.to_vec();
            if poly[0] != poly[poly.len() - 1] {
                poly.push(poly[0]);
            }
            return Some(area_signed_reference(&poly).abs() / cell_area);
        };
        let (h, w) = (b.height(), b.width());
        let corners = [
            (c(b.xmin, b.ymin), 0.0),
            (c(b.xmin, b.ymax), h),
            (c(b.xmax, b.ymax), h + w),
            (c(b.xmax, b.ymin), 2.0 * h + w),
        ];
        let mut in_arc: Vec<(Coord, f64)> = Vec::new();
        for (cc, pd) in corners {
            let mut d = exit_pd - pd;
            if d < 0.0 {
                d += perim;
            }
            if d > 1e-12 && d < arc - 1e-12 {
                in_arc.push((cc, d));
            }
        }
        for i in 1..in_arc.len() {
            let key = in_arc[i];
            let mut j = i;
            while j > 0 && in_arc[j - 1].1 > key.1 {
                in_arc[j] = in_arc[j - 1];
                j -= 1;
            }
            in_arc[j] = key;
        }
        let mut poly = coords.to_vec();
        poly.extend(in_arc.iter().map(|t| t.0));
        poly.push(poly[0]);
        Some(area_signed_reference(&poly).abs() / cell_area)
    }

    // exactextract test/test_traversal_areas.cpp "Exit from same side as entry"
    #[test]
    fn lha_exit_same_side() {
        let t = [c(7.0, 0.0), c(7.0, 1.0), c(6.0, 1.0), c(6.0, 0.0)];
        assert_eq!(lha(&[&t]), Some(1.0));
        let r: Vec<Coord> = t.iter().rev().copied().collect();
        assert_eq!(lha(&[&r]), Some(99.0));
    }

    // "Enter bottom, exit left" and "Enter bottom, exit top"
    #[test]
    fn lha_single_traversals() {
        assert_eq!(lha(&[&[c(5.0, 0.0), c(5.0, 5.0), c(0.0, 5.0)]]), Some(25.0));
        assert_eq!(lha(&[&[c(4.0, 0.0), c(4.0, 10.0)]]), Some(40.0));
    }

    // "Multiple traversals (basic)" and "Multiple traversals"
    #[test]
    fn lha_multiple_traversals() {
        assert_eq!(lha(&[&[c(2.0, 10.0), c(2.0, 0.0)], &[c(4.0, 0.0), c(4.0, 10.0)]]), Some(20.0));
        let t1 = [c(2.0, 0.0), c(2.0, 2.0), c(0.0, 2.0)];
        let t2 = [c(3.0, 10.0), c(3.0, 0.0)];
        let t3 = [c(5.0, 0.0), c(5.0, 10.0)];
        let t4 = [c(8.0, 10.0), c(10.0, 8.0)];
        let t5 = [c(10.0, 6.0), c(8.0, 6.0), c(8.0, 3.0), c(10.0, 3.0)];
        let t6 = [c(10.0, 4.0), c(9.0, 4.0), c(9.0, 5.0), c(10.0, 5.0)];
        let t7 = [c(10.0, 3.0), c(8.0, 3.0), c(8.0, 0.0)];
        assert_eq!(
            lha(&[&t1, &t2, &t3, &t4, &t5, &t6, &t7]),
            Some(4.0 + 20.0 + 2.0 + 6.0 - 1.0 + 6.0)
        );
    }

    #[test]
    fn analytical_matches_lha_on_single_traversals() {
        let b = b10();
        let cases: [&[Coord]; 3] = [
            &[c(5.0, 0.0), c(5.0, 5.0), c(0.0, 5.0)],
            &[c(4.0, 0.0), c(4.0, 10.0)],
            &[c(7.0, 0.0), c(7.0, 1.0), c(6.0, 1.0), c(6.0, 0.0)],
        ];
        for t in cases {
            let a = analytical_fraction(&b, t).expect("on boundary");
            let l = lha(&[t]).expect("lha") / b.area();
            assert!((a - l).abs() < 1e-12, "{t:?}: {a} vs {l}");
        }
    }

    #[test]
    fn closed_ring_inside_cell() {
        let ring = [c(1.0, 1.0), c(2.0, 1.0), c(2.0, 2.0), c(1.0, 2.0), c(1.0, 1.0)];
        assert_eq!(closed_ring_fraction(&b10(), &ring), 0.01);
    }

    #[test]
    fn to_f32_narrows() {
        assert_eq!(to_f32(0.25), 0.25_f32);
    }

    #[test]
    fn cell_class_thresholds() {
        assert_eq!(CellClass::of(0.0), CellClass::Empty);
        assert_eq!(CellClass::of(1e-7), CellClass::Empty);
        assert_eq!(CellClass::of(0.5), CellClass::Edge(0.5));
        assert_eq!(CellClass::of(1.0 - 1e-7), CellClass::Full);
        assert_eq!(CellClass::of(1.0), CellClass::Full);
    }

    /// A point on the boundary of `b` at perimeter position `t` in [0, 4).
    fn on_boundary(b: &CellBox, t: f64) -> Coord {
        let (w, h) = (b.width(), b.height());
        match t {
            t if t < 1.0 => c(b.xmin, b.ymin + t * h),
            t if t < 2.0 => c(b.xmin + (t - 1.0) * w, b.ymax),
            t if t < 3.0 => c(b.xmax, b.ymax - (t - 2.0) * h),
            t => c(b.xmax - (t - 3.0) * w, b.ymin),
        }
    }

    proptest! {
        // O3: the buffer-free version is bit-equal to the C++-shaped reference.
        #[test]
        fn analytical_is_bit_equal_to_reference(
            x0 in -1e3..1e3_f64, y0 in -1e3..1e3_f64, w in 1e-3..1e2_f64, h in 1e-3..1e2_f64,
            t_in in 0.0..4.0_f64, t_out in 0.0..4.0_f64,
            inner in proptest::collection::vec((0.01..0.99_f64, 0.01..0.99_f64), 0..6),
        ) {
            let b = CellBox { xmin: x0, ymin: y0, xmax: x0 + w, ymax: y0 + h };
            let mut t = vec![on_boundary(&b, t_in)];
            t.extend(inner.iter().map(|&(u, v)| c(x0 + u * w, y0 + v * h)));
            t.push(on_boundary(&b, t_out));
            let got = analytical_fraction(&b, &t);
            let want = analytical_fraction_reference(&b, &t);
            prop_assert_eq!(got.map(f64::to_bits), want.map(f64::to_bits));
        }
    }
}
