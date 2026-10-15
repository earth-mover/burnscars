//! Ring areas with the two upstream sign conventions.
//!
//! `area_signed_ee`: exactextract `area_signed` (Apache-2.0, Copyright 2018-2020 ISciences,
//! LLC), commit 67969d99d1dfe1115be955087e2c48e1aab133e2, `src/measures.cpp`; clockwise is
//! positive. `signed_area_ccw`: controlledburn `signed_area` (Apache-2.0, Copyright 2025
//! Michael Sumner), commit 1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558, `cpp/src/scanline.cpp`;
//! counter-clockwise is positive.

use super::coord::Coord;

/// Port of exactextract `area_signed`
/// ([src/measures.cpp](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/measures.cpp)).
/// A left fold in input order, so the bits equal the C++ loop.
pub(crate) fn area_signed_ee<I: IntoIterator<Item = Coord>>(ring: I) -> f64 {
    let mut it = ring.into_iter();
    let (Some(p0), Some(mut cur)) = (it.next(), it.next()) else {
        return 0.0;
    };
    let mut prev = p0;
    let sum = it.fold(0.0, |sum, next| {
        let s = sum + (cur.x - p0.x) * (prev.y - next.y);
        prev = cur;
        cur = next;
        s
    });
    sum / 2.0
}

/// Port of controlledburn `signed_area`
/// ([cpp/src/scanline.cpp](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp)).
pub(crate) fn signed_area_ccw(ring: &[Coord]) -> f64 {
    let Some(p0) = ring.first() else {
        return 0.0;
    };
    let sum = ring.windows(3).fold(0.0, |sum, w| sum + (w[1].x - p0.x) * (w[2].y - w[0].y));
    sum / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square_ccw() -> Vec<Coord> {
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]
            .into_iter()
            .map(|(x, y)| Coord::new(x, y))
            .collect()
    }

    #[test]
    fn sign_conventions_are_opposite() {
        let sq = square_ccw();
        assert_eq!(area_signed_ee(sq.iter().copied()), -1.0);
        assert_eq!(signed_area_ccw(&sq), 1.0);
        let rev: Vec<Coord> = sq.iter().rev().copied().collect();
        assert_eq!(area_signed_ee(rev.iter().copied()), 1.0);
        assert_eq!(signed_area_ccw(&rev), -1.0);
    }

    #[test]
    fn short_rings_have_zero_area() {
        let sq = square_ccw();
        assert_eq!(area_signed_ee(sq[..2].iter().copied()), 0.0);
        assert_eq!(signed_area_ccw(&sq[..2]), 0.0);
    }
}
