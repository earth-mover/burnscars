//! Cell box geometry.
//!
//! Ported from exactextract (Apache-2.0, Copyright 2018-2022 ISciences, LLC),
//! commit 67969d99d1dfe1115be955087e2c48e1aab133e2: `src/box.h`, `src/box.cpp`.

use super::{coord::Coord, crossing::Crossing, side::Side};

/// Port of exactextract `Box`
/// ([src/box.h](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/box.h)).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct CellBox {
    pub(crate) xmin: f64,
    pub(crate) ymin: f64,
    pub(crate) xmax: f64,
    pub(crate) ymax: f64,
}

/// controlledburn `Location` ([cpp/src/scanline.cpp#L99](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp#L99)).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Location {
    Inside,
    Boundary,
    Outside,
}

impl CellBox {
    pub(crate) fn width(&self) -> f64 {
        self.xmax - self.xmin
    }

    pub(crate) fn height(&self) -> f64 {
        self.ymax - self.ymin
    }

    pub(crate) fn area(&self) -> f64 {
        self.width() * self.height()
    }

    /// Port of exactextract `Box::perimeter` ([src/box.h#L69](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/box.h#L69)):
    /// `2 * width() + 2 * height()`.
    pub(crate) fn perimeter(&self) -> f64 {
        2.0 * self.width() + 2.0 * self.height()
    }

    /// Port of exactextract `Box::side` ([src/box.cpp#L28](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/box.cpp#L28)).
    pub(crate) fn side(&self, c: Coord) -> Option<Side> {
        if c.x == self.xmin {
            Some(Side::Left)
        } else if c.x == self.xmax {
            Some(Side::Right)
        } else if c.y == self.ymin {
            Some(Side::Bottom)
        } else if c.y == self.ymax {
            Some(Side::Top)
        } else {
            None
        }
    }

    /// Port of exactextract `Box::contains(const Coordinate&)`
    /// ([src/box.cpp#L127](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/box.cpp#L127)).
    pub(crate) fn contains(&self, c: Coord) -> bool {
        c.x >= self.xmin && c.x <= self.xmax && c.y >= self.ymin && c.y <= self.ymax
    }

    /// Port of exactextract `Box::strictly_contains` ([src/box.cpp#L133](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/box.cpp#L133)).
    pub(crate) fn strictly_contains(&self, c: Coord) -> bool {
        c.x > self.xmin && c.x < self.xmax && c.y > self.ymin && c.y < self.ymax
    }

    /// controlledburn `point_location` ([cpp/src/scanline.cpp#L102](https://github.com/hypertidy/controlledburn/blob/1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558/cpp/src/scanline.cpp#L102)).
    pub(crate) fn location(&self, c: Coord) -> Location {
        if self.strictly_contains(c) {
            Location::Inside
        } else if self.contains(c) {
            Location::Boundary
        } else {
            Location::Outside
        }
    }

    /// Port of exactextract `Box::crossing`
    /// ([src/box.cpp](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/box.cpp)).
    /// `None` replaces the "Never get here" throw.
    pub(crate) fn crossing(&self, c1: Coord, c2: Coord) -> Option<Crossing> {
        let cr = |side, x, y| Some(Crossing { side, coord: Coord::new(x, y) });
        // vertical line
        if c1.x == c2.x {
            return if c2.y >= self.ymax {
                cr(Side::Top, c1.x, self.ymax)
            } else if c2.y <= self.ymin {
                cr(Side::Bottom, c1.x, self.ymin)
            } else {
                None
            };
        }
        // horizontal line
        if c1.y == c2.y {
            return if c2.x >= self.xmax {
                cr(Side::Right, self.xmax, c1.y)
            } else if c2.x <= self.xmin {
                cr(Side::Left, self.xmin, c1.y)
            } else {
                None
            };
        }
        let m = ((c2.y - c1.y) / (c2.x - c1.x)).abs();
        let up = c2.y > c1.y;
        let right = c2.x > c1.x;
        match (up, right) {
            (true, true) => {
                let y2 = c1.y + m * (self.xmax - c1.x);
                if y2 < self.ymax {
                    cr(Side::Right, self.xmax, clamp(y2, self.ymin, self.ymax))
                } else {
                    let x2 = c1.x + (self.ymax - c1.y) / m;
                    cr(Side::Top, clamp(x2, self.xmin, self.xmax), self.ymax)
                }
            }
            (true, false) => {
                let y2 = c1.y + m * (c1.x - self.xmin);
                if y2 < self.ymax {
                    cr(Side::Left, self.xmin, clamp(y2, self.ymin, self.ymax))
                } else {
                    let x2 = c1.x - (self.ymax - c1.y) / m;
                    cr(Side::Top, clamp(x2, self.xmin, self.xmax), self.ymax)
                }
            }
            (false, true) => {
                let y2 = c1.y - m * (self.xmax - c1.x);
                if y2 > self.ymin {
                    cr(Side::Right, self.xmax, clamp(y2, self.ymin, self.ymax))
                } else {
                    let x2 = c1.x + (c1.y - self.ymin) / m;
                    cr(Side::Bottom, clamp(x2, self.xmin, self.xmax), self.ymin)
                }
            }
            (false, false) => {
                let y2 = c1.y - m * (c1.x - self.xmin);
                if y2 > self.ymin {
                    cr(Side::Left, self.xmin, clamp(y2, self.ymin, self.ymax))
                } else {
                    let x2 = c1.x - (c1.y - self.ymin) / m;
                    cr(Side::Bottom, clamp(x2, self.xmin, self.xmax), self.ymin)
                }
            }
        }
    }
}

/// exactextract `clamp` ([src/box.cpp#L22](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/box.cpp#L22)): `min(max(x, low), high)`, not `f64::clamp`.
/// Uses the `std::max` / `std::min` comparisons, so signed zero and NaN match C++.
pub(crate) fn clamp(x: f64, low: f64, high: f64) -> f64 {
    let m = if x < low { low } else { x };
    if high < m { high } else { m }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bx(xmin: f64, ymin: f64, xmax: f64, ymax: f64) -> CellBox {
        CellBox { xmin, ymin, xmax, ymax }
    }
    fn near(a: Coord, x: f64, y: f64) -> bool {
        (a.x - x).abs() < 1e-14 && (a.y - y).abs() < 1e-14
    }

    // exactextract test/test_box.cpp "Box dimensions are calculated correctly"
    #[test]
    fn dimensions() {
        let b = bx(5.0, 7.0, 7.0, 10.0);
        assert_eq!((b.width(), b.height(), b.area(), b.perimeter()), (2.0, 3.0, 6.0, 10.0));
    }

    // "Coordinate sides are correctly identified"
    #[test]
    fn sides() {
        let b = bx(0.0, 0.0, 13.0, 17.0);
        assert_eq!(b.side(Coord::new(0.0, 5.0)), Some(Side::Left));
        assert_eq!(b.side(Coord::new(13.0, 5.0)), Some(Side::Right));
        assert_eq!(b.side(Coord::new(5.0, 17.0)), Some(Side::Top));
        assert_eq!(b.side(Coord::new(5.0, 0.0)), Some(Side::Bottom));
        assert_eq!(b.side(Coord::new(5.0, 5.0)), None);
    }

    #[test]
    fn location() {
        let b = bx(0.0, 0.0, 1.0, 1.0);
        assert_eq!(b.location(Coord::new(0.5, 0.5)), Location::Inside);
        assert_eq!(b.location(Coord::new(1.0, 0.5)), Location::Boundary);
        assert_eq!(b.location(Coord::new(1.5, 0.5)), Location::Outside);
    }

    // "Crossing is correctly determined"
    #[test]
    fn crossings() {
        let b = bx(0.0, 0.0, 1.0, 1.0);
        let o = Coord::new(0.4, 0.6);
        let cases = [
            ((0.4, 11.0), Side::Top, (0.4, 1.0)),
            ((0.4, -11.0), Side::Bottom, (0.4, 0.0)),
            ((-11.0, 0.6), Side::Left, (0.0, 0.6)),
            ((11.0, 0.6), Side::Right, (1.0, 0.6)),
            ((1.2, 0.8), Side::Right, (1.0, 0.75)),
            ((0.7, 1.2), Side::Top, (0.6, 1.0)),
            ((0.1, 1.2), Side::Top, (0.2, 1.0)),
            ((-0.4, 1.2), Side::Left, (0.0, 0.9)),
            ((-0.4, 0.2), Side::Left, (0.0, 0.4)),
            ((0.1, -0.3), Side::Bottom, (0.2, 0.0)),
            ((0.7, -0.3), Side::Bottom, (0.6, 0.0)),
            ((1.6, 0.0), Side::Right, (1.0, 0.3)),
        ];
        for ((x2, y2), side, (ex, ey)) in cases {
            let c = b.crossing(o, Coord::new(x2, y2)).expect("crossing");
            assert_eq!(c.side, side, "to ({x2}, {y2})");
            assert!(near(c.coord, ex, ey), "to ({x2}, {y2}): {:?}", c.coord);
        }
    }

    // C++ `std::min(std::max(x, low), high)` keeps x for -0.0 and NaN.
    #[test]
    fn clamp_matches_cpp_for_signed_zero_and_nan() {
        assert_eq!(clamp(-0.0, 0.0, 1.0).to_bits(), (-0.0_f64).to_bits());
        assert!(clamp(f64::NAN, 0.0, 1.0).is_nan());
        assert_eq!(clamp(2.0, 0.0, 1.0), 1.0);
        assert_eq!(clamp(-2.0, 0.0, 1.0), 0.0);
    }

    #[test]
    fn crossing_that_does_not_leave_is_none() {
        let b = bx(0.0, 0.0, 1.0, 1.0);
        assert_eq!(b.crossing(Coord::new(0.4, 0.2), Coord::new(0.4, 0.8)), None);
    }
}
