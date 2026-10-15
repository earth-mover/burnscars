//! Perimeter distance along a cell boundary.
//!
//! Ported from exactextract (Apache-2.0, Copyright 2018-2022 ISciences, LLC),
//! commit 67969d99d1dfe1115be955087e2c48e1aab133e2: `src/perimeter_distance.cpp`.

use super::{box_::CellBox, coord::Coord};

/// Port of exactextract `perimeter_distance`
/// ([src/perimeter_distance.cpp](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/perimeter_distance.cpp)).
/// `None` replaces the "Never get here" throw.
pub(crate) fn perimeter_distance(b: &CellBox, c: Coord) -> Option<f64> {
    if c.x == b.xmin {
        Some(c.y - b.ymin) // left
    } else if c.y == b.ymax {
        Some((b.ymax - b.ymin) + c.x - b.xmin) // top
    } else if c.x == b.xmax {
        Some((b.xmax - b.xmin) + (b.ymax - b.ymin) + b.ymax - c.y) // right
    } else if c.y == b.ymin {
        Some((b.xmax - b.xmin) + 2.0 * (b.ymax - b.ymin) + (b.xmax - c.x)) // bottom
    } else {
        None
    }
}

/// Port of exactextract `perimeter_distance_ccw`.
pub(crate) fn perimeter_distance_ccw(measure1: f64, measure2: f64, perimeter: f64) -> f64 {
    if measure2 <= measure1 { measure1 - measure2 } else { perimeter + measure1 - measure2 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perimeter_distance_from_origin() {
        let b = CellBox { xmin: 0.0, ymin: 0.0, xmax: 13.0, ymax: 17.0 };
        let cases = [
            ((0.0, 0.0), 0.0),
            ((0.0, 5.0), 5.0),
            ((0.0, 17.0), 17.0),
            ((6.0, 17.0), 23.0),
            ((13.0, 17.0), 30.0),
            ((13.0, 14.0), 33.0),
            ((13.0, 0.0), 47.0),
            ((5.0, 0.0), 55.0),
        ];
        for ((x, y), d) in cases {
            assert_eq!(perimeter_distance(&b, Coord::new(x, y)), Some(d), "({x}, {y})");
        }
        assert_eq!(perimeter_distance(&b, Coord::new(5.0, 5.0)), None);
    }

    #[test]
    fn ccw_distance_between_measures() {
        assert_eq!(perimeter_distance_ccw(3.0, 2.0, 100.0), 1.0);
        assert_eq!(perimeter_distance_ccw(0.5, 3.5, 4.0), 1.0);
        assert_eq!(perimeter_distance_ccw(1.0, 1.0, 100.0), 0.0);
    }
}
