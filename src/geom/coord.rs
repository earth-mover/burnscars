//! Planar coordinate. Port of exactextract `Coordinate` (`src/coordinate.h`,
//! commit 67969d99d1dfe1115be955087e2c48e1aab133e2, Apache-2.0).

/// Port of exactextract `Coordinate`
/// ([src/coordinate.h](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/coordinate.h)).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Coord {
    pub(crate) x: f64,
    pub(crate) y: f64,
}

impl Coord {
    pub(crate) const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}
