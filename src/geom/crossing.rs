//! Port of exactextract `Crossing` (`src/crossing.h`, commit
//! 67969d99d1dfe1115be955087e2c48e1aab133e2, Apache-2.0).

use super::{coord::Coord, side::Side};

/// Port of exactextract `Crossing`
/// ([src/crossing.h](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/crossing.h)).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Crossing {
    pub(crate) side: Side,
    pub(crate) coord: Coord,
}
