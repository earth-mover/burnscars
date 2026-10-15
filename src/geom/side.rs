//! Cell side. Port of exactextract `Side` (`src/side.h`, commit
//! 67969d99d1dfe1115be955087e2c48e1aab133e2, Apache-2.0). `Side::NONE` becomes `Option<Side>`.

/// Port of exactextract `Side`
/// ([src/side.h](https://github.com/isciences/exactextract/blob/67969d99d1dfe1115be955087e2c48e1aab133e2/src/side.h)).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Side {
    Left,
    Right,
    Top,
    Bottom,
}
