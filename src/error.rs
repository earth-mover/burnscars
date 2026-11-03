//! Errors for one burn.

use std::fmt;

use crate::grid::GeomId;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum GeomKind {
    Point,
    LineString,
    Polygon,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
    GeometryCollection,
    Rect,
    Triangle,
    Line,
}

/// Where in a geometry an error occurred.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct At {
    pub id: GeomId,
    pub part: usize,
    pub ring: usize,
    pub vertex: usize,
}

impl fmt::Display for At {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "geometry {}, part {}, ring {}, vertex {}",
            self.id, self.part, self.ring, self.vertex
        )
    }
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum BurnError {
    #[error("geometry {id}: unsupported type {kind:?}; split collections into one type each")]
    Unsupported { id: GeomId, kind: GeomKind },
    #[error("{at}: non-finite coordinate ({x}, {y})")]
    NonFinite { at: At, x: f64, y: f64 },
    #[error("{at}: {n} coordinates; a closed ring needs at least 4")]
    ShortRing { at: At, n: usize },
    #[error("{at}: 1 coordinate; a line needs at least 2")]
    ShortLine { at: At },
    #[error("{at}: walk did not converge after {coords} coordinates")]
    WalkDiverged { at: At, coords: usize },
    #[error("geometry {id}: x extent {width} is wider than the lattice period {period}")]
    WiderThanPeriod { id: GeomId, width: f64, period: f64 },
    #[error("geometry {id}: internal error: {detail}")]
    Internal { id: GeomId, detail: &'static str },
}
