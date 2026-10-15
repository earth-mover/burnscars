//! Per-worker buffers, reused across geometries (spec: Storage, O6).

use crate::approx::ApproxScratch;
use crate::coverage::{CoverageScratch, walk::Arena};
use crate::geom::coord::Coord;
use crate::input::InputBuf;

#[derive(Debug, Default)]
pub(crate) struct WalkScratch {
    pub(crate) ring: Vec<Coord>,
    pub(crate) arena: Arena,
}

#[derive(Debug, Default)]
pub(crate) struct Scratch {
    pub(crate) input: InputBuf,
    pub(crate) walk: WalkScratch,
    pub(crate) cov: CoverageScratch,
    pub(crate) approx: ApproxScratch,
}
