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

impl Scratch {
    /// The capacity of every buffer, in bytes.
    pub(crate) fn bytes(&self) -> usize {
        let Self { input, walk, cov, approx } = self;
        input.bytes() + vec_bytes(&walk.ring) + walk.arena.bytes() + cov.bytes() + approx.bytes()
    }
}

/// The capacity of `v`, in bytes.
pub(crate) fn vec_bytes<T>(v: &Vec<T>) -> usize {
    v.capacity() * size_of::<T>()
}
