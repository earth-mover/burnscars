//! Public entry point.

use geo_traits::GeometryTrait;

use crate::approx;
use crate::coverage;
use crate::error::BurnError;
use crate::grid::{GeomId, Grid};
use crate::input::{Loaded, for_each_xy};
use crate::output::BurnResult;
use crate::scratch::Scratch;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Mode {
    /// Exact coverage fractions (controlledburn `BurnMode::Coverage`).
    Coverage,
    /// Cell-centre rule (controlledburn `BurnMode::Approx`, fasterize semantics).
    Approx,
}

/// One per worker. `Send`, not `Sync`. Owns the scratch buffers and reuses them.
#[derive(Debug)]
pub struct Burner {
    mode: Mode,
    scratch: Scratch,
}

impl Burner {
    #[must_use]
    pub fn new(mode: Mode) -> Self {
        Self { mode, scratch: Scratch::default() }
    }

    #[must_use]
    pub fn mode(&self) -> Mode {
        self.mode
    }

    /// Burns each `(id, geometry)` onto the window of `grid` and appends to `out`.
    ///
    /// Stops at the first error. There is no rollback. Each geometry is validated before
    /// it writes, and its `GeomIndex` is pushed only when it is complete. After
    /// `WalkDiverged` or `Internal`, records of the failing geometry can follow the last
    /// `GeomIndex`; truncate each table to the end of the last `GeomIndex` to drop them.
    pub fn burn_into<G, I>(
        &mut self,
        grid: &Grid,
        geoms: I,
        out: &mut BurnResult,
    ) -> Result<(), BurnError>
    where
        G: GeometryTrait<T = f64>,
        I: IntoIterator<Item = (GeomId, G)>,
    {
        geoms.into_iter().try_for_each(|(id, g)| self.burn_one(grid, id, &g, out))
    }

    pub fn burn<G, I>(&mut self, grid: &Grid, geoms: I) -> Result<BurnResult, BurnError>
    where
        G: GeometryTrait<T = f64>,
        I: IntoIterator<Item = (GeomId, G)>,
    {
        let mut out = BurnResult::default();
        self.burn_into(grid, geoms, &mut out)?;
        Ok(out)
    }

    fn burn_one<G: GeometryTrait<T = f64>>(
        &mut self,
        grid: &Grid,
        id: GeomId,
        g: &G,
        out: &mut BurnResult,
    ) -> Result<(), BurnError> {
        let start = out.start_geom();
        let s = &mut self.scratch;
        match s.input.load(g, id)? {
            Loaded::Empty => {}
            Loaded::Polygons => match self.mode {
                Mode::Coverage => coverage::burn(&s.input, grid, id, &mut s.walk, &mut s.cov, out)?,
                Mode::Approx => approx::burn(&s.input, grid, id, &mut s.icpt, out)?,
            },
            Loaded::Lines => {
                return Err(BurnError::Internal { id, detail: "line path lands in Task 10" });
            }
            Loaded::Points => {
                return Err(BurnError::Internal { id, detail: "point path lands in Task 10" });
            }
        }
        out.end_geom(id, start);
        Ok(())
    }
}

/// Envelope of a geometry (finite coordinates only), as `[xmin, ymin, xmax, ymax]`.
/// Lets an orchestrator select the geometries for a chunk.
pub fn bbox<G: GeometryTrait<T = f64>>(geom: &G) -> Option<[f64; 4]> {
    let mut env: Option<[f64; 4]> = None;
    for_each_xy(geom, &mut |x, y| {
        if x.is_finite() && y.is_finite() {
            env = Some(match env {
                None => [x, y, x, y],
                Some([a, b, c, d]) => [a.min(x), b.min(y), c.max(x), d.max(y)],
            });
        }
    });
    env
}
