//! Public entry point.

use geo_traits::GeometryTrait;

use crate::approx;
use crate::coverage;
use crate::error::BurnError;
use crate::grid::{GeomId, Grid};
use crate::input::{Loaded, for_each_xy};
use crate::output::{BurnResult, GeomStart};
use crate::scratch::Scratch;
use crate::{line, point};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Mode {
    /// Exact coverage fractions (controlledburn `BurnMode::Coverage`).
    Coverage,
    /// Cell-centre rule (controlledburn `BurnMode::Approx`, fasterize semantics).
    Approx,
}

/// `Send` and `Sync`; use one per worker, because it owns reused scratch buffers.
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
        let loaded = self.scratch.input.load(g, id)?;
        if grid.lattice().is_north_up() {
            self.burn_loaded(loaded, grid, id, out)?;
        } else {
            self.burn_mapped(loaded, grid, id, start, out)?;
        }
        out.end_geom(id, start);
        Ok(())
    }

    /// The ported algorithms, on a north-up grid.
    fn burn_loaded(
        &mut self,
        loaded: Loaded,
        grid: &Grid,
        id: GeomId,
        out: &mut BurnResult,
    ) -> Result<(), BurnError> {
        let s = &mut self.scratch;
        match loaded {
            Loaded::Empty => {}
            Loaded::Polygons => match self.mode {
                Mode::Coverage => coverage::burn(&s.input, grid, id, &mut s.walk, &mut s.cov, out)?,
                Mode::Approx => approx::burn(&s.input, grid, id, &mut s.approx, out)?,
            },
            Loaded::Lines => line::burn(&s.input, grid, id, &mut s.walk, out)?,
            Loaded::Points => point::burn(&s.input, grid, id, out),
        }
        Ok(())
    }

    /// Signed steps and a period: burns on the north-up grid, one copy per period that
    /// meets the geometry, then maps the records to the caller's indices.
    fn burn_mapped(
        &mut self,
        loaded: Loaded,
        grid: &Grid,
        id: GeomId,
        start: GeomStart,
        out: &mut BurnResult,
    ) -> Result<(), BurnError> {
        let n = grid.north_up();
        let (period, ks) = match (n.period, self.scratch.input.envelope()) {
            (Some(period), Some(env)) => {
                let width = env.xmax - env.xmin;
                if width > period {
                    return Err(BurnError::WiderThanPeriod { id, width, period });
                }
                let l = n.grid.lattice();
                let west = l.x_edge(i64::from(n.grid.cols().start().0));
                let east = l.x_edge(i64::from(n.grid.cols().end().0));
                // Copy k covers [west + k p, east + k p]; keep each copy that meets the envelope.
                let lo = to_i64(((env.xmin - east) / period).ceil());
                let hi = to_i64(((env.xmax - west) / period).floor());
                (period, lo..=hi)
            }
            _ => (0.0, 0..=0),
        };
        let copies = ks.map(|k| n.shifted(k, period));
        if loaded == Loaded::Points {
            point::burn_shifted(&self.scratch.input, &copies, id, out);
        } else {
            let mut copies_hit = 0;
            for g in copies {
                let before = out.len();
                self.burn_loaded(loaded, &g, id, out)?;
                copies_hit += usize::from(out.len() > before);
            }
            if copies_hit > 1 {
                out.merge_tail(id, start);
            }
        }
        out.map_tail(start, n.rows, n.cols);
        Ok(())
    }
}

/// Saturating: a copy that far away does not meet the window.
#[expect(clippy::cast_possible_truncation, reason = "saturating float to int cast")]
fn to_i64(v: f64) -> i64 {
    v as i64
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
