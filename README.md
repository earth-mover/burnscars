# burnscars

Sparse rasterization of polygons, lines and points onto a regular grid, in Rust.
Polygons get exact per-cell coverage fractions (coverage mode) or the cell-centre
rule (approx mode). Time and memory are O(perimeter): the output is sparse tables,
never a dense raster.

burnscars is a port of the controlledburn C++ core. It has no internal parallelism.
Its grid windows use global lattice indices, so a caller can split the work by
geometry or by grid chunk and run it in parallel.

## Inspiration and credits

- [controlledburn](https://github.com/hypertidy/controlledburn) (Michael Sumner): the
  code that this crate ports; sparse output, the winding sweep, the approx mode.
- [exactextract](https://github.com/isciences/exactextract) (Daniel Baston, ISciences):
  the exact coverage algorithm (cell traversal, crossing, `left_hand_area`); the
  reference for correctness.
- [fasterize](https://github.com/hypertidy/fasterize) (moved from
  `ecohealthalliance/fasterize` in March 2025): the cell-centre rule and the pixel
  functions that approx mode and a later materialize follow.
- [rusterize](https://github.com/ttrotto/rusterize): prior art for scanline
  rasterization in Rust.
- GDAL rasterize: the cell-centre oracle for approx mode (through rasterio).

## Upstream and rules

- Ported from controlledburn `1f8a25fd030ac6ec1aa2fbf2c69ee35df8c17558` and
  exactextract `67969d99d1dfe1115be955087e2c48e1aab133e2`.
- If controlledburn and exactextract disagree, exactextract wins. The oracle suite
  decides.
- No code comes from GEOS or JTS. GEOS is used only as a black-box test oracle
  (through shapely, in `oracle/`).

## Output

Four tables in struct-of-arrays form: `runs` (interior intervals), `edges` (boundary
cells with a coverage fraction), `lines` (length in cell), `points`. Indices are
global, 0-based lattice indices; `col_end` is exclusive. See the crate docs for the
relation to CSR and CSC sparse matrices.

## Development

`just` lists the tasks. `just ci` runs format, lint, tests and doctests.
`just oracle` regenerates the oracle fixtures (needs `uv`).

## TODO

- R1: rotate each ring to start on a cell boundary, which removes the reprocess
  branch. Do this only after the oracle suite passes.
- R2: snap each vertex when the walker reads it, so the walker borrows the input with
  no copy (needs R1).
- R3: prepared geometry: row-bucketed edges in a `Sync` struct, shared with `Arc`
  across brix chunks, so each chunk scans only its row band.
- R4: pre-clip each ring to the window, for very large polygons on small chunks.
- Exact predicates: use `robust` (Shewchuk `orient2d`) for the ring orientation and the
  side tests, in place of float arithmetic. exactextract does not do this, so the
  oracle fixtures can change on near-degenerate input. Decide which result is correct
  before the change, and record each changed fixture as a deliberate difference.
- Merge adjacent one-cell runs and interior runs in a row.
- A performance regression gate on the benchmarks.
- Antimeridian, as an opt-in `Wrap::Longitude` lattice option. Conventions: A, split
  at ±180 (works already as a multipolygon); B, continuous past 180 (burn shifted
  copies at k·360 and wrap the columns); C, jumps at ±180 (unwrap where |Δx| > 180,
  then B); D, rings that enclose a pole (close the ring through the pole). Test: burn
  the split form A and compare.
- DE-9IM per cell from the compressed form: fraction, a one-byte contact mask per
  edge cell, a `touches` table, and one bit per geometry (see the design spec).
- Materialize and pixel functions, as a brix adapter.
- Benchmarks against controlledburn C++ and rusterize.
- Python bindings.
