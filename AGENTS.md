# AGENTS.md

Rust library for geometry, mesh generation/IO, and numerical integration for FEM/FEA.

## Build backend (read this first)

`russell_lab` links BLAS/LAPACK through the C shim in
`russell_lab/c_code/interface_blas.c`. Two backends are available:

- **OpenBLAS** (default) — resolved via `pkg-config`, with a Homebrew fallback.
  This requires the OpenBLAS/LAPACK development headers (e.g. `lapack.h`).
- **Intel MKL** — enabled with the `intel_mkl` feature. MKL is expected under
  `/opt/intel/oneapi/mkl/<version>`; set `MKL_VERSION` to override the default
  (`latest`).

On machines without a complete OpenBLAS/LAPACK install the default backend fails
to compile. **Use `--features intel_mkl` for local builds and tests.**

The `intel_mkl` feature is forwarded from every dependent crate
(`russell_tensor`, `russell_sparse`, `russell_stat`, `russell_ode`,
`russell_pde`, `russell_nonlin`).

## Code intelligence (CodeGraph)

This repository is indexed by CodeGraph (a `.codegraph/` directory exists at the
root; see the global AGENTS.md for the full guidance). Reach for
`codegraph_explore` **before** grep/read when you need to understand or locate
code — one call returns the relevant symbols' verbatim source plus the call
paths, which is cheaper and more accurate than a search/read loop.

The index keeps itself fresh: a file watcher with a debounced auto-sync
(~2 s, `CODEGRAPH_WATCH_DEBOUNCE_MS`), a per-file staleness banner on tool
responses, and a connect-time catch-up. A path can therefore look stale for a
few seconds right after a rename/create — this once surfaced `eigen2_values.rs`
moments after it had been renamed to `eigen_values.rs` (the next check already
showed the correct name). When that happens, or when a tool response carries the
staleness banner, verify the filename on disk (`ls` / glob) and `Read` the
specific file for line-level edits.

If a path is still wrong after the debounce window — or the watcher is disabled
(e.g. a sandbox, or `CODEGRAPH_NO_DAEMON=1`) — force a refresh:
`codegraph status` (reports a `### Pending sync:` list), `codegraph sync`
(incremental) or `codegraph index` (full rebuild); `codegraph unlock` clears a
stale lock file.

## Extra features:

- `local_sparse` (and `cudss`) on `russell_sparse`, `russell_ode`,
  `russell_pde`, `russell_nonlin` — require locally compiled MUMPS/SuiteSparse
  (see `zscripts/*-compile-mumps.bash` and `zscripts/*-compile-suitesparse.bash`).

## Architecture
- Modules: `geometry` (entities, triangle/tetrahedron), `graph` (directed/undirected), `integ`
  (element matrices `mat_*`/vectors `vec_*`, Gauss rules), `mesh` (+ `mesh/algorithms`), `recovery`
  (extrapolation/interpolation), `shapes` (element shape functions + `Scratchpad`), `util`.
- Common API is re-exported through `gemlab::prelude::*` (Mesh, Cell, Features, Draw, Scratchpad, GeoKind...).
- Errors are `gemlab::StrError = &'static str`; fallible functions return `Result<T, StrError>`.

## Build prerequisites (system libraries)
The crate links BLAS/LAPACK + SuiteSparse via `russell`:
- Arch: `pacman -Syu blas-openblas python-matplotlib suitesparse`
- Debian/Ubuntu: `sudo apt-get install -y liblapacke-dev libopenblas-dev libsuitesparse-dev python3-matplotlib`
- Features (re-exported to `russell_*`): default uses OpenBLAS+SuiteSparse; `intel_mkl`, `local_sparse`,
  `cudss` require the matching libraries installed locally.

## Commands
- Full tests (same as CI): `cargo test -- --nocapture`
- Single test: `cargo test test_name`
- Run an example: `cargo run --example <name>` (see `examples/README.md`)
- All examples: `bash zscripts/run-examples.bash` (skips `check_grid_search_performance`, ~7 min)
- Format only: `cargo fmt` (`rustfmt.toml` sets `max_width = 120`); no clippy/lint gate in CI.

## Gotchas
- `src/lib.rs` compiles `README.md` as a doctest (`#[cfg(doctest)]`): keep README Rust blocks valid.
- Tests/examples save plots to `/tmp/gemlab` via `plotpy`; `data/figures/` holds generated SVGs.
- Unit tests are inline `#[cfg(test)]` modules under `src/**`; integration tests are `tests/test_*.rs`.
- Mesh fixtures: `data/meshes/*.msh` (`bad_*.msh` are intentionally invalid); JSON inputs in `data/input/`.
- `Cargo.lock` is gitignored (library crate) — don't try to commit it.
- CLI binaries in `src/bin/` (`drawmsh`, `hex2msh`, `msh2tet`, `msh2tri`, `qua2msh`) use `structopt`.
