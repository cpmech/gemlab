# AGENTS.md

Rust library for geometry, mesh generation/IO, and numerical integration for FEM/FEA.

## Build backend (read this first)

`russell_lab` 4 links **nothing external by default**: the BLAS/LAPACK kernels are
pure Rust (`algebra_linear`), so a plain `cargo build` and `cargo test` work on a
machine with no BLAS installed. Two opt-in backends exist:

- **OpenBLAS** — `open_blas`, resolved via `pkg-config`, with a Homebrew fallback.
  Requires the OpenBLAS/LAPACK development headers (e.g. `lapack.h`).
- **Intel MKL** — `intel_mkl`, built from `/opt/intel/oneapi/mkl/<version>`; set
  `MKL_VERSION` to override the default (`latest`).

Both are forwarded from every dependent crate (`russell_tensor`, `russell_sparse`,
`russell_stat`, `russell_ode`, `russell_pde`, `russell_nonlin`). Reach for
`--features intel_mkl` when you want to *test the accelerated path*, not because the
default build needs it.

The one native library the repository does need is SuiteSparse, for the `umfpack`
feature.


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

## Test / lint commands

```bash
# Run tests (no test is feature-gated, so the default build covers the suite)
cargo test

# Formatting (CI runs the check)
cargo fmt --all
cargo fmt --all -- --check

# Lint
cargo clippy --all-targets --all-features

# Coverage: one instrumented run, which also runs the tests.
# CI does this as a single step -- `cpmech/rust-cov@v1`, min 90 -- because
# `--fail-under-lines` counts a line once per monomorphization and reports this
# crate lower than it is. To reproduce locally:
cargo llvm-cov --workspace --features umfpack \
  --ignore-filename-regex 'build.rs|src/bin/' --lcov --output-path lcov.info
curl -sL https://raw.githubusercontent.com/cpmech/rust-cov/v1/check-coverage.bash \
  | bash -s -- lcov.info 90
```

`src/bin/` is excluded from the gate on purpose: the five mesh converters
(`msh2tri`, `msh2tet`, `hex2msh`, `qua2msh`, `drawmsh`) have no tests, and the gate
measures library code. Coverage of the library is 91.94% on that basis; the bare
`cargo llvm-cov` summary reads 91.29% because of same-line re-counting.


## Extra features:

- `umfpack` and `mumps` forward to `russell_sparse`. `mumps` requires
  `RUSSELL_LOCAL_PREFIX` with `include/mumps` and `lib/mumps` (see
  `zscripts/arch-compile-mumps.bash`); `umfpack` prefers
  `$RUSSELL_LOCAL_PREFIX/include/umfpack` and otherwise uses the system library.

## Architecture
- Modules: `geometry` (entities, triangle/tetrahedron), `graph` (directed/undirected), `integ`
  (element matrices `mat_*`/vectors `vec_*`, Gauss rules), `mesh` (+ `mesh/algorithms`), `recovery`
  (extrapolation/interpolation), `shapes` (element shape functions + `Scratchpad`), `util`.
- Common API is re-exported through `gemlab::prelude::*` (Mesh, Cell, Features, Draw, Scratchpad, GeoKind...).
- Errors are `gemlab::StrError = &'static str`; fallible functions return `Result<T, StrError>`.

## Build prerequisites (system libraries)
Nothing, for a default build: `russell_lab` 4 links no external BLAS/LAPACK. The
optional backends need the matching library:
- `open_blas`: OpenBLAS with its LAPACK headers
- `umfpack`: SuiteSparse, or a local install under `RUSSELL_LOCAL_PREFIX`
- `mumps`: locally compiled MUMPS at `RUSSELL_LOCAL_PREFIX`
- Arch: `pacman -Syu blas-openblas suitesparse python-matplotlib`
- Debian/Ubuntu: `sudo apt-get install -y libopenblas-dev libsuitesparse-dev python3-matplotlib`

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
