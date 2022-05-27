# Geometry, meshes, and integration for finite element analyses

This repository contains structures and functions to perform geometry computations, generate meshes, and perform numerical integration for finite element analyses (FEM/FEA).

We use `Vector` and `Matrix` from [Russell Lab](https://github.com/cpmech/russell), thus some Debian packages are required.

Documentation:

- [API reference (docs.rs)](https://docs.rs/gemlab)

## Installation

Install some libraries:

```bash
sudo apt-get install \
    liblapacke-dev \
    libopenblas-dev
```

Add this to your Cargo.toml:

```toml
[dependencies]
gemlab = "0.2"
```

## Todo

- [x] Implement read/write mesh functions
- [x] Add tests for the numerical integrations
- [ ] Implement triangle and tetrahedron generators
- [ ] Implement drawing functions

## Examples

TODO
