# Geometry laboratory for finite element analyses

This repository contains structures and functions to perform geometry computations and generate meshes for finite element analyses.

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

## Examples

TODO
