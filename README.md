# Geometry, meshes, and numerical integration for finite element analyses <!-- omit from toc -->

[![codecov](https://codecov.io/gh/cpmech/gemlab/graph/badge.svg?token=OjQKQ0PrNF)](https://codecov.io/gh/cpmech/gemlab)
[![Test & Coverage](https://github.com/cpmech/gemlab/actions/workflows/test_and_coverage.yml/badge.svg)](https://github.com/cpmech/gemlab/actions/workflows/test_and_coverage.yml)

[![documentation](https://img.shields.io/badge/gemlab-documentation-blue)](https://docs.rs/gemlab)

## Contents <!-- omit from toc -->

- [Introduction](#introduction)
  - [Documentation](#documentation)
- [Installation](#installation)
  - [TL;DR (Debian/Ubuntu/Linux)](#tldr-debianubuntulinux)
  - [Details](#details)
  - [Setting Cargo.toml](#setting-cargotoml)
- [Examples](#examples)
- [Roadmap](#roadmap)
- [Appendix A - Shapes and local numbering of nodes](#appendix-a---shapes-and-local-numbering-of-nodes)
  - [Lines (Lin)](#lines-lin)
  - [Triangles (Tri)](#triangles-tri)
  - [Quadrilaterals (Qua)](#quadrilaterals-qua)
  - [Tetrahedra (Tet)](#tetrahedra-tet)
  - [Hexahedra (Hex)](#hexahedra-hex)
- [Appendix B - Geometry versus space dimensions](#appendix-b---geometry-versus-space-dimensions)



## Introduction

This crate contains structures and functions for geometry computations, generate meshes, and perform numerical integration for finite element analyses (FEM/FEA).

### Documentation

[![documentation](https://img.shields.io/badge/gemlab-documentation-blue)](https://docs.rs/gemlab)


## Installation

At this moment, Gemlab works on **Linux** (Debian/Ubuntu; and maybe Arch).

### TL;DR (Debian/Ubuntu/Linux)

First:

```bash
sudo apt-get install -y --no-install-recommends \
    g++ \
    gdb \
    gfortran \
    liblapacke-dev \
    libmumps-seq-dev \
    libopenblas-dev \
    libsuitesparse-dev
```

Then:

```bash
cargo add gemlab
```

### Details

This crates depends on `russell_lab` and, hence, needs some external libraries. See the [installation of required dependencies](https://github.com/cpmech/russell) on `russell_lab`.

### Setting Cargo.toml

[![Crates.io](https://img.shields.io/crates/v/gemlab.svg)](https://crates.io/crates/gemlab)

👆 Check the crate version and update your Cargo.toml accordingly:

```toml
[dependencies]
gemlab = "*"
```



## Examples

```rust
use gemlab::integ;
use gemlab::mesh::{At, Features, Mesh};
use gemlab::shapes::Scratchpad;
use gemlab::StrError;
use std::collections::HashSet;

fn main() -> Result<(), StrError> {
    // Input the raw mesh data using a text file
    //
    // 1.0  5------,6.------7
    //      | [3],'   `.[4] |
    //      |  ,'       `.  |
    //      |,'           `.|
    // 0.5  3      [2]      4
    //      |`.           .'|
    //      |  `.       .'  |
    //      | [0]`.   .'[1] |
    // 0.0  0------`1'------2
    //     0.0     0.5     1.0
    let path = "./data/meshes/four_tri3_one_qua4.msh";
    let mesh = Mesh::from_text_file(path)?;

    // Extract features such boundary edges and faces.
    // Search entities along the boundary of the mesh given coordinates.
    // The `At` enum provides an easy way to define the type of the
    // constraint such as line, plane, circle, etc.
    let feat = Features::new(&mesh, false);
    assert_eq!(feat.search_point_ids(At::Y(0.5), |_| true)?, &[3, 4]);
    assert_eq!(feat.search_edge_keys(At::X(1.0), |_| true)?, &[(2, 4), (4, 7)]);

    // Perform numerical integration to compute
    // the area of cell # 2
    let ndim = 2;
    let cell_2 = &mesh.cells[2];
    let mut pad = Scratchpad::new(ndim, cell_2.kind)?;
    mesh.set_pad(&mut pad, &cell_2.points);
    let ips = integ::default_points(cell_2.kind);
    let mut area = 0.0;
    for p in 0..ips.len() {
        let iota = &ips[p];
        let weight = ips[p][3];
        let det_jac = pad.calc_jacobian(iota)?;
        area += weight * det_jac;
    }
    assert_eq!(area, 0.5);
    Ok(())
}
```



## Roadmap

- [x] Implement read/write mesh functions
- [x] Add tests for the numerical integrations
- [x] Implement triangle and tetrahedron generators
- [x] Implement drawing functions



## Appendix A - Shapes and local numbering of nodes

### Lines (Lin)

![lin_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_1_lin.svg)

### Triangles (Tri)

![tri_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_2_tri.svg)

### Quadrilaterals (Qua)

![qua_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_3_qua.svg)

### Tetrahedra (Tet)

![tet_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_4_tet.svg)

### Hexahedra (Hex)

![hex_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_5_hex.svg)



## Appendix B - Geometry versus space dimensions

The following table shows what combinations of geometry-number-of-dimensions (`geo_ndim`) and
space-number-of-dimensions (`space_ndim`) are possible. There are three cases:

1. Case `CABLE` -- `geo_ndim = 1` and `space_ndim = 2 or 3`; e.g., line in 2D or 3D (cables and rods)
2. Case `SHELL` -- `geo_ndim = 2` and `space_ndim = 3`; e.g. Tri or Qua in 3D (shells and surfaces)
3. Case `SOLID` -- `geo_ndim = space_ndim`; e.g., Tri and Qua in 2D or Tet and Hex in 3D

| `geo_ndim` | `space_ndim = 2` | `space_ndim = 3` |
| :--------: | :--------------: | :--------------: |
|     1      |     `CABLE`      |     `CABLE`      |
|     2      |     `SOLID`      |     `SHELL`      |
|     3      |    impossible    |     `SOLID`      |
