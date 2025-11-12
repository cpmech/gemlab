//! Mesh definitions and generation, including tools to find features
//!
//! This module defines the [Mesh] structure and auxiliary functions for mesh generation,
//! search features such as [Point], edges, and faces, and other algorithms such as
//! merging meshes and drawing.
//!
//! A [Mesh] is defined by [Point]s and [Cell]s with the associated [Features] being the edges
//! and faces.
//!
//! Below are some example of [Cell]s, classified according to [super::shapes::GeoClass].
//! The numbers are the local numbers of the cell points (nodes).
//!
//! # Linear cells -- Lin
//!
//! ![lin_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_1_lin.svg)
//!
//! # Triangles -- Tri
//!
//! ![tri_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_2_tri.svg)
//!
//! # Quadrilaterals -- Qua
//!
//! ![qua_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_3_qua.svg)
//!
//! # Tetrahedra -- Tet
//!
//! ![tet_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_4_tet.svg)
//!
//! # Hexahedra -- Hex
//!
//! ![hex_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_5_hex.svg)
//!

mod algorithms;
mod block;
mod blocks2d;
mod check;
mod convert_2d;
mod draw_cell;
mod edges;
mod enums;
mod features;
mod figure;
mod generators_qua_hex;
mod generators_tri_tet;
mod grid_cells;
mod join_meshes;
mod mesh;
mod paraview;
mod read_text_mesh;
mod samples;
mod write_text_file;

pub use block::*;
pub use blocks2d::*;
pub use enums::*;
pub use features::*;
pub use figure::*;
pub use generators_qua_hex::*;
pub use generators_tri_tet::*;
pub use grid_cells::*;
pub use join_meshes::*;
pub use mesh::*;
pub use samples::*;

// re-export GeoKind
pub use crate::shapes::GeoCase;
pub use crate::shapes::GeoClass;
pub use crate::shapes::GeoKind;
