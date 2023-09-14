//! Mesh definitions and generation, including tools to find features
//!
//! This module defines the [Mesh] structure and auxiliary functions for mesh generation,
//! finding features such as [Point], edges, and faces, and other algorithms such as
//! merging meshes and drawing.
//!
//! A [Mesh] is composed of [Point]s and [Cell]s whereas the secondary features are edges
//! and faces. The structure [Features] holds the (secondary) features.
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
mod check;
mod draw;
mod features;
mod find;
mod generators;
mod generators_tri;
mod grid_cells;
mod helpers;
mod join_meshes;
mod mesh;
mod paraview;
mod read_text_mesh;
mod samples;
mod upgrade_mesh_2d;
mod write_text_file;
pub use crate::mesh::block::*;
pub use crate::mesh::check::*;
pub use crate::mesh::draw::*;
pub use crate::mesh::features::*;
pub use crate::mesh::find::*;
pub use crate::mesh::generators::*;
pub use crate::mesh::generators_tri::*;
pub use crate::mesh::grid_cells::*;
pub use crate::mesh::helpers::*;
pub use crate::mesh::join_meshes::*;
pub use crate::mesh::mesh::*;
pub use crate::mesh::paraview::*;
pub use crate::mesh::samples::*;
pub use crate::mesh::upgrade_mesh_2d::*;
pub use crate::mesh::write_text_file::*;
