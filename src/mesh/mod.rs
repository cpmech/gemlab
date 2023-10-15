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
mod check;
mod convert_2d;
mod draw_cell;
mod drawing;
mod enums;
mod features;
mod generators_qua_hex;
mod generators_tri_tet;
mod grid_cells;
mod join_meshes;
mod mesh;
mod paraview;
mod read_text_mesh;
mod samples;
mod write_text_file;
pub use crate::mesh::block::*;
pub use crate::mesh::check::*;
pub use crate::mesh::convert_2d::*;
pub use crate::mesh::draw_cell::*;
pub use crate::mesh::drawing::*;
pub use crate::mesh::enums::*;
pub use crate::mesh::features::*;
pub use crate::mesh::generators_qua_hex::*;
pub use crate::mesh::generators_tri_tet::*;
pub use crate::mesh::grid_cells::*;
pub use crate::mesh::join_meshes::*;
pub use crate::mesh::mesh::*;
pub use crate::mesh::paraview::*;
pub use crate::mesh::samples::*;
pub use crate::mesh::write_text_file::*;
