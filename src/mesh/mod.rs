//! Mesh definitions, read/write, boundary, find, and generators for FEA
//!
//! It is recommended to look at [Region] first, since it is the struct
//! of highest level in this module.
//!
//! Then, it is worth looking at [Mesh] which contains the raw mesh data
//! and [Features] which contains post-processed information given mesh data.
//!

mod algorithms;
mod block;
mod check;
mod draw;
mod features;
mod find;
mod generators;
mod generators_tri;
mod helpers;
mod join_meshes;
mod mesh;
mod read_text_mesh;
mod region;
mod samples;
pub use crate::mesh::block::*;
pub use crate::mesh::check::*;
pub use crate::mesh::draw::*;
pub use crate::mesh::features::*;
pub use crate::mesh::find::*;
pub use crate::mesh::generators::*;
pub use crate::mesh::generators_tri::*;
pub use crate::mesh::helpers::*;
pub use crate::mesh::join_meshes::*;
pub use crate::mesh::mesh::*;
pub use crate::mesh::region::*;
pub use crate::mesh::samples::*;
