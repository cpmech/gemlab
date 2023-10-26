//! Internal algorithms for use in the mesh module

mod extract_all_2d_edges;
mod extract_all_faces;
mod extract_features_2d;
mod extract_features_3d;
mod graph;
pub(crate) use crate::mesh::algorithms::extract_all_2d_edges::*;
pub(crate) use crate::mesh::algorithms::extract_all_faces::*;
pub(crate) use crate::mesh::algorithms::extract_features_2d::*;
pub(crate) use crate::mesh::algorithms::extract_features_3d::*;
pub use crate::mesh::algorithms::graph::*;
