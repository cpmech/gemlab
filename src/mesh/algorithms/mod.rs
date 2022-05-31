//! Internal algorithms for use in the mesh module

mod extract_all_edges_2d;
mod extract_all_faces_3d;
mod extract_features_2d;
mod extract_features_3d;
mod extract_rods_and_shells;
pub(crate) use crate::mesh::algorithms::extract_all_edges_2d::*;
pub(crate) use crate::mesh::algorithms::extract_all_faces_3d::*;
pub(crate) use crate::mesh::algorithms::extract_features_2d::*;
pub(crate) use crate::mesh::algorithms::extract_features_3d::*;
pub(crate) use crate::mesh::algorithms::extract_rods_and_shells::*;
