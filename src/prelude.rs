//! Makes available common structures needed to generate or use meshes
//!
//! You may write `use gemlab::prelude::*` in your code and obtain
//! access to commonly used functionality.

pub use crate::mesh::{
    check_all, display_features, draw_mesh, join_meshes, At, Block, Cell, CellAttribute, CellId, Draw, Extract,
    Feature, Features, Find, Mesh, Point, PointId, Structured, Unstructured,
};
pub use crate::shapes::{GeoClass, GeoKind};
pub use crate::util::any_x;
