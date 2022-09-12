//! Makes available common structures needed to generate or use meshes
//!
//! You may write `use gemlab::prelude::*` in your code and obtain
//! access to commonly used functionality.

pub use crate::mesh::{
    check_all, draw_mesh, join_meshes, At, Block, Cell, CellAttributeId, CellId, Draw, Extract, Feature, Features,
    Find, Mesh, Point, PointId, Structured, Unstructured,
};
pub use crate::shapes::{GeoClass, GeoKind};
pub use crate::util::any;
