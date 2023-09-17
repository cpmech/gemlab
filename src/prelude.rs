//! Makes available common structures needed to generate or use meshes
//!
//! You may write `use gemlab::prelude::*` in your code and obtain
//! access to commonly used functionality.

pub use crate::mesh::{
    join_meshes, At, Block, Cell, CellAttribute, CellId, Extract, Feature, Features, Figure, Find, Mesh, Point,
    PointId, Structured, Unstructured,
};
pub use crate::shapes::{GeoClass, GeoKind, Scratchpad};
pub use crate::util::any_x;
