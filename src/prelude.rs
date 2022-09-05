//! Makes available common structures needed to generate or use meshes
//!
//! You may write `use gemlab::prelude::*` in your code and obtain
//! access to commonly used functionality.

pub use crate::mesh::{
    check_all, draw_mesh, At, Block, Cell, CellAttributeId, CellId, Draw, Extract, Feature, Find, Mesh, Point, PointId,
};
pub use crate::shapes::{GeoClass, GeoKind};
