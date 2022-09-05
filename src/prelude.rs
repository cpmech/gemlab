//! Makes available common structures needed to generate or use meshes
//!
//! You may write `use gemlab::prelude::*` in your code and obtain
//! access to commonly used functionality.

pub use crate::mesh::{
    check_all, draw_mesh, Block, Cell, CellAttributeId, CellId, Draw, Feature, Find, Mesh, Point, PointId,
};
pub use crate::shapes::{GeoClass, GeoKind};
pub use crate::StrError;
