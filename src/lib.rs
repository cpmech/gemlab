//! Gemlab

pub type StrError = &'static str;

mod as_array;
mod constants;
mod fem;
mod geometry;
mod grid_search;
mod mesh;
mod shapes;
mod structured;
pub use crate::as_array::*;
pub use crate::constants::*;
pub use crate::fem::*;
pub use crate::geometry::*;
pub use crate::grid_search::*;
pub use crate::mesh::*;
pub use crate::shapes::*;
pub use crate::structured::*;
