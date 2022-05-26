//! Gemlab -- Geometry, meshes, and integration for finite element analyses

/// Defines a type alias for the error type as a static string
pub type StrError = &'static str;

pub mod geometry;
pub mod integ;
pub mod mesh;
pub mod shapes;
pub mod util;
