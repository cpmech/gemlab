//! Gemlab -- Geometry And Meshes laboratory for Finite Element Analyses

/// Defines a type alias for the error type as a static string
pub type StrError = &'static str;

pub mod geometry;
pub mod integ;
pub mod mesh;
pub mod shapes;
pub mod util;
