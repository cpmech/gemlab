//! Gemlab -- Geometry And Meshes laboratory for Finite Element Analyses

/// Defines a type alias for the error type as a static string
pub type StrError = &'static str;

pub mod fem;
pub mod geometry;
pub mod integration;
pub mod mesh;
pub mod shapes;
pub mod util;
