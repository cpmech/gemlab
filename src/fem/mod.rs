//! Implements a simple linear finite element solver

mod definitions;
mod element;
mod element_solid;
mod simulation;
pub use crate::fem::definitions::*;
pub use crate::fem::element::*;
pub use crate::fem::element_solid::*;
pub use crate::fem::simulation::*;
