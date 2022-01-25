//! Implements a simple linear finite element solver

mod attribute;
mod boundary_condition;
mod element;
mod element_solid;
mod equation_numbers;
mod simulation;
pub use crate::fem::attribute::*;
pub use crate::fem::boundary_condition::*;
pub use crate::fem::element::*;
pub use crate::fem::element_solid::*;
pub use crate::fem::equation_numbers::*;
pub use crate::fem::simulation::*;
