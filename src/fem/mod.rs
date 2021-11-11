//! Implements a simple linear finite element solver

mod attribute;
mod boundary_condition;
mod degree_of_freedom;
mod element;
mod element_solid;
mod simulation;
pub use crate::fem::attribute::*;
pub use crate::fem::boundary_condition::*;
pub use crate::fem::degree_of_freedom::*;
pub use crate::fem::element::*;
pub use crate::fem::element_solid::*;
pub use crate::fem::simulation::*;
