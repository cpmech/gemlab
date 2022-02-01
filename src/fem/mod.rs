//! Implements a simple linear finite element solver

mod boundary_condition;
mod element;
mod element_porous;
mod element_solid;
mod equation_numbers;
mod model_seepage;
mod model_seepage_standard;
mod model_solid;
mod model_solid_linear_elastic;
// mod simulation;
pub use crate::fem::boundary_condition::*;
pub use crate::fem::element::*;
pub use crate::fem::element_porous::*;
pub use crate::fem::element_solid::*;
pub use crate::fem::equation_numbers::*;
pub use crate::fem::model_seepage::*;
pub use crate::fem::model_seepage_standard::*;
pub use crate::fem::model_solid::*;
pub use crate::fem::model_solid_linear_elastic::*;
// pub use crate::fem::simulation::*;
