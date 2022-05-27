//! Functions to perform numerical integration using Shapes
//!
//! # Definitions
//!
//! * `n_integ_point` -- number of integration (Gauss) points
//!
//! When performing numerical integrations, we use the following notation:
//! |J| is the determinant of the Jacobian, ||J|| is the norm of the Jacobian vector
//! for line in multi-dimensions, `n_integ_point` is the number of integration points,
//! `ιp := ξp` is the reference coordinate of the integration point,
//! and `wp` is the weight of the p-th integration point.

mod analytical_tet4;
mod analytical_tri3;
mod integ_points;
mod matrix_cases;
mod vector_cases;
pub use crate::integ::analytical_tet4::*;
pub use crate::integ::analytical_tri3::*;
pub use crate::integ::integ_points::*;
pub use crate::integ::matrix_cases::*;
pub use crate::integ::vector_cases::*;
