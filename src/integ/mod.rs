//! Functions to perform numerical integration using Shapes
//!
//! # Definitions
//!
//! * `n_integ_point` -- number of integration (Gauss) points
//! * `|J|` --is the determinant of the Jacobian
//! * `||J||` -- is the norm of the Jacobian vector for lines in multi-dimensions
//! * `ιp := ξp` -- iota-p gets ksi-p (or xi-p) -- is the reference coordinate of the integration point
//! * `wp` -- is the weight of the p-th integration point

mod analytical_tet4;
mod analytical_tri3;
mod calc_ips_coords;
mod integ_points;
mod matrix_cases;
mod scalar_field;
mod vector_cases;
pub use crate::integ::analytical_tet4::*;
pub use crate::integ::analytical_tri3::*;
pub use crate::integ::calc_ips_coords::*;
pub use crate::integ::integ_points::*;
pub use crate::integ::matrix_cases::*;
pub use crate::integ::scalar_field::*;
pub use crate::integ::vector_cases::*;
