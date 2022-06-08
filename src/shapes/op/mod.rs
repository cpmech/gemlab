//! Functions (operators) to perform numerical integration and more

mod approximate_ksi;
mod calc_coords;
mod calc_gradient;
mod calc_jacobian;
mod calc_normal_vector;
mod draw_shape;
mod testing;
pub use crate::shapes::op::approximate_ksi::*;
pub use crate::shapes::op::calc_coords::*;
pub use crate::shapes::op::calc_gradient::*;
pub use crate::shapes::op::calc_jacobian::*;
pub use crate::shapes::op::calc_normal_vector::*;
pub use crate::shapes::op::draw_shape::*;
