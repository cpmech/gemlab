//! Defines shapes and interpolation functions
//!
//! # Example
//!
//! ```
//! // import
//! use gemlab::*;
//! use russell_lab::Vector;
//!
//! // shape object
//! let mut shape = new_shape(Kind::Hex8);
//! let ndim = shape.get_ndim();
//!
//! // compute interp fn and deriv @ ksi=(0,0,0) for all points
//! let ksi = Vector::new(ndim);
//! shape.calc_interp(&ksi);
//! shape.calc_deriv(&ksi);
//!
//! // get interp fn and deriv for point m
//! let m = 0;
//! let sm = shape.get_interp(m);
//! let mut dsm_dksi = Vector::new(ndim);
//! shape.get_deriv(&mut dsm_dksi, m);
//!
//! // check interp fn and deriv @ ksi for point m
//! assert_eq!(sm, 0.125);
//! assert_eq!(
//!     format!("{}", dsm_dksi),
//!     "┌        ┐\n\
//!      │ -0.125 │\n\
//!      │ -0.125 │\n\
//!      │ -0.125 │\n\
//!      └        ┘"
//! );
//! ```

mod enums;
mod hex20;
mod hex8;
mod qua4;
mod qua8;
mod shape;
pub use crate::shapes::enums::*;
pub use crate::shapes::hex20::*;
pub use crate::shapes::hex8::*;
pub use crate::shapes::qua4::*;
pub use crate::shapes::qua8::*;
pub use crate::shapes::shape::*;
