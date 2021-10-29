//! Defines shapes and interpolation functions
//! # Example
//!
//! ```
//! // import
//! use gemlab::shape;
//! use russell_lab::Vector;
//!
//! // shape object
//! let mut shape = shape::new(shape::Kind::Hex8);
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

mod hex8;
mod qua4;
mod shape;
pub use crate::shape::hex8::*;
pub use crate::shape::qua4::*;
pub use crate::shape::shape::*;
