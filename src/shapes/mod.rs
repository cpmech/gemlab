//! Defines shapes and interpolation functions
//!
//! # Example
//!
//! ```
//! // import
//! use gemlab::shapes::*;
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

mod hex20;
mod hex8;
mod lin2;
mod lin3;
mod lin4;
mod lin5;
mod qua12;
mod qua16;
mod qua17;
mod qua4;
mod qua8;
mod qua9;
mod shape;
pub use crate::shapes::hex20::*;
pub use crate::shapes::hex8::*;
pub use crate::shapes::lin2::*;
pub use crate::shapes::lin3::*;
pub use crate::shapes::lin4::*;
pub use crate::shapes::lin5::*;
pub use crate::shapes::qua12::*;
pub use crate::shapes::qua16::*;
pub use crate::shapes::qua17::*;
pub use crate::shapes::qua4::*;
pub use crate::shapes::qua8::*;
pub use crate::shapes::qua9::*;
pub use crate::shapes::shape::*;
