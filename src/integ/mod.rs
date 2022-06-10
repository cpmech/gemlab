//! Functions to perform numerical integration using Shapes
//!
//! # Definitions
//!
//! * `n_integ_point` -- number of integration (Gauss) points
//! * `|J|` -- the determinant of the Jacobian dx/dξ
//! * `||J||` -- the norm of the Jacobian vector for lines in multi-dimensions
//! * `ιᵖ := ξᵖ` -- (**iota**-p gets ksi-p (or xi-p)) --
//!    the coordinate of the integration point on the reference (natural) space
//! * `wᵖ` -- the weight of the p-th integration point
//!
//! # Expressions
//!
//! See also [crate::shapes]
//!
//! ## Xᵀ: Transposed matrix of coordinates
//!
//! ```text
//!      ┌                              ┐  superscript = node
//!      | x⁰₀  x¹₀  x²₀  x³₀       xᴹ₀ |  subscript = space dimension
//! Xᵀ = | x⁰₁  x¹₁  x²₁  x³₁  ...  xᴹ₁ |
//!      | x⁰₂  x¹₂  x²₂  x³₂       xᴹ₂ |
//!      └                              ┘_(space_ndim,nnode)
//! where `M = nnode - 1`
//! ```
//!
//! ## N: Shape functions
//!
//! Shape function of node m at ξ (ksi; i.e., xi):
//!
//! ```text
//! Nᵐ(ξ)
//! matrix notation: N is an (nnode) vector
//! ```
//!
//! ## L: Derivatives of shape functions
//!
//! Derivatives of shape functions with respect to the natural coordinates:
//!
//! ```text
//!             →
//! →  →    dNᵐ(ξ)
//! Lᵐ(ξ) = ——————
//!            →
//!           dξ
//! matrix notation: L is an (nnode,geo_ndim) matrix
//! ```
//!
//! ## Isoparametric formula
//!
//! ```text
//! → →         →  →
//! x(ξ) = Σ Nᵐ(ξ) xᵐ
//!        m
//! matrix notation: x = Xᵀ ⋅ N
//! ```
//!
//! ## J: Jacobian tensor (geo_ndim = space_ndim = 2 or 3)
//!
//! ```text
//!         →
//!   →    dx     →    →
//! J(ξ) = —— = Σ xᵐ ⊗ Lᵐ
//!         →   m
//!        dξ
//! matrix notation: J = Xᵀ · L
//! J is a (space_ndim,geo_ndim) matrix
//! ```
//!
//! ## J: Jacobian vector (geo_ndim = 1 and space_ndim = 2 or 3)
//!
//! ```text
//!                          →
//! →    →     →    →  →    dx
//! J := Jline(ξ) = g₁(ξ) = ——
//!                         dξ
//! matrix notation: J = Jline = Xᵀ · L
//! J is a (space_ndim,1) matrix; i.e., a vector
//! ```
//!
//! ## G: Gradients of shape functions (derivatives w.r.t real coordinates x)
//!
//! ```text
//!             →
//! →  →    dNᵐ(ξ)
//! Gᵐ(ξ) = ——————
//!            →
//!           dx
//! matrix notation: G = L · J⁻¹
//! G is an (nnode,space_ndim) matrix
//! ```

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
