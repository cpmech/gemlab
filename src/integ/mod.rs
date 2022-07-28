//! Functions to perform numerical integration using Shapes
//!
//! # Definitions
//!
//! * `n_integ_point` -- number of integration (Gauss) points
//! * `|J|` -- the determinant of the Jacobian dx/dξ
//! * `||J||` -- the norm of the Jacobian vector for lines in multi-dimensions
//! * `ξ` -- ksi (or xi) -- coordinates in the reference space
//! * `ιᵖ := ξᵖ` -- (**iota**-p gets ksi-p (or xi-p)) --
//!    the coordinate of the integration point on the reference (natural) space
//! * `wᵖ` -- the weight of the p-th integration point
//! * Xᵀ -- Transposed matrix of coordinates
//! * N -- Shape functions
//! * L -- Derivatives of shape functions
//! * J -- Jacobian tensor
//! * J⁻¹ -- Inverse Jacobian matrix (only available if geo_ndim = space_ndim)
//! * G -- Gradients of shape functions (derivatives w.r.t real coordinates x)
//! * see also the subsection named **Expressions** below
//!
//! # Functionality
//!
//! ## Integration of scalar field over a geometric shape
//!
//! Function [scalar_field()]
//!
//! ```text
//!     ⌠   → →
//! I = │ s(x(ξ)) dΩ
//!     ⌡
//!     Ωₑ
//! ```
//!
//! ## Integration of some combinations involving N and G resulting in vectors
//!
//! Interpolation functions times scalar field [vec_a()]:
//!
//! ```text
//!      ⌠    → →     →
//! aᵐ = │ Nᵐ(x(ξ)) s(x) dΩ
//!      ⌡
//!      Ωₑ
//! ```
//!
//! Interpolation functions times vector field [vec_b()]:
//!
//! ```text
//! →    ⌠    → →   → →
//! bᵐ = │ Nᵐ(x(ξ)) v(x) dΩ
//!      ⌡
//!      Ωₑ
//! ```
//!
//! Vector dot gradient [vec_c()]:
//!
//! ```text
//!      ⌠ → →    →  → →
//! cᵐ = │ w(x) · Gᵐ(x(ξ)) dΩ
//!      ⌡
//!      Ωₑ
//! ```
//!
//! Tensor dot gradient [vec_d()]:
//!
//! ```text
//! →    ⌠   →    →  → →
//! dᵐ = │ σ(x) · Gᵐ(x(ξ)) dΩ
//!      ⌡ ▔
//!      Ωₑ
//! ```
//!
//! ## Integration of some combinations involving N, tensors, and G, resulting in matrices
//!
//! Gradient(G) dot 4th-tensor(D) dot gradient(G) integration (stiffness matrix) [mat_gdg()]:
//!
//! ```text
//!       ⌠               →    →
//! Kᵐⁿ = │ Gᵐₖ Dᵢₖⱼₗ Gⁿₗ eᵢ ⊗ eⱼ dΩ
//! ▔     ⌡
//!       Ωₑ
//! ```
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
//! ## J: Jacobian tensor (SOLID case with geo_ndim = space_ndim = 2 or 3)
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
//! ## J: Jacobian vector (CABLE case with geo_ndim = 1 and space_ndim = 2 or 3)
//!
//! ```text
//!                           →
//! →    →      →    →  →    dx
//! J := Jcable(ξ) = g₁(ξ) = ——
//!                          dξ
//! matrix notation: J = Jcable = Xᵀ · L
//! J is a (space_ndim,1) matrix; i.e., a vector
//! ```
//!
//! ## J: Jacobian matrix (SHELL case with geo_ndim = 2 and space_ndim = 3)
//!
//! ```text
//!                 dx
//! J(ξ) = Jshell = ——
//!                 dξ
//! matrix notation: J = Jshell = Xᵀ · L
//! J is a (3,2) matrix
//! ```
//!
//! ## G: Gradients of shape functions (derivatives w.r.t real coordinates x) (only for SOLID case)
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

mod analytical_qua4;
mod analytical_qua8;
mod analytical_tet4;
mod analytical_tri3;
mod integ_points;
mod matrix_cases;
mod matrix_cases_coupled;
mod matrix_cases_nonsym;
mod point_coords;
mod scalar_field;
mod testing;
mod vector_cases;
mod vector_cases_boundary;
pub use crate::integ::analytical_qua4::*;
pub use crate::integ::analytical_qua8::*;
pub use crate::integ::analytical_tet4::*;
pub use crate::integ::analytical_tri3::*;
pub use crate::integ::integ_points::*;
pub use crate::integ::matrix_cases::*;
pub use crate::integ::matrix_cases_coupled::*;
pub use crate::integ::matrix_cases_nonsym::*;
pub use crate::integ::point_coords::*;
pub use crate::integ::scalar_field::*;
pub use crate::integ::vector_cases::*;
pub use crate::integ::vector_cases_boundary::*;
