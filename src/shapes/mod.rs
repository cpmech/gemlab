//! Interpolation functions and derivatives for geometric shapes
//!
//! # Definitions
//!
//! Here, we consider the following dimensions:
//!
//! * `space_ndim` -- is the number of dimensions of the space under study (2 or 3)
//! * `geo_ndim` -- is the number of dimensions of the geometry element (shape),
//!                  for instance, a line in the 2D space has `geo_ndim = 1` and
//!                  `space_ndim = 2`. Another example is a triangle in the 3D space
//!                  which has `geo_ndim = 2` and `space_ndim = 3`.
//! * `local` refers to a numbering scheme for the nodes of the shape (or element)
//! * `global` refers to a numbering scheme applied for the whole mesh
//!
//! We also consider the following counting variables:
//!
//! * `nnode` -- (local) number of points that define the shape (number of nodes)
//! * `npoint` -- (global) number of points in the whole mesh; not used in this module but important to remember
//! * `nedge` -- number of edges
//! * `nface` -- number of faces
//! * `edge_nnode` -- number of points that define the edge (number of nodes)
//! * `face_nnode` -- number of points that define the face (number of nodes)
//! * `face_nedge` -- face's number of edges
//! * `n_integ_point` -- number of integration (Gauss) points
//!
//! When performing numerical integrations, we use the following notation:
//! |J| is the determinant of the Jacobian, ||J|| is the norm of the Jacobian vector
//! for line in multi-dimensions, `n_integ_point` is the number of integration points,
//! `ιp := ξp` is the reference coordinate of the integration point,
//! and `wp` is the weight of the p-th integration point.
//!
//! # Isoparametric formulation
//!
//! The isoparametric formulation establishes that
//!
//! ```text
//! → →         →  →
//! x(ξ) = Σ Nᵐ(ξ) xᵐ
//!        m         
//! ```
//!
//! where `x` is the (space_ndim) vector of real coordinates, `ξ` is the (geo_ndim)
//! vector of reference coordinates, `Nm` are the (nnode) interpolation functions,
//! and `xm` are the (nnode) coordinates of each m-node of the geometric shape.
//!
//! Given an (nnode,space_ndim) matrix of coordinates X, we can calculate the
//! (space_ndim) vector of coordinates x by means of
//!
//! ```text
//! x = Xᵀ ⋅ N
//! ```
//!
//! where `N` is an (nnode) array formed with all `Nm`.
//!
//! # General case (geo_ndim == space_ndim)
//!
//! If `geo_ndim == space_ndim`, we define the Jacobian tensor as
//!
//! ```text
//!         →
//!   →    dx     →    →
//! J(ξ) = —— = Σ xᵐ ⊗ Lᵐ
//!         →   m
//!        dξ
//! ```
//!
//! where
//!
//! ```text
//!             →
//! →  →    dNᵐ(ξ)
//! Lᵐ(ξ) = ——————
//!            →
//!           dξ
//! ```
//!
//! are the derivatives of each interpolation function `Nm` with respect to the
//! reference coordinate. `Lm` are (geo_ndim) vectors and can be organized in
//! an (nnode,geo_ndim) matrix `L` of "local" derivatives.
//!
//! We can write the Jacobian in matrix notation as follows
//!
//! ```text
//! J = Xᵀ · L
//! ```
//!
//! where X is the (nnode,space_ndim) matrix of coordinates and L is the (nnode,geo_ndim) matrix.
//!
//! We define the gradient of interpolation functions (i.e., derivatives of interpolation
//! functions w.r.t real coordinates) by means of
//!
//! ```text
//!             →
//! →  →    dNᵐ(ξ)
//! Gᵐ(ξ) = ——————
//!            →
//!           dx
//! ```
//!
//! which can be organized in an (nnode,space_ndim) matrix `G`.
//!
//! The inverse Jacobian allows us to determine the gradient vectors G as follows
//!
//! ```text
//! →       →  →        →
//! Gᵐ(ξ) = Lᵐ(ξ) · J⁻¹(ξ)
//! ```
//!
//! Or, in matrix notation,
//!
//! ```text
//! G = L · J⁻¹
//! ```
//!
//! where G is an (nnode,space_ndim) matrix.
//!
//! # Line in multi-dimensions (geo_ndim == 1 and space_ndim > 1)
//!
//! In this case, the Jacobian equals the (space_ndim,1) base vector `g1` tangent
//! to the line element, i.e.,
//!
//! ```text
//!                          →
//! →    →     →            dx
//! J := Jline(ξ) = g₁(ξ) = —— = Xᵀ · L
//!                         dξ
//! ```
//!
//! We also consider a parametric coordinate `ℓ` which varies
//! from 0 to `ℓ_max` (the length of the line) according to
//!
//! ```text
//!                ℓ_max
//! ℓ(ξ) = (1 + ξ) —————
//!                  2
//!
//!        2 · ℓ
//! ξ(ℓ) = ————— - 1
//!        ℓ_max
//! ```
//!
//! ```text
//! 0 ≤ ℓ ≤ ℓ_max
//!
//! -1 ≤ ξ ≤ +1
//! ```

mod enums;
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
mod state_of_shape;
mod tet10;
mod tet4;
mod tri10;
mod tri15;
mod tri3;
mod tri6;
pub use crate::shapes::enums::*;
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
pub use crate::shapes::state_of_shape::*;
pub use crate::shapes::tet10::*;
pub use crate::shapes::tet4::*;
pub use crate::shapes::tri10::*;
pub use crate::shapes::tri15::*;
pub use crate::shapes::tri3::*;
pub use crate::shapes::tri6::*;
