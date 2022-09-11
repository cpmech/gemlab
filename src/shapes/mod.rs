//! Interpolation functions and derivatives for geometric shapes (elements)
//!
//! # Definitions
//!
//! Here, we consider the following definitions:
//!
//! * `space_ndim` -- is the number of dimensions of the space under study (2 or 3)
//! * `geo_ndim` -- is the number of dimensions of the geometry element (shape).
//!                  For instance, a line in the 2D space has `geo_ndim = 1` and
//!                  `space_ndim = 2`. Another example: a triangle in the 3D space
//!                  has `geo_ndim = 2` and `space_ndim = 3`.
//! * `local` -- refers to a numbering scheme for the nodes of the shape (or element)
//! * `global` -- refers to a numbering scheme applied for the whole mesh
//! * `spatial (real) space` -- is the "real" space mapped by the x₀,x₁,x₂ coordinates (see figure below)
//! * `reference (natural) space` -- is the "virtual" space mapped by the ξ₀,ξ₁,ξ₂ coordinates (see figure below)
//!
//! **Note:** In the code, we use (r,s,t) as (ξ₀,ξ₁,ξ₂) and `ksi` to represent the vector `ξ`.
//!
//! ![Real and reference space mapping](https://github.com/cpmech/gemlab/raw/main/data/figures/mapping-real-to-reference.png)
//!
//! We also consider the following counting variables:
//!
//! * `nnode` -- (local) number of points (aka nodes) that define the shape/element.
//! * `npoint` -- (global) number of points in the whole mesh; not used in this module
//!               but important to remember
//! * `nedge` -- number of edges on the shape (2D or 3D)
//! * `nface` -- number of faces on the shape (3D only)
//! * `edge_nnode` -- number of points/nodes that define the edge
//! * `face_nnode` -- number of points/nodes that define the face
//! * `face_nedge` -- number of edges on the face
//!
//! Geometry cases regarding the number of dimensions (geo vs space)
//!
//! 1. Case `CABLE` -- `geo_ndim = 1` and `space_ndim = 2 or 3`; e.g., line in 2D or 3D (cables and rods)
//! 2. Case `SHELL` -- `geo_ndim = 2` and `space_ndim = 3`; e.g. Tri or Qua in 3D (shells and surfaces)
//! 3. Case `SOLID` -- `geo_ndim = space_ndim`; e.g., Tri and Qua in 2D or Tet and Hex in 3D
//!
//! | `geo_ndim` | `space_ndim = 2` | `space_ndim = 3` |
//! |:----------:|:----------------:|:----------------:|
//! |     1      |     `CABLE`      |     `CABLE`      |
//! |     2      |     `SOLID`      |     `SHELL`      |
//! |     3      |    impossible    |     `SOLID`      |
//!
//! # Isoparametric formulation
//!
//! The isoparametric formulation establishes that we can calculate the coordinates `x(ξ)`
//! within the shape/element from the shape functions `Nᵐ(ξ)` and the coordinates at each
//! node by using the formula:
//!
//! ```text
//! → →         →  →
//! x(ξ) = Σ Nᵐ(ξ) xᵐ
//!        m
//! ```
//!
//! where `x` is the (space_ndim) vector of real coordinates, `ξ` is the (geo_ndim)
//! vector of reference coordinates, `Nᵐ` are the (nnode) interpolation functions,
//! and `xᵐ` are the (nnode) coordinates of each m-node of the geometric shape.
//!
//! Given an (nnode,space_ndim) **matrix** of coordinates `X`, we can calculate the
//! (space_ndim) **vector** of coordinates `x` by means of
//!
//! ```text
//! x = Xᵀ ⋅ N
//! ```
//!
//! where `N` is an (nnode) **vector** formed with all `Nᵐ` values.
//!
//! # Derivatives on the reference space and gradients on the real space
//!
//! Here, we consider two cases:
//!
//! * General case with geo_ndim = space_ndim; and
//! * Line in multi-dimensions with geo_ndim = 1 and space_ndim > 1.
//!
//! ## SOLID case with geo_ndim = space_ndim
//!
//! If `SOLID` (`geo_ndim = space_ndim = 2 or 3`), we define the Jacobian tensor as
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
//! are the derivatives of each interpolation function `Nᵐ` with respect to the
//! reference coordinate. `Lᵐ` are (geo_ndim) vectors and can be organized in
//! an (nnode,geo_ndim) matrix `L` of **local** derivatives.
//!
//! We can write the Jacobian in matrix (space_ndim,geo_ndim) notation as follows
//!
//! ```text
//! J = Xᵀ · L
//! ```
//!
//! where `X` is the (nnode,space_ndim) matrix of coordinates and `L` is the (nnode,geo_ndim) matrix.
//!
//! Next, we define the gradient of interpolation functions (i.e., derivatives of interpolation
//! functions with respect to real coordinates) by means of
//!
//! ```text
//!             →
//! →  →    dNᵐ(ξ)
//! Bᵐ(ξ) = ——————
//!            →
//!           dx
//! ```
//!
//! which can be organized in an (nnode,space_ndim) matrix `G`.
//!
//! The inverse Jacobian allows us to determine the gradient vectors `G` as follows
//!
//! ```text
//! →       →  →        →
//! Bᵐ(ξ) = Lᵐ(ξ) · J⁻¹(ξ)
//! ```
//!
//! Or, in matrix notation,
//!
//! ```text
//! G = L · J⁻¹
//! ```
//!
//! where `G` is an (nnode,space_ndim) matrix.
//!
//! ## SHELL case with geo_ndim = 2 and space_ndim = 3
//!
//! In this case, the Jacobian matrix is (3,2) and can also be computed by the following matrix
//! multiplication
//!
//! ```text
//!        dx
//! J(ξ) = ——
//!        dξ
//! ```
//!
//! Or, in matrix notation,
//!
//! ```text
//! J = Jshell = Xᵀ · L
//! ```
//!
//! However, the inverse Jacobian and gradients are not available in this case.
//!
//! ## CABLE case with geo_ndim = 1 and space_ndim = 2 or 3
//!
//! In this case, the Jacobian equals the (space_ndim,1) base vector `g₁` which
//! is tangent to the line element, i.e.,
//!
//! ```text
//!                           →
//! →    →      →    →  →    dx
//! J := Jcable(ξ) = g₁(ξ) = ——
//!                          dξ
//! matrix notation: Jcable = Xᵀ · L
//! ```
//!
//! We also consider a parametric coordinate `ℓ` which varies
//! from `0` to `ℓ_max` (the length of the line) according to
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
//! Note that:
//!
//! ```text
//! 0 ≤ ℓ ≤ ℓ_max
//!
//! -1 ≤ ξ ≤ +1
//! ```
//!
//! # Normal vectors
//!
//! ## CABLE case with geo_ndim = 1 and space_ndim = 2 or 3
//!
//! Base vector tangent to the line:
//!
//! ```text
//!          →
//!         dx
//! g₁(ξ) = —— = Xᵀ · L = first_column(J)
//!         dξ
//! ```
//!
//! Normal vector:
//!
//! ```text
//! →   →    →
//! n = e₃ × g₁
//!
//!   →       →
//! ||n|| = ||g₁||
//! ```
//!
//! Thus
//!
//! ```text
//!        →           →
//! dℓ = ||g₁|| dξ = ||n|| dξ
//! ```
//!
//! For a straight line (segment):
//!
//! ```text
//!   →
//! ||n|| = Δℓ / Δξ = L / 2
//! ```
//!
//! because all [GeoClass::Lin] have `Δξ = 2`.
//!
//! ## SHELL case with geo_ndim = 2 and space_ndim = 3
//!
//! Base vectors tangent to the surface:
//!
//! ```text
//!          →
//! →  →    dx
//! g₁(ξ) = ——— = first_column(J)
//!         dξ₁
//!
//!          →
//! →  →    dx
//! g₂(ξ) = ——— = second_column(J)
//!         dξ₂
//! ```
//!
//! Normal vector:
//!
//! ```text
//! →   →    →
//! n = g₁ × g₂
//! ```
//!
//! Thus
//!
//! ```text
//!         →
//! dA := ||n|| dξ₁ dξ₂
//! ```
//!
//! For flat quadrilateral faces with sides perpendicular one with another
//!
//! ```text
//!   →
//! ||n|| = A / (Δξ₁ Δξ₂) = A / 4
//! ```
//!
//! because all [GeoClass::Qua] have `Δξᵢ = 2`.
//!
//! # GeoKind and GeoClass
//!
//! [GeoKind] is perhaps the most important *enum* in this module.
//! It defines the actual *shape* used in finite element analyses.
//!
//! Below are some example of shapes, classified according to [GeoClass].
//! The numbers are the local numbers (nodes).
//!
//! # Lines -- Lin
//!
//! ![lin_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_1_lin.svg)
//!
//! # Triangles -- Tri
//!
//! ![tri_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_2_tri.svg)
//!
//! # Quadrilaterals -- Qua
//!
//! ![qua_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_3_qua.svg)
//!
//! # Tetrahedra -- Tet
//!
//! ![tet_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_4_tet.svg)
//!
//! # Hexahedra -- Hex
//!
//! ![hex_cells](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_cells_and_points_work_5_hex.svg)
//!

mod enums;
mod hex20;
mod hex32;
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
mod scratchpad;
mod scratchpad_approximate_ksi;
mod scratchpad_calc_coords;
mod scratchpad_calc_gradient;
mod scratchpad_calc_jacobian;
mod scratchpad_calc_normal_vector;
mod scratchpad_draw_shape;
mod scratchpad_testing;
mod tet10;
mod tet20;
mod tet4;
mod tri10;
mod tri15;
mod tri3;
mod tri6;
pub use crate::shapes::enums::*;
pub use crate::shapes::hex20::*;
pub use crate::shapes::hex32::*;
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
pub use crate::shapes::scratchpad::*;
pub use crate::shapes::scratchpad_approximate_ksi::*;
pub use crate::shapes::scratchpad_calc_coords::*;
pub use crate::shapes::scratchpad_calc_gradient::*;
pub use crate::shapes::scratchpad_calc_jacobian::*;
pub use crate::shapes::scratchpad_calc_normal_vector::*;
pub use crate::shapes::scratchpad_draw_shape::*;
pub use crate::shapes::tet10::*;
pub use crate::shapes::tet20::*;
pub use crate::shapes::tet4::*;
pub use crate::shapes::tri10::*;
pub use crate::shapes::tri15::*;
pub use crate::shapes::tri3::*;
pub use crate::shapes::tri6::*;
