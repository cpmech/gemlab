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
//! # Integration of scalar field over a geometric shape
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
//! # Vector results: Integration of some combinations involving N and G resulting in vectors
//!
//! ## VEC 01: Interpolation functions times scalar field [vec_01_ns()]:
//!
//! ```text
//!      ⌠    → →     →
//! aᵐ = │ Nᵐ(x(ξ)) s(x) dΩ
//!      ⌡
//!      Ωₑ
//! ```
//!
//! ## VEC 02: Interpolation functions times vector field [vec_02_nv()]:
//!
//! ```text
//! →    ⌠    → →   → →
//! bᵐ = │ Nᵐ(x(ξ)) v(x) dΩ
//!      ⌡
//!      Ωₑ
//! ```
//!
//! ## VEC 03: Vector dot gradient [vec_03_vg()]:
//!
//! ```text
//!      ⌠ → →    →  → →
//! cᵐ = │ w(x) · Gᵐ(x(ξ)) dΩ
//!      ⌡
//!      Ωₑ
//! ```
//!
//! ## VEC 04: Tensor dot gradient [vec_04_tg()]:
//!
//! ```text
//! →    ⌠   →    →  → →
//! dᵐ = │ σ(x) · Gᵐ(x(ξ)) dΩ
//!      ⌡ ▔
//!      Ωₑ
//! ```
//!
//! # Matrix results: Integration of some combinations involving N, tensors, and G, resulting in matrices
//!
//! ![cases](https://github.com/cpmech/gemlab/raw/main/data/figures/integ-matrix-cases.png)
//!
//! ## MAT 01: Shape(N) times scalar(S) times shape(N) (e.g., diffusion matrix)
//!
//! [mat_01_nsn()]
//!
//! ```text
//!       ⌠
//! Kᵐⁿ = │ Nᵐ s Nⁿ dΩ
//!       ⌡
//!       Ωₑ
//! ```
//!
//! ## MAT 02: Gradient(G) dot vector(V) times shape(N) (e.g., compressibility matrix)
//!
//! [mat_02_gvn()]
//!
//! ```text
//!       ⌠ →    →
//! Kᵐⁿ = │ Gᵐ ⋅ v Nⁿ dΩ
//!       ⌡
//!       Ωₑ
//! ```
//!
//! ## MAT 03: Gradient(G) dot tensor(T) dot gradient(G) (e.g., conductivity matrix)
//!
//! [mat_03_gtg()]
//!
//! ```text
//!       ⌠ →        →
//! Kᵐⁿ = │ Gᵐ ⋅ T ⋅ Gⁿ dΩ
//!       ⌡      ▔
//!       Ωₑ
//! ```
//!
//! ## MAT 04: shape(Nb) time scalar(S) time gradient(G) (e.g., coupling matrix)
//!
//! [mat_04_nsg()]
//!
//! ```text
//! →     ⌠       →
//! Kᵐⁿ = │ Nbᵐ s Gⁿ dΩ
//!       ⌡
//!       Ωₑ
//! ```
//!
//! ## MAT 05: Gradient(Gb) time tensor(T) time shape(N) (e.g., coupling matrix)
//!
//! [mat_05_gtn()]
//!
//! ```text
//! →     ⌠ →
//! Kᵐⁿ = │ Gbᵐ ⋅ T Nⁿ dΩ
//!       ⌡       ▔
//!       Ωₑ
//! ```
//!
//! ## MAT 06: Shape(N) time vector(V) time shape(Nb) (e.g., coupling matrix)
//!
//! [mat_06_nvn()]
//!
//! ```text
//! →     ⌠    →
//! Kᵐⁿ = │ Nᵐ v Nbⁿ dΩ
//!       ⌡
//!       Ωₑ
//! ```
//!
//! ## MAT 07: Gradient(G) time scalar(S) time shape(Nb) (e.g., coupling matrix)
//!
//! [mat_07_gsn()]
//!
//! ```text
//! →     ⌠ →
//! Kᵐⁿ = │ Gᵐ s Nbⁿ dΩ
//!       ⌡
//!       Ωₑ
//! ```
//!
//! ## MAT 08: Shape(N) times tensor(T) times shape(N) (e.g., mass matrix)
//!
//! [mat_08_ntn()]
//!
//! ```text
//!       ⌠
//! Kᵐⁿ = │ Nᵐ T Nⁿ dΩ
//! ▔     ⌡    ▔
//!       Ωₑ
//! ```
//!
//! ## MAT 09: Shape(N) times vector(V) dot gradient(G) (e.g., variable density matrix)
//!
//! [mat_09_nvg()]
//!
//! ```text
//!       ⌠    →   →
//! Kᵐⁿ = │ Nᵐ v ⊗ Gⁿ dΩ
//! ▔     ⌡
//!       Ωₑ
//! ```
//!
//! ## MAT 10: Gradient(G) dot 4th-tensor(D) dot gradient(G) (e.g., stiffness matrix)
//!
//! [mat_10_gdg()]
//!
//! ```text
//!       ⌠                       →    →
//! Kᵐⁿ = │ Σ Σ Σ Σ Gᵐₖ Dᵢₖⱼₗ Gⁿₗ eᵢ ⊗ eⱼ dΩ
//! ▔     ⌡ i j k l
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
mod mat_01_nsn;
mod mat_02_gvn;
mod mat_03_gtg;
mod mat_04_nsg;
mod mat_05_gtn;
mod mat_06_nvn;
mod mat_07_gsn;
mod mat_08_ntn;
mod mat_09_nvg;
mod mat_10_gdg;
mod point_coords;
mod scalar_field;
mod testing;
mod vec_01_ns;
mod vec_02_nv;
mod vec_03_vg;
mod vec_04_tg;
mod vector_cases_boundary;
pub use crate::integ::analytical_qua4::*;
pub use crate::integ::analytical_qua8::*;
pub use crate::integ::analytical_tet4::*;
pub use crate::integ::analytical_tri3::*;
pub use crate::integ::integ_points::*;
pub use crate::integ::mat_01_nsn::*;
pub use crate::integ::mat_02_gvn::*;
pub use crate::integ::mat_03_gtg::*;
pub use crate::integ::mat_04_nsg::*;
pub use crate::integ::mat_05_gtn::*;
pub use crate::integ::mat_06_nvn::*;
pub use crate::integ::mat_07_gsn::*;
pub use crate::integ::mat_08_ntn::*;
pub use crate::integ::mat_09_nvg::*;
pub use crate::integ::mat_10_gdg::*;
pub use crate::integ::point_coords::*;
pub use crate::integ::scalar_field::*;
pub use crate::integ::vec_01_ns::*;
pub use crate::integ::vec_02_nv::*;
pub use crate::integ::vec_03_vg::*;
pub use crate::integ::vec_04_tg::*;
pub use crate::integ::vector_cases_boundary::*;
