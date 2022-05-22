use super::{geo_class_and_kind, ref_domain_limits, GeoClass, GeoKind, IntegPointData, StateOfShape};
use super::{
    Hex20, Hex8, Lin2, Lin3, Lin4, Lin5, Qua12, Qua16, Qua17, Qua4, Qua8, Qua9, Tet10, Tet4, Tri10, Tri15, Tri3, Tri6,
};
use crate::StrError;
use russell_lab::{inverse, mat_mat_mul, mat_vec_mul, vector_norm, Matrix, NormVec, Vector};
use std::fmt;

/// Defines an alias for interpolation functions
#[derive(Clone)]
struct FnInterp(fn(&mut Vector, &[f64]));

/// Defines an alias for derivative of interpolation functions
#[derive(Clone)]
struct FnDeriv(fn(&mut Matrix, &[f64]));

impl fmt::Debug for FnInterp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "FnInterp")
    }
}

impl fmt::Debug for FnDeriv {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "FnDeriv")
    }
}

/// Implements an isoparametric geometric shape for numerical integration and more
///
/// # Definitions
///
/// Here, we consider the following dimensions:
///
/// * `space_ndim` -- is the number of dimensions of the space under study (2 or 3)
/// * `geo_ndim` -- is the number of dimensions of the geometry element (shape),
///                  for instance, a line in the 2D space has `geo_ndim = 1` and
///                  `space_ndim = 2`. Another example is a triangle in the 3D space
///                  which has `geo_ndim = 2` and `space_ndim = 3`.
/// * `local` refers to a numbering scheme for the nodes of the shape (or element)
/// * `global` refers to a numbering scheme applied for the whole mesh
///
/// We also consider the following counting variables:
///
/// * `nnode` -- (local) number of points that define the shape (number of nodes)
/// * `npoint` -- (global) number of points in the whole mesh; not used here but important to remember
/// * `nedge` -- number of edges
/// * `nface` -- number of faces
/// * `edge_nnode` -- number of points that define the edge (number of nodes)
/// * `face_nnode` -- number of points that define the face (number of nodes)
/// * `face_nedge` -- face's number of edges
/// * `n_integ_point` -- number of integration (Gauss) points
///
/// When performing numerical integrations, we use the following notation:
/// |J| is the determinant of the Jacobian, ||J|| is the norm of the Jacobian vector
/// for line in multi-dimensions, `n_integ_point` is the number of integration points,
/// `ιp := ξp` is the reference coordinate of the integration point,
/// and `wp` is the weight of the p-th integration point.
///
/// # Isoparametric formulation
///
/// The isoparametric formulation establishes that
///
/// ```text
/// → →         →  →
/// x(ξ) = Σ Nᵐ(ξ) xᵐ
///        m         
/// ```
///
/// where `x` is the (space_ndim) vector of real coordinates, `ξ` is the (geo_ndim)
/// vector of reference coordinates, `Nm` are the (nnode) interpolation functions,
/// and `xm` are the (nnode) coordinates of each m-node of the geometric shape.
///
/// Given an (nnode,space_ndim) matrix of coordinates X, we can calculate the
/// (space_ndim) vector of coordinates x by means of
///
/// ```text
/// x = Xᵀ ⋅ N
/// ```
///
/// where `N` is an (nnode) array formed with all `Nm`.
///
/// # General case (geo_ndim == space_ndim)
///
/// If `geo_ndim == space_ndim`, we define the Jacobian tensor as
///
/// ```text
///         →
///   →    dx     →    →
/// J(ξ) = —— = Σ xᵐ ⊗ Lᵐ
///         →   m
///        dξ
/// ```
///
/// where
///
/// ```text
///             →
/// →  →    dNᵐ(ξ)
/// Lᵐ(ξ) = ——————
///            →
///           dξ
/// ```
///
/// are the derivatives of each interpolation function `Nm` with respect to the
/// reference coordinate. `Lm` are (geo_ndim) vectors and can be organized in
/// an (nnode,geo_ndim) matrix `L` of "local" derivatives.
///
/// We can write the Jacobian in matrix notation as follows
///
/// ```text
/// J = Xᵀ · L
/// ```
///
/// where X is the (nnode,space_ndim) matrix of coordinates and L is the (nnode,geo_ndim) matrix.
///
/// We define the gradient of interpolation functions (i.e., derivatives of interpolation
/// functions w.r.t real coordinates) by means of
///
/// ```text
///             →
/// →  →    dNᵐ(ξ)
/// Gᵐ(ξ) = ——————
///            →
///           dx
/// ```
///
/// which can be organized in an (nnode,space_ndim) matrix `G`.
///
/// The inverse Jacobian allows us to determine the gradient vectors G as follows
///
/// ```text
/// →       →  →        →
/// Gᵐ(ξ) = Lᵐ(ξ) · J⁻¹(ξ)
/// ```
///
/// Or, in matrix notation,
///
/// ```text
/// G = L · J⁻¹
/// ```
///
/// where G is an (nnode,space_ndim) matrix.
///
/// # Line in multi-dimensions (geo_ndim == 1 and space_ndim > 1)
///
/// In this case, the Jacobian equals the (space_ndim,1) base vector `g1` tangent
/// to the line element, i.e.,
///
/// ```text
///                          →
/// →    →     →            dx
/// J := Jline(ξ) = g₁(ξ) = —— = Xᵀ · L
///                         dξ
/// ```
///
/// We also consider a parametric coordinate `ℓ` which varies
/// from 0 to `ℓ_max` (the length of the line) according to
///
/// ```text
///                ℓ_max
/// ℓ(ξ) = (1 + ξ) —————
///                  2
///
///        2 · ℓ
/// ξ(ℓ) = ————— - 1
///        ℓ_max
/// ```
///
/// ```text
/// 0 ≤ ℓ ≤ ℓ_max
///
/// -1 ≤ ξ ≤ +1
/// ```
///
/// # Warning
///
/// All public properties are **readonly** and should **not** be modified externally.
#[derive(Clone, Debug)]
pub struct Shape {
    /// Geometry class
    pub class: GeoClass,

    /// Geometry kind
    pub kind: GeoKind,

    /// Space ndim
    pub space_ndim: usize,

    /// Geometry ndim
    pub geo_ndim: usize,

    /// Number of points that define the shape (shape's number of nodes)
    pub nnode: usize,

    /// Number of edges
    pub nedge: usize,

    /// Number of faces
    pub nface: usize,

    /// Number of points that define the edge (edge's number of nodes)
    pub edge_nnode: usize,

    /// Number of points that define the face (face's number of nodes)
    pub face_nnode: usize,

    /// Face's number of edges
    pub face_nedge: usize,

    /// Callback to evaluate interpolation functions
    fn_interp: FnInterp,

    /// Callback to evaluate local derivatives (w.r.t. ksi) of interpolation functions
    fn_deriv: FnDeriv,
}

impl Shape {
    /// Creates a new geometric shape
    ///
    /// # Input
    ///
    /// * `space_ndim` -- the space dimension (2 or 3)
    /// * `geo_ndim` -- the dimension of the shape; e.g. a 1D line in a 3D space
    /// * `nnode` -- the number of points defining the shape (number of nodes)
    ///
    /// # Note
    ///
    /// Some methods require that the coordinates matrix be set first.
    /// This can be accomplished by calling the `set_node` method.
    pub fn new(space_ndim: usize, geo_ndim: usize, nnode: usize) -> Result<Self, StrError> {
        if space_ndim < 2 || space_ndim > 3 {
            return Err("space_ndim must be 2 or 3");
        }
        // collect geometry data
        let (class, kind) = geo_class_and_kind(geo_ndim, nnode)?;
        let (nedge, nface, edge_nnode, face_nnode, face_nedge, fn_interp, fn_deriv): (
            usize,
            usize,
            usize,
            usize,
            usize,
            fn(&mut Vector, &[f64]),
            fn(&mut Matrix, &[f64]),
        ) = match kind {
            // Lin
            GeoKind::Lin2 => (
                Lin2::NEDGE,
                Lin2::NFACE,
                Lin2::EDGE_NNODE,
                Lin2::FACE_NNODE,
                Lin2::FACE_NEDGE,
                Lin2::calc_interp,
                Lin2::calc_deriv,
            ),
            GeoKind::Lin3 => (
                Lin3::NEDGE,
                Lin3::NFACE,
                Lin3::EDGE_NNODE,
                Lin3::FACE_NNODE,
                Lin3::FACE_NEDGE,
                Lin3::calc_interp,
                Lin3::calc_deriv,
            ),
            GeoKind::Lin4 => (
                Lin4::NEDGE,
                Lin4::NFACE,
                Lin4::EDGE_NNODE,
                Lin4::FACE_NNODE,
                Lin4::FACE_NEDGE,
                Lin4::calc_interp,
                Lin4::calc_deriv,
            ),
            GeoKind::Lin5 => (
                Lin5::NEDGE,
                Lin5::NFACE,
                Lin5::EDGE_NNODE,
                Lin5::FACE_NNODE,
                Lin5::FACE_NEDGE,
                Lin5::calc_interp,
                Lin5::calc_deriv,
            ),

            // Tri
            GeoKind::Tri3 => (
                Tri3::NEDGE,
                Tri3::NFACE,
                Tri3::EDGE_NNODE,
                Tri3::FACE_NNODE,
                Tri3::FACE_NEDGE,
                Tri3::calc_interp,
                Tri3::calc_deriv,
            ),
            GeoKind::Tri6 => (
                Tri6::NEDGE,
                Tri6::NFACE,
                Tri6::EDGE_NNODE,
                Tri6::FACE_NNODE,
                Tri6::FACE_NEDGE,
                Tri6::calc_interp,
                Tri6::calc_deriv,
            ),
            GeoKind::Tri10 => (
                Tri10::NEDGE,
                Tri10::NFACE,
                Tri10::EDGE_NNODE,
                Tri10::FACE_NNODE,
                Tri10::FACE_NEDGE,
                Tri10::calc_interp,
                Tri10::calc_deriv,
            ),
            GeoKind::Tri15 => (
                Tri15::NEDGE,
                Tri15::NFACE,
                Tri15::EDGE_NNODE,
                Tri15::FACE_NNODE,
                Tri15::FACE_NEDGE,
                Tri15::calc_interp,
                Tri15::calc_deriv,
            ),

            // Qua
            GeoKind::Qua4 => (
                Qua4::NEDGE,
                Qua4::NFACE,
                Qua4::EDGE_NNODE,
                Qua4::FACE_NNODE,
                Qua4::FACE_NEDGE,
                Qua4::calc_interp,
                Qua4::calc_deriv,
            ),
            GeoKind::Qua8 => (
                Qua8::NEDGE,
                Qua8::NFACE,
                Qua8::EDGE_NNODE,
                Qua8::FACE_NNODE,
                Qua8::FACE_NEDGE,
                Qua8::calc_interp,
                Qua8::calc_deriv,
            ),
            GeoKind::Qua9 => (
                Qua9::NEDGE,
                Qua9::NFACE,
                Qua9::EDGE_NNODE,
                Qua9::FACE_NNODE,
                Qua9::FACE_NEDGE,
                Qua9::calc_interp,
                Qua9::calc_deriv,
            ),
            GeoKind::Qua12 => (
                Qua12::NEDGE,
                Qua12::NFACE,
                Qua12::EDGE_NNODE,
                Qua12::FACE_NNODE,
                Qua12::FACE_NEDGE,
                Qua12::calc_interp,
                Qua12::calc_deriv,
            ),
            GeoKind::Qua16 => (
                Qua16::NEDGE,
                Qua16::NFACE,
                Qua16::EDGE_NNODE,
                Qua16::FACE_NNODE,
                Qua16::FACE_NEDGE,
                Qua16::calc_interp,
                Qua16::calc_deriv,
            ),
            GeoKind::Qua17 => (
                Qua17::NEDGE,
                Qua17::NFACE,
                Qua17::EDGE_NNODE,
                Qua17::FACE_NNODE,
                Qua17::FACE_NEDGE,
                Qua17::calc_interp,
                Qua17::calc_deriv,
            ),

            // Tet
            GeoKind::Tet4 => (
                Tet4::NEDGE,
                Tet4::NFACE,
                Tet4::EDGE_NNODE,
                Tet4::FACE_NNODE,
                Tet4::FACE_NEDGE,
                Tet4::calc_interp,
                Tet4::calc_deriv,
            ),
            GeoKind::Tet10 => (
                Tet10::NEDGE,
                Tet10::NFACE,
                Tet10::EDGE_NNODE,
                Tet10::FACE_NNODE,
                Tet10::FACE_NEDGE,
                Tet10::calc_interp,
                Tet10::calc_deriv,
            ),

            // Hex
            GeoKind::Hex8 => (
                Hex8::NEDGE,
                Hex8::NFACE,
                Hex8::EDGE_NNODE,
                Hex8::FACE_NNODE,
                Hex8::FACE_NEDGE,
                Hex8::calc_interp,
                Hex8::calc_deriv,
            ),
            GeoKind::Hex20 => (
                Hex20::NEDGE,
                Hex20::NFACE,
                Hex20::EDGE_NNODE,
                Hex20::FACE_NNODE,
                Hex20::FACE_NEDGE,
                Hex20::calc_interp,
                Hex20::calc_deriv,
            ),
        };

        // return new Shape
        Ok(Shape {
            class,
            kind,
            space_ndim,
            geo_ndim,
            nnode,
            nedge,
            nface,
            edge_nnode,
            face_nnode,
            face_nedge,
            fn_interp: FnInterp(fn_interp),
            fn_deriv: FnDeriv(fn_deriv),
        })
    }

    /// Calculates the interpolation functions
    ///
    /// Computes Nᵐ from:
    ///
    /// ```text
    /// → →         →  →
    /// u(ξ) = Σ Nᵐ(ξ) uᵐ
    ///        m         
    /// ```
    ///
    /// # Output
    ///
    /// * `state.interp` -- interpolation functions (nnode)
    ///
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate. The length of ξ must be equal to geo_ndim at least,
    ///            while lengths greater than geo_ndim are allowed (and ignored). In this way,
    ///            we can pass a slice with integration point data such as `[f64; 4]`.
    ///
    /// # Panics
    ///
    /// * This function does NOT check for sizes in `state`;
    ///   thus, make sure that `state` is compatible with this `shape`.
    #[inline]
    pub fn calc_interp(&self, state: &mut StateOfShape, ksi: &[f64]) -> Result<(), StrError> {
        if ksi.len() < self.geo_ndim {
            return Err("ksi.len() must be ≥ geo_ndim");
        }
        (self.fn_interp.0)(&mut state.interp, ksi);
        Ok(())
    }

    /// Calculates the derivatives of interpolation functions w.r.t reference coordinate
    ///
    /// Computes Lᵐ from:
    ///
    /// ```text
    ///             →
    /// →  →    dNᵐ(ξ)
    /// Lᵐ(ξ) = ——————
    ///            →
    ///           dξ
    /// ```
    ///
    /// # Output
    ///
    /// * `state.deriv` -- interpolation functions (nnode,geo_ndim)
    ///
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate. The length of ξ must be equal to geo_ndim at least,
    ///            while lengths greater than geo_ndim are allowed (and ignored). In this way,
    ///            we can pass a slice with integration point data such as `[f64; 4]`.
    ///
    /// # Panics
    ///
    /// * This function does NOT check for sizes in `state`;
    ///   thus, make sure that `state` is compatible with this `shape`.
    #[inline]
    pub fn calc_deriv(&self, state: &mut StateOfShape, ksi: &[f64]) -> Result<(), StrError> {
        if ksi.len() < self.geo_ndim {
            return Err("ksi.len() must be ≥ geo_ndim");
        }
        (self.fn_deriv.0)(&mut state.deriv, ksi);
        Ok(())
    }

    /// Calculates the real coordinates x from reference coordinates ξ
    ///
    /// The isoparametric formulation establishes:
    ///
    /// ```text
    /// → →         →  →
    /// x(ξ) = Σ Nᵐ(ξ) xᵐ
    ///        m         
    /// ```
    ///
    /// or
    ///
    /// ```text
    /// x := Xᵀ ⋅ N
    /// ```
    ///
    /// # Output
    ///
    /// * `x` -- real coordinate (space_ndim)
    /// * `state.interp` -- interpolation functions (nnode)
    ///
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate. The length of ξ must be equal to geo_ndim at least,
    ///            while lengths greater than geo_ndim are allowed (and ignored). In this way,
    ///            we can pass a slice with integration point data such as `[f64; 4]`.
    ///
    /// # Panics
    ///
    /// * This function does NOT check for sizes in `state`;
    ///   thus, make sure that `state` is compatible with this `shape`.
    pub fn calc_coords(&self, x: &mut Vector, state: &mut StateOfShape, ksi: &[f64]) -> Result<(), StrError> {
        if x.dim() != self.space_ndim {
            return Err("x.dim() must equal space_ndim");
        }
        self.calc_interp(state, ksi)?;
        mat_vec_mul(x, 1.0, &state.coords_transp, &state.interp)
    }

    /// Calculates the Jacobian of the mapping from general to reference space
    ///
    /// The components of the Jacobian matrix are
    ///
    /// ```text
    ///        ∂xᵢ
    /// Jᵢⱼ := ——— = Σ X[m][i] * L[m][j]
    ///        ∂ξⱼ   m
    /// ```
    ///
    /// Thus
    ///
    /// ```text
    /// jacobian := J = Xᵀ · L
    ///
    /// or (line in multi-dimensions, geom_ndim < space_ndim)
    ///
    /// jacobian := Jline = Xᵀ · L
    ///
    /// or (3D surface, geo_ndim = 2 and space_ndim = 3)
    ///
    /// jacobian := Jsurf = Xᵀ · L
    /// ```
    ///
    /// If `geo_ndim = space_ndim`, we also compute the inverse Jacobian
    ///
    /// ```text
    /// inv_jacobian := J⁻¹
    /// ```
    ///
    /// # Output
    ///
    /// * `state.deriv` -- interpolation functions (nnode)
    /// * `state.jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    /// * `state.inv_jacobian` -- If `geo_ndim = space_ndim`: inverse Jacobian matrix (space_ndim,space_ndim)
    ///
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate. The length of ξ must be equal to geo_ndim at least,
    ///            while lengths greater than geo_ndim are allowed (and ignored). In this way,
    ///            we can pass a slice with integration point data such as `[f64; 4]`.
    ///
    /// # Result
    ///
    /// * If `geo_ndim = space_ndim`, returns the determinant of the Jacobian
    /// * If `geo_ndim = 1` and `space_ndim > 1`, returns the norm of the Jacobian vector
    /// * Otherwise, returns zero
    ///
    /// # Panics
    ///
    /// * This function does NOT check for sizes in `state`;
    ///   thus, make sure that `state` is compatible with this `shape`.
    pub fn calc_jacobian(&self, state: &mut StateOfShape, ksi: &[f64]) -> Result<f64, StrError> {
        self.calc_deriv(state, ksi)?;
        mat_mat_mul(&mut state.jacobian, 1.0, &state.coords_transp, &state.deriv)?;
        if self.geo_ndim == self.space_ndim {
            inverse(&mut state.inv_jacobian, &state.jacobian)
        } else {
            if self.geo_ndim == 1 {
                let mut norm_jac = 0.0;
                for i in 0..self.space_ndim {
                    norm_jac += state.jacobian[i][0] * state.jacobian[i][0];
                }
                return Ok(f64::sqrt(norm_jac));
            }
            Ok(0.0)
        }
    }

    /// Computes the boundary normal vector
    ///
    /// **Note:** This function works with `geo_ndim < space_ndim` only. In particular we must have:
    ///
    /// * `geo_ndim = 1` and `space_ndim = 2` -- line in 2D, or
    /// * `geo_ndim = 2` and `space_ndim = 3` -- surface in 3D.
    ///
    /// # Line in multi-dimensions (geo_ndim == 1 and space_ndim > 1)
    ///
    /// Base vector tangent with the line:
    ///
    /// ```text
    ///          →
    ///         dx
    /// g₁(ξ) = —— = Xᵀ · L = first_column(J)
    ///         dξ
    /// ```
    ///
    /// Normal vector:
    ///
    /// ```text
    /// →   →    →
    /// n = e₃ × g₁
    ///
    ///   →       →
    /// ||n|| = ||g₁||
    /// ```
    ///
    /// Thus
    ///
    /// ```text
    ///        →           →
    /// dℓ = ||g₁|| dξ = ||n|| dξ
    /// ```
    ///
    /// # Boundary surface (geo_ndim == 2 and space_ndim == 3)
    ///
    /// Base vectors tangent to the surface:
    ///
    /// ```text
    ///          →
    /// →  →    dx
    /// g₁(ξ) = ——— = first_column(J)
    ///         dξ₁
    ///
    ///          →
    /// →  →    dx
    /// g₂(ξ) = ——— = second_column(J)
    ///         dξ₂
    /// ```
    ///
    /// Normal vector:
    ///
    /// ```text
    /// →   →    →
    /// n = g₁ × g₂
    /// ```
    ///
    /// Thus
    ///
    /// ```text
    ///         →
    /// dA := ||n|| dξ₁ dξ₂
    /// ```
    ///
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate. The length of ξ must be equal to geo_ndim at least,
    ///            while lengths greater than geo_ndim are allowed (and ignored). In this way,
    ///            we can pass a slice with integration point data such as `[f64; 4]`.
    ///
    /// # Output
    ///
    /// * `normal` -- (space_ndim) the boundary normal vector; not necessarily unitary
    /// * `state.deriv` -- interpolation functions (nnode)
    /// * `state.jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    ///
    /// # Panics
    ///
    /// * This function does NOT check for sizes in `state`;
    ///   thus, make sure that `state` is compatible with this `shape`.
    pub fn calc_boundary_normal(
        &self,
        normal: &mut Vector,
        state: &mut StateOfShape,
        ksi: &[f64],
    ) -> Result<(), StrError> {
        // check
        if self.geo_ndim >= self.space_ndim {
            return Err("geo_ndim must be smaller than space_ndim");
        }
        if normal.dim() != self.space_ndim {
            return Err("normal.dim() must equal space_ndim");
        }

        // compute Jacobian
        self.calc_deriv(state, ksi)?;
        mat_mat_mul(&mut state.jacobian, 1.0, &state.coords_transp, &state.deriv)?;

        // line in 2D (geo_ndim == 1 && self.space_ndim == 2)
        if self.space_ndim == 2 {
            normal[0] = -state.jacobian[1][0];
            normal[1] = state.jacobian[0][0];
            return Ok(());
        }

        // surface in 3D (geo_ndim == 2 && space_ndim == 3)
        let jj = &state.jacobian;
        normal[0] = jj[1][0] * jj[2][1] - jj[2][0] * jj[1][1];
        normal[1] = jj[2][0] * jj[0][1] - jj[0][0] * jj[2][1];
        normal[2] = jj[0][0] * jj[1][1] - jj[1][0] * jj[0][1];
        Ok(())
    }

    /// Approximates the reference coordinates from given real coordinates (inverse mapping)
    ///
    /// **Note:** This function works with `geo_ndim == space_ndim` only.
    ///
    /// We use Newton iterations with the inverse of the Jacobian to compute ξ(x).
    ///
    /// # Output
    ///
    /// * `ksi` -- ξ reference coordinates (geo_ndim = space_ndim)
    /// * `state.interp` -- interpolation functions (nnode)
    /// * `state.deriv` -- interpolation functions (nnode)
    /// * `state.jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    /// * `state.inv_jacobian` -- If `geo_ndim = space_ndim`: inverse Jacobian matrix (space_ndim,space_ndim)
    ///
    /// # Input
    ///
    /// * `x` -- real coordinates (space_ndim = geo_ndim)
    /// * `nit_max` -- maximum number of iterations (e.g., 10)
    /// * `tol` -- tolerance for the norm of the difference x - x(ξ) (e.g., 1e-14)
    ///
    /// # Results
    ///
    /// * Returns the number of iterations
    ///
    /// # Panics
    ///
    /// * This function does NOT check for sizes in `state`;
    ///   thus, make sure that `state` is compatible with this `shape`.
    pub fn approximate_ksi(
        &self,
        ksi: &mut [f64],
        state: &mut StateOfShape,
        x: &Vector,
        nit_max: usize,
        tol: f64,
    ) -> Result<usize, StrError> {
        // check
        if self.geo_ndim != self.space_ndim {
            return Err("geo_ndim must equal space_ndim");
        }
        if x.dim() != self.space_ndim {
            return Err("x.dim() must equal space_ndim");
        }
        if ksi.len() != self.geo_ndim {
            return Err("ksi.len() must equal geo_ndim");
        }

        // use linear scale to guess ksi
        let (min_ksi, _, del_ksi) = ref_domain_limits(self.class);
        for j in 0..self.geo_ndim {
            ksi[j] = (x[j] - state.coords_min[j]) / (state.coords_max[j] - state.coords_min[j]) * del_ksi + min_ksi;
        }

        // perform iterations
        let mut residual = Vector::new(self.space_ndim);
        let mut x_at_ksi = Vector::new(self.space_ndim);
        let mut delta_ksi = Vector::new(self.geo_ndim);
        for it in 0..nit_max {
            self.calc_coords(&mut x_at_ksi, state, &ksi)?;
            for i in 0..self.space_ndim {
                residual[i] = x[i] - x_at_ksi[i];
            }
            if vector_norm(&residual, NormVec::Euc) <= tol {
                return Ok(it);
            }
            self.calc_jacobian(state, ksi)?;
            mat_vec_mul(&mut delta_ksi, 1.0, &state.inv_jacobian, &residual).unwrap(); // cannot fail because all dims have been checked
            for j in 0..self.geo_ndim {
                ksi[j] += delta_ksi[j];
            }
        }

        Err("approximate_ksi failed to converge")
    }

    /// Calculates the gradient of the interpolation functions
    ///
    /// **Note:** This function works with `geo_ndim == space_ndim` only
    ///
    /// ```text
    ///             →
    /// →  →    dNᵐ(ξ)
    /// Gᵐ(ξ) = ——————
    ///            →
    ///           dx
    /// ```
    ///
    /// which can be organized in an (nnode,space_ndim) matrix `G` as follows
    ///
    /// ```text
    /// G = L · J⁻¹
    /// ```
    ///
    /// # Output
    ///
    /// * `state.deriv` -- interpolation functions (nnode)
    /// * `state.jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    /// * `state.inv_jacobian` -- inverse Jacobian matrix (space_ndim,space_ndim)
    /// * `state.gradient` -- gradient matrix (nnode,space_ndim)
    ///
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate. The length of ξ must be equal to geo_ndim at least,
    ///            while lengths greater than geo_ndim are allowed (and ignored). In this way,
    ///            we can pass a slice with integration point data such as `[f64; 4]`.
    ///
    /// # Results
    ///
    /// * Returns the determinant of the Jacobian.
    ///
    /// # Panics
    ///
    /// * This function does NOT check for sizes in `state`;
    ///   thus, make sure that `state` is compatible with this `shape`.
    pub fn calc_gradient(&self, state: &mut StateOfShape, ksi: &[f64]) -> Result<f64, StrError> {
        if self.geo_ndim != self.space_ndim {
            return Err("geo_ndim must equal space_ndim");
        }
        let det_jac = self.calc_jacobian(state, ksi)?;
        mat_mat_mul(&mut state.gradient, 1.0, &state.deriv, &state.inv_jacobian).unwrap(); // cannot fail because the dims are properly defined
        Ok(det_jac)
    }

    /// Calculates the real coordinates of all integration points
    ///
    /// We use the `calc_coords` method to apply the isoparametric formula
    /// to all p-th integration points `ιᵖ`. For example,
    ///
    /// ```text
    /// → →          → →   →
    /// x(ιᵖ) = Σ Nᵐ(ξ=ιᵖ) xᵐ
    ///         m         
    /// ```
    ///
    /// # Output
    ///
    /// * `state.interp` -- interpolation functions (nnode)
    ///
    /// # Input
    ///
    /// * `integ_points` -- Integration points' constants (n_integ_point)
    ///
    /// # Results
    ///
    /// * Returns an array with `n_integ_point` (number of integration points) vectors, where
    ///   each vector has a dimension equal to `space_ndim`.
    ///
    /// # Panics
    ///
    /// * This function does NOT check for sizes in `state`;
    ///   thus, make sure that `state` is compatible with this `shape`.
    pub fn calc_integ_points_coords(
        &self,
        state: &mut StateOfShape,
        integ_points: IntegPointData,
    ) -> Result<Vec<Vector>, StrError> {
        let mut all_coords = Vec::new();
        for iota in integ_points {
            let mut x = Vector::new(self.space_ndim);
            self.calc_coords(&mut x, state, iota)?;
            all_coords.push(x);
        }
        Ok(all_coords)
    }

    // --- getters ------------------------------------------------------------------------------

    /// Returns the local id of node on edge
    ///
    /// # Input
    ///
    /// * `e` -- index of edge in [0, nedge-1]
    /// * `i` -- index of local node [0, edge_nnode-1]
    pub fn edge_node_id(&self, e: usize, i: usize) -> usize {
        match self.kind {
            GeoKind::Lin2 => 0,
            GeoKind::Lin3 => 0,
            GeoKind::Lin4 => 0,
            GeoKind::Lin5 => 0,
            GeoKind::Tri3 => Tri3::EDGE_NODE_IDS[e][i],
            GeoKind::Tri6 => Tri6::EDGE_NODE_IDS[e][i],
            GeoKind::Tri10 => Tri10::EDGE_NODE_IDS[e][i],
            GeoKind::Tri15 => Tri15::EDGE_NODE_IDS[e][i],
            GeoKind::Qua4 => Qua4::EDGE_NODE_IDS[e][i],
            GeoKind::Qua8 => Qua8::EDGE_NODE_IDS[e][i],
            GeoKind::Qua9 => Qua9::EDGE_NODE_IDS[e][i],
            GeoKind::Qua12 => Qua12::EDGE_NODE_IDS[e][i],
            GeoKind::Qua16 => Qua16::EDGE_NODE_IDS[e][i],
            GeoKind::Qua17 => Qua17::EDGE_NODE_IDS[e][i],
            GeoKind::Tet4 => Tet4::EDGE_NODE_IDS[e][i],
            GeoKind::Tet10 => Tet10::EDGE_NODE_IDS[e][i],
            GeoKind::Hex8 => Hex8::EDGE_NODE_IDS[e][i],
            GeoKind::Hex20 => Hex20::EDGE_NODE_IDS[e][i],
        }
    }

    /// Returns the local id of node on face
    ///
    /// # Input
    ///
    /// * `f` -- index of face in [0, nface-1]
    /// * `i` -- index of local node [0, face_nnode-1]
    pub fn face_node_id(&self, f: usize, i: usize) -> usize {
        match self.kind {
            GeoKind::Lin2 => 0,
            GeoKind::Lin3 => 0,
            GeoKind::Lin4 => 0,
            GeoKind::Lin5 => 0,
            GeoKind::Tri3 => 0,
            GeoKind::Tri6 => 0,
            GeoKind::Tri10 => 0,
            GeoKind::Tri15 => 0,
            GeoKind::Qua4 => 0,
            GeoKind::Qua8 => 0,
            GeoKind::Qua9 => 0,
            GeoKind::Qua12 => 0,
            GeoKind::Qua16 => 0,
            GeoKind::Qua17 => 0,
            GeoKind::Tet4 => Tet4::FACE_NODE_IDS[f][i],
            GeoKind::Tet10 => Tet10::FACE_NODE_IDS[f][i],
            GeoKind::Hex8 => Hex8::FACE_NODE_IDS[f][i],
            GeoKind::Hex20 => Hex20::FACE_NODE_IDS[f][i],
        }
    }

    /// Returns the local node id on an edge on the face
    ///
    /// # Input
    ///
    /// * `f` -- index of face in [0, nface-1]
    /// * `k` -- index of face's edge (not the index of cell's edge) in [0, face_nedge-1]
    /// * `i` -- index of local node [0, edge_nnode-1]
    pub fn face_edge_node_id(&self, f: usize, k: usize, i: usize) -> usize {
        match self.kind {
            GeoKind::Lin2 => 0,
            GeoKind::Lin3 => 0,
            GeoKind::Lin4 => 0,
            GeoKind::Lin5 => 0,
            GeoKind::Tri3 => 0,
            GeoKind::Tri6 => 0,
            GeoKind::Tri10 => 0,
            GeoKind::Tri15 => 0,
            GeoKind::Qua4 => 0,
            GeoKind::Qua8 => 0,
            GeoKind::Qua9 => 0,
            GeoKind::Qua12 => 0,
            GeoKind::Qua16 => 0,
            GeoKind::Qua17 => 0,
            GeoKind::Tet4 => Tet4::FACE_EDGE_NODE_IDS[f][k][i],
            GeoKind::Tet10 => Tet10::FACE_EDGE_NODE_IDS[f][k][i],
            GeoKind::Hex8 => Hex8::FACE_EDGE_NODE_IDS[f][k][i],
            GeoKind::Hex20 => Hex20::FACE_EDGE_NODE_IDS[f][k][i],
        }
    }

    /// Returns the reference coordinates at node m
    ///
    /// # Output
    ///
    /// * `ksi` -- (geo_ndim) reference coordinates `ξᵐ` at node m
    pub fn reference_coords(&self, m: usize) -> &'static [f64] {
        match self.kind {
            GeoKind::Lin2 => &Lin2::NODE_REFERENCE_COORDS[m],
            GeoKind::Lin3 => &Lin3::NODE_REFERENCE_COORDS[m],
            GeoKind::Lin4 => &Lin4::NODE_REFERENCE_COORDS[m],
            GeoKind::Lin5 => &Lin5::NODE_REFERENCE_COORDS[m],
            GeoKind::Tri3 => &Tri3::NODE_REFERENCE_COORDS[m],
            GeoKind::Tri6 => &Tri6::NODE_REFERENCE_COORDS[m],
            GeoKind::Tri10 => &Tri10::NODE_REFERENCE_COORDS[m],
            GeoKind::Tri15 => &Tri15::NODE_REFERENCE_COORDS[m],
            GeoKind::Qua4 => &Qua4::NODE_REFERENCE_COORDS[m],
            GeoKind::Qua8 => &Qua8::NODE_REFERENCE_COORDS[m],
            GeoKind::Qua9 => &Qua9::NODE_REFERENCE_COORDS[m],
            GeoKind::Qua12 => &Qua12::NODE_REFERENCE_COORDS[m],
            GeoKind::Qua16 => &Qua16::NODE_REFERENCE_COORDS[m],
            GeoKind::Qua17 => &Qua17::NODE_REFERENCE_COORDS[m],
            GeoKind::Tet4 => &Tet4::NODE_REFERENCE_COORDS[m],
            GeoKind::Tet10 => &Tet10::NODE_REFERENCE_COORDS[m],
            GeoKind::Hex8 => &Hex8::NODE_REFERENCE_COORDS[m],
            GeoKind::Hex20 => &Hex20::NODE_REFERENCE_COORDS[m],
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{GeoClass, GeoKind, Shape};
    use crate::shapes::{ref_domain_limits, StateOfShape, IP_LIN_LEGENDRE_2};
    use crate::util::{PI, SQRT_3};
    use crate::StrError;
    use russell_chk::{assert_approx_eq, assert_deriv_approx_eq, assert_vec_approx_eq};
    use russell_lab::{copy_vector, scale_vector, vector_norm, Matrix, NormVec, Vector};
    use std::collections::HashMap;

    #[test]
    fn derive_works() {
        let shape = Shape::new(2, 1, 2).unwrap();
        let shape_clone = shape.clone();
        assert_eq!(format!("{:?}", shape), "Shape { class: Lin, kind: Lin2, space_ndim: 2, geo_ndim: 1, nnode: 2, nedge: 0, nface: 0, edge_nnode: 0, face_nnode: 0, face_nedge: 0, fn_interp: FnInterp, fn_deriv: FnDeriv }");
        assert_eq!(shape_clone.class, shape.class);
        assert_eq!(shape_clone.kind, shape.kind);
        assert_eq!(shape_clone.space_ndim, shape.space_ndim);
        assert_eq!(shape_clone.geo_ndim, shape.geo_ndim);
        assert_eq!(shape_clone.nnode, shape.nnode);
    }

    const RMIN: f64 = 1.0;
    const RMAX: f64 = 10.0;
    const AMIN: f64 = 30.0 * PI / 180.0;
    const AMAX: f64 = 60.0 * PI / 180.0;

    /// Generate coordinates
    ///
    /// The shape is the area indicated with "?" or the edge with "%".
    /// If class == Tri, the shape is half of the highlighted wedge.
    /// In 3D, an extrusion is applied along the out-of-plane direction.
    ///
    /// ```text
    ///   |            /
    ///   |           / αmax
    ///   ***=---__  /
    ///   |         % _
    ///   |        % ? *._          ,
    ///   |       % ????? *.     ,-'
    ///   ***=-_ % ???????? *.,-' αmin
    ///   |     % - ?????? ,-'*
    ///   |    /    *.? ,-'    *
    ///   |   /      ,*'        *
    ///   |  /    ,-'  *         *
    ///   | /  ,-'      *         *
    ///   |/.-'         #         #
    ///   o ----------- # ------- # --> r
    ///               rmin       rmax
    /// ```
    ///
    /// Intermediary mapping:
    ///
    /// r(ξ₀,ξ₁,ξ₂) = rmin + (ξ₀ - ξ₀min) * Δr / Δξ₀
    /// α(ξ₀,ξ₁,ξ₂) = αmin + (ξ₁ - ξ₁min) * Δα / Δξ₁
    /// z(ξ₀,ξ₁,ξ₂) = ξ₂
    ///
    /// Cylindrical coordinates:
    ///
    /// x₀ := r * cos(α)
    /// x₁ := r * sin(α)
    /// x₂ := z
    fn gen_coords(x: &mut Vector, ksi: &[f64], class: GeoClass) {
        assert_eq!(x.dim(), ksi.len());
        let (min_ksi, _, del_ksi) = ref_domain_limits(class);
        let r = RMIN + (ksi[0] - min_ksi) * (RMAX - RMIN) / del_ksi;
        let a = AMIN + (ksi[1] - min_ksi) * (AMAX - AMIN) / del_ksi;
        x[0] = r * f64::cos(a);
        x[1] = r * f64::sin(a);
        if x.dim() == 3 {
            x[2] = ksi[2];
        }
    }

    // Generates state with a pre-set matrix of coordinates
    fn gen_state_with_coords_matrix(shape: &Shape) -> Result<StateOfShape, StrError> {
        let mut x = Vector::new(shape.space_ndim);
        let mut ksi_aux = vec![0.0; shape.space_ndim];
        let mut coords: Vec<Vec<f64>> = Vec::new();
        for m in 0..shape.nnode {
            let ksi = shape.reference_coords(m);
            if shape.geo_ndim == shape.space_ndim {
                gen_coords(&mut x, ksi, shape.class);
            } else if shape.geo_ndim == 1 && shape.space_ndim == 2 {
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = 1.0;
                gen_coords(&mut x, &ksi_aux, shape.class);
            } else {
                // shape.geo_ndim == 2 && shape.space_ndim == 3
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = ksi[1];
                ksi_aux[2] = 1.0;
                gen_coords(&mut x, &ksi_aux, shape.class);
            }
            coords.push(x.as_data().clone());
        }
        StateOfShape::new(shape.geo_ndim, &coords)
    }

    // Generates state with a pre-set matrix of coordinates such that the
    // shape is parallel to the x,y,z axes.
    // For triangles and tetrahedra, one edge/face will cut the axes equally
    // * Qua and Hex will have real coordinates equal to the natural coordinates
    // * Tri and Tet will be scaled using the natural coordinates
    // * All edges parallel to the x,y,z axes will have lengths equal to 2.0
    fn gen_state_with_coords_matrix_parallel(shape: &Shape) -> Result<StateOfShape, StrError> {
        let mut scale = 1.0;
        if shape.class == GeoClass::Tri || shape.class == GeoClass::Tet {
            scale = 2.0;
        }
        let mut coords: Vec<Vec<f64>> = Vec::new();
        for m in 0..shape.nnode {
            let ksi = shape.reference_coords(m);
            coords.push(ksi.iter().map(|value| scale * value).collect());
        }
        StateOfShape::new(shape.geo_ndim, &coords)
    }

    #[test]
    fn capture_some_wrong_input() -> Result<(), StrError> {
        // new
        assert_eq!(Shape::new(1, 1, 1).err(), Some("space_ndim must be 2 or 3"));
        assert_eq!(
            Shape::new(2, 1, 1).err(),
            Some("(geo_ndim,nnode) combination is invalid")
        );

        // shape and state
        let shape = Shape::new(2, 1, 2)?;
        let mut state = StateOfShape::new(shape.geo_ndim, &[[0.0, 0.0], [1.0, 1.0]])?;

        // calc_interp
        assert_eq!(
            shape.calc_interp(&mut state, &[]).err(),
            Some("ksi.len() must be ≥ geo_ndim")
        );

        // calc_deriv
        assert_eq!(
            shape.calc_deriv(&mut state, &[]).err(),
            Some("ksi.len() must be ≥ geo_ndim")
        );

        // calc_coords
        let mut x = Vector::new(3);
        assert_eq!(
            shape.calc_coords(&mut x, &mut state, &[0.0, 0.0]).err(),
            Some("x.dim() must equal space_ndim")
        );
        let mut x = Vector::new(2);
        assert_eq!(
            shape.calc_coords(&mut x, &mut state, &[]).err(),
            Some("ksi.len() must be ≥ geo_ndim")
        );

        // calc_jacobian
        assert_eq!(
            shape.calc_jacobian(&mut state, &[]).err(),
            Some("ksi.len() must be ≥ geo_ndim")
        );
        let mut impossible = state.clone(); // this would never happen, unless there is a bug in StateOfShape
        impossible.jacobian = Matrix::new(0, 0);
        assert_eq!(
            shape.calc_jacobian(&mut impossible, &[0.0, 0.0]).err(),
            Some("matrices are incompatible")
        );

        // calc_boundary_normal
        let shape = Shape::new(2, 3, 4)?;
        let mut normal = Vector::new(2);
        let ksi = &[0.0, 0.0, 0.0];
        assert_eq!(
            shape.calc_boundary_normal(&mut normal, &mut state, ksi).err(),
            Some("geo_ndim must be smaller than space_ndim")
        );
        let shape = Shape::new(2, 1, 2)?;
        let mut normal = Vector::new(1);
        let mut state = StateOfShape::new(shape.geo_ndim, &[[0.0, 0.0], [1.0, 1.0]])?;
        assert_eq!(
            shape.calc_boundary_normal(&mut normal, &mut state, &[0.0]).err(),
            Some("normal.dim() must equal space_ndim")
        );
        let mut normal = Vector::new(2);
        assert_eq!(
            shape.calc_boundary_normal(&mut normal, &mut state, &[]).err(),
            Some("ksi.len() must be ≥ geo_ndim")
        );
        assert_eq!(
            shape
                .calc_boundary_normal(&mut normal, &mut impossible, &[0.0, 0.0])
                .err(),
            Some("matrices are incompatible")
        );

        // approximate_ksi
        let shape = Shape::new(2, 1, 2)?;
        let mut state = StateOfShape::new(shape.geo_ndim, &[[0.0, 0.0], [1.0, 1.0]])?;
        let mut ksi = vec![0.0; 1];
        let x = Vector::new(2);
        assert_eq!(
            shape.approximate_ksi(&mut ksi, &mut state, &x, 2, 1e-5).err(),
            Some("geo_ndim must equal space_ndim")
        );
        let shape = Shape::new(2, 2, 3)?;
        let mut state = StateOfShape::new(shape.geo_ndim, &[[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]])?;
        let x = Vector::new(1);
        assert_eq!(
            shape.approximate_ksi(&mut ksi, &mut state, &x, 2, 1e-5).err(),
            Some("x.dim() must equal space_ndim")
        );
        let mut ksi = vec![0.0; 3];
        let x = Vector::new(2);
        assert_eq!(
            shape.approximate_ksi(&mut ksi, &mut state, &x, 2, 1e-5).err(),
            Some("ksi.len() must equal geo_ndim")
        );
        let mut ksi = vec![0.0; 2];
        assert_eq!(
            shape.approximate_ksi(&mut ksi, &mut state, &x, 0, 1e-5).err(),
            Some("approximate_ksi failed to converge")
        );
        let mut bug_for_calc_coords = state.clone(); // this would never happen, unless there is a bug in StateOfShape
        bug_for_calc_coords.coords_transp = Matrix::new(0, 0);
        assert_eq!(
            shape
                .approximate_ksi(&mut ksi, &mut bug_for_calc_coords, &x, 2, 1e-5)
                .err(),
            Some("matrix and vectors are incompatible")
        );
        let mut bug_for_calc_jacobian = state.clone(); // this would never happen, unless there is a bug in StateOfShape
        bug_for_calc_jacobian.jacobian = Matrix::new(0, 0);
        let x = Vector::from(&[1000.0, 1000.0]);
        assert_eq!(
            shape
                .approximate_ksi(&mut ksi, &mut bug_for_calc_jacobian, &x, 2, 1e-5)
                .err(),
            Some("matrices are incompatible")
        );

        // calc_gradient
        assert_eq!(
            shape.calc_gradient(&mut bug_for_calc_jacobian, &[0.0, 0.0]).err(),
            Some("matrices are incompatible")
        );
        let shape = Shape::new(2, 1, 2)?;
        let mut state = StateOfShape::new(shape.geo_ndim, &[[0.0, 0.0], [1.0, 1.0]])?;
        assert_eq!(
            shape.calc_gradient(&mut state, &[0.0]).err(),
            Some("geo_ndim must equal space_ndim")
        );

        // calc_integ_points_coords
        assert_eq!(
            shape
                .calc_integ_points_coords(&mut bug_for_calc_coords, &IP_LIN_LEGENDRE_2)
                .err(),
            Some("matrix and vectors are incompatible")
        );

        Ok(())
    }

    #[test]
    fn getters_work() -> Result<(), StrError> {
        for (geo_ndim, nnode) in GeoKind::PAIRS {
            let space_ndim = usize::max(2, geo_ndim);
            let shape = &mut Shape::new(space_ndim, geo_ndim, nnode)?;
            match shape.class {
                GeoClass::Lin => assert_eq!(shape.edge_node_id(0, 0), 0),
                GeoClass::Tri => assert_eq!(shape.edge_node_id(0, 0), 1),
                GeoClass::Qua => assert_eq!(shape.edge_node_id(0, 0), 1),
                GeoClass::Tet => assert_eq!(shape.edge_node_id(0, 0), 0),
                GeoClass::Hex => assert_eq!(shape.edge_node_id(0, 0), 0),
            }
            match shape.class {
                GeoClass::Lin => assert_eq!(shape.face_node_id(0, 0), 0),
                GeoClass::Tri => assert_eq!(shape.face_node_id(0, 0), 0),
                GeoClass::Qua => assert_eq!(shape.face_node_id(0, 0), 0),
                GeoClass::Tet => assert_eq!(shape.face_node_id(0, 0), 0),
                GeoClass::Hex => assert_eq!(shape.face_node_id(0, 0), 0),
            }
            match shape.class {
                GeoClass::Lin => assert_eq!(shape.face_edge_node_id(0, 0, 0), 0),
                GeoClass::Tri => assert_eq!(shape.face_edge_node_id(0, 0, 0), 0),
                GeoClass::Qua => assert_eq!(shape.face_edge_node_id(0, 0, 0), 0),
                GeoClass::Tet => assert_eq!(shape.face_edge_node_id(0, 0, 0), 0),
                GeoClass::Hex => assert_eq!(shape.face_edge_node_id(0, 0, 0), 0),
            }
        }
        Ok(())
    }

    #[test]
    fn calc_interp_works() -> Result<(), StrError> {
        // define tolerances
        let mut tols = HashMap::new();
        tols.insert(GeoKind::Lin2, 1e-15);
        tols.insert(GeoKind::Lin3, 1e-15);
        tols.insert(GeoKind::Lin4, 1e-15);
        tols.insert(GeoKind::Lin5, 1e-15);
        tols.insert(GeoKind::Tri3, 1e-15);
        tols.insert(GeoKind::Tri6, 1e-15);
        tols.insert(GeoKind::Tri10, 1e-15);
        tols.insert(GeoKind::Tri15, 1e-15);
        tols.insert(GeoKind::Qua4, 1e-15);
        tols.insert(GeoKind::Qua8, 1e-15);
        tols.insert(GeoKind::Qua9, 1e-15);
        tols.insert(GeoKind::Qua12, 1e-15);
        tols.insert(GeoKind::Qua16, 1e-15);
        tols.insert(GeoKind::Qua17, 1e-15);
        tols.insert(GeoKind::Tet4, 1e-15);
        tols.insert(GeoKind::Tet10, 1e-15);
        tols.insert(GeoKind::Hex8, 1e-15);
        tols.insert(GeoKind::Hex20, 1e-15);

        // loop over shapes
        for (geo_ndim, nnode) in GeoKind::PAIRS {
            // allocate shape and state
            let space_ndim = usize::max(2, geo_ndim);
            let shape = Shape::new(space_ndim, geo_ndim, nnode)?;
            let mut state = StateOfShape::new(
                shape.geo_ndim,
                &(0..shape.nnode)
                    .map(|_| vec![0.0; shape.space_ndim])
                    .collect::<Vec<_>>(),
            )?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // loop over nodes of shape
            for m in 0..shape.nnode {
                // get ξᵐ corresponding to node m
                let ksi = shape.reference_coords(m);

                // compute interpolation function Nⁿ(ξᵐ)
                shape.calc_interp(&mut state, ksi)?;

                // check: Nⁿ(ξᵐ) = 1 if m==n; 0 otherwise
                for n in 0..shape.nnode {
                    if m == n {
                        assert_approx_eq!(state.interp[n], 1.0, tol);
                    } else {
                        assert_approx_eq!(state.interp[n], 0.0, tol);
                    }
                }
            }
        }
        Ok(())
    }

    // Holds arguments for numerical differentiation of N w.r.t ξ => L (deriv) matrix
    struct ArgsNumL<'a> {
        shape: &'a Shape,    // shape
        state: StateOfShape, // auxiliary (copy) state
        at_ksi: Vec<f64>,    // at reference coord value
        ksi: Vec<f64>,       // temporary reference coord
        m: usize,            // node index from 0 to nnode
        j: usize,            // dimension index from 0 to geom_ndim
    }

    // Computes Nᵐ(ξ) with variable v := ξⱼ
    fn aux_deriv(v: f64, args: &mut ArgsNumL) -> f64 {
        args.ksi.copy_from_slice(&args.at_ksi);
        args.ksi[args.j] = v;
        args.shape.calc_interp(&mut args.state, &args.ksi).unwrap();
        args.state.interp[args.m]
    }

    #[test]
    fn calc_deriv_works() -> Result<(), StrError> {
        // define tolerances
        let mut tols = HashMap::new();
        tols.insert(GeoKind::Lin2, 1e-13);
        tols.insert(GeoKind::Lin3, 1e-13);
        tols.insert(GeoKind::Lin4, 1e-10);
        tols.insert(GeoKind::Lin5, 1e-10);
        tols.insert(GeoKind::Tri3, 1e-12);
        tols.insert(GeoKind::Tri6, 1e-12);
        tols.insert(GeoKind::Tri10, 1e-10);
        tols.insert(GeoKind::Tri15, 1e-9);
        tols.insert(GeoKind::Qua4, 1e-13);
        tols.insert(GeoKind::Qua8, 1e-12);
        tols.insert(GeoKind::Qua9, 1e-13);
        tols.insert(GeoKind::Qua12, 1e-10);
        tols.insert(GeoKind::Qua16, 1e-10);
        tols.insert(GeoKind::Qua17, 1e-10);
        tols.insert(GeoKind::Hex8, 1e-13);
        tols.insert(GeoKind::Tet4, 1e-12);
        tols.insert(GeoKind::Tet10, 1e-12);
        tols.insert(GeoKind::Hex20, 1e-12);

        // loop over shapes
        for (geo_ndim, nnode) in GeoKind::PAIRS {
            // allocate shape and state
            let space_ndim = usize::max(2, geo_ndim);
            let shape = Shape::new(space_ndim, geo_ndim, nnode)?;
            let mut state = StateOfShape::new(
                shape.geo_ndim,
                &(0..shape.nnode)
                    .map(|_| vec![0.0; shape.space_ndim])
                    .collect::<Vec<_>>(),
            )?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // set ξ within reference space
            let at_ksi = vec![0.25; shape.geo_ndim];

            // compute all derivatives of interpolation functions w.r.t ξ
            shape.calc_deriv(&mut state, &at_ksi)?;

            // set arguments for numerical integration
            let args = &mut ArgsNumL {
                shape: &shape,
                state: state.clone(),
                at_ksi,
                ksi: vec![0.0; shape.geo_ndim],
                m: 0,
                j: 0,
            };

            // check Lᵐ(ξ) = dNᵐ(ξ)/dξ
            for m in 0..shape.nnode {
                args.m = m;
                for j in 0..shape.geo_ndim {
                    args.j = j;
                    // Lᵐⱼ := dNᵐ/dξⱼ
                    assert_deriv_approx_eq!(state.deriv[m][j], args.at_ksi[j], aux_deriv, args, tol);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn calc_coords_works() -> Result<(), StrError> {
        // define dims and number of nodes
        let pairs = vec![
            (2, 3),
            (2, 6),
            (2, 10),
            (2, 15),
            (2, 4),
            (2, 8),
            (2, 17),
            (3, 4),
            (3, 10),
            (3, 8),
            (3, 20),
        ];

        // define tolerances
        let mut tols = HashMap::new();
        tols.insert(GeoKind::Qua4, 1e-15);
        tols.insert(GeoKind::Tri3, 1e-15);
        tols.insert(GeoKind::Tri6, 1e-15);
        tols.insert(GeoKind::Tri10, 1e-14);
        tols.insert(GeoKind::Tri15, 1e-15);
        tols.insert(GeoKind::Qua8, 1e-15);
        tols.insert(GeoKind::Qua17, 1e-15);
        tols.insert(GeoKind::Hex8, 1e-15);
        tols.insert(GeoKind::Tet4, 1e-15);
        tols.insert(GeoKind::Tet10, 1e-15);
        tols.insert(GeoKind::Hex20, 1e-15);

        // define tolerances for node in the reference domain
        let mut tols_in = HashMap::new();
        tols_in.insert(GeoKind::Tri3, 0.45); // linear maps are inaccurate for the circular wedge
        tols_in.insert(GeoKind::Tri6, 0.02); // << quadratic mapping is inaccurate as well
        tols_in.insert(GeoKind::Tri10, 1e-14);
        tols_in.insert(GeoKind::Tri15, 1e-4); // << this triangle is inaccurate as well here
        tols_in.insert(GeoKind::Qua4, 0.14); // linear maps are inaccurate for the circular wedge
        tols_in.insert(GeoKind::Qua8, 1e-14);
        tols_in.insert(GeoKind::Qua17, 1e-14);
        tols_in.insert(GeoKind::Tet4, 0.45); // linear tetrahedron is also inaccurate here
        tols_in.insert(GeoKind::Tet10, 0.02); // quadratic tetrahedron is also inaccurate here
        tols_in.insert(GeoKind::Hex8, 0.14); // bi-linear maps are inaccurate for the circular wedge
        tols_in.insert(GeoKind::Hex20, 1e-14);

        // loop over shapes
        for (geo_ndim, nnode) in pairs {
            // allocate shape
            let space_ndim = usize::max(2, geo_ndim);
            let shape = Shape::new(space_ndim, geo_ndim, nnode)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();
            let tol_in = *tols_in.get(&shape.kind).unwrap();

            // generate state with coordinates matrix
            let mut state = gen_state_with_coords_matrix(&shape)?;

            // loop over nodes of shape
            let mut x = Vector::new(shape.space_ndim);
            let mut x_correct = Vector::new(shape.space_ndim);
            for m in 0..shape.nnode {
                // get ξᵐ corresponding to node m
                let ksi = shape.reference_coords(m);

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                shape.calc_coords(&mut x, &mut state, ksi)?;

                // compare xᵐ with generated coordinates
                gen_coords(&mut x_correct, ksi, shape.class);
                assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol);
            }

            // test again inside the reference domain
            let ksi_in = if shape.class == GeoClass::Tri || shape.class == GeoClass::Tet {
                vec![1.0 / 3.0; shape.geo_ndim]
            } else {
                vec![0.0; shape.geo_ndim]
            };
            shape.calc_coords(&mut x, &mut state, &ksi_in)?;
            gen_coords(&mut x_correct, &ksi_in, shape.class);
            assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol_in);
        }
        Ok(())
    }

    // Holds arguments for numerical differentiation of x w.r.t ξ => Jacobian
    struct ArgsNumJ<'a> {
        shape: &'a Shape,    // shape
        state: StateOfShape, // auxiliary (copy) state
        at_ksi: Vec<f64>,    // at reference coord value
        ksi: Vec<f64>,       // temporary reference coord
        x: Vector,           // (space_ndim) coordinates at ξ
        i: usize,            // dimension index from 0 to space_ndim
        j: usize,            // dimension index from 0 to geo_ndim
    }

    // Computes xᵢ(ξ) with variable v := ξⱼ
    fn aux_jacobian(v: f64, args: &mut ArgsNumJ) -> f64 {
        args.ksi.copy_from_slice(&args.at_ksi);
        args.ksi[args.j] = v;
        args.shape.calc_coords(&mut args.x, &mut args.state, &args.ksi).unwrap();
        args.x[args.i]
    }

    #[test]
    fn calc_jacobian_works() -> Result<(), StrError> {
        // define tolerances
        let mut tols = HashMap::new();
        tols.insert(GeoKind::Lin2, 1e-11);
        tols.insert(GeoKind::Lin3, 1e-11);
        tols.insert(GeoKind::Lin4, 1e-11);
        tols.insert(GeoKind::Lin5, 1e-11);
        tols.insert(GeoKind::Tri3, 1e-11);
        tols.insert(GeoKind::Tri6, 1e-11);
        tols.insert(GeoKind::Tri10, 1e-11);
        tols.insert(GeoKind::Tri15, 1e-10);
        tols.insert(GeoKind::Qua4, 1e-11);
        tols.insert(GeoKind::Qua8, 1e-11);
        tols.insert(GeoKind::Qua9, 1e-11);
        tols.insert(GeoKind::Qua12, 1e-10);
        tols.insert(GeoKind::Qua16, 1e-10);
        tols.insert(GeoKind::Qua17, 1e-10);
        tols.insert(GeoKind::Hex8, 1e-11);
        tols.insert(GeoKind::Hex20, 1e-11);
        tols.insert(GeoKind::Tet4, 1e-12);
        tols.insert(GeoKind::Tet10, 1e-12);

        // loop over shapes
        for (geo_ndim, nnode) in GeoKind::PAIRS {
            // allocate shape
            let space_ndim = usize::max(2, geo_ndim);
            let shape = Shape::new(space_ndim, geo_ndim, nnode)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // generate state with coordinates matrix
            let mut state = gen_state_with_coords_matrix(&shape)?;

            // set ξ within reference space
            let at_ksi = vec![0.25; shape.geo_ndim];

            // compute Jacobian, its inverse, and determinant
            let det_jac = shape.calc_jacobian(&mut state, &at_ksi)?;
            assert!(det_jac > 0.0);

            // set arguments for numerical integration
            let args = &mut ArgsNumJ {
                shape: &shape,
                state: state.clone(),
                at_ksi,
                ksi: vec![0.0; shape.geo_ndim],
                x: Vector::new(shape.space_ndim),
                i: 0,
                j: 0,
            };

            // check J(ξ) = dx(ξ)/dξ
            for i in 0..shape.space_ndim {
                args.i = i;
                for j in 0..shape.geo_ndim {
                    args.j = j;
                    // Jᵢⱼ := dxᵢ/dξⱼ
                    assert_deriv_approx_eq!(state.jacobian[i][j], args.at_ksi[j], aux_jacobian, args, tol);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn calc_jacobian_special_cases_work() -> Result<(), StrError> {
        let shape = Shape::new(2, 1, 2)?;
        let l = 3.5;
        let mut state = StateOfShape::new(shape.geo_ndim, &[[0.0, 0.0], [l, 0.0]])?;
        let norm_jac_vec = shape.calc_jacobian(&mut state, &[0.0])?;
        assert_eq!(norm_jac_vec, l / 2.0); // 2.0 = length of shape in the reference space

        let shape = Shape::new(3, 2, 3)?;
        let mut state = StateOfShape::new(shape.geo_ndim, &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 1.0]])?;
        let norm_jac_vec = shape.calc_jacobian(&mut state, &[0.0, 0.0])?;
        assert_eq!(norm_jac_vec, 0.0);
        Ok(())
    }

    #[test]
    fn calc_boundary_normal_works() -> Result<(), StrError> {
        // allocate surface shape and state
        let shape = Shape::new(3, 2, 17)?;

        // generate state with coordinates matrix
        let mut state = gen_state_with_coords_matrix(&shape)?;

        // compute boundary normal vector
        let at_ksi = vec![0.0; shape.geo_ndim];
        let mut normal = Vector::new(shape.space_ndim);
        shape.calc_boundary_normal(&mut normal, &mut state, &at_ksi)?;

        // check magnitude of normal vector
        let mag_normal = vector_norm(&normal, NormVec::Euc);
        let area = PI * (RMAX * RMAX - RMIN * RMIN) / 12.0;
        let ref_area = 4.0;
        let area_ratio = area / ref_area;
        assert_approx_eq!(mag_normal, area_ratio, 1e-4);

        // check direction of normal vector
        let mut unit_normal = Vector::from(normal.as_data());
        scale_vector(&mut unit_normal, 1.0 / mag_normal);
        assert_vec_approx_eq!(unit_normal.as_data(), &[0.0, 0.0, 1.0], 1e-15);
        Ok(())
    }

    #[test]
    fn calc_boundary_normal_edge_works() -> Result<(), StrError> {
        //  The edge is indicated by "%" in the figure below
        // |           / αmax
        // ***=---__  /
        // |        `%._
        // |        %   *._          ,
        // |       %       *.     ,-'
        // ***=-_ %          *.,-' αmin
        // |     % -        ,-'*
        // |    /    *.  ,-'    *
        // |   /      ,*'        *
        // |  /    ,-'  *         *
        // | /  ,-'      *         *
        // |/.-'         #         #
        // o ----------- # ------- # --> r
        //             rmin       rmax

        // allocate boundary edge
        let shape = Shape::new(2, 1, 5)?;

        // generate state with coordinates matrix
        let mut state = gen_state_with_coords_matrix(&shape)?;

        // compute boundary normal vector
        let at_ksi = vec![0.0; shape.geo_ndim];
        let mut normal = Vector::new(shape.space_ndim);
        shape.calc_boundary_normal(&mut normal, &mut state, &at_ksi)?;

        // check magnitude of normal vector
        let mag_normal = vector_norm(&normal, NormVec::Euc);
        let length = RMAX - RMIN;
        let ref_length = 2.0;
        let length_ratio = length / ref_length;
        assert_approx_eq!(mag_normal, length_ratio, 1e-15);

        // check direction of normal vector
        let mut unit_normal = Vector::from(normal.as_data());
        scale_vector(&mut unit_normal, 1.0 / mag_normal);
        assert_vec_approx_eq!(unit_normal.as_data(), &[-f64::sin(AMAX), f64::cos(AMAX)], 1e-15);
        Ok(())
    }

    #[test]
    fn normals_are_outward_2d() -> Result<(), StrError> {
        // define dims and number of nodes
        let pairs = vec![
            // Tri
            (2, 3),
            (2, 6),
            (2, 10),
            (2, 15),
            // Qua
            (2, 4),
            (2, 8),
            (2, 9),
            (2, 12),
            (2, 16),
            (2, 17),
        ];

        // solution
        //
        //   →     Δℓ_edge   Δℓ_edge
        // ||n|| = ——————— = ———————
        //         Δξ_lin       2
        //                                                 →     2 √2
        // For the diagonal of Tri: Δℓ_edge = 2 √2, thus ||n|| = ————— = √2
        //                                                         2
        let tri_correct = vec![
            &[0.0, -1.0], // bottom
            &[1.0, 1.0],  // diagonal
            &[-1.0, 0.0], // left
        ];
        let qua_correct = vec![
            &[0.0, -1.0], // bottom
            &[1.0, 0.0],  // right
            &[0.0, 1.0],  // top
            &[-1.0, 0.0], // left
        ];

        // auxiliary
        let mut normal = Vector::new(2);
        let ksi = &[0.0, 0.0, 0.0];

        // loop over shapes
        for (geo_ndim, nnode) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let shape = Shape::new(space_ndim, geo_ndim, nnode)?;

            // generate state with coordinates matrix
            let state = gen_state_with_coords_matrix_parallel(&shape)?;

            // loop over edges
            for e in 0..shape.nedge {
                let edge_shape = Shape::new(space_ndim, 1, shape.edge_nnode)?;
                let mut edge_state = StateOfShape::new(
                    edge_shape.geo_ndim,
                    &(0..edge_shape.nnode)
                        .map(|i| {
                            let m = shape.edge_node_id(e, i);
                            (0..shape.space_ndim).map(|j| state.coords_transp[j][m]).collect()
                        })
                        .collect::<Vec<_>>(),
                )?;

                // calc normal vector
                edge_shape.calc_boundary_normal(&mut normal, &mut edge_state, ksi)?;
                if shape.class == GeoClass::Tri {
                    // check triangle
                    assert_vec_approx_eq!(normal.as_data(), tri_correct[e], 1e-15);
                } else {
                    // check quadrilateral
                    assert_vec_approx_eq!(normal.as_data(), qua_correct[e], 1e-15);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn normals_are_outward_3d() -> Result<(), StrError> {
        // define dims and number of nodes
        let pairs = vec![
            // Tet
            (3, 4),
            (3, 10),
            // Hex
            (3, 8),
            (3, 20),
        ];

        // solution
        // Hex with sides equal to h=2:
        //
        //   →         ΔA_face       h²
        // ||n|| = ——————————————— = —— = 1
        //         Δξ₁_qua Δξ₂_qua    4
        //
        // Tet with sides equal to h=2:
        //
        //   Faces orthogonal to the x,y,z axes:
        //
        //      →     ΔA_face   h²/2
        //    ||n|| = ——————— = ———— = 4
        //             ΔA_tri    1/2
        //
        //   Face orthogonal to the the diagonal:
        //
        //      →     ΔA_face   √3 h²/2
        //    ||n|| = ——————— = ——————— = 4 √3
        //             ΔA_tri     1/2
        let tet_correct = vec![
            &[-4.0, 0.0, 0.0], // negative-x face
            &[0.0, -4.0, 0.0], // negative-y face
            &[0.0, 0.0, -4.0], // negative-z face
            &[4.0, 4.0, 4.0],  // face orthogonal to the diagonal
        ];
        let hex_correct = vec![
            &[-1.0, 0.0, 0.0], // behind
            &[1.0, 0.0, 0.0],  // front
            &[0.0, -1.0, 0.0], // left
            &[0.0, 1.0, 0.0],  // right
            &[0.0, 0.0, -1.0], // bottom
            &[0.0, 0.0, 1.0],  // top
        ];

        // auxiliary
        let mut normal = Vector::new(3);
        let ksi = &[0.0, 0.0, 0.0];

        // loop over shapes
        for (geo_ndim, nnode) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let shape = Shape::new(space_ndim, geo_ndim, nnode)?;

            // generate state with coordinates matrix
            let state = gen_state_with_coords_matrix_parallel(&shape)?;

            // loop over faces
            for f in 0..shape.nface {
                let face_shape = Shape::new(space_ndim, 2, shape.face_nnode)?;
                let mut face_state = StateOfShape::new(
                    face_shape.geo_ndim,
                    &(0..face_shape.nnode)
                        .map(|i| {
                            let m = shape.face_node_id(f, i);
                            (0..shape.space_ndim).map(|j| state.coords_transp[j][m]).collect()
                        })
                        .collect::<Vec<_>>(),
                )?;

                // calc normal vector
                face_shape.calc_boundary_normal(&mut normal, &mut face_state, ksi)?;
                if shape.class == GeoClass::Tet {
                    // check tetrahedron
                    assert_vec_approx_eq!(normal.as_data(), tet_correct[f], 1e-15);
                } else {
                    // check hexahedron
                    assert_vec_approx_eq!(normal.as_data(), hex_correct[f], 1e-15);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn approximate_ksi_works() -> Result<(), StrError> {
        // define dims and number of nodes
        let pairs = vec![(2, 3), (2, 6), (2, 4), (2, 8), (3, 4), (3, 10), (3, 8), (3, 20)];

        // define tolerances
        let mut tols = HashMap::new();
        tols.insert(GeoKind::Tri3, 1e-14);
        tols.insert(GeoKind::Tri6, 1e-15);
        tols.insert(GeoKind::Qua4, 1e-15);
        tols.insert(GeoKind::Qua8, 1e-15);
        tols.insert(GeoKind::Tet4, 1e-15);
        tols.insert(GeoKind::Tet10, 1e-15);
        tols.insert(GeoKind::Hex8, 1e-15);
        tols.insert(GeoKind::Hex20, 1e-15);

        // loop over shapes
        for (geo_ndim, nnode) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let shape = Shape::new(space_ndim, geo_ndim, nnode)?;

            // generate state with coordinates matrix
            let mut state = gen_state_with_coords_matrix(&shape)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // loop over nodes of shape
            let mut x = Vector::new(shape.space_ndim);
            let mut ksi = vec![0.0; shape.geo_ndim];
            for m in 0..shape.nnode {
                // get ξᵐ corresponding to node m
                let ksi_ref = shape.reference_coords(m);

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                shape.calc_coords(&mut x, &mut state, ksi_ref)?;

                // compute approximation of the inverse mapping ξᵐ(xᵐ)
                let nit = shape.approximate_ksi(&mut ksi, &mut state, &x, 10, 1e-14)?;

                // check (linear and bi-linear shapes converge with nit = 1)
                if shape.kind == GeoKind::Tri3 || shape.kind == GeoKind::Qua4 || shape.kind == GeoKind::Hex8 {
                    assert_eq!(nit, 1);
                }
                assert_vec_approx_eq!(ksi, ksi_ref, tol);
            }

            // test again inside the reference domain
            let ksi_in = if shape.class == GeoClass::Tri || shape.class == GeoClass::Tet {
                vec![1.0 / 3.0; shape.geo_ndim]
            } else {
                vec![0.0; shape.geo_ndim]
            };
            shape.calc_coords(&mut x, &mut state, &ksi_in)?;
            shape.approximate_ksi(&mut ksi, &mut state, &x, 10, 1e-14)?;
            assert_vec_approx_eq!(ksi, ksi_in, tol);
        }
        Ok(())
    }

    /// Returns Shape and StateOfShape for an equilateral triangle
    ///
    /// ```text
    /// Equilateral triangle
    ///
    ///           /
    ///        2   \
    ///       / \   \
    ///      / ↑ \   l
    ///     5  h  4   \
    ///    /   ↓   \   \
    ///   /         \   /
    ///  0-----3-----1
    ///
    ///  |--s--|--s--|
    ///
    ///  |-----l-----|
    /// ```
    fn equilateral_triangle() -> Result<(Shape, StateOfShape), StrError> {
        let l = 5.0;
        let s = l / 2.0;
        let h = l * SQRT_3 / 2.0;
        let (x0, y0) = (3.0, 4.0);
        let (x1, y1) = (x0 + l, y0);
        let (x2, y2) = (x0 + s, y0 + h);
        let (x3, y3) = (x0 + s, y0);
        let (x4, y4) = (x0 + 1.5 * s, y0 + 0.5 * h);
        let (x5, y5) = (x0 + 0.5 * s, y0 + 0.5 * h);
        let shape = Shape::new(2, 2, 6)?;
        let state = StateOfShape::new(
            shape.geo_ndim,
            &[[x0, y0], [x1, y1], [x2, y2], [x3, y3], [x4, y4], [x5, y5]],
        )?;
        Ok((shape, state))
    }

    #[test]
    fn approximate_ksi_works_outside() -> Result<(), StrError> {
        let (shape, mut state) = equilateral_triangle()?;
        assert_eq!(
            format!("{:.2}", state.coords_transp),
            "┌                               ┐\n\
             │ 3.00 8.00 5.50 5.50 6.75 4.25 │\n\
             │ 4.00 4.00 8.33 4.00 6.17 6.17 │\n\
             └                               ┘"
        );
        let mut ksi = vec![0.0; shape.geo_ndim];
        for (nit_correct, x_data, tol) in &[
            (0, [3.0, 4.0], 1e-15),
            (0, [8.0, 4.0], 1e-15),
            (1, [5.5, 8.33], 1e-14),
            (0, [5.5, 4.0], 1e-15),
            (1, [6.75, 6.17], 1e-14),
            (1, [4.25, 6.17], 1e-14),
            (2, [10.0, 10.0], 1e-14),
            (2, [-10.0, -10.0], 1e-13),
            (6, [100.0, 100.0], 1e-12),
        ] {
            let x = Vector::from(x_data);
            let nit = shape.approximate_ksi(&mut ksi, &mut state, &x, 30, *tol)?;
            let mut x_out = Vector::new(2);
            shape.calc_coords(&mut x_out, &mut state, &ksi)?;
            assert_vec_approx_eq!(x.as_data(), x_out.as_data(), *tol);
            assert_eq!(nit, *nit_correct);
        }
        Ok(())
    }

    // Holds arguments for numerical differentiation of N w.r.t x => G (gradient) matrix
    struct ArgsNumG<'a> {
        shape: &'a Shape,    // shape
        state: StateOfShape, // auxiliary (copy) state
        at_x: Vector,        // at x coord value
        x: Vector,           // temporary x coord
        ksi: Vec<f64>,       // temporary reference coord
        m: usize,            // node index from 0 to nnode
        j: usize,            // dimension index from 0 to space_ndim
    }

    // Computes Nᵐ(ξ(x)) with variable v := xⱼ
    fn aux_grad(v: f64, args: &mut ArgsNumG) -> f64 {
        copy_vector(&mut args.x, &args.at_x).unwrap();
        args.x[args.j] = v;
        args.shape
            .approximate_ksi(&mut args.ksi, &mut args.state, &args.x, 10, 1e-14)
            .unwrap();
        args.shape.calc_interp(&mut args.state, &args.ksi).unwrap();
        args.state.interp[args.m]
    }

    #[test]
    fn fn_gradient_works() -> Result<(), StrError> {
        // define dims and number of nodes
        let pairs = vec![(2, 6), (2, 4), (2, 8), (3, 10), (3, 8), (3, 20)];

        // define tolerances
        let mut tols = HashMap::new();
        tols.insert(GeoKind::Tri6, 1e-9);
        tols.insert(GeoKind::Qua4, 1e-11);
        tols.insert(GeoKind::Qua8, 1e-10);
        tols.insert(GeoKind::Tet10, 1e-9);
        tols.insert(GeoKind::Hex8, 1e-11);
        tols.insert(GeoKind::Hex20, 1e-10);

        // loop over shapes
        for (geo_ndim, nnode) in pairs {
            // allocate shape
            let space_ndim = usize::max(2, geo_ndim);
            let shape = Shape::new(space_ndim, geo_ndim, nnode)?;

            // generate state with coordinates matrix
            let mut state = gen_state_with_coords_matrix(&shape)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // set ξ within reference space
            let at_ksi = vec![0.25; shape.geo_ndim];

            // compute x corresponding to ξ using the isoparametric formula
            let mut at_x = Vector::new(shape.space_ndim);
            shape.calc_coords(&mut at_x, &mut state, &at_ksi)?;

            // compute gradient
            let det_jac = shape.calc_gradient(&mut state, &at_ksi)?;
            assert!(det_jac > 0.0);

            // set arguments for numerical integration
            let args = &mut ArgsNumG {
                shape: &shape,
                state: state.clone(),
                at_x,
                x: Vector::new(shape.space_ndim),
                ksi: vec![0.0; shape.geo_ndim],
                m: 0,
                j: 0,
            };

            // check Gᵐ(ξ(x)) = dNᵐ(ξ(x))/dx
            for m in 0..shape.nnode {
                args.m = m;
                for j in 0..shape.geo_ndim {
                    args.j = j;
                    // Gᵐⱼ := dNᵐ/dxⱼ
                    assert_deriv_approx_eq!(state.gradient[m][j], args.at_x[j], aux_grad, args, tol);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn calc_integ_points_coords() -> Result<(), StrError> {
        let shape = Shape::new(2, 1, 2)?;
        let (xa, xb) = (2.0, 5.0);
        let mut state = StateOfShape::new(shape.geo_ndim, &[[xa, 8.0], [xb, 8.0]])?;
        let integ_points = shape.calc_integ_points_coords(&mut state, &IP_LIN_LEGENDRE_2)?;
        assert_eq!(integ_points.len(), 2);
        let ksi_a = -1.0 / SQRT_3;
        let ksi_b = 1.0 / SQRT_3;
        let x_ksi_a = xa + (ksi_a + 1.0) * (xb - xa) / 2.0;
        let x_ksi_b = xa + (ksi_b + 1.0) * (xb - xa) / 2.0;
        assert_eq!(integ_points[0].dim(), 2);
        assert_eq!(integ_points[1].dim(), 2);
        assert_vec_approx_eq!(integ_points[0].as_data(), &[x_ksi_a, 8.0], 1e-15);
        assert_vec_approx_eq!(integ_points[1].as_data(), &[x_ksi_b, 8.0], 1e-15);
        Ok(())
    }
}
