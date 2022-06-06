use super::{GeoKind, StateOfShape};
use super::{
    Hex20, Hex32, Hex8, Lin2, Lin3, Lin4, Lin5, Qua12, Qua16, Qua17, Qua4, Qua8, Qua9, Tet10, Tet20, Tet4, Tri10,
    Tri15, Tri3, Tri6,
};
use crate::StrError;
use russell_lab::{inverse, mat_mat_mul, mat_vec_mul, vector_norm, Matrix, NormVec, Vector};
use std::fmt;

/// Defines an alias for interpolation functions
#[derive(Clone)]
struct ToDeleteFnInterp(fn(&mut Vector, &[f64]));

/// Defines an alias for derivative of interpolation functions
#[derive(Clone)]
struct ToDeleteFnDeriv(fn(&mut Matrix, &[f64]));

impl fmt::Debug for ToDeleteFnInterp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "FnInterp")
    }
}

impl fmt::Debug for ToDeleteFnDeriv {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "FnDeriv")
    }
}

/// Collects functions for interpolations and computing the derivatives
/// and gradients related to geometric shapes (elements)
///
/// # Warning
///
/// All public properties here are **readonly** and must **not** be modified externally.
#[derive(Clone, Debug)]
pub struct Shape {
    /// Geometry kind
    pub kind: GeoKind,

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
    fn_interp: ToDeleteFnInterp,

    /// Callback to evaluate local derivatives (with respect to ksi) of interpolation functions
    fn_deriv: ToDeleteFnDeriv,
}

impl Shape {
    /// Creates a new geometric shape
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::{GeoClass, GeoKind, Shape};
    /// use gemlab::StrError;
    ///
    /// fn main() {
    ///     //          .4--------------7
    ///     //        ,' |            ,'|         ξ₀   ξ₁   ξ₂
    ///     //      ,'              ,'  |  node    r    s    t
    ///     //    ,'     |        ,'    |     0 -1.0 -1.0 -1.0
    ///     //  5'==============6'      |     1  1.0 -1.0 -1.0
    ///     //  |               |       |     2  1.0  1.0 -1.0
    ///     //  |        |      |       |     3 -1.0  1.0 -1.0
    ///     //  |       ,0- - - | - - - 3     4 -1.0 -1.0  1.0
    ///     //  |     ,'        |     ,'      5  1.0 -1.0  1.0
    ///     //  |   ,'          |   ,'        6  1.0  1.0  1.0
    ///     //  | ,'            | ,'          7 -1.0  1.0  1.0
    ///     //  1'--------------2'
    ///     let shape = Shape::new(GeoKind::Hex8);
    ///     assert_eq!(shape.kind, GeoKind::Hex8);
    ///     assert_eq!(shape.geo_ndim, 3);
    ///     assert_eq!(shape.nnode, 8);
    ///     assert_eq!(shape.nedge, 12);
    ///     assert_eq!(shape.nface, 6);
    ///     assert_eq!(shape.edge_nnode, 2);
    ///     assert_eq!(shape.face_nnode, 4);
    ///     assert_eq!(shape.face_nedge, 4);
    /// }
    /// ```
    pub fn new(kind: GeoKind) -> Self {
        // collect geometry data
        let (geo_ndim, nnode, nedge, nface, edge_nnode, face_nnode, face_nedge, fn_interp, fn_deriv): (
            usize,
            usize,
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
                Lin2::GEO_NDIM,
                Lin2::NNODE,
                Lin2::NEDGE,
                Lin2::NFACE,
                Lin2::EDGE_NNODE,
                Lin2::FACE_NNODE,
                Lin2::FACE_NEDGE,
                Lin2::calc_interp,
                Lin2::calc_deriv,
            ),
            GeoKind::Lin3 => (
                Lin3::GEO_NDIM,
                Lin3::NNODE,
                Lin3::NEDGE,
                Lin3::NFACE,
                Lin3::EDGE_NNODE,
                Lin3::FACE_NNODE,
                Lin3::FACE_NEDGE,
                Lin3::calc_interp,
                Lin3::calc_deriv,
            ),
            GeoKind::Lin4 => (
                Lin4::GEO_NDIM,
                Lin4::NNODE,
                Lin4::NEDGE,
                Lin4::NFACE,
                Lin4::EDGE_NNODE,
                Lin4::FACE_NNODE,
                Lin4::FACE_NEDGE,
                Lin4::calc_interp,
                Lin4::calc_deriv,
            ),
            GeoKind::Lin5 => (
                Lin5::GEO_NDIM,
                Lin5::NNODE,
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
                Tri3::GEO_NDIM,
                Tri3::NNODE,
                Tri3::NEDGE,
                Tri3::NFACE,
                Tri3::EDGE_NNODE,
                Tri3::FACE_NNODE,
                Tri3::FACE_NEDGE,
                Tri3::calc_interp,
                Tri3::calc_deriv,
            ),
            GeoKind::Tri6 => (
                Tri6::GEO_NDIM,
                Tri6::NNODE,
                Tri6::NEDGE,
                Tri6::NFACE,
                Tri6::EDGE_NNODE,
                Tri6::FACE_NNODE,
                Tri6::FACE_NEDGE,
                Tri6::calc_interp,
                Tri6::calc_deriv,
            ),
            GeoKind::Tri10 => (
                Tri10::GEO_NDIM,
                Tri10::NNODE,
                Tri10::NEDGE,
                Tri10::NFACE,
                Tri10::EDGE_NNODE,
                Tri10::FACE_NNODE,
                Tri10::FACE_NEDGE,
                Tri10::calc_interp,
                Tri10::calc_deriv,
            ),
            GeoKind::Tri15 => (
                Tri15::GEO_NDIM,
                Tri15::NNODE,
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
                Qua4::GEO_NDIM,
                Qua4::NNODE,
                Qua4::NEDGE,
                Qua4::NFACE,
                Qua4::EDGE_NNODE,
                Qua4::FACE_NNODE,
                Qua4::FACE_NEDGE,
                Qua4::calc_interp,
                Qua4::calc_deriv,
            ),
            GeoKind::Qua8 => (
                Qua8::GEO_NDIM,
                Qua8::NNODE,
                Qua8::NEDGE,
                Qua8::NFACE,
                Qua8::EDGE_NNODE,
                Qua8::FACE_NNODE,
                Qua8::FACE_NEDGE,
                Qua8::calc_interp,
                Qua8::calc_deriv,
            ),
            GeoKind::Qua9 => (
                Qua9::GEO_NDIM,
                Qua9::NNODE,
                Qua9::NEDGE,
                Qua9::NFACE,
                Qua9::EDGE_NNODE,
                Qua9::FACE_NNODE,
                Qua9::FACE_NEDGE,
                Qua9::calc_interp,
                Qua9::calc_deriv,
            ),
            GeoKind::Qua12 => (
                Qua12::GEO_NDIM,
                Qua12::NNODE,
                Qua12::NEDGE,
                Qua12::NFACE,
                Qua12::EDGE_NNODE,
                Qua12::FACE_NNODE,
                Qua12::FACE_NEDGE,
                Qua12::calc_interp,
                Qua12::calc_deriv,
            ),
            GeoKind::Qua16 => (
                Qua16::GEO_NDIM,
                Qua16::NNODE,
                Qua16::NEDGE,
                Qua16::NFACE,
                Qua16::EDGE_NNODE,
                Qua16::FACE_NNODE,
                Qua16::FACE_NEDGE,
                Qua16::calc_interp,
                Qua16::calc_deriv,
            ),
            GeoKind::Qua17 => (
                Qua17::GEO_NDIM,
                Qua17::NNODE,
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
                Tet4::GEO_NDIM,
                Tet4::NNODE,
                Tet4::NEDGE,
                Tet4::NFACE,
                Tet4::EDGE_NNODE,
                Tet4::FACE_NNODE,
                Tet4::FACE_NEDGE,
                Tet4::calc_interp,
                Tet4::calc_deriv,
            ),
            GeoKind::Tet10 => (
                Tet10::GEO_NDIM,
                Tet10::NNODE,
                Tet10::NEDGE,
                Tet10::NFACE,
                Tet10::EDGE_NNODE,
                Tet10::FACE_NNODE,
                Tet10::FACE_NEDGE,
                Tet10::calc_interp,
                Tet10::calc_deriv,
            ),
            GeoKind::Tet20 => (
                Tet20::GEO_NDIM,
                Tet20::NNODE,
                Tet20::NEDGE,
                Tet20::NFACE,
                Tet20::EDGE_NNODE,
                Tet20::FACE_NNODE,
                Tet20::FACE_NEDGE,
                Tet20::calc_interp,
                Tet20::calc_deriv,
            ),

            // Hex
            GeoKind::Hex8 => (
                Hex8::GEO_NDIM,
                Hex8::NNODE,
                Hex8::NEDGE,
                Hex8::NFACE,
                Hex8::EDGE_NNODE,
                Hex8::FACE_NNODE,
                Hex8::FACE_NEDGE,
                Hex8::calc_interp,
                Hex8::calc_deriv,
            ),
            GeoKind::Hex20 => (
                Hex20::GEO_NDIM,
                Hex20::NNODE,
                Hex20::NEDGE,
                Hex20::NFACE,
                Hex20::EDGE_NNODE,
                Hex20::FACE_NNODE,
                Hex20::FACE_NEDGE,
                Hex20::calc_interp,
                Hex20::calc_deriv,
            ),
            GeoKind::Hex32 => (
                Hex32::GEO_NDIM,
                Hex32::NNODE,
                Hex32::NEDGE,
                Hex32::NFACE,
                Hex32::EDGE_NNODE,
                Hex32::FACE_NNODE,
                Hex32::FACE_NEDGE,
                Hex32::calc_interp,
                Hex32::calc_deriv,
            ),
        };

        // return new Shape
        Shape {
            kind,
            geo_ndim,
            nnode,
            nedge,
            nface,
            edge_nnode,
            face_nnode,
            face_nedge,
            fn_interp: ToDeleteFnInterp(fn_interp),
            fn_deriv: ToDeleteFnDeriv(fn_deriv),
        }
    }

    /// Calculates the interpolation functions
    ///
    /// Computes Nᵐ used in interpolations such as:
    ///
    /// ```text
    /// → →         →  →
    /// u(ξ) = Σ Nᵐ(ξ) uᵐ
    ///        m         
    /// ```
    ///
    /// # Updated
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
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     //  3-------------2         ξ₀   ξ₁
    ///     //  |      ξ₁     |  node    r    s
    ///     //  |      |      |     0 -1.0 -1.0
    ///     //  |      +--ξ₀  |     1  1.0 -1.0
    ///     //  |             |     2  1.0  1.0
    ///     //  |             |     3 -1.0  1.0
    ///     //  0-------------1
    ///     let coords = &[
    ///         [0.0, 0.0],
    ///         [1.0, 0.0],
    ///         [1.0, 1.0],
    ///         [0.0, 1.0],
    ///     ];
    ///     let shape = Shape::new(GeoKind::Qua4);
    ///     let mut state = StateOfShape::new(shape.kind, coords)?;
    ///
    ///     shape.calc_interp(&mut state, &[-1.0, -1.0])?;
    ///     assert_eq!(state.interp.as_data(), &[1.0, 0.0, 0.0, 0.0]);
    ///
    ///     shape.calc_interp(&mut state, &[1.0, -1.0])?;
    ///     assert_eq!(state.interp.as_data(), &[0.0, 1.0, 0.0, 0.0]);
    ///
    ///     shape.calc_interp(&mut state, &[1.0, 1.0])?;
    ///     assert_eq!(state.interp.as_data(), &[0.0, 0.0, 1.0, 0.0]);
    ///
    ///     shape.calc_interp(&mut state, &[-1.0, 1.0])?;
    ///     assert_eq!(state.interp.as_data(), &[0.0, 0.0, 0.0, 1.0]);
    ///
    ///     shape.calc_interp(&mut state, &[0.0, 0.0])?;
    ///     assert_eq!(state.interp.as_data(), &[0.25, 0.25, 0.25, 0.25]);
    ///     Ok(())
    /// }
    /// ```
    #[inline]
    pub fn calc_interp(&self, state: &mut StateOfShape, ksi: &[f64]) -> Result<(), StrError> {
        if ksi.len() < self.geo_ndim {
            return Err("ksi.len() must be ≥ geo_ndim");
        }
        (self.fn_interp.0)(&mut state.interp, ksi);
        Ok(())
    }

    /// Calculates the derivatives of interpolation functions with respect to reference coordinate
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
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     //  3-------------2         ξ₀   ξ₁
    ///     //  |      ξ₁     |  node    r    s
    ///     //  |      |      |     0 -1.0 -1.0
    ///     //  |      +--ξ₀  |     1  1.0 -1.0
    ///     //  |             |     2  1.0  1.0
    ///     //  |             |     3 -1.0  1.0
    ///     //  0-------------1
    ///     let coords = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    ///     let shape = Shape::new(GeoKind::Qua4);
    ///     let mut state = StateOfShape::new(shape.kind, coords)?;
    ///
    ///     shape.calc_deriv(&mut state, &[-1.0, -1.0])?;
    ///     assert_eq!(
    ///         format!("{}", state.deriv),
    ///         "┌           ┐\n\
    ///          │ -0.5 -0.5 │\n\
    ///          │  0.5    0 │\n\
    ///          │    0    0 │\n\
    ///          │    0  0.5 │\n\
    ///          └           ┘"
    ///     );
    ///
    ///     shape.calc_deriv(&mut state, &[1.0, 1.0])?;
    ///     assert_eq!(
    ///         format!("{}", state.deriv),
    ///         "┌           ┐\n\
    ///          │    0    0 │\n\
    ///          │    0 -0.5 │\n\
    ///          │  0.5  0.5 │\n\
    ///          │ -0.5    0 │\n\
    ///          └           ┘"
    ///     );
    ///     Ok(())
    /// }
    /// ```
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
    ///
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
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
    /// use gemlab::StrError;
    /// use russell_chk::assert_vec_approx_eq;
    /// use russell_lab::Vector;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     //  3-------------2         ξ₀   ξ₁
    ///     //  |      ξ₁     |  node    r    s
    ///     //  |      |      |     0 -1.0 -1.0
    ///     //  |      +--ξ₀  |     1  1.0 -1.0
    ///     //  |             |     2  1.0  1.0
    ///     //  |             |     3 -1.0  1.0
    ///     //  0-------------1
    ///     let (x0, y0) = (10.0, 10.0);
    ///     let (w, h) = (10.0, 5.0);
    ///     let coords = &[
    ///         [10.0,  5.0],
    ///         [20.0,  5.0],
    ///         [20.0, 10.0],
    ///         [10.0, 10.0],
    ///     ];
    ///     let shape = Shape::new(GeoKind::Qua4);
    ///     let mut state = StateOfShape::new(shape.kind, coords)?;
    ///
    ///     let mut x = Vector::new(2);
    ///     shape.calc_coords(&mut x, &mut state, &[0.0, 0.0])?;
    ///     assert_vec_approx_eq!(x.as_data(), &[15.0, 7.5], 1e-15);
    ///     Ok(())
    /// }
    /// ```
    pub fn calc_coords(&self, x: &mut Vector, state: &mut StateOfShape, ksi: &[f64]) -> Result<(), StrError> {
        if x.dim() != state.coords_min.len() {
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
    /// Thus, in matrix notation
    ///
    /// ```text
    /// jacobian := J = Xᵀ · L
    /// jacobian := Jline = Xᵀ · L
    /// jacobian := Jsurf = Xᵀ · L
    /// ```
    ///
    /// where:
    ///
    /// * `Jline`` -- Jacobian for line in multi-dimensions (geom_ndim < space_ndim)
    /// * `Jsurf`` -- Jacobian for 3D surfaces (geo_ndim = 2 and space_ndim = 3)
    ///
    /// If `geo_ndim = space_ndim`, we also compute the inverse Jacobian
    ///
    /// ```text
    /// inv_jacobian := J⁻¹
    /// ```
    ///
    /// # Output
    ///
    /// * `state.deriv` -- derivatives of the interpolation functions (nnode); `L` matrix
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
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
    /// use gemlab::StrError;
    /// use russell_chk::assert_approx_eq;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     //  3-------------2         ξ₀   ξ₁
    ///     //  |      ξ₁     |  node    r    s
    ///     //  |      |      |     0 -1.0 -1.0
    ///     //  |      +--ξ₀  |     1  1.0 -1.0
    ///     //  |             |     2  1.0  1.0
    ///     //  |             |     3 -1.0  1.0
    ///     //  0-------------1
    ///     let a = 3.0;
    ///     let coords = &[[0.0, 0.0], [2.0 * a, 0.0], [2.0 * a, a], [0.0, a]];
    ///     let shape = Shape::new(GeoKind::Qua4);
    ///     let mut state = StateOfShape::new(shape.kind, coords)?;
    ///
    ///     let det_jj = shape.calc_jacobian(&mut state, &[0.0, 0.0])?;
    ///     assert_approx_eq!(det_jj, a * a / 2.0, 1e-15);
    ///
    ///     // the solution is
    ///     //  ┌         ┐
    ///     //  │  a   0  │
    ///     //  │  0  a/2 │
    ///     //  └         ┘
    ///     assert_eq!(
    ///         format!("{}", state.jacobian),
    ///         "┌         ┐\n\
    ///          │   3   0 │\n\
    ///          │   0 1.5 │\n\
    ///          └         ┘"
    ///     );
    ///     Ok(())
    /// }
    /// ```
    pub fn calc_jacobian(&self, state: &mut StateOfShape, ksi: &[f64]) -> Result<f64, StrError> {
        self.calc_deriv(state, ksi)?;
        mat_mat_mul(&mut state.jacobian, 1.0, &state.coords_transp, &state.deriv)?;
        let space_ndim = state.coords_min.len();
        if self.geo_ndim == space_ndim {
            inverse(&mut state.inv_jacobian, &state.jacobian)
        } else {
            if self.geo_ndim == 1 {
                let mut norm_jac = 0.0;
                for i in 0..space_ndim {
                    norm_jac += state.jacobian[i][0] * state.jacobian[i][0];
                }
                return Ok(f64::sqrt(norm_jac));
            }
            Ok(0.0)
        }
    }

    /// Computes the boundary normal vector
    ///
    /// **Important:** This function only works with:
    ///
    /// * `geo_ndim = 1` and `space_ndim = 2` -- line in 2D, or
    /// * `geo_ndim = 2` and `space_ndim = 3` -- surface in 3D.
    ///
    /// i.e., `geo_ndim < space_ndim`.
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
    ///
    /// # Examples
    ///
    /// ## Line in multi-dimensions (geo_ndim = 1 and space_ndim > 1)
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
    /// use gemlab::StrError;
    /// use russell_lab::Vector;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     //                 .
    ///     //                /|\  →
    ///     //                 |   n
    ///     //                 |
    ///     //  0----+----2----+----1
    ///     const L: f64 = 5.0;
    ///     let coords = &[[0.0, 0.0], [L, 0.0], [L / 2.0, 0.0]];
    ///     let shape = Shape::new(GeoKind::Lin3);
    ///     let mut state = StateOfShape::new(shape.kind, coords)?;
    ///     let mut normal = Vector::new(2);
    ///     shape.calc_boundary_normal(&mut normal, &mut state, &[0.5, 0.0])?;
    ///     assert_eq!(normal.as_data(), &[0.0, L / 2.0]);
    ///     Ok(())
    /// }
    /// ```
    ///
    /// ## Boundary surface (geo_ndim = 2 and space_ndim = 3)
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
    /// use gemlab::StrError;
    /// use russell_lab::Vector;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     //           .   .  .   . ,.2|
    ///     //         ' .           ,,'||
    ///     //       '   .         ,,'  ||
    ///     //     '     .       .,'    ||  →
    ///     //  .  .   . .   .  3'      ||  n
    ///     //           z     ||   ==========)
    ///     //  .        |     ||       ||
    ///     //          ,*---y || .  . ,1
    ///     //  .      x       ||    ,,'
    ///     //      ,'         ||  ,,'
    ///     //  . ,'           ||,,'
    ///     //  . . .   .   .  |0'
    ///     let coords = &[
    ///         [1.0, 1.0, 0.0],
    ///         [0.0, 1.0, 0.0],
    ///         [0.0, 1.0, 1.0],
    ///         [1.0, 1.0, 1.0],
    ///     ];
    ///     let shape = Shape::new(GeoKind::Qua4);
    ///     let mut state = StateOfShape::new(shape.kind, coords)?;
    ///     let mut normal = Vector::new(3);
    ///     shape.calc_boundary_normal(&mut normal, &mut state, &[0.0, 0.0, 0.0])?;
    ///     const A: f64 = 1.0;
    ///     assert_eq!(normal.as_data(), &[0.0, A / 4.0, 0.0]);
    ///     Ok(())
    /// }
    /// ```
    pub fn calc_boundary_normal(
        &self,
        normal: &mut Vector,
        state: &mut StateOfShape,
        ksi: &[f64],
    ) -> Result<(), StrError> {
        // check
        let space_ndim = state.coords_min.len();
        if self.geo_ndim >= space_ndim {
            return Err("geo_ndim must be smaller than space_ndim");
        }
        if normal.dim() != space_ndim {
            return Err("normal.dim() must equal space_ndim");
        }

        // compute Jacobian
        self.calc_deriv(state, ksi)?;
        mat_mat_mul(&mut state.jacobian, 1.0, &state.coords_transp, &state.deriv)?;

        // line in 2D (geo_ndim = 1 and self.space_ndim = 2)
        if space_ndim == 2 {
            normal[0] = -state.jacobian[1][0];
            normal[1] = state.jacobian[0][0];
            return Ok(());
        }

        // surface in 3D (geo_ndim = 2 and space_ndim = 3)
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
    /// We use Newton iterations with the inverse of the Jacobian to compute `ξ(x)`.
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
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
    /// use gemlab::StrError;
    /// use russell_chk::assert_vec_approx_eq;
    /// use russell_lab::Vector;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     // 7.0        2                ξ₀   ξ₁
    ///     //           / `.       node    r    s
    ///     //          /    `.        0  0.0  0.0
    ///     //     (3.5,6.0)   `.      1  1.0  0.0
    ///     //        /          `.    2  0.0  1.0
    ///     //       /             `.
    ///     // 5.0  0-----------------1
    ///     //     3.0   4.0   5.0   6.0
    ///     #[rustfmt::skip]
    ///     let coords = &[
    ///         [3.0, 5.0],
    ///         [6.0, 5.0],
    ///         [4.0, 7.0],
    ///     ];
    ///     let shape = Shape::new(GeoKind::Tri3);
    ///     let mut state = StateOfShape::new(shape.kind, coords)?;
    ///
    ///     // x @ middle of edge (0,2)
    ///     let x = Vector::from(&[3.5, 6.0]);
    ///
    ///     // find ξ corresponding to x @ middle of edge (0,2)
    ///     let mut ksi = vec![0.0; 2];
    ///     shape.approximate_ksi(&mut ksi, &mut state, &x, 10, 1e-8)?;
    ///     assert_vec_approx_eq!(ksi, &[0.0, 0.5], 1e-8);
    ///     Ok(())
    /// }
    /// ```
    pub fn approximate_ksi(
        &self,
        ksi: &mut [f64],
        state: &mut StateOfShape,
        x: &Vector,
        nit_max: usize,
        tol: f64,
    ) -> Result<usize, StrError> {
        // check
        let space_ndim = state.coords_min.len();
        if self.geo_ndim != space_ndim {
            return Err("geo_ndim must equal space_ndim");
        }
        if x.dim() != space_ndim {
            return Err("x.dim() must equal space_ndim");
        }
        if ksi.len() != self.geo_ndim {
            return Err("ksi.len() must equal geo_ndim");
        }

        // use linear scale to guess ksi
        let (min_ksi, _, del_ksi) = self.kind.reference_limits();
        for j in 0..self.geo_ndim {
            ksi[j] = (x[j] - state.coords_min[j]) / (state.coords_max[j] - state.coords_min[j]) * del_ksi + min_ksi;
        }

        // perform iterations
        let mut residual = Vector::new(space_ndim);
        let mut x_at_ksi = Vector::new(space_ndim);
        let mut delta_ksi = Vector::new(self.geo_ndim);
        for it in 0..nit_max {
            self.calc_coords(&mut x_at_ksi, state, &ksi)?;
            for i in 0..space_ndim {
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
    /// **Note:** This function works with `geo_ndim == space_ndim` only.
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
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
    /// use gemlab::StrError;
    /// use russell_chk::assert_vec_approx_eq;
    /// use russell_lab::Matrix;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     //  3-------------2         ξ₀   ξ₁
    ///     //  |      ξ₁     |  node    r    s
    ///     //  |      |      |     0 -1.0 -1.0
    ///     //  |      +--ξ₀  |     1  1.0 -1.0
    ///     //  |             |     2  1.0  1.0
    ///     //  |             |     3 -1.0  1.0
    ///     //  0-------------1
    ///     let a = 3.0;
    ///     let coords = &[[0.0, 0.0], [2.0 * a, 0.0], [2.0 * a, a], [0.0, a]];
    ///     let shape = Shape::new(GeoKind::Qua4);
    ///     let mut state = StateOfShape::new(shape.kind, coords)?;
    ///
    ///     shape.calc_gradient(&mut state, &[0.0, 0.0])?;
    ///
    ///     let correct_gg = Matrix::from(&[
    ///         [-1.0/(4.0*a), -1.0/(2.0*a)],
    ///         [ 1.0/(4.0*a), -1.0/(2.0*a)],
    ///         [ 1.0/(4.0*a),  1.0/(2.0*a)],
    ///         [-1.0/(4.0*a),  1.0/(2.0*a)],
    ///     ]);
    ///     assert_vec_approx_eq!(state.gradient.as_data(), correct_gg.as_data(), 1e-15);
    ///     Ok(())
    /// }
    /// ```
    pub fn calc_gradient(&self, state: &mut StateOfShape, ksi: &[f64]) -> Result<f64, StrError> {
        let space_ndim = state.coords_min.len();
        if self.geo_ndim != space_ndim {
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
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     //  6 2
    ///     //  5 | `.    * indicates the
    ///     //  4 | * `.    location of ips
    ///     //  3 |     `.
    ///     //  2 |       `.
    ///     //  1 | *     * `.
    ///     //  0 0-----------1
    ///     //    0 1 2 3 4 5 6
    ///     #[rustfmt::skip]
    ///     let coords = &[
    ///         [0.0, 0.0],
    ///         [6.0, 0.0],
    ///         [0.0, 6.0],
    ///     ];
    ///     let shape = Shape::new(GeoKind::Tri3);
    ///     let mut state = StateOfShape::new(shape.kind, coords)?;
    ///
    ///     const IP_TRI_INTERNAL_3: [[f64; 4]; 3] = [
    ///         [1.0/6.0, 1.0/6.0, 0.0, 1.0/6.0], // last column
    ///         [2.0/3.0, 1.0/6.0, 0.0, 1.0/6.0], // is the weight
    ///         [1.0/6.0, 2.0/3.0, 0.0, 1.0/6.0],
    ///     ];
    ///
    ///     let x_ips = shape.calc_integ_points_coords(&mut state, &IP_TRI_INTERNAL_3)?;
    ///     assert_eq!(x_ips[0].as_data(), &[1.0, 1.0]);
    ///     assert_eq!(x_ips[1].as_data(), &[4.0, 1.0]);
    ///     assert_eq!(x_ips[2].as_data(), &[1.0, 4.0]);
    ///     Ok(())
    /// }
    /// ```
    pub fn calc_integ_points_coords(
        &self,
        state: &mut StateOfShape,
        integ_points: &[[f64; 4]],
    ) -> Result<Vec<Vector>, StrError> {
        let space_ndim = state.coords_min.len();
        let mut all_coords = Vec::new();
        for iota in integ_points {
            let mut x = Vector::new(space_ndim);
            self.calc_coords(&mut x, state, iota)?;
            all_coords.push(x);
        }
        Ok(all_coords)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Shape;
    use crate::shapes::{GeoKind, StateOfShape};
    use crate::util::{PI, SQRT_3};
    use crate::StrError;
    use russell_chk::{assert_approx_eq, assert_deriv_approx_eq, assert_vec_approx_eq};
    use russell_lab::{copy_vector, scale_vector, vector_norm, Matrix, NormVec, Vector};
    use std::collections::HashMap;

    #[test]
    fn derive_works() {
        let shape = Shape::new(GeoKind::Lin2);
        let shape_clone = shape.clone();
        assert_eq!(format!("{:?}", shape), "Shape { kind: Lin2, geo_ndim: 1, nnode: 2, nedge: 0, nface: 0, edge_nnode: 0, face_nnode: 0, face_nedge: 0, fn_interp: FnInterp, fn_deriv: FnDeriv }");
        assert_eq!(shape_clone.kind, shape.kind);
        assert_eq!(shape_clone.geo_ndim, shape.geo_ndim);
        assert_eq!(shape_clone.nnode, shape.nnode);
    }

    #[test]
    fn new_works() -> Result<(), StrError> {
        let nnode = [
            2, 3, 4, 5, // Lin
            3, 6, 10, 15, // Tri
            4, 8, 9, 12, 16, 17, // Qua
            4, 10, 20, // Tet
            8, 20, 32, // Hex
        ];
        let geo_ndim = [
            1, 1, 1, 1, // Lin
            2, 2, 2, 2, // Tri
            2, 2, 2, 2, 2, 2, // Qua
            3, 3, 3, // Tet
            3, 3, 3, // Hex
        ];
        for i in 0..GeoKind::VALUES.len() {
            let shape = Shape::new(GeoKind::VALUES[i]);
            assert_eq!(shape.nnode, nnode[i]);
            assert_eq!(shape.geo_ndim, geo_ndim[i]);
        }
        Ok(())
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
    fn gen_coords(x: &mut Vector, ksi: &[f64], kind: GeoKind) {
        assert_eq!(x.dim(), ksi.len());
        let (min_ksi, _, del_ksi) = kind.reference_limits();
        let r = RMIN + (ksi[0] - min_ksi) * (RMAX - RMIN) / del_ksi;
        let a = AMIN + (ksi[1] - min_ksi) * (AMAX - AMIN) / del_ksi;
        x[0] = r * f64::cos(a);
        x[1] = r * f64::sin(a);
        if x.dim() == 3 {
            x[2] = ksi[2];
        }
    }

    // Generates state with a pre-set matrix of coordinates
    fn gen_state_with_coords_matrix(space_ndim: usize, shape: &Shape) -> Result<StateOfShape, StrError> {
        let mut x = Vector::new(space_ndim);
        let mut ksi_aux = vec![0.0; space_ndim];
        let mut coords: Vec<Vec<f64>> = Vec::new();
        for m in 0..shape.nnode {
            let ksi = shape.kind.reference_coords(m);
            if shape.geo_ndim == space_ndim {
                gen_coords(&mut x, ksi, shape.kind);
            } else if shape.geo_ndim == 1 && space_ndim == 2 {
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = 1.0;
                gen_coords(&mut x, &ksi_aux, shape.kind);
            } else {
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = ksi[1];
                ksi_aux[2] = 1.0;
                gen_coords(&mut x, &ksi_aux, shape.kind);
            }
            coords.push(x.as_data().clone());
        }
        StateOfShape::new(shape.kind, &coords)
    }

    // Generates state with a pre-set matrix of coordinates such that the
    // shape is parallel to the x,y,z axes.
    // For triangles and tetrahedra, one edge/face will cut the axes equally
    // * Qua and Hex will have real coordinates equal to the natural coordinates
    // * Tri and Tet will be scaled using the natural coordinates
    // * All edges parallel to the x,y,z axes will have lengths equal to 2.0
    fn gen_state_with_coords_matrix_parallel(shape: &Shape) -> Result<StateOfShape, StrError> {
        let scale = if shape.kind.is_tri_or_tet() { 2.0 } else { 1.0 };
        let mut coords: Vec<Vec<f64>> = Vec::new();
        for m in 0..shape.nnode {
            let ksi = shape.kind.reference_coords(m);
            coords.push(ksi.iter().map(|value| scale * value).collect());
        }
        StateOfShape::new(shape.kind, &coords)
    }

    #[test]
    fn capture_some_wrong_input() {
        // shape and state
        let shape = Shape::new(GeoKind::Lin2);
        let mut state = StateOfShape::new(shape.kind, &[[0.0, 0.0], [1.0, 1.0]]).unwrap();

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
        let shape = Shape::new(GeoKind::Tet4);
        let mut normal = Vector::new(2);
        let ksi = &[0.0, 0.0, 0.0];
        assert_eq!(
            shape.calc_boundary_normal(&mut normal, &mut state, ksi).err(),
            Some("geo_ndim must be smaller than space_ndim")
        );
        let shape = Shape::new(GeoKind::Lin2);
        let mut normal = Vector::new(1);
        let mut state = StateOfShape::new(shape.kind, &[[0.0, 0.0], [1.0, 1.0]]).unwrap();
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
        let shape = Shape::new(GeoKind::Lin2);
        let mut state = StateOfShape::new(shape.kind, &[[0.0, 0.0], [1.0, 1.0]]).unwrap();
        let mut ksi = vec![0.0; 1];
        let x = Vector::new(2);
        assert_eq!(
            shape.approximate_ksi(&mut ksi, &mut state, &x, 2, 1e-5).err(),
            Some("geo_ndim must equal space_ndim")
        );
        let shape = Shape::new(GeoKind::Tri3);
        let mut state = StateOfShape::new(shape.kind, &[[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]]).unwrap();
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
        let shape = Shape::new(GeoKind::Lin2);
        let mut state = StateOfShape::new(shape.kind, &[[0.0, 0.0], [1.0, 1.0]]).unwrap();
        assert_eq!(
            shape.calc_gradient(&mut state, &[0.0]).err(),
            Some("geo_ndim must equal space_ndim")
        );

        // calc_integ_points_coords
        assert_eq!(
            shape
                .calc_integ_points_coords(&mut bug_for_calc_coords, &[[0.0, 0.0, 0.0, 2.0]])
                .err(),
            Some("matrix and vectors are incompatible")
        );
    }

    #[test]
    fn calc_interp_works() -> Result<(), StrError> {
        // define tolerances
        let tols = HashMap::from([
            (GeoKind::Lin2, 1e-15),
            (GeoKind::Lin3, 1e-15),
            (GeoKind::Lin4, 1e-15),
            (GeoKind::Lin5, 1e-15),
            (GeoKind::Tri3, 1e-15),
            (GeoKind::Tri6, 1e-15),
            (GeoKind::Tri10, 1e-15),
            (GeoKind::Tri15, 1e-15),
            (GeoKind::Qua4, 1e-15),
            (GeoKind::Qua8, 1e-15),
            (GeoKind::Qua9, 1e-15),
            (GeoKind::Qua12, 1e-15),
            (GeoKind::Qua16, 1e-15),
            (GeoKind::Qua17, 1e-15),
            (GeoKind::Tet4, 1e-15),
            (GeoKind::Tet10, 1e-15),
            (GeoKind::Tet20, 1e-15),
            (GeoKind::Hex8, 1e-15),
            (GeoKind::Hex20, 1e-15),
            (GeoKind::Hex32, 1e-15),
        ]);

        // loop over shapes
        for kind in GeoKind::VALUES {
            // allocate shape and state
            let shape = Shape::new(kind);
            let space_ndim = usize::max(2, shape.geo_ndim);
            let mut state = StateOfShape::new(
                shape.kind,
                &(0..shape.nnode).map(|_| vec![0.0; space_ndim]).collect::<Vec<_>>(),
            )?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // loop over nodes of shape
            for m in 0..shape.nnode {
                // get ξᵐ corresponding to node m
                let ksi = shape.kind.reference_coords(m);

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

    // Holds arguments for numerical differentiation of N with respect to ξ => L (deriv) matrix
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
        let tols = HashMap::from([
            (GeoKind::Lin2, 1e-13),
            (GeoKind::Lin3, 1e-13),
            (GeoKind::Lin4, 1e-10),
            (GeoKind::Lin5, 1e-10),
            (GeoKind::Tri3, 1e-12),
            (GeoKind::Tri6, 1e-12),
            (GeoKind::Tri10, 1e-10),
            (GeoKind::Tri15, 1e-9),
            (GeoKind::Qua4, 1e-13),
            (GeoKind::Qua8, 1e-12),
            (GeoKind::Qua9, 1e-13),
            (GeoKind::Qua12, 1e-10),
            (GeoKind::Qua16, 1e-10),
            (GeoKind::Qua17, 1e-10),
            (GeoKind::Hex8, 1e-13),
            (GeoKind::Tet4, 1e-12),
            (GeoKind::Tet10, 1e-12),
            (GeoKind::Tet20, 1e-10),
            (GeoKind::Hex20, 1e-12),
            (GeoKind::Hex32, 1e-10),
        ]);

        // loop over shapes
        for kind in GeoKind::VALUES {
            // allocate shape and state
            let shape = Shape::new(kind);
            let space_ndim = usize::max(2, shape.geo_ndim);
            let mut state = StateOfShape::new(
                shape.kind,
                &(0..shape.nnode).map(|_| vec![0.0; space_ndim]).collect::<Vec<_>>(),
            )?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // set ξ within reference space
            let at_ksi = vec![0.25; shape.geo_ndim];

            // compute all derivatives of interpolation functions with respect to ξ
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
        // define kinds
        let kinds = vec![
            GeoKind::Tri3,
            GeoKind::Tri6,
            GeoKind::Tri10,
            GeoKind::Tri15,
            GeoKind::Qua4,
            GeoKind::Qua8,
            GeoKind::Qua17,
            GeoKind::Tet4,
            GeoKind::Tet10,
            GeoKind::Tet20,
            GeoKind::Hex8,
            GeoKind::Hex20,
            GeoKind::Hex32,
        ];

        // define tolerances
        let tols = HashMap::from([
            (GeoKind::Qua4, 1e-15),
            (GeoKind::Tri3, 1e-15),
            (GeoKind::Tri6, 1e-15),
            (GeoKind::Tri10, 1e-14),
            (GeoKind::Tri15, 1e-15),
            (GeoKind::Qua8, 1e-15),
            (GeoKind::Qua17, 1e-15),
            (GeoKind::Hex8, 1e-15),
            (GeoKind::Tet4, 1e-15),
            (GeoKind::Tet10, 1e-15),
            (GeoKind::Tet20, 1e-14),
            (GeoKind::Hex20, 1e-15),
            (GeoKind::Hex32, 1e-15),
        ]);

        // define tolerances for node in the reference domain
        let tols_in = HashMap::from([
            (GeoKind::Tri3, 0.45), // linear maps are inaccurate for the circular wedge
            (GeoKind::Tri6, 0.02), // << quadratic mapping is inaccurate as well
            (GeoKind::Tri10, 1e-14),
            (GeoKind::Tri15, 1e-4), // << this triangle is inaccurate as well here
            (GeoKind::Qua4, 0.14),  // linear maps are inaccurate for the circular wedge
            (GeoKind::Qua8, 1e-14),
            (GeoKind::Qua17, 1e-14),
            (GeoKind::Tet4, 0.45),   // linear tetrahedron is also inaccurate here
            (GeoKind::Tet10, 0.02),  // quadratic tetrahedron is also inaccurate here
            (GeoKind::Tet20, 1e-14), // cubic tetrahedron
            (GeoKind::Hex8, 0.14),   // bi-linear maps are inaccurate for the circular wedge
            (GeoKind::Hex20, 1e-14),
            (GeoKind::Hex32, 1e-4), // TODO: check why this tolerance is high
        ]);

        // loop over shapes
        for kind in kinds {
            // allocate shape
            let shape = Shape::new(kind);

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();
            let tol_in = *tols_in.get(&shape.kind).unwrap();

            // generate state with coordinates matrix
            let space_ndim = usize::max(2, shape.geo_ndim);
            let mut state = gen_state_with_coords_matrix(space_ndim, &shape)?;

            // loop over nodes of shape
            let mut x = Vector::new(space_ndim);
            let mut x_correct = Vector::new(space_ndim);
            for m in 0..shape.nnode {
                // get ξᵐ corresponding to node m
                let ksi = shape.kind.reference_coords(m);

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                shape.calc_coords(&mut x, &mut state, ksi)?;

                // compare xᵐ with generated coordinates
                gen_coords(&mut x_correct, ksi, shape.kind);
                assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol);
            }

            // test again inside the reference domain
            let ksi_in = if shape.kind.is_tri_or_tet() {
                vec![1.0 / 3.0; shape.geo_ndim]
            } else {
                vec![0.0; shape.geo_ndim]
            };
            shape.calc_coords(&mut x, &mut state, &ksi_in)?;
            gen_coords(&mut x_correct, &ksi_in, shape.kind);
            assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol_in);
        }
        Ok(())
    }

    // Holds arguments for numerical differentiation of x with respect to ξ => Jacobian
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
        let tols = HashMap::from([
            (GeoKind::Lin2, 1e-11),
            (GeoKind::Lin3, 1e-11),
            (GeoKind::Lin4, 1e-11),
            (GeoKind::Lin5, 1e-11),
            (GeoKind::Tri3, 1e-11),
            (GeoKind::Tri6, 1e-11),
            (GeoKind::Tri10, 1e-11),
            (GeoKind::Tri15, 1e-10),
            (GeoKind::Qua4, 1e-11),
            (GeoKind::Qua8, 1e-11),
            (GeoKind::Qua9, 1e-11),
            (GeoKind::Qua12, 1e-10),
            (GeoKind::Qua16, 1e-10),
            (GeoKind::Qua17, 1e-10),
            (GeoKind::Hex8, 1e-11),
            (GeoKind::Hex20, 1e-11),
            (GeoKind::Hex32, 1e-10),
            (GeoKind::Tet4, 1e-12),
            (GeoKind::Tet10, 1e-12),
            (GeoKind::Tet20, 1e-10),
        ]);

        // loop over shapes
        for kind in GeoKind::VALUES {
            // allocate shape
            let shape = Shape::new(kind);

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // generate state with coordinates matrix
            let space_ndim = usize::max(2, shape.geo_ndim);
            let mut state = gen_state_with_coords_matrix(space_ndim, &shape)?;

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
                x: Vector::new(space_ndim),
                i: 0,
                j: 0,
            };

            // check J(ξ) = dx(ξ)/dξ
            for i in 0..space_ndim {
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
        let shape = Shape::new(GeoKind::Lin2);
        let l = 3.5;
        let mut state = StateOfShape::new(shape.kind, &[[0.0, 0.0], [l, 0.0]])?;
        let norm_jac_vec = shape.calc_jacobian(&mut state, &[0.0])?;
        assert_eq!(norm_jac_vec, l / 2.0); // 2.0 = length of shape in the reference space

        let shape = Shape::new(GeoKind::Tri3);
        let mut state = StateOfShape::new(shape.kind, &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 1.0]])?;
        let norm_jac_vec = shape.calc_jacobian(&mut state, &[0.0, 0.0])?;
        assert_eq!(norm_jac_vec, 0.0);
        Ok(())
    }

    #[test]
    fn calc_boundary_normal_works() -> Result<(), StrError> {
        // allocate surface shape and state
        let shape = Shape::new(GeoKind::Qua17);

        // generate state with coordinates matrix
        let space_ndim = 3;
        let mut state = gen_state_with_coords_matrix(space_ndim, &shape)?;

        // compute boundary normal vector
        let at_ksi = vec![0.0; shape.geo_ndim];
        let mut normal = Vector::new(space_ndim);
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
        let shape = Shape::new(GeoKind::Lin5);

        // generate state with coordinates matrix
        let space_ndim = 2;
        let mut state = gen_state_with_coords_matrix(space_ndim, &shape)?;

        // compute boundary normal vector
        let at_ksi = vec![0.0; shape.geo_ndim];
        let mut normal = Vector::new(space_ndim);
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
        // select Tri and Qua
        let kinds = vec![
            // Tri
            GeoKind::Tri3,
            GeoKind::Tri6,
            GeoKind::Tri10,
            GeoKind::Tri15,
            // Qua
            GeoKind::Qua4,
            GeoKind::Qua8,
            GeoKind::Qua9,
            GeoKind::Qua12,
            GeoKind::Qua16,
            GeoKind::Qua17,
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
        let space_ndim = 2;
        for kind in kinds {
            // allocate shape
            let shape = Shape::new(kind);

            // generate state with coordinates matrix
            let state = gen_state_with_coords_matrix_parallel(&shape)?;

            // loop over edges
            for e in 0..shape.nedge {
                let edge_shape = Shape::new(shape.kind.edge_kind().unwrap());
                let mut edge_state = StateOfShape::new(
                    edge_shape.kind,
                    &(0..edge_shape.nnode)
                        .map(|i| {
                            let m = shape.kind.edge_node_id(e, i);
                            (0..space_ndim).map(|j| state.coords_transp[j][m]).collect()
                        })
                        .collect::<Vec<_>>(),
                )?;

                // calc normal vector
                edge_shape.calc_boundary_normal(&mut normal, &mut edge_state, ksi)?;
                if shape.kind.is_tri_or_tet() {
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
        // select Tet and Hex
        let kinds = vec![
            // Tet
            GeoKind::Tet4,
            GeoKind::Tet10,
            GeoKind::Tet20,
            // Hex
            GeoKind::Hex8,
            GeoKind::Hex20,
            GeoKind::Hex32,
        ];

        // define tolerances
        let tols = HashMap::from([
            (GeoKind::Tet4, 1e-15),
            (GeoKind::Tet10, 1e-15),
            (GeoKind::Tet20, 1e-14),
            (GeoKind::Hex8, 1e-15),
            (GeoKind::Hex20, 1e-15),
            (GeoKind::Hex32, 1e-15),
        ]);

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
        let space_ndim = 3;
        for kind in kinds {
            // allocate shape
            let shape = Shape::new(kind);

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // generate state with coordinates matrix
            let state = gen_state_with_coords_matrix_parallel(&shape)?;

            // loop over faces
            for f in 0..shape.nface {
                let face_shape = Shape::new(shape.kind.face_kind().unwrap());
                let mut face_state = StateOfShape::new(
                    face_shape.kind,
                    &(0..face_shape.nnode)
                        .map(|i| {
                            let m = shape.kind.face_node_id(f, i);
                            (0..space_ndim).map(|j| state.coords_transp[j][m]).collect()
                        })
                        .collect::<Vec<_>>(),
                )?;

                // calc normal vector
                face_shape.calc_boundary_normal(&mut normal, &mut face_state, ksi)?;
                if shape.kind.is_tri_or_tet() {
                    // check tetrahedron
                    assert_vec_approx_eq!(normal.as_data(), tet_correct[f], tol);
                } else {
                    // check hexahedron
                    assert_vec_approx_eq!(normal.as_data(), hex_correct[f], tol);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn approximate_ksi_works() -> Result<(), StrError> {
        // select all kinds, except Lin
        let kinds = vec![
            // Tri
            GeoKind::Tri3,
            GeoKind::Tri6,
            GeoKind::Tri10,
            GeoKind::Tri15,
            // Qua
            GeoKind::Qua4,
            GeoKind::Qua8,
            GeoKind::Qua9,
            GeoKind::Qua12,
            GeoKind::Qua16,
            GeoKind::Qua17,
            // Tet
            GeoKind::Tet4,
            GeoKind::Tet10,
            GeoKind::Tet20,
            // Hex
            GeoKind::Hex8,
            GeoKind::Hex20,
            GeoKind::Hex32,
        ];

        // define tolerances
        let tols = HashMap::from([
            // Tri
            (GeoKind::Tri3, 1e-14),
            (GeoKind::Tri6, 1e-15),
            (GeoKind::Tri10, 1e-15),
            (GeoKind::Tri15, 1e-14),
            // Qua
            (GeoKind::Qua4, 1e-15),
            (GeoKind::Qua8, 1e-15),
            (GeoKind::Qua9, 1e-15),
            (GeoKind::Qua12, 1e-15),
            (GeoKind::Qua16, 1e-15),
            (GeoKind::Qua17, 1e-13),
            // Tet
            (GeoKind::Tet4, 1e-15),
            (GeoKind::Tet10, 1e-15),
            (GeoKind::Tet20, 1e-15),
            // Hex
            (GeoKind::Hex8, 1e-15),
            (GeoKind::Hex20, 1e-15),
            (GeoKind::Hex32, 1e-14),
        ]);

        // loop over shapes
        for kind in kinds {
            // allocate shape
            let shape = Shape::new(kind);

            // generate state with coordinates matrix
            let space_ndim = shape.geo_ndim;
            let mut state = gen_state_with_coords_matrix(space_ndim, &shape)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // loop over nodes of shape
            let mut x = Vector::new(space_ndim);
            let mut ksi = vec![0.0; shape.geo_ndim];
            for m in 0..shape.nnode {
                // get ξᵐ corresponding to node m
                let ksi_ref = shape.kind.reference_coords(m);

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
            let ksi_in = if shape.kind.is_tri_or_tet() {
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

    #[test]
    fn approximate_ksi_works_outside() -> Result<(), StrError> {
        // Equilateral triangle
        //
        //           /
        //        2   \
        //       / \   \
        //      / ↑ \   l
        //     5  h  4   \
        //    /   ↓   \   \
        //   /         \   /
        //  0-----3-----1
        //
        //  |--s--|--s--|
        //
        //  |-----l-----|
        //
        // area = l * h / 2.0;
        let l = 5.0;
        let s = l / 2.0;
        let h = l * SQRT_3 / 2.0;
        let (x0, y0) = (3.0, 4.0);
        let (x1, y1) = (x0 + l, y0);
        let (x2, y2) = (x0 + s, y0 + h);
        let (x3, y3) = (x0 + s, y0);
        let (x4, y4) = (x0 + 1.5 * s, y0 + 0.5 * h);
        let (x5, y5) = (x0 + 0.5 * s, y0 + 0.5 * h);
        let shape = Shape::new(GeoKind::Tri6);
        let mut state = StateOfShape::new(
            shape.kind,
            &[[x0, y0], [x1, y1], [x2, y2], [x3, y3], [x4, y4], [x5, y5]],
        )
        .unwrap();
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

    // Holds arguments for numerical differentiation of N with respect to x => G (gradient) matrix
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
        // select all kinds, except Lin
        let kinds = vec![
            // Tri
            GeoKind::Tri3,
            GeoKind::Tri6,
            GeoKind::Tri10,
            GeoKind::Tri15,
            // Qua
            GeoKind::Qua4,
            GeoKind::Qua8,
            GeoKind::Qua9,
            GeoKind::Qua12,
            GeoKind::Qua16,
            GeoKind::Qua17,
            // Tet
            GeoKind::Tet4,
            GeoKind::Tet10,
            GeoKind::Tet20,
            // Hex
            GeoKind::Hex8,
            GeoKind::Hex20,
            GeoKind::Hex32,
        ];

        // define tolerances
        let tols = HashMap::from([
            // Tri
            (GeoKind::Tri3, 1e-11),
            (GeoKind::Tri6, 1e-10),
            (GeoKind::Tri10, 1e-10),
            (GeoKind::Tri15, 1e-9),
            // Qua
            (GeoKind::Qua4, 1e-11),
            (GeoKind::Qua8, 1e-10),
            (GeoKind::Qua9, 1e-10),
            (GeoKind::Qua12, 1e-10),
            (GeoKind::Qua16, 1e-9),
            (GeoKind::Qua17, 1e-10),
            // Tet
            (GeoKind::Tet4, 1e-11),
            (GeoKind::Tet10, 1e-10),
            (GeoKind::Tet20, 1e-9),
            // Hex
            (GeoKind::Hex8, 1e-11),
            (GeoKind::Hex20, 1e-10),
            (GeoKind::Hex32, 1e-10),
        ]);

        // loop over shapes
        for kind in kinds {
            // allocate shape
            let shape = Shape::new(kind);

            // generate state with coordinates matrix
            let space_ndim = usize::max(2, shape.geo_ndim);
            let mut state = gen_state_with_coords_matrix(space_ndim, &shape)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // set ξ within reference space
            let at_ksi = vec![0.25; shape.geo_ndim];

            // compute x corresponding to ξ using the isoparametric formula
            let mut at_x = Vector::new(space_ndim);
            shape.calc_coords(&mut at_x, &mut state, &at_ksi)?;

            // compute gradient
            let det_jac = shape.calc_gradient(&mut state, &at_ksi)?;
            assert!(det_jac > 0.0);

            // set arguments for numerical integration
            let args = &mut ArgsNumG {
                shape: &shape,
                state: state.clone(),
                at_x,
                x: Vector::new(space_ndim),
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
        let shape = Shape::new(GeoKind::Lin2);
        let (xa, xb) = (2.0, 5.0);
        let mut state = StateOfShape::new(shape.kind, &[[xa, 8.0], [xb, 8.0]])?;
        let ip_lin_legendre_2: [[f64; 4]; 2] = [
            [-0.5773502691896257, 0.0, 0.0, 1.0],
            [0.5773502691896257, 0.0, 0.0, 1.0],
        ];
        let integ_points = shape.calc_integ_points_coords(&mut state, &ip_lin_legendre_2)?;
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
