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
