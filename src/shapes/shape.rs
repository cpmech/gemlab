use super::*;
use crate::StrError;
use russell_lab::{inverse, mat_mat_mul, mat_vec_mul, vector_norm, Matrix, NormVec, Vector};
use serde::{Deserialize, Serialize};

/// Implements an isoparametric geometric shape for numerical integration and more
///
/// # Definitions
///
/// Here, we consider the following dimensions:
///
/// * `space_ndim` -- is the number of dimensions of the space under study
/// * `geo_ndim` -- is the number of dimensions of the geometry element (shape),
///                  for instance, a line in the 2D space has `geo_ndim = 1` and
///                  `space_ndim = 2`. Another example is a triangle in the 3D space
///                  which has `geo_ndim = 2` and `space_ndim = 3`.
///
/// We also consider the following counting variables:
///
/// * `nnode` -- number of points that define the shape (number of nodes)
/// * `nedge` -- number of edges
/// * `nface` -- number of faces
/// * `edge_nnode` -- number of points that define the edge (number of nodes)
/// * `face_nnode` -- number of points that define the face (number of nodes)
/// * `face_nedge` -- face's number of edges
///
/// When performing numerical integrations, we use the following notation:
/// |J| is the determinant of the Jacobian, ||J|| is the norm of the Jacobian vector
/// for line in multi-dimensions, `nip` is the number of integration points,
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
/// # Note
///
/// All public members are **readonly** and should not be modified externally.
#[derive(Clone, Debug, Deserialize, Serialize)]
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

    /// Function to calculate interpolation functions
    #[serde(
        serialize_with = "shape_serialize_fn_interp",
        deserialize_with = "shape_deserialize_fn_interp"
    )]
    fn_interp: FnInterp,

    /// Function to calculate local derivatives (w.r.t. ksi) of interpolation functions
    #[serde(
        serialize_with = "shape_serialize_fn_deriv",
        deserialize_with = "shape_deserialize_fn_deriv"
    )]
    fn_deriv: FnDeriv,

    /// Matrix Xᵀ: (space_ndim,nnode) transposed coordinates matrix (real space)
    pub coords_transp: Matrix,

    /// Tells if at least the last entry of the coordinates matrix has been input
    ok_last_coord: bool,

    /// Minimum (space_ndim) coordinates from the X matrix
    pub min_coords: Vec<f64>,

    /// Maximum (space_ndim) coordinates from the X matrix
    pub max_coords: Vec<f64>,
}

impl Shape {
    /// Creates a new geometric shape
    ///
    /// # Input
    ///
    /// * `space_ndim` -- the space dimension (1, 2, or 3)
    /// * `geo_ndim` -- the dimension of the shape; e.g. a 1D line in a 3D space
    /// * `nnode` -- the number of points defining the shape (number of nodes)
    ///
    /// # Note
    ///
    /// Some methods require that the coordinates matrix be set first.
    /// This can be accomplished by calling the `set_node` method.
    pub fn new(space_ndim: usize, geo_ndim: usize, nnode: usize) -> Result<Self, StrError> {
        // collect geometry data
        let (class, kind, nedge, nface, edge_nnode, face_nnode, face_nedge): (
            GeoClass,
            GeoKind,
            usize,
            usize,
            usize,
            usize,
            usize,
        ) = match (geo_ndim, nnode) {
            // Lin
            (1, 2) => (
                GeoClass::Lin,
                GeoKind::Lin2,
                Lin2::NEDGE,
                Lin2::NFACE,
                Lin2::EDGE_NNODE,
                Lin2::FACE_NNODE,
                Lin2::FACE_NEDGE,
            ),
            (1, 3) => (
                GeoClass::Lin,
                GeoKind::Lin3,
                Lin3::NEDGE,
                Lin3::NFACE,
                Lin3::EDGE_NNODE,
                Lin3::FACE_NNODE,
                Lin3::FACE_NEDGE,
            ),
            (1, 4) => (
                GeoClass::Lin,
                GeoKind::Lin4,
                Lin4::NEDGE,
                Lin4::NFACE,
                Lin4::EDGE_NNODE,
                Lin4::FACE_NNODE,
                Lin4::FACE_NEDGE,
            ),
            (1, 5) => (
                GeoClass::Lin,
                GeoKind::Lin5,
                Lin5::NEDGE,
                Lin5::NFACE,
                Lin5::EDGE_NNODE,
                Lin5::FACE_NNODE,
                Lin5::FACE_NEDGE,
            ),

            // Tri
            (2, 3) => (
                GeoClass::Tri,
                GeoKind::Tri3,
                Tri3::NEDGE,
                Tri3::NFACE,
                Tri3::EDGE_NNODE,
                Tri3::FACE_NNODE,
                Tri3::FACE_NEDGE,
            ),
            (2, 6) => (
                GeoClass::Tri,
                GeoKind::Tri6,
                Tri6::NEDGE,
                Tri6::NFACE,
                Tri6::EDGE_NNODE,
                Tri6::FACE_NNODE,
                Tri6::FACE_NEDGE,
            ),
            (2, 10) => (
                GeoClass::Tri,
                GeoKind::Tri10,
                Tri10::NEDGE,
                Tri10::NFACE,
                Tri10::EDGE_NNODE,
                Tri10::FACE_NNODE,
                Tri10::FACE_NEDGE,
            ),
            (2, 15) => (
                GeoClass::Tri,
                GeoKind::Tri15,
                Tri15::NEDGE,
                Tri15::NFACE,
                Tri15::EDGE_NNODE,
                Tri15::FACE_NNODE,
                Tri15::FACE_NEDGE,
            ),

            // Qua
            (2, 4) => (
                GeoClass::Qua,
                GeoKind::Qua4,
                Qua4::NEDGE,
                Qua4::NFACE,
                Qua4::EDGE_NNODE,
                Qua4::FACE_NNODE,
                Qua4::FACE_NEDGE,
            ),
            (2, 8) => (
                GeoClass::Qua,
                GeoKind::Qua8,
                Qua8::NEDGE,
                Qua8::NFACE,
                Qua8::EDGE_NNODE,
                Qua8::FACE_NNODE,
                Qua8::FACE_NEDGE,
            ),
            (2, 9) => (
                GeoClass::Qua,
                GeoKind::Qua9,
                Qua9::NEDGE,
                Qua9::NFACE,
                Qua9::EDGE_NNODE,
                Qua9::FACE_NNODE,
                Qua9::FACE_NEDGE,
            ),
            (2, 12) => (
                GeoClass::Qua,
                GeoKind::Qua12,
                Qua12::NEDGE,
                Qua12::NFACE,
                Qua12::EDGE_NNODE,
                Qua12::FACE_NNODE,
                Qua12::FACE_NEDGE,
            ),
            (2, 16) => (
                GeoClass::Qua,
                GeoKind::Qua16,
                Qua16::NEDGE,
                Qua16::NFACE,
                Qua16::EDGE_NNODE,
                Qua16::FACE_NNODE,
                Qua16::FACE_NEDGE,
            ),
            (2, 17) => (
                GeoClass::Qua,
                GeoKind::Qua17,
                Qua17::NEDGE,
                Qua17::NFACE,
                Qua17::EDGE_NNODE,
                Qua17::FACE_NNODE,
                Qua17::FACE_NEDGE,
            ),

            // Tet
            (3, 4) => (
                GeoClass::Tet,
                GeoKind::Tet4,
                Tet4::NEDGE,
                Tet4::NFACE,
                Tet4::EDGE_NNODE,
                Tet4::FACE_NNODE,
                Tet4::FACE_NEDGE,
            ),
            (3, 10) => (
                GeoClass::Tet,
                GeoKind::Tet10,
                Tet10::NEDGE,
                Tet10::NFACE,
                Tet10::EDGE_NNODE,
                Tet10::FACE_NNODE,
                Tet10::FACE_NEDGE,
            ),

            // Hex
            (3, 8) => (
                GeoClass::Hex,
                GeoKind::Hex8,
                Hex8::NEDGE,
                Hex8::NFACE,
                Hex8::EDGE_NNODE,
                Hex8::FACE_NNODE,
                Hex8::FACE_NEDGE,
            ),
            (3, 20) => (
                GeoClass::Hex,
                GeoKind::Hex20,
                Hex20::NEDGE,
                Hex20::NFACE,
                Hex20::EDGE_NNODE,
                Hex20::FACE_NNODE,
                Hex20::FACE_NEDGE,
            ),
            _ => return Err("(geo_ndim,nnode) combination is invalid"),
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
            fn_interp: i32_to_fn_interp(kind as i32),
            fn_deriv: i32_to_fn_deriv(kind as i32),
            coords_transp: Matrix::new(space_ndim, nnode),
            ok_last_coord: false,
            min_coords: vec![f64::MAX; space_ndim],
            max_coords: vec![f64::MIN; space_ndim],
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
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate (geo_ndim)
    ///
    /// # Updated variables
    ///
    /// * `interp` -- interpolation functions (nnode)
    #[inline]
    pub fn calc_interp(&self, state: &mut ShapeState, ksi: &[f64]) {
        (self.fn_interp.1)(&mut state.interp, ksi);
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
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate (geo_ndim)
    ///
    /// # Updated variables
    ///
    /// * `deriv` -- interpolation functions (nnode)
    #[inline]
    pub fn calc_deriv(&self, state: &mut ShapeState, ksi: &[f64]) {
        (self.fn_deriv.1)(&mut state.deriv, ksi);
    }

    /// Sets a component of the coordinates matrix
    ///
    /// ```text
    ///     ┌               ┐
    ///     | x⁰₀  x⁰₁  x⁰₂ |
    ///     | x¹₀  x¹₁  x¹₂ |
    /// X = | x²₀  x²₁  x²₂ |
    ///     |      ...      |
    ///     | xᴹ₀  xᴹ₁  xᴹ₂ |
    ///     └               ┘_(nnode,space_ndim)
    /// ```
    ///
    /// where `M = nnode - 1`
    ///
    /// # Input
    ///
    /// * `m` -- node index from 0 to nnode - 1
    /// * `j` -- dimension index from 0 to space_ndim - 1
    /// * `value` -- the X(m,j) component
    pub fn set_node(&mut self, m: usize, j: usize, value: f64) -> Result<(), StrError> {
        if m >= self.nnode {
            return Err("index of node is invalid");
        }
        if j >= self.space_ndim {
            return Err("index of space dimension is invalid");
        }
        self.coords_transp[j][m] = value;
        if value < self.min_coords[j] {
            self.min_coords[j] = value;
        }
        if value > self.max_coords[j] {
            self.max_coords[j] = value;
        }
        if m == self.nnode - 1 && j == self.space_ndim - 1 {
            self.ok_last_coord = true;
        } else {
            self.ok_last_coord = false;
        }
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
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate. The length of ξ must be equal to geo_ndim at least,
    ///            while lengths greater than geo_ndim are allowed (and ignored). In this way,
    ///            we can pass a slice with integration point data such as `[f64; 4]`.
    ///
    /// # Output
    ///
    /// * `x` -- real coordinate (space_ndim)
    ///
    /// # Warning
    ///
    /// You must set the coordinates matrix first (via the `set_node` method),
    /// otherwise the computations will generate wrong results.
    pub fn calc_coords(&self, state: &mut ShapeState, x: &mut Vector, ksi: &[f64]) -> Result<(), StrError> {
        if x.dim() != self.space_ndim {
            return Err("x.dim() must equal space_ndim");
        }
        if ksi.len() < self.geo_ndim {
            return Err("ksi.len() must equal geo_ndim at least");
        }
        if !self.ok_last_coord {
            return Err("the last node coordinate has not been input yet");
        }
        self.calc_interp(state, ksi);
        mat_vec_mul(x, 1.0, &self.coords_transp, &state.interp)
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
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate. The length of ξ must be equal to geo_ndim at least,
    ///            while lengths greater than geo_ndim are allowed (and ignored). In this way,
    ///            we can pass a slice with integration point data such as `[f64; 4]`.
    ///
    /// # Output
    ///
    /// * If `geo_ndim = space_ndim`, returns the determinant of the Jacobian
    /// * If `geo_ndim = 1` and `space_ndim > 1`, returns the norm of the Jacobian vector
    /// * Otherwise, returns zero
    ///
    /// # Updated variables
    ///
    /// * `jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    /// * `inv_jacobian` -- If `geo_ndim = space_ndim`: inverse Jacobian matrix (space_ndim,space_ndim)
    ///
    /// # Warning
    ///
    /// You must set the coordinates matrix first (via the `set_node` method),
    /// otherwise the computations will generate wrong results.
    pub fn calc_jacobian(&self, state: &mut ShapeState, ksi: &[f64]) -> Result<f64, StrError> {
        if ksi.len() < self.geo_ndim {
            return Err("ksi.len() must equal geo_ndim at least");
        }
        if !self.ok_last_coord {
            return Err("the last node coordinate has not been input yet");
        }
        self.calc_deriv(state, ksi);
        mat_mat_mul(&mut state.jacobian, 1.0, &self.coords_transp, &state.deriv)?;
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
    ///
    /// # Updated variables
    ///
    /// * `jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    ///
    /// # Warning
    ///
    /// You must set the coordinates matrix first (via the `set_node` method),
    /// otherwise the computations will generate wrong results.
    pub fn calc_boundary_normal(
        &self,
        state: &mut ShapeState,
        normal: &mut Vector,
        ksi: &[f64],
    ) -> Result<(), StrError> {
        // check
        if self.geo_ndim == self.space_ndim {
            return Err("geo_ndim must be smaller than space_ndim");
        }
        if normal.dim() != self.space_ndim {
            return Err("normal.dim() must equal space_ndim");
        }
        if ksi.len() < self.geo_ndim {
            return Err("ksi.len() must equal geo_ndim at least");
        }
        if !self.ok_last_coord {
            return Err("the last node coordinate has not been input yet");
        }

        // compute Jacobian
        self.calc_deriv(state, ksi);
        mat_mat_mul(&mut state.jacobian, 1.0, &self.coords_transp, &state.deriv)?;

        // line in 2D
        if self.geo_ndim == 1 && self.space_ndim == 2 {
            normal[0] = -state.jacobian[1][0];
            normal[1] = state.jacobian[0][0];
            return Ok(());
        }

        // surface in 3D
        if self.geo_ndim == 2 && self.space_ndim == 3 {
            normal[0] = state.jacobian[1][0] * state.jacobian[2][1] - state.jacobian[2][0] * state.jacobian[1][1];
            normal[1] = state.jacobian[2][0] * state.jacobian[0][1] - state.jacobian[0][0] * state.jacobian[2][1];
            normal[2] = state.jacobian[0][0] * state.jacobian[1][1] - state.jacobian[1][0] * state.jacobian[0][1];
            return Ok(());
        }

        Err("geo_ndim and space_ndim combination is invalid for boundary normal")
    }

    /// Approximates the reference coordinates from given real coordinates (inverse mapping)
    ///
    /// **Note:** This function works with `geo_ndim == space_ndim` only.
    ///
    /// We use Newton iterations with the inverse of the Jacobian to compute ξ(x).
    ///
    /// # Input
    ///
    /// * `x` -- real coordinates (space_ndim = geo_ndim)
    /// * `nit_max` -- maximum number of iterations (e.g., 10)
    /// * `tol` -- tolerance for the norm of the difference x - x(ξ) (e.g., 1e-14)
    ///
    /// # Output
    ///
    /// * `ksi` -- ξ reference coordinates (geo_ndim = space_ndim)
    /// * returns the number of iterations
    ///
    /// # Warning
    ///
    /// You must set the coordinates matrix first (via the `set_node` method),
    /// otherwise the computations will generate wrong results.
    pub fn approximate_ksi(
        &self,
        state: &mut ShapeState,
        ksi: &mut [f64],
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
        if !self.ok_last_coord {
            return Err("the last node coordinate has not been input yet");
        }

        // use linear scale to guess ksi
        let (min_ksi, max_ksi, del_ksi) = ref_domain_limits(self.class);
        for j in 0..self.geo_ndim {
            ksi[j] = (x[j] - self.min_coords[j]) / (self.max_coords[j] - self.min_coords[j]) * del_ksi + min_ksi;
            if ksi[j] < min_ksi {
                ksi[j] = min_ksi;
            }
            if ksi[j] > max_ksi {
                ksi[j] = max_ksi;
            }
        }

        // perform iterations
        let mut residual = Vector::new(self.space_ndim);
        let mut x_at_ksi = Vector::new(self.space_ndim);
        let mut delta_ksi = Vector::new(self.geo_ndim);
        for it in 0..nit_max {
            self.calc_coords(state, &mut x_at_ksi, &ksi)?;
            for i in 0..self.space_ndim {
                residual[i] = x[i] - x_at_ksi[i];
            }
            if vector_norm(&residual, NormVec::Euc) <= tol {
                return Ok(it);
            }
            self.calc_jacobian(state, ksi)?;
            mat_vec_mul(&mut delta_ksi, 1.0, &state.inv_jacobian, &residual)?;
            for j in 0..self.geo_ndim {
                ksi[j] += delta_ksi[j];
            }
        }

        Err("inverse mapping failed to converge")
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
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate. The length of ξ must be equal to geo_ndim at least,
    ///            while lengths greater than geo_ndim are allowed (and ignored). In this way,
    ///            we can pass a slice with integration point data such as `[f64; 4]`.
    ///
    /// # Output
    ///
    /// Returns the determinant of the Jacobian.
    ///
    /// # Updated variables
    ///
    /// * `jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    /// * `inv_jacobian` -- inverse Jacobian matrix (space_ndim,space_ndim)
    /// * `gradient` -- gradient matrix (nnode,space_ndim)
    ///
    /// # Warning
    ///
    /// You must set the coordinates matrix first (via the `set_node` method),
    /// otherwise the computations will generate wrong results.
    pub fn calc_gradient(&self, state: &mut ShapeState, ksi: &[f64]) -> Result<f64, StrError> {
        if self.geo_ndim != self.space_ndim {
            return Err("geo_ndim must equal space_ndim");
        }
        if ksi.len() < self.geo_ndim {
            return Err("ksi.len() must equal geo_ndim at least");
        }
        if !self.ok_last_coord {
            return Err("the last node coordinate has not been input yet");
        }
        let det_jac = self.calc_jacobian(state, ksi)?;
        mat_mat_mul(&mut state.gradient, 1.0, &state.deriv, &state.inv_jacobian)?;
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
    /// Returns an array with `nip` (number of integration point) vectors, where
    /// each vector has a dimension equal to `space_ndim`.
    ///
    /// # Warning
    ///
    /// You must set the coordinates matrix first (via the `set_node` method),
    /// otherwise the computations will generate wrong results.
    pub fn calc_int_points_coords(&self, state: &mut ShapeState) -> Result<Vec<Vector>, StrError> {
        if !self.ok_last_coord {
            return Err("the last node coordinate has not been input yet");
        }
        let mut all_coords = Vec::new();
        for iota in state.ip_data {
            let mut x = Vector::new(self.space_ndim);
            self.calc_coords(state, &mut x, iota).unwrap();
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
    pub fn get_edge_node_id(&self, e: usize, i: usize) -> usize {
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
    pub fn get_face_node_id(&self, f: usize, i: usize) -> usize {
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
    pub fn get_face_edge_node_id(&self, f: usize, k: usize, i: usize) -> usize {
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
    pub fn get_reference_coords(&self, m: usize) -> &'static [f64] {
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

/// Returns the (min,max,delta) limits on the reference domain
///
/// # Tri and Tet
///
/// * `ξ_min` = 0.0
/// * `ξ_max` = 1.0
/// * `Δξ` = 1.0
///
/// # Lin, Qua, Hex
///
/// * `ξ_min` = -1.0
/// * `ξ_max` = +1.0
/// * `Δξ` = 2.0
#[inline]
pub fn ref_domain_limits(class: GeoClass) -> (f64, f64, f64) {
    if class == GeoClass::Tri || class == GeoClass::Tet {
        (0.0, 1.0, 1.0)
    } else {
        (-1.0, 1.0, 2.0)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{ref_domain_limits, GeoClass, GeoKind, Shape, ShapeState};
    use crate::util::{PI, SQRT_3};
    use crate::StrError;
    use russell_chk::{assert_approx_eq, assert_deriv_approx_eq, assert_vec_approx_eq};
    use russell_lab::{copy_vector, scale_vector, vector_norm, NormVec, Vector};
    use serde::{Deserialize, Serialize};
    use std::collections::HashMap;

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
        if x.dim() == 1 {
            x[0] = r;
            return;
        }
        let a = AMIN + (ksi[1] - min_ksi) * (AMAX - AMIN) / del_ksi;
        x[0] = r * f64::cos(a);
        x[1] = r * f64::sin(a);
        if x.dim() == 3 {
            x[2] = ksi[2];
        }
    }

    // Generates matrix of coordinates
    fn set_coords_matrix(shape: &mut Shape) {
        let mut x = Vector::new(shape.space_ndim);
        let mut ksi_aux = vec![0.0; shape.space_ndim];
        for m in 0..shape.nnode {
            let ksi = shape.get_reference_coords(m);
            if shape.geo_ndim == shape.space_ndim {
                gen_coords(&mut x, ksi, shape.class);
            } else if shape.geo_ndim == 1 && shape.space_ndim == 2 {
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = 1.0;
                gen_coords(&mut x, &ksi_aux, shape.class);
            } else if shape.geo_ndim == 2 && shape.space_ndim == 3 {
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = ksi[1];
                ksi_aux[2] = 1.0;
                gen_coords(&mut x, &ksi_aux, shape.class);
            } else {
                panic!("(geo_ndim,space_ndim) pair is invalid");
            }
            for j in 0..shape.space_ndim {
                shape.set_node(m, j, x[j]).unwrap();
            }
        }
    }

    // Sets the coordinates such that the shape is parallel to the x,y,z axes
    // For triangles and tetrahedra, one edge/face will cut the axes equally
    // * Qua and Hex will have real coordinates equal to the natural coordinates
    // * Tri and Tet will be scaled using the natural coordinates
    // * All edges parallel to the x,y,z axes will have lengths equal to 2.0
    fn set_coords_matrix_parallel(shape: &mut Shape) -> Result<(), StrError> {
        let mut scale = 1.0;
        if shape.class == GeoClass::Tri || shape.class == GeoClass::Tet {
            scale = 2.0;
        }
        for m in 0..shape.nnode {
            let ksi = shape.get_reference_coords(m);
            for j in 0..shape.geo_ndim {
                shape.set_node(m, j, scale * ksi[j])?;
            }
        }
        Ok(())
    }

    #[test]
    fn capture_some_wrong_input() -> Result<(), StrError> {
        // new
        assert_eq!(
            Shape::new(1, 1, 1).err(),
            Some("(geo_ndim,nnode) combination is invalid")
        );

        // set_node
        let mut shape = Shape::new(1, 1, 2)?;
        let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
        assert_eq!(shape.set_node(3, 0, 0.0).err(), Some("index of node is invalid"));
        assert_eq!(
            shape.set_node(0, 2, 0.0).err(),
            Some("index of space dimension is invalid")
        );

        // calc_jacobian
        assert_eq!(
            shape.calc_jacobian(&mut state, &[]).err(),
            Some("ksi.len() must equal geo_ndim at least")
        );
        assert_eq!(
            shape.calc_jacobian(&mut state, &[0.0]).err(),
            Some("the last node coordinate has not been input yet")
        );

        // calc_coords
        let mut x = Vector::new(1);
        assert_eq!(
            shape.calc_coords(&mut state, &mut x, &[0.0]).err(),
            Some("the last node coordinate has not been input yet")
        );
        assert_eq!(
            shape.calc_coords(&mut state, &mut x, &[]).err(),
            Some("ksi.len() must equal geo_ndim at least")
        );

        // calc_coords
        shape.set_node(0, 0, 0.0)?;
        shape.set_node(1, 0, 1.0)?;
        let mut x = Vector::new(2);
        assert_eq!(
            shape.calc_coords(&mut state, &mut x, &[0.0]).err(),
            Some("x.dim() must equal space_ndim")
        );

        // calc_boundary_normal
        let mut normal = Vector::new(1);
        assert_eq!(
            shape.calc_boundary_normal(&mut state, &mut normal, &[0.0]).err(),
            Some("geo_ndim must be smaller than space_ndim")
        );
        let shape = Shape::new(2, 1, 2)?;
        let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
        assert_eq!(
            shape.calc_boundary_normal(&mut state, &mut normal, &[0.0]).err(),
            Some("normal.dim() must equal space_ndim")
        );
        let mut normal = Vector::new(2);
        assert_eq!(
            shape.calc_boundary_normal(&mut state, &mut normal, &[]).err(),
            Some("ksi.len() must equal geo_ndim at least")
        );
        assert_eq!(
            shape.calc_boundary_normal(&mut state, &mut normal, &[0.0]).err(),
            Some("the last node coordinate has not been input yet")
        );

        // approximate_ksi
        let shape = Shape::new(2, 1, 2)?;
        let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
        let mut ksi = vec![0.0; 1];
        let x = Vector::new(2);
        assert_eq!(
            shape.approximate_ksi(&mut state, &mut ksi, &x, 2, 1e-5).err(),
            Some("geo_ndim must equal space_ndim")
        );
        let shape = Shape::new(1, 1, 2)?;
        let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
        let x = Vector::new(2);
        assert_eq!(
            shape.approximate_ksi(&mut state, &mut ksi, &x, 2, 1e-5).err(),
            Some("x.dim() must equal space_ndim")
        );
        let x = Vector::new(1);
        let mut ksi = vec![0.0; 0];
        assert_eq!(
            shape.approximate_ksi(&mut state, &mut ksi, &x, 2, 1e-5).err(),
            Some("ksi.len() must equal geo_ndim")
        );
        let mut ksi = vec![0.0; 1];
        assert_eq!(
            shape.approximate_ksi(&mut state, &mut ksi, &x, 2, 1e-5).err(),
            Some("the last node coordinate has not been input yet")
        );

        // calc_gradient
        let shape = Shape::new(2, 1, 2)?;
        let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
        assert_eq!(
            shape.calc_gradient(&mut state, &[0.0]).err(),
            Some("geo_ndim must equal space_ndim")
        );
        let shape = Shape::new(1, 1, 2)?;
        let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
        assert_eq!(
            shape.calc_gradient(&mut state, &[]).err(),
            Some("ksi.len() must equal geo_ndim at least")
        );
        assert_eq!(
            shape.calc_gradient(&mut state, &[0.0]).err(),
            Some("the last node coordinate has not been input yet")
        );

        // calc_int_points_coords
        assert_eq!(
            shape.calc_int_points_coords(&mut state).err(),
            Some("the last node coordinate has not been input yet")
        );
        Ok(())
    }

    #[test]
    fn set_node_resets_ok_last_coord() -> Result<(), StrError> {
        let mut shape = Shape::new(1, 1, 2)?;
        assert_eq!(shape.ok_last_coord, false);
        shape.set_node(0, 0, 0.0)?;
        assert_eq!(shape.ok_last_coord, false);
        shape.set_node(1, 0, 0.0)?;
        assert_eq!(shape.ok_last_coord, true);
        shape.set_node(0, 0, 0.0)?; // must reset
        assert_eq!(shape.ok_last_coord, false);
        Ok(())
    }

    #[test]
    fn getters_work() -> Result<(), StrError> {
        let pairs = vec![
            // Lin
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5),
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
            // Tet
            (3, 4),
            (3, 10),
            // Hex
            (3, 8),
            (3, 20),
        ];
        for (geo_ndim, nnode) in pairs {
            let space_ndim = geo_ndim;
            let shape = &mut Shape::new(space_ndim, geo_ndim, nnode)?;
            match shape.class {
                GeoClass::Lin => assert_eq!(shape.get_edge_node_id(0, 0), 0),
                GeoClass::Tri => assert_eq!(shape.get_edge_node_id(0, 0), 1),
                GeoClass::Qua => assert_eq!(shape.get_edge_node_id(0, 0), 1),
                GeoClass::Tet => assert_eq!(shape.get_edge_node_id(0, 0), 0),
                GeoClass::Hex => assert_eq!(shape.get_edge_node_id(0, 0), 0),
            }
            match shape.class {
                GeoClass::Lin => assert_eq!(shape.get_face_node_id(0, 0), 0),
                GeoClass::Tri => assert_eq!(shape.get_face_node_id(0, 0), 0),
                GeoClass::Qua => assert_eq!(shape.get_face_node_id(0, 0), 0),
                GeoClass::Tet => assert_eq!(shape.get_face_node_id(0, 0), 0),
                GeoClass::Hex => assert_eq!(shape.get_face_node_id(0, 0), 0),
            }
            match shape.class {
                GeoClass::Lin => assert_eq!(shape.get_face_edge_node_id(0, 0, 0), 0),
                GeoClass::Tri => assert_eq!(shape.get_face_edge_node_id(0, 0, 0), 0),
                GeoClass::Qua => assert_eq!(shape.get_face_edge_node_id(0, 0, 0), 0),
                GeoClass::Tet => assert_eq!(shape.get_face_edge_node_id(0, 0, 0), 0),
                GeoClass::Hex => assert_eq!(shape.get_face_edge_node_id(0, 0, 0), 0),
            }
        }
        Ok(())
    }

    #[test]
    fn calc_interp_works() -> Result<(), StrError> {
        // define dims and number of nodes
        let pairs = vec![
            // Lin
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5),
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
            // Tet
            (3, 4),
            (3, 10),
            // Hex
            (3, 8),
            (3, 20),
        ];

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
        for (geo_ndim, nnode) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let shape = Shape::new(space_ndim, geo_ndim, nnode)?;
            let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // loop over nodes of shape
            for m in 0..shape.nnode {
                // get ξᵐ corresponding to node m
                let ksi = shape.get_reference_coords(m);

                // compute interpolation function Nⁿ(ξᵐ)
                shape.calc_interp(&mut state, ksi);

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
    struct ArgsNumL {
        shape: Shape,      // auxiliary (copy) shape
        state: ShapeState, // auxiliary state
        at_ksi: Vec<f64>,  // at reference coord value
        ksi: Vec<f64>,     // temporary reference coord
        m: usize,          // node index from 0 to nnode
        j: usize,          // dimension index from 0 to geom_ndim
    }

    // Computes Nᵐ(ξ) with variable v := ξⱼ
    fn aux_deriv(v: f64, args: &mut ArgsNumL) -> f64 {
        args.ksi.copy_from_slice(&args.at_ksi);
        args.ksi[args.j] = v;
        args.shape.calc_interp(&mut args.state, &args.ksi);
        args.state.interp[args.m]
    }

    #[test]
    fn calc_deriv_works() -> Result<(), StrError> {
        // define dims and number of nodes
        let pairs = vec![
            // Lin
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5),
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
            // Tet
            (3, 4),
            (3, 10),
            // Hex
            (3, 8),
            (3, 20),
        ];

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
        for (geo_ndim, nnode) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let shape = Shape::new(space_ndim, geo_ndim, nnode)?;
            let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // set ξ within reference space
            let at_ksi = vec![0.25; shape.geo_ndim];

            // compute all derivatives of interpolation functions w.r.t ξ
            shape.calc_deriv(&mut state, &at_ksi);

            // set arguments for numerical integration
            let aux_shape = Shape::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
            let aux_state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
            let args = &mut ArgsNumL {
                shape: aux_shape,
                state: aux_state,
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
            let space_ndim = geo_ndim;
            let mut shape = Shape::new(space_ndim, geo_ndim, nnode)?;
            let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();
            let tol_in = *tols_in.get(&shape.kind).unwrap();

            // set coordinates matrix
            set_coords_matrix(&mut shape);

            // loop over nodes of shape
            let mut x = Vector::new(shape.space_ndim);
            let mut x_correct = Vector::new(shape.space_ndim);
            for m in 0..shape.nnode {
                // get ξᵐ corresponding to node m
                let ksi = shape.get_reference_coords(m);

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                shape.calc_coords(&mut state, &mut x, ksi)?;

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
            shape.calc_coords(&mut state, &mut x, &ksi_in)?;
            gen_coords(&mut x_correct, &ksi_in, shape.class);
            assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol_in);
        }
        Ok(())
    }

    // Holds arguments for numerical differentiation of x w.r.t ξ => Jacobian
    struct ArgsNumJ {
        shape: Shape,      // auxiliary (copy) shape
        state: ShapeState, // auxiliary state
        at_ksi: Vec<f64>,  // at reference coord value
        ksi: Vec<f64>,     // temporary reference coord
        x: Vector,         // (space_ndim) coordinates at ξ
        i: usize,          // dimension index from 0 to space_ndim
        j: usize,          // dimension index from 0 to geo_ndim
    }

    // Computes xᵢ(ξ) with variable v := ξⱼ
    fn aux_jacobian(v: f64, args: &mut ArgsNumJ) -> f64 {
        args.ksi.copy_from_slice(&args.at_ksi);
        args.ksi[args.j] = v;
        args.shape.calc_coords(&mut args.state, &mut args.x, &args.ksi).unwrap();
        args.x[args.i]
    }

    #[test]
    fn calc_jacobian_works() -> Result<(), StrError> {
        // define dims and number of nodes
        let pairs = vec![
            // Lin
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5),
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
            // Tet
            (3, 4),
            (3, 10),
            // Hex
            (3, 8),
            (3, 20),
        ];

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
        for (geo_ndim, nnode) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let mut shape = Shape::new(space_ndim, geo_ndim, nnode)?;
            let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // set coordinates matrix
            set_coords_matrix(&mut shape);

            // set ξ within reference space
            let at_ksi = vec![0.25; shape.geo_ndim];

            // compute Jacobian, its inverse, and determinant
            let det_jac = shape.calc_jacobian(&mut state, &at_ksi)?;
            assert!(det_jac > 0.0);

            // set arguments for numerical integration
            let aux_shape = Shape::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
            let aux_state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
            let args = &mut ArgsNumJ {
                shape: aux_shape,
                state: aux_state,
                at_ksi,
                ksi: vec![0.0; shape.geo_ndim],
                x: Vector::new(shape.space_ndim),
                i: 0,
                j: 0,
            };
            set_coords_matrix(&mut args.shape);

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
        let mut shape = Shape::new(2, 1, 2)?;
        let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
        let l = 3.5;
        shape.set_node(0, 0, 0.0)?;
        shape.set_node(0, 1, 0.0)?;
        shape.set_node(1, 0, l)?;
        shape.set_node(1, 1, 0.0)?;
        let norm_jac_vec = shape.calc_jacobian(&mut state, &[0.0])?;
        assert_eq!(norm_jac_vec, l / 2.0); // 2.0 = length of shape in the reference space

        let mut shape = Shape::new(3, 2, 3)?;
        let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
        shape.set_node(0, 0, 0.0)?;
        shape.set_node(0, 1, 0.0)?;
        shape.set_node(0, 2, 0.0)?;
        shape.set_node(1, 0, 1.0)?;
        shape.set_node(1, 1, 0.0)?;
        shape.set_node(1, 2, 0.0)?;
        shape.set_node(2, 0, 0.5)?;
        shape.set_node(2, 1, 1.0)?;
        shape.set_node(2, 2, 1.0)?;
        let norm_jac_vec = shape.calc_jacobian(&mut state, &[0.0, 0.0])?;
        assert_eq!(norm_jac_vec, 0.0);
        Ok(())
    }

    #[test]
    fn calc_boundary_normal_works() -> Result<(), StrError> {
        // allocate surface shape
        let mut shape = Shape::new(3, 2, 17)?;
        let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;

        // set coordinates matrix
        set_coords_matrix(&mut shape);

        // compute boundary normal vector
        let at_ksi = vec![0.0; shape.geo_ndim];
        let mut normal = Vector::new(shape.space_ndim);
        shape.calc_boundary_normal(&mut state, &mut normal, &at_ksi)?;

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
        let mut shape = Shape::new(2, 1, 5)?;
        let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;

        // set coordinates matrix
        set_coords_matrix(&mut shape);

        // compute boundary normal vector
        let at_ksi = vec![0.0; shape.geo_ndim];
        let mut normal = Vector::new(shape.space_ndim);
        shape.calc_boundary_normal(&mut state, &mut normal, &at_ksi)?;

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
            let mut shape = Shape::new(space_ndim, geo_ndim, nnode)?;

            // set coordinates matrix
            set_coords_matrix_parallel(&mut shape)?;

            // loop over edges
            for e in 0..shape.nedge {
                let edge_nnode = shape.edge_nnode;
                let mut edge_shape = Shape::new(space_ndim, 1, shape.edge_nnode)?;
                let mut edge_state = ShapeState::new(edge_shape.space_ndim, edge_shape.geo_ndim, edge_shape.nnode)?;

                // set edge coordinates
                for i in 0..edge_nnode {
                    for j in 0..space_ndim {
                        let m = shape.get_edge_node_id(e, i);
                        edge_shape.set_node(i, j, shape.coords_transp[j][m])?;
                    }
                }

                // calc normal vector
                edge_shape.calc_boundary_normal(&mut edge_state, &mut normal, ksi)?;
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
            let mut shape = Shape::new(space_ndim, geo_ndim, nnode)?;

            // set coordinates matrix
            set_coords_matrix_parallel(&mut shape)?;

            // loop over faces
            for f in 0..shape.nface {
                let face_nnode = shape.face_nnode;
                let mut face_shape = Shape::new(space_ndim, 2, shape.face_nnode)?;
                let mut face_state = ShapeState::new(face_shape.space_ndim, face_shape.geo_ndim, face_shape.nnode)?;

                // set face coordinates
                for i in 0..face_nnode {
                    for j in 0..space_ndim {
                        let m = shape.get_face_node_id(f, i);
                        face_shape.set_node(i, j, shape.coords_transp[j][m])?;
                    }
                }

                // calc normal vector
                face_shape.calc_boundary_normal(&mut face_state, &mut normal, ksi)?;
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
            let mut shape = Shape::new(space_ndim, geo_ndim, nnode)?;
            let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // set coordinates matrix
            set_coords_matrix(&mut shape);

            // loop over nodes of shape
            let mut x = Vector::new(shape.space_ndim);
            let mut ksi = vec![0.0; shape.geo_ndim];
            for m in 0..shape.nnode {
                // get ξᵐ corresponding to node m
                let ksi_ref = shape.get_reference_coords(m);

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                shape.calc_coords(&mut state, &mut x, ksi_ref)?;

                // compute approximation of the inverse mapping ξᵐ(xᵐ)
                let nit = shape.approximate_ksi(&mut state, &mut ksi, &x, 10, 1e-14)?;

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
            shape.calc_coords(&mut state, &mut x, &ksi_in)?;
            shape.approximate_ksi(&mut state, &mut ksi, &x, 10, 1e-14)?;
            assert_vec_approx_eq!(ksi, ksi_in, tol);
        }
        Ok(())
    }

    // Holds arguments for numerical differentiation of N w.r.t x => G (gradient) matrix
    struct ArgsNumG {
        shape: Shape,      // auxiliary (copy) shape
        state: ShapeState, // auxiliary state
        at_x: Vector,      // at x coord value
        x: Vector,         // temporary x coord
        ksi: Vec<f64>,     // temporary reference coord
        m: usize,          // node index from 0 to nnode
        j: usize,          // dimension index from 0 to space_ndim
    }

    // Computes Nᵐ(ξ(x)) with variable v := xⱼ
    fn aux_grad(v: f64, args: &mut ArgsNumG) -> f64 {
        copy_vector(&mut args.x, &args.at_x).unwrap();
        args.x[args.j] = v;
        args.shape
            .approximate_ksi(&mut args.state, &mut args.ksi, &args.x, 10, 1e-14)
            .unwrap();
        args.shape.calc_interp(&mut args.state, &args.ksi);
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
            let space_ndim = geo_ndim;
            let mut shape = Shape::new(space_ndim, geo_ndim, nnode)?;
            let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();

            // set coordinates matrix
            set_coords_matrix(&mut shape);

            // set ξ within reference space
            let at_ksi = vec![0.25; shape.geo_ndim];

            // compute x corresponding to ξ using the isoparametric formula
            let mut at_x = Vector::new(shape.space_ndim);
            shape.calc_coords(&mut state, &mut at_x, &at_ksi)?;

            // compute gradient
            let det_jac = shape.calc_gradient(&mut state, &at_ksi)?;
            assert!(det_jac > 0.0);

            // set arguments for numerical integration
            let aux_shape = Shape::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
            let aux_state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
            let args = &mut ArgsNumG {
                shape: aux_shape,
                state: aux_state,
                at_x,
                x: Vector::new(shape.space_ndim),
                ksi: vec![0.0; shape.geo_ndim],
                m: 0,
                j: 0,
            };
            set_coords_matrix(&mut args.shape);

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
    fn calc_int_points_coords() -> Result<(), StrError> {
        let mut shape = Shape::new(1, 1, 2)?;
        let mut state = ShapeState::new(shape.space_ndim, shape.geo_ndim, shape.nnode)?;
        let (xa, xb) = (2.0, 5.0);
        shape.set_node(0, 0, xa)?;
        shape.set_node(1, 0, xb)?;
        let int_points = shape.calc_int_points_coords(&mut state)?;
        assert_eq!(int_points.len(), 2);
        let ksi_a = -1.0 / SQRT_3;
        let ksi_b = 1.0 / SQRT_3;
        let x_ksi_a = xa + (ksi_a + 1.0) * (xb - xa) / 2.0;
        let x_ksi_b = xa + (ksi_b + 1.0) * (xb - xa) / 2.0;
        assert_eq!(int_points[0].dim(), 1);
        assert_eq!(int_points[1].dim(), 1);
        assert_approx_eq!(int_points[0][0], x_ksi_a, 1e-15);
        assert_approx_eq!(int_points[1][0], x_ksi_b, 1e-15);
        Ok(())
    }

    #[test]
    fn clone_and_serialize_work() -> Result<(), StrError> {
        // original
        let mut shape = Shape::new(1, 1, 2)?;
        shape.set_node(0, 0, 1.0)?;
        shape.set_node(1, 0, 2.0)?;
        let before = format!("{:?}", shape);
        assert_eq!(
            format!("{}", shape.coords_transp),
            "┌     ┐\n\
             │ 1 2 │\n\
             └     ┘"
        );

        // cloned
        let mut cloned = shape.clone();
        assert_eq!(format!("{:?}", cloned), before);
        cloned.set_node(0, 0, 3.0)?;
        cloned.set_node(1, 0, 4.0)?;
        assert_eq!(
            format!("{}", cloned.coords_transp),
            "┌     ┐\n\
             │ 3 4 │\n\
             └     ┘"
        );
        assert_eq!(format!("{:?}", shape), before);

        // serialize to JSON
        let json = serde_json::to_string_pretty(&shape).map_err(|_| "json encode failed")?;
        assert_eq!(
            format!("{}", json),
            "{\n\
             \x20\x20\"class\": \"Lin\",\n\
             \x20\x20\"kind\": \"Lin2\",\n\
             \x20\x20\"space_ndim\": 1,\n\
             \x20\x20\"geo_ndim\": 1,\n\
             \x20\x20\"nnode\": 2,\n\
             \x20\x20\"nedge\": 0,\n\
             \x20\x20\"nface\": 0,\n\
             \x20\x20\"edge_nnode\": 0,\n\
             \x20\x20\"face_nnode\": 0,\n\
             \x20\x20\"face_nedge\": 0,\n\
             \x20\x20\"fn_interp\": 1002,\n\
             \x20\x20\"fn_deriv\": 1002,\n\
             \x20\x20\"coords_transp\": {\n\
             \x20\x20\x20\x20\"nrow\": 1,\n\
             \x20\x20\x20\x20\"ncol\": 2,\n\
             \x20\x20\x20\x20\"data\": [\n\
             \x20\x20\x20\x20\x20\x201.0,\n\
             \x20\x20\x20\x20\x20\x202.0\n\
             \x20\x20\x20\x20]\n\
             \x20\x20},\n\
             \x20\x20\"ok_last_coord\": true,\n\
             \x20\x20\"min_coords\": [\n\
             \x20\x20\x20\x201.0\n\
             \x20\x20],\n\
             \x20\x20\"max_coords\": [\n\
             \x20\x20\x20\x202.0\n\
             \x20\x20]\n\
            }"
        );

        // deserialize from JSON
        let shape_json: Shape = serde_json::from_str(&json).map_err(|_| "json decode failed")?;
        assert_eq!(format!("{:?}", shape_json), before);

        // serialize to BIN
        let mut bin = Vec::new();
        let mut ser = rmp_serde::Serializer::new(&mut bin);
        shape.serialize(&mut ser).map_err(|_| "bin encode failed")?;
        assert!(bin.len() > 0);

        // deserialize from BIN
        let mut des = rmp_serde::Deserializer::new(&bin[..]);
        let shape_bin: Shape = Deserialize::deserialize(&mut des).map_err(|_| "bin decode failed")?;
        assert_eq!(format!("{:?}", shape_bin), before);
        Ok(())
    }
}
