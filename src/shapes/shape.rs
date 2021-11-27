use super::{Hex20, Hex8, Qua4, Qua8};
use crate::StrError;
use russell_lab::{inverse, mat_mat_mul, mat_vec_mul, Matrix, Vector};

/// Defines the kind of shape
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum GeoKind {
    Lin2,
    Lin3,
    Lin4,
    Lin5,
    Tri3,
    Tri6,
    Tri10,
    Tri15,
    Qua4,
    Qua8,
    Qua9,
    Qua12,
    Qua16,
    Qua17,
    Tet4,
    Tet10,
    Hex8,
    Hex20,
}

// alias for interpolation functions
type FnInterp = fn(&mut Vector, &Vector);

// alias for derivative of interpolation functions
type FnDeriv = fn(&mut Matrix, &Vector);

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
/// * `npoint` -- number of points
/// * `nedge` -- number of edges
/// * `nface` -- number of faces
/// * `edge_npoint` -- edge's number of points
/// * `face_npoint` -- face's number of points
/// * `face_nedge` -- face's number of edges
///
/// # Formulae
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
/// vector of reference coordinates, `Nm` are the (npoint) interpolation functions,
/// and `xm` are the (npoint) coordinates of each m-point of the geometric shape.
///
/// Given an (npoint,space_ndim) matrix of coordinates X, we can calculate the
/// (space_ndim) vector of coordinates x by means of
///
/// ```text
/// x = Xᵀ ⋅ N
/// ```
///
/// where `N` is an (npoint) array formed with all `Nm`.
///
/// ## First case: geo_ndim = space_ndim
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
/// →       dNᵐ(ξ)
/// Lᵐ(ξ) = ——————
///            →
///           dξ
/// ```
///
/// are the derivatives of each interpolation function `Nm` with respect to the
/// reference coordinate. `Lm` are (geo_ndim) vectors and can be organized in
/// an (npoint,geo_ndim) matrix `L` of "local" derivatives.
///
/// We can write the Jacobian in matrix notation as follows
///
/// ```text
/// J = Xᵀ · L
/// ```
///
/// where X is the (npoint,space_ndim) matrix of coordinates and L is the (npoint,geo_ndim) matrix.
///
/// We define the gradient of interpolation functions (i.e., derivatives of interpolation
/// functions w.r.t real coordinates) by means of
///
/// ```text
///             →
/// →       dNᵐ(ξ)
/// Gᵐ(ξ) = ——————
///            →
///           dx
/// ```
///
/// which can be organized in an (npoint,space_ndim) matrix `G`.
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
/// where G is an (npoint,space_ndim) matrix.
///
/// We can then perform integrations in the reference space as follows:
///
/// ```text
///       ⌠   →       ⌠   → →            →
/// res = │ f(x) dΩ = │ f(x(ξ)) ⋅ det(J)(ξ) dΩ
///       ⌡           ⌡
///       Ω           Ωref
/// ```
///
/// which is replaced by numerical integration according to:
///
/// ```text
///       nip-1  →             →
/// res ≈  Σ   f(ιp)) ⋅ det(J)(ιp) ⋅ wp
///       p=0
/// ```
///
/// where `nip` is the number of integration points, `ιp := ξp` is the reference
/// coordinate of the integration point, and `wp` is the weight attached to the
/// p-th integration point.
///
/// ## Second case: geo_ndim < space_ndim
///
/// This case corresponds to a line or surface in a multi-dimensional space.
/// For instance, a line in 2D or 3D or a triangle in 3D.
///
/// If `geo_ndim < space_ndim`, we must consider some particular cases. For now, the only
/// cases we can handle are:
///
/// * `geo_ndim = 1` and `space_ndim = 2 or 3` -- Line in multi-dimensions (boundary or not)
/// * `geo_ndim = 2` and `space_ndim = 3` -- 3D surface (boundary only)
///
/// ### Line in multi-dimensions (geo_ndim = 1 and space_ndim > 1)
///
/// In this case, the Jacobian equals the (space_ndim,1) base vector `g1` aligned
/// with the line element, i.e.,
///
/// ```text
///                     →
/// →          →       dx
/// Jline(ξ) = g1(ξ) = —— = Xᵀ · L
///                    dξ
/// ```
///
/// We further consider a parametric coordinate `ell` which varies from 0 to `ell_max`
/// (the length of the line) according to
///
/// ```text
///                  ell_max
/// ell(ξ) = (1 + ξ) ———————
///                     2
///
///          2 · ell
/// ξ(ell) = ——————— - 1
///          ell_max
/// ```
///
/// ```text
/// 0 ≤ ell ≤ ell_max
///
/// -1 ≤ ξ ≤ +1
/// ```
///
/// Then, we replace the line integrals over the real space with line integrals over
/// the mapped space as follows:
///
/// ```text
///       ⌠               ⌠
/// res = │ f(ell) dell = │ f(ξ(ell)) ⋅ ||Jline||(ξ) dξ
///       ⌡               ⌡
///       Ω               Ωref
/// ```
///
/// where `||Jline||` is the Euclidean norm of `Jline`.
///
/// The above integral is replaced by numerical integration according to:
///
/// ```text
///       nip-1      →               →
/// res ≈  Σ   f(ell(ιp)) ⋅ ||Jline||(ιp) ⋅ wp
///       p=0
/// ```
///
/// where `nip` is the number of integration points, `ιp := ξp` is the reference coordinate
/// of the integration point, and `wp` is the weight attached to the p-th integration point.
/// Furthermore, we can use the definition of `ell` to compute `ell(ξ:=ιp)`.
///
/// **Boundary line in 2D (geo_ndim = 1 and space_ndim = 2)**
///
/// If the line defines a boundary in 2D, we compute a normal vector by means of
///
/// ```text
/// →   →    →
/// n = e3 × g1
/// ```
///
/// Then, we can compute the following boundary integral
///
/// ```text
///       ⌠             →        ⌠           →    →
/// res = │ q(ell) unit_n dell = │ q(ell) ⋅ (e3 × g1) dξ
///       ⌡                      ⌡
///       Γ                      Γref
/// ```
///
/// where `unit_n` is the unit normal vector.
///
/// The above integral is replaced by numerical integration according to:
///
/// ```text
///       nip-1      →       →  →     →  →
/// res ≈  Σ   q(ell(ιp)) ⋅ (e3(ιp) × g1(ιp)) ⋅ wp
///       p=0
/// ```
///
/// where we can use the definition of `ell` to compute `ell(ξ:=ιp)`.
///
/// ### 3D surface (boundary only) (geo_ndim = 2 and space_ndim = 3)
///
/// In this case, we use the normal vector to the surface to replace the surface
/// integrations in the real space with integrations over the mapped space.
/// Instead of using the Jacobian, we use the normal vector.
///
/// Considering convective coordinates (ξ1,ξ2) on the surface, we compute the
/// following base vectors
///
/// ```text
///          →
/// →  →    dx
/// g1(ξ) = ——— = first_column(Jsurf)
///         dξ1
///
///          →
/// →  →    dx
/// g2(ξ) = ——— = second_column(Jsurf)
///         dξ1
/// ```
///
/// where we defined
///
/// ```text
/// Jsurf = Xᵀ · L
/// ```
///
/// We can then perform integrations in the reference space as follows:
///
/// ```text
///       ⌠   →       →      ⌠   → →      →    →
/// res = │ f(x) unit_n dA = │ f(x(ξ)) ⋅ (g1 × g2) dξ1 dξ2
///       ⌡                  ⌡
///       Γ                  Γref
/// ```
///
/// where `unit_n` is the unit normal vector.
///
/// The above integral is replaced by numerical integration according to:
///
/// ```text
///       nip-1   →      →   →    →   →
/// res ≈  Σ   f(ιp)) ⋅ (g1(ιp) × g2(ιp)) ⋅ wp
///       p=0
/// ```
///
/// ## Third case: geo_ndim > space_ndim
///
/// This case (e.g., a cube in the 1D space) is not allowed.
pub struct Shape {
    /// Geometry kind
    pub kind: GeoKind,

    /// Space ndim
    pub space_ndim: usize,

    /// Geometry ndim
    pub geo_ndim: usize,

    /// Number of points
    pub npoint: usize,

    /// Number of edges
    pub nedge: usize,

    /// Number of faces
    pub nface: usize,

    /// Edge's number of points
    pub edge_npoint: usize,

    /// Face's number of points
    pub face_npoint: usize,

    /// Face's number of edges
    pub face_nedge: usize,

    /// Function to calculate interpolation functions
    pub fn_interp: FnInterp,

    /// Function to calculate local derivatives (w.r.t. ksi) of interpolation functions
    pub fn_deriv: FnDeriv,

    /// Array N: (npoint) interpolation functions at reference coordinate ksi
    pub interp: Vector,

    /// Matrix L: (npoint,geo_ndim) derivatives of interpolation functions w.r.t reference coordinate ksi
    pub deriv: Matrix,

    /// Matrix J: (space_ndim,geo_ndim) Jacobian matrix
    pub jacobian: Matrix,

    /// Matrix inv(J): (space_ndim,space_ndim) Inverse Jacobian matrix (only if geo_ndim == space_ndim)
    pub inv_jacobian: Matrix,

    /// Matrix G: (npoint,space_ndim) Gradient of shape functions (only if geo_ndim == space_ndim)
    pub gradient: Matrix,

    /// Matrix X: (space_ndim,npoint) transposed coordinates matrix (real space)
    coords_transp: Matrix,
}

impl Shape {
    /// Creates a new geometric shape
    ///
    /// # Input
    ///
    /// * `space_ndim` -- the space dimension (1, 2, or 3)
    /// * `geo_ndim` -- the dimension of the shape; e.g. a 1D line in a 3D space
    /// * `npoint` -- the number of points defining the shape
    ///
    /// # Note
    ///
    /// Some methods require that the coordinates matrix be set first.
    /// This can be accomplished by calling the `set_coords` method.
    pub fn new(space_ndim: usize, geo_ndim: usize, npoint: usize) -> Result<Self, StrError> {
        let interp = Vector::new(npoint);
        let deriv = Matrix::new(npoint, geo_ndim);
        let jacobian = Matrix::new(space_ndim, geo_ndim);
        let (inv_jacobian, gradient) = if geo_ndim == space_ndim {
            (Matrix::new(space_ndim, space_ndim), Matrix::new(npoint, space_ndim))
        } else {
            (Matrix::new(0, 0), Matrix::new(0, 0))
        };
        let coords_transp = Matrix::new(space_ndim, npoint);
        match (geo_ndim, npoint) {
            (1, 2) => Err("Lin2 is not available yet"),
            (1, 3) => Err("Lin3 is not available yet"),
            (1, 4) => Err("Lin4 is not available yet"),
            (1, 5) => Err("Lin5 is not available yet"),
            (2, 3) => Err("Tri3 is not available yet"),
            (2, 6) => Err("Tri6 is not available yet"),
            (2, 10) => Err("Tri10 is not available yet"),
            (2, 15) => Err("Tri15 is not available yet"),
            (2, 4) => Ok(Shape {
                kind: GeoKind::Qua4,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Qua4::NEDGE,
                nface: Qua4::NFACE,
                edge_npoint: Qua4::EDGE_NPOINT,
                face_npoint: Qua4::FACE_NPOINT,
                face_nedge: Qua4::FACE_NEDGE,
                fn_interp: Qua4::calc_interp,
                fn_deriv: Qua4::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                coords_transp,
            }),
            (2, 8) => Ok(Shape {
                kind: GeoKind::Qua8,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Qua8::NEDGE,
                nface: Qua8::NFACE,
                edge_npoint: Qua8::EDGE_NPOINT,
                face_npoint: Qua8::FACE_NPOINT,
                face_nedge: Qua8::FACE_NEDGE,
                fn_interp: Qua8::calc_interp,
                fn_deriv: Qua8::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                coords_transp,
            }),
            (2, 9) => Err("Qua9 is not available yet"),
            (2, 12) => Err("Qua12 is not available yet"),
            (2, 16) => Err("Qua16 is not available yet"),
            (2, 17) => Err("Qua17 is not available yet"),
            (3, 4) => Err("Tet4 is not available yet"),
            (3, 10) => Err("Tet10 is not available yet"),
            (3, 8) => Ok(Shape {
                kind: GeoKind::Hex8,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Hex8::NEDGE,
                nface: Hex8::NFACE,
                edge_npoint: Hex8::EDGE_NPOINT,
                face_npoint: Hex8::FACE_NPOINT,
                face_nedge: Hex8::FACE_NEDGE,
                fn_interp: Hex8::calc_interp,
                fn_deriv: Hex8::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                coords_transp,
            }),
            (3, 20) => Ok(Shape {
                kind: GeoKind::Hex20,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Hex20::NEDGE,
                nface: Hex20::NFACE,
                edge_npoint: Hex20::EDGE_NPOINT,
                face_npoint: Hex20::FACE_NPOINT,
                face_nedge: Hex20::FACE_NEDGE,
                fn_interp: Hex20::calc_interp,
                fn_deriv: Hex20::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                coords_transp,
            }),
            _ => Err("(space_ndim, geo_ndim, npoint) combination is invalid"),
        }
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
    ///     └               ┘_(npoint,space_ndim)
    /// ```
    ///
    /// where `M = npoint - 1`
    ///
    /// # Input
    ///
    /// * `m` -- point index in form 0 to npoint-1
    /// * `j` -- dimension index from 0 to space_ndim-1
    /// * `value` -- the X(m,j) component
    pub fn set_point_coord(&mut self, m: usize, j: usize, value: f64) -> Result<(), StrError> {
        if m >= self.npoint {
            return Err("index of point is invalid");
        }
        if j >= self.space_ndim {
            return Err("index of space dimension is invalid");
        }
        self.coords_transp[j][m] = value;
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
    /// # Warning
    ///
    /// You must set the coordinates matrix first, otherwise the computations
    /// will generate wrong results.
    pub fn calc_coords(&mut self, x: &mut Vector, ksi: &Vector) -> Result<(), StrError> {
        (self.fn_interp)(&mut self.interp, ksi);
        mat_vec_mul(x, 1.0, &self.coords_transp, &self.interp)
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
    /// * `ksi` -- ξ reference coordinate
    ///
    /// # Output
    ///
    /// If `geo_ndim = space_ndim`, returns the determinant of the Jacobian.
    /// Otherwise, returns zero.
    ///
    /// # Updated variables
    ///
    /// * `jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    /// * `inv_jacobian` -- If `geo_ndim = space_ndim`: inverse Jacobian matrix (space_ndim,space_ndim)
    ///
    /// # Warning
    ///
    /// You must set the coordinates matrix first, otherwise the computations
    /// will generate wrong results.
    pub fn calc_jacobian(&mut self, ksi: &Vector) -> Result<f64, StrError> {
        (self.fn_deriv)(&mut self.deriv, ksi);
        mat_mat_mul(&mut self.jacobian, 1.0, &self.coords_transp, &self.deriv)?;
        if self.geo_ndim == self.space_ndim {
            inverse(&mut self.inv_jacobian, &self.jacobian)
        } else {
            Ok(0.0)
        }
    }

    /// Calculates the gradient of the interpolation functions
    ///
    /// **Note:** This function works with `geo_ndim == space_ndim` only
    ///
    /// ```text
    ///             →
    /// →       dNᵐ(ξ)
    /// Gᵐ(ξ) = ——————
    ///            →
    ///           dx
    /// ```
    ///
    /// which can be organized in an (npoint,space_ndim) matrix `G` as follows
    ///
    /// ```text
    /// G = L · J⁻¹
    /// ```
    ///
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate
    ///
    /// # Output
    ///
    /// Returns the determinant of the Jacobian.
    ///
    /// # Updated variables
    ///
    /// * `jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    /// * `inv_jacobian` -- inverse Jacobian matrix (space_ndim,space_ndim)
    /// * `gradient` -- gradient matrix (npoint,space_ndim)
    ///
    /// # Warning
    ///
    /// You must set the coordinates matrix first, otherwise the computations
    /// will generate wrong results.
    pub fn calc_gradient(&mut self, ksi: &Vector) -> Result<f64, StrError> {
        if self.geo_ndim != self.space_ndim {
            return Err("geo_ndim must equal space_ndim");
        }
        let det_jac = self.calc_jacobian(ksi)?;
        mat_mat_mul(&mut self.gradient, 1.0, &self.deriv, &self.inv_jacobian)?;
        Ok(det_jac)
    }

    // --- getters ------------------------------------------------------------------------------

    /// Returns the local id of point on edge
    ///
    /// # Input
    ///
    /// * `e` -- index of edge in [0, nedge-1]
    /// * `i` -- index of local point [0, edge_npoint-1]
    pub fn get_edge_point_id(&self, e: usize, i: usize) -> usize {
        match self.kind {
            GeoKind::Qua4 => Qua4::EDGE_POINT_IDS[e][i],
            GeoKind::Qua8 => Qua8::EDGE_POINT_IDS[e][i],
            GeoKind::Hex8 => Hex8::EDGE_POINT_IDS[e][i],
            GeoKind::Hex20 => Hex20::EDGE_POINT_IDS[e][i],
            _ => panic!("ShapeKind is not available yet"),
        }
    }

    /// Returns the local id of point on face
    ///
    /// # Input
    ///
    /// * `f` -- index of face in [0, nface-1]
    /// * `i` -- index of local point [0, face_npoint-1]
    pub fn get_face_point_id(&self, f: usize, i: usize) -> usize {
        match self.kind {
            GeoKind::Qua4 => 0,
            GeoKind::Qua8 => 0,
            GeoKind::Hex8 => Hex8::FACE_POINT_IDS[f][i],
            GeoKind::Hex20 => Hex20::FACE_POINT_IDS[f][i],
            _ => panic!("ShapeKind is not available yet"),
        }
    }

    /// Returns the local point id on an edge on the face
    ///
    /// # Input
    ///
    /// * `f` -- index of face in [0, nface-1]
    /// * `k` -- index of face's edge (not the index of cell's edge) in [0, face_nedge-1]
    /// * `i` -- index of local point [0, edge_npoint-1]
    pub fn get_face_edge_point_id(&self, f: usize, k: usize, i: usize) -> usize {
        match self.kind {
            GeoKind::Qua4 => 0,
            GeoKind::Qua8 => 0,
            GeoKind::Hex8 => Hex8::FACE_EDGE_POINT_IDS[f][k][i],
            GeoKind::Hex20 => Hex20::FACE_EDGE_POINT_IDS[f][k][i],
            _ => panic!("ShapeKind is not available yet"),
        }
    }

    /// Returns the reference coordinates at point m
    ///
    /// # Output
    ///
    /// * `ksi` -- (geo_ndim) reference coordinates `ξᵐ` at point m
    pub fn get_reference_coords(&self, ksi: &mut Vector, m: usize) {
        match self.kind {
            GeoKind::Qua4 => {
                ksi[0] = Qua4::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Qua4::POINT_REFERENCE_COORDS[m][1];
            }
            GeoKind::Qua8 => {
                ksi[0] = Qua8::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Qua8::POINT_REFERENCE_COORDS[m][1];
            }
            GeoKind::Hex8 => {
                ksi[0] = Hex8::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Hex8::POINT_REFERENCE_COORDS[m][1];
                ksi[2] = Hex8::POINT_REFERENCE_COORDS[m][2];
            }
            GeoKind::Hex20 => {
                ksi[0] = Hex20::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Hex20::POINT_REFERENCE_COORDS[m][1];
                ksi[2] = Hex20::POINT_REFERENCE_COORDS[m][2];
            }
            _ => panic!("ShapeKind is not available yet"),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use russell_chk::*;
    use std::collections::HashMap;

    struct TestData {
        shapes: Vec<Shape>,
        tolerances_fn_interp: HashMap<GeoKind, f64>,
        tolerances_fn_deriv: HashMap<GeoKind, f64>,
        tolerances_calc_coords: HashMap<GeoKind, f64>,
        tolerances_calc_coords_high: HashMap<GeoKind, f64>, // high tol for center
        tolerances_calc_jacobian: HashMap<GeoKind, f64>,
        tolerances_calc_gradient: HashMap<GeoKind, f64>,
    }

    // Generates test data
    fn gen_test_data() -> TestData {
        let mut data = TestData {
            shapes: Vec::new(),
            tolerances_fn_interp: HashMap::new(),
            tolerances_fn_deriv: HashMap::new(),
            tolerances_calc_coords: HashMap::new(),
            tolerances_calc_coords_high: HashMap::new(),
            tolerances_calc_jacobian: HashMap::new(),
            tolerances_calc_gradient: HashMap::new(),
        };
        let pairs = vec![(2, 4), (2, 8), (3, 8), (3, 20)];
        for (geo_ndim, npoint) in pairs {
            let space_ndim = geo_ndim;
            data.shapes.push(Shape::new(space_ndim, geo_ndim, npoint).unwrap());
            let shape = data.shapes.last().unwrap();
            data.tolerances_fn_interp.insert(shape.kind, 1e-15);
            data.tolerances_fn_deriv.insert(shape.kind, 1e-12);
            data.tolerances_calc_coords.insert(shape.kind, 1e-15);
        }
        data.tolerances_calc_coords_high.insert(GeoKind::Qua4, 0.15);
        data.tolerances_calc_coords_high.insert(GeoKind::Qua8, 1e-14);
        data.tolerances_calc_coords_high.insert(GeoKind::Hex8, 0.15);
        data.tolerances_calc_coords_high.insert(GeoKind::Hex20, 1e-14);
        data.tolerances_calc_jacobian.insert(GeoKind::Qua4, 1e-11);
        data.tolerances_calc_jacobian.insert(GeoKind::Qua8, 1e-11);
        data.tolerances_calc_jacobian.insert(GeoKind::Hex8, 1e-11);
        data.tolerances_calc_jacobian.insert(GeoKind::Hex20, 1e-11);
        data
    }

    const RMIN: f64 = 1.0;
    const RMAX: f64 = 10.0;
    const AMIN: f64 = 30.0 * std::f64::consts::PI / 180.0;
    const AMAX: f64 = 60.0 * std::f64::consts::PI / 180.0;

    // Generate coordinates
    //
    //                         |            /
    //                                     /
    //                 ,,****#"|""#****,, /
    //            ,***""'              `""**,             r
    //         ,**#"           |        /??"#**,       ,~'
    //       ,*#"                      /??????"#*,  ,~'
    //     ,*#"         _,,**#"|"#**,,_?????????"#*,
    //    ,#"         ,*#"'         `"#*,??????,~ "#,
    //   ,#'        ,*"        |    /   "*,?,~'    `#,
    //  ,#'        *"              / αmax "*        `#,
    //  *'        *'           |  /    ,~' `*        `*
    //  #         #              /  ,~'     #         #
    // ,#        ,#            | ,~' αmin    #,        #,
    // #         #             o - - - - - - #  - - -  # - - -
    // `#        '#                    rmin #'   rmax #'
    //  #         #,                       ,#         #
    //  #,         #*                     *#         ,#
    //  `#,         "#*                 *#"         ,#'
    //   `#,          "#*,_         _,*#"          ,#'
    //    `#*           `""#*******#""'           *#'
    //     `#**                                 **#'
    //       "#**                             **#"
    //         `"#**,                     ,**#"'
    //            `"##**,             ,**##"'
    //                 ``""##*****##""''
    //
    // if 3D => cylindrical coordinates
    //
    // r(ξ₀,ξ₁,ξ₂) = rmin + (1+ξ₀) * Δr/2
    // α(ξ₀,ξ₁,ξ₂) = αmin + (1+ξ₁) * Δα/2
    // z(ξ₀,ξ₁,ξ₂) = ξ₂
    //
    // x₀ := r * cos(α)
    // x₁ := r * sin(α)
    // x₂ := z
    fn gen_coords(x: &mut Vector, ksi: &Vector) {
        assert_eq!(x.dim(), ksi.dim());
        let r = RMIN + (1.0 + ksi[0]) * (RMAX - RMIN) / 2.0;
        if x.dim() == 1 {
            x[0] = r;
            return;
        }
        let a = AMIN + (1.0 + ksi[1]) * (AMAX - AMIN) / 2.0;
        x[0] = r * f64::cos(a);
        x[1] = r * f64::sin(a);
        if x.dim() == 3 {
            x[2] = ksi[2];
        }
    }

    // Generate reference coordinates
    //
    // if 3D => cylindrical coordinates
    //
    // ξ₀(r,α,z) = (2/Δr) * (r-rmin) - 1
    // ξ₁(r,α,z) = (2/Δα) * (α-αmin) - 1
    // ξ₂(r,α,z) = z
    //
    // r(x₀,x₁,x₂) = sqrt(x₀ x₀ + x₁ x₁)
    // α(x₀,x₁,x₂) = atan(x₁ / x₀)
    // z(x₀,x₁,x₂) = x₂
    fn gen_ref_coords(ksi: &mut Vector, x: &Vector) {
        assert_eq!(x.dim(), ksi.dim());
        let r = f64::sqrt(x[0] * x[0] + x[1] * x[1]);
        ksi[0] = (2.0 / (RMAX - RMIN)) * (r - RMIN) - 1.0;
        if x.dim() == 1 {
            return;
        }
        let a = f64::atan2(x[1], x[0]);
        ksi[1] = (2.0 / (AMAX - AMIN)) * (a - AMIN) - 1.0;
        if x.dim() == 3 {
            ksi[2] = x[2];
        }
    }

    // Generates matrix of coordinates for geo_ndim = space_ndim
    fn set_coords_matrix(shape: &mut Shape) {
        assert_eq!(shape.geo_ndim, shape.space_ndim);
        let mut x = Vector::new(shape.space_ndim);
        let mut ksi = Vector::new(shape.geo_ndim);
        for m in 0..shape.npoint {
            shape.get_reference_coords(&mut ksi, m);
            gen_coords(&mut x, &ksi);
            for j in 0..shape.geo_ndim {
                shape.set_point_coord(m, j, x[j]).unwrap();
            }
        }
    }

    // Holds arguments for numerical differentiation of N w.r.t ξ => L (deriv) matrix
    struct ArgsInterpVarKsi {
        shape: Shape,   // auxiliary (copy) shape
        at_ksi: Vector, // at reference coord value
        ksi: Vector,    // temporary reference coord
        m: usize,       // point index from 0 to npoint
        j: usize,       // dimension index from 0 to geom_ndim
    }

    // Holds arguments for numerical differentiation of x w.r.t ξ => Jacobian
    struct ArgsCoordsVarKsi {
        shape: Shape,   // auxiliary (copy) shape
        at_ksi: Vector, // at reference coord value
        ksi: Vector,    // temporary reference coord
        x: Vector,      // (space_ndim) coordinates at ξ
        i: usize,       // dimension index from 0 to space_ndim
        j: usize,       // dimension index from 0 to geo_ndim
    }

    // Holds arguments for numerical differentiation of N w.r.t x => G (gradient) matrix
    struct ArgsInterpVarX {
        shape: Shape, // auxiliary (copy) shape
        at_x: Vector, // at x coord value
        x: Vector,    // temporary x coord
        ksi: Vector,  // temporary reference coord
        m: usize,     // point index from 0 to npoint
        j: usize,     // dimension index from 0 to space_ndim
    }

    // Computes Nᵐ(ξ) with variable v := ξⱼ
    fn calc_interp_m_var_ksi(v: f64, args: &mut ArgsInterpVarKsi) -> f64 {
        for j in 0..args.shape.geo_ndim {
            args.ksi[j] = args.at_ksi[j];
        }
        args.ksi[args.j] = v;
        (args.shape.fn_interp)(&mut args.shape.interp, &args.ksi);
        args.shape.interp[args.m]
    }

    // Computes xᵢ(ξ) with variable v := ξⱼ
    fn calc_coord_i_var_ksi(v: f64, args: &mut ArgsCoordsVarKsi) -> f64 {
        for j in 0..args.shape.geo_ndim {
            args.ksi[j] = args.at_ksi[j];
        }
        args.ksi[args.j] = v;
        args.shape.calc_coords(&mut args.x, &args.ksi).unwrap();
        args.x[args.i]
    }

    // Computes Nᵐ(ξ(x)) with variable v := xⱼ
    fn calc_interp_m_var_x(v: f64, args: &mut ArgsInterpVarX) -> f64 {
        for j in 0..args.shape.space_ndim {
            args.x[j] = args.at_x[j];
        }
        args.x[args.j] = v;
        gen_ref_coords(&mut args.ksi, &args.x);
        (args.shape.fn_interp)(&mut args.shape.interp, &args.ksi);
        args.shape.interp[args.m]
    }

    #[test]
    fn fn_interp_works() -> Result<(), StrError> {
        let mut data = gen_test_data();
        for shape in &mut data.shapes {
            let tol = *data.tolerances_fn_interp.get(&shape.kind).unwrap();
            println!("{:?}: tol={:e}", shape.kind, tol);
            let mut ksi = Vector::new(shape.geo_ndim);
            for m in 0..shape.npoint {
                shape.get_reference_coords(&mut ksi, m);
                (shape.fn_interp)(&mut shape.interp, &ksi);
                for n in 0..shape.npoint {
                    let interp_mn = shape.interp[n];
                    if m == n {
                        assert_approx_eq!(interp_mn, 1.0, tol);
                    } else {
                        assert_approx_eq!(interp_mn, 0.0, tol);
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn fn_deriv_works() -> Result<(), StrError> {
        let mut data = gen_test_data();
        for shape in &mut data.shapes {
            let tol = *data.tolerances_fn_deriv.get(&shape.kind).unwrap();
            println!("{:?}: tol={:e}", shape.kind, tol);
            let mut at_ksi = Vector::new(shape.geo_ndim);
            for j in 0..shape.geo_ndim {
                at_ksi[j] = 0.25;
            }
            let args = &mut ArgsInterpVarKsi {
                shape: Shape::new(shape.space_ndim, shape.geo_ndim, shape.npoint)?,
                at_ksi,
                ksi: Vector::new(shape.geo_ndim),
                m: 0,
                j: 0,
            };
            (shape.fn_deriv)(&mut shape.deriv, &args.at_ksi);
            for m in 0..shape.npoint {
                args.m = m;
                for j in 0..shape.geo_ndim {
                    args.j = j;
                    // Lᵐⱼ := dNᵐ/dξⱼ
                    assert_deriv_approx_eq!(shape.deriv[m][j], args.at_ksi[j], calc_interp_m_var_ksi, args, tol);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn calc_coords_works() -> Result<(), StrError> {
        let mut data = gen_test_data();
        for shape in &mut data.shapes {
            let tol = *data.tolerances_calc_coords.get(&shape.kind).unwrap();
            let tol_high = *data.tolerances_calc_coords_high.get(&shape.kind).unwrap();
            println!("{:?}: tol={:e}, tol_high={:e}", shape.kind, tol, tol_high);
            set_coords_matrix(shape);
            let mut ksi = Vector::new(shape.geo_ndim);
            let mut x = Vector::new(shape.space_ndim);
            let mut x_correct = Vector::new(shape.space_ndim);
            for m in 0..shape.npoint {
                shape.get_reference_coords(&mut ksi, m);
                shape.calc_coords(&mut x, &ksi)?;
                gen_coords(&mut x_correct, &ksi);
                assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol);
            }
            // for ksi=(0,0,0), use higher tolerance because the cylinder mapping is inaccurate
            let ksi_0 = Vector::new(shape.geo_ndim);
            shape.calc_coords(&mut x, &ksi_0)?;
            gen_coords(&mut x_correct, &ksi_0);
            assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol_high);
        }
        Ok(())
    }

    #[test]
    fn calc_jacobian_works() -> Result<(), StrError> {
        let mut data = gen_test_data();
        for shape in &mut data.shapes {
            let tol = *data.tolerances_calc_jacobian.get(&shape.kind).unwrap();
            println!("{:?}: tol={:e}", shape.kind, tol);
            set_coords_matrix(shape);
            let mut at_ksi = Vector::new(shape.geo_ndim);
            for i in 0..shape.geo_ndim {
                at_ksi[i] = 0.25;
            }
            let args = &mut ArgsCoordsVarKsi {
                shape: Shape::new(shape.space_ndim, shape.geo_ndim, shape.npoint)?,
                at_ksi,
                ksi: Vector::new(shape.geo_ndim),
                x: Vector::new(shape.space_ndim),
                i: 0,
                j: 0,
            };
            set_coords_matrix(&mut args.shape);
            shape.calc_jacobian(&args.at_ksi)?;
            for i in 0..shape.space_ndim {
                args.i = i;
                for j in 0..shape.geo_ndim {
                    args.j = j;
                    // Jᵢⱼ := dxᵢ/dξⱼ
                    assert_deriv_approx_eq!(shape.jacobian[i][j], args.at_ksi[j], calc_coord_i_var_ksi, args, tol);
                }
            }
        }
        Ok(())
    }
}
