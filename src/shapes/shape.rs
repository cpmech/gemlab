use super::{
    Hex20, Hex8, Lin2, Lin3, Lin4, Lin5, Qua12, Qua16, Qua17, Qua4, Qua8, Qua9, Tet10, Tet4, Tri10, Tri15, Tri3, Tri6,
};
use crate::StrError;
use russell_lab::{inverse, mat_mat_mul, mat_vec_mul, Matrix, NormVec, Vector};

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

    /// Minimum (space_ndim) coordinates from the X matrix
    pub min_coords: Vec<f64>,

    /// Maximum (space_ndim) coordinates from the X matrix
    pub max_coords: Vec<f64>,

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
        let min_coords = vec![f64::MAX; space_ndim];
        let max_coords = vec![f64::MIN; space_ndim];
        let coords_transp = Matrix::new(space_ndim, npoint);
        match (geo_ndim, npoint) {
            (1, 2) => Ok(Shape {
                kind: GeoKind::Lin2,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Lin2::NEDGE,
                nface: Lin2::NFACE,
                edge_npoint: Lin2::EDGE_NPOINT,
                face_npoint: Lin2::FACE_NPOINT,
                face_nedge: Lin2::FACE_NEDGE,
                fn_interp: Lin2::calc_interp,
                fn_deriv: Lin2::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
            (1, 3) => Ok(Shape {
                kind: GeoKind::Lin3,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Lin3::NEDGE,
                nface: Lin3::NFACE,
                edge_npoint: Lin3::EDGE_NPOINT,
                face_npoint: Lin3::FACE_NPOINT,
                face_nedge: Lin3::FACE_NEDGE,
                fn_interp: Lin3::calc_interp,
                fn_deriv: Lin3::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
            (1, 4) => Ok(Shape {
                kind: GeoKind::Lin4,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Lin4::NEDGE,
                nface: Lin4::NFACE,
                edge_npoint: Lin4::EDGE_NPOINT,
                face_npoint: Lin4::FACE_NPOINT,
                face_nedge: Lin4::FACE_NEDGE,
                fn_interp: Lin4::calc_interp,
                fn_deriv: Lin4::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
            (1, 5) => Ok(Shape {
                kind: GeoKind::Lin5,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Lin5::NEDGE,
                nface: Lin5::NFACE,
                edge_npoint: Lin5::EDGE_NPOINT,
                face_npoint: Lin5::FACE_NPOINT,
                face_nedge: Lin5::FACE_NEDGE,
                fn_interp: Lin5::calc_interp,
                fn_deriv: Lin5::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
            (2, 3) => Ok(Shape {
                kind: GeoKind::Tri3,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Tri3::NEDGE,
                nface: Tri3::NFACE,
                edge_npoint: Tri3::EDGE_NPOINT,
                face_npoint: Tri3::FACE_NPOINT,
                face_nedge: Tri3::FACE_NEDGE,
                fn_interp: Tri3::calc_interp,
                fn_deriv: Tri3::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
            (2, 6) => Ok(Shape {
                kind: GeoKind::Tri6,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Tri6::NEDGE,
                nface: Tri6::NFACE,
                edge_npoint: Tri6::EDGE_NPOINT,
                face_npoint: Tri6::FACE_NPOINT,
                face_nedge: Tri6::FACE_NEDGE,
                fn_interp: Tri6::calc_interp,
                fn_deriv: Tri6::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
            (2, 10) => Ok(Shape {
                kind: GeoKind::Tri10,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Tri10::NEDGE,
                nface: Tri10::NFACE,
                edge_npoint: Tri10::EDGE_NPOINT,
                face_npoint: Tri10::FACE_NPOINT,
                face_nedge: Tri10::FACE_NEDGE,
                fn_interp: Tri10::calc_interp,
                fn_deriv: Tri10::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
            (2, 15) => Ok(Shape {
                kind: GeoKind::Tri15,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Tri15::NEDGE,
                nface: Tri15::NFACE,
                edge_npoint: Tri15::EDGE_NPOINT,
                face_npoint: Tri15::FACE_NPOINT,
                face_nedge: Tri15::FACE_NEDGE,
                fn_interp: Tri15::calc_interp,
                fn_deriv: Tri15::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
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
                min_coords,
                max_coords,
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
                min_coords,
                max_coords,
                coords_transp,
            }),
            (2, 9) => Ok(Shape {
                kind: GeoKind::Qua9,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Qua9::NEDGE,
                nface: Qua9::NFACE,
                edge_npoint: Qua9::EDGE_NPOINT,
                face_npoint: Qua9::FACE_NPOINT,
                face_nedge: Qua9::FACE_NEDGE,
                fn_interp: Qua9::calc_interp,
                fn_deriv: Qua9::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
            (2, 12) => Ok(Shape {
                kind: GeoKind::Qua12,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Qua12::NEDGE,
                nface: Qua12::NFACE,
                edge_npoint: Qua12::EDGE_NPOINT,
                face_npoint: Qua12::FACE_NPOINT,
                face_nedge: Qua12::FACE_NEDGE,
                fn_interp: Qua12::calc_interp,
                fn_deriv: Qua12::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
            (2, 16) => Ok(Shape {
                kind: GeoKind::Qua16,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Qua16::NEDGE,
                nface: Qua16::NFACE,
                edge_npoint: Qua16::EDGE_NPOINT,
                face_npoint: Qua16::FACE_NPOINT,
                face_nedge: Qua16::FACE_NEDGE,
                fn_interp: Qua16::calc_interp,
                fn_deriv: Qua16::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
            (2, 17) => Ok(Shape {
                kind: GeoKind::Qua17,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Qua17::NEDGE,
                nface: Qua17::NFACE,
                edge_npoint: Qua17::EDGE_NPOINT,
                face_npoint: Qua17::FACE_NPOINT,
                face_nedge: Qua17::FACE_NEDGE,
                fn_interp: Qua17::calc_interp,
                fn_deriv: Qua17::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
            (3, 4) => Ok(Shape {
                kind: GeoKind::Tet4,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Tet4::NEDGE,
                nface: Tet4::NFACE,
                edge_npoint: Tet4::EDGE_NPOINT,
                face_npoint: Tet4::FACE_NPOINT,
                face_nedge: Tet4::FACE_NEDGE,
                fn_interp: Tet4::calc_interp,
                fn_deriv: Tet4::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
            (3, 10) => Ok(Shape {
                kind: GeoKind::Tet10,
                space_ndim,
                geo_ndim,
                npoint,
                nedge: Tet10::NEDGE,
                nface: Tet10::NFACE,
                edge_npoint: Tet10::EDGE_NPOINT,
                face_npoint: Tet10::FACE_NPOINT,
                face_nedge: Tet10::FACE_NEDGE,
                fn_interp: Tet10::calc_interp,
                fn_deriv: Tet10::calc_deriv,
                interp,
                deriv,
                jacobian,
                inv_jacobian,
                gradient,
                min_coords,
                max_coords,
                coords_transp,
            }),
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
                min_coords,
                max_coords,
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
                min_coords,
                max_coords,
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
    /// * `m` -- point index in form 0 to npoint - 1
    /// * `j` -- dimension index from 0 to space_ndim - 1
    /// * `value` -- the X(m,j) component
    pub fn set_point_coord(&mut self, m: usize, j: usize, value: f64) -> Result<(), StrError> {
        if m >= self.npoint {
            return Err("index of point is invalid");
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
    /// * `ksi` -- ξ reference coordinate (geo_ndim)
    ///
    /// # Output
    ///
    /// * `x` -- real coordinate (space_ndim)
    ///
    /// # Warning
    ///
    /// You must set the coordinates matrix first, otherwise the computations
    /// will generate wrong results.
    pub fn calc_coords(&mut self, x: &mut Vector, ksi: &Vector) -> Result<(), StrError> {
        if x.dim() != self.space_ndim {
            return Err("x.dim() must equal space_ndim");
        }
        if ksi.dim() != self.geo_ndim {
            return Err("ksi.dim() must equal geo_ndim");
        }
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
    /// * `ksi` -- ξ reference coordinate (geo_ndim)
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
        if ksi.dim() != self.geo_ndim {
            return Err("ksi.dim() must equal geo_ndim");
        }
        (self.fn_deriv)(&mut self.deriv, ksi);
        mat_mat_mul(&mut self.jacobian, 1.0, &self.coords_transp, &self.deriv)?;
        if self.geo_ndim == self.space_ndim {
            inverse(&mut self.inv_jacobian, &self.jacobian)
        } else {
            Ok(0.0)
        }
    }

    /// Computes the boundary normal vector
    ///
    /// **Note:** This function works with `geo_ndim < space_ndim` only. In particular:
    ///
    /// * `geo_ndim = 1` and `space_ndim = 2`: line in 2D, or
    /// * `geo_ndim = 2` and `space_ndim = 3`: surface in 3D.
    ///
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate at the boundary (geo_ndim)
    ///
    /// # Output
    ///
    /// * `normal` -- (space_ndim) the boundary normal vector; not necessarily unitary
    ///
    /// # Updated variables
    ///
    /// * `jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    pub fn calc_boundary_normal(&mut self, normal: &mut Vector, ksi: &Vector) -> Result<(), StrError> {
        // check
        if self.geo_ndim == self.space_ndim {
            return Err("geo_ndim must be smaller than space_ndim");
        }
        if normal.dim() != self.space_ndim {
            return Err("normal.dim() must equal space_ndim");
        }
        if ksi.dim() != self.geo_ndim {
            return Err("ksi.dim() must equal geo_ndim");
        }

        // compute Jacobian
        (self.fn_deriv)(&mut self.deriv, ksi);
        mat_mat_mul(&mut self.jacobian, 1.0, &self.coords_transp, &self.deriv)?;

        // line in 2D
        if self.geo_ndim == 1 && self.space_ndim == 2 {
            normal[0] = -self.jacobian[1][0];
            normal[1] = self.jacobian[0][0];
            return Ok(());
        }

        // surface in 3D
        if self.geo_ndim == 2 && self.space_ndim == 3 {
            normal[0] = self.jacobian[1][0] * self.jacobian[2][1] - self.jacobian[2][0] * self.jacobian[1][1];
            normal[1] = self.jacobian[2][0] * self.jacobian[0][1] - self.jacobian[0][0] * self.jacobian[2][1];
            normal[2] = self.jacobian[0][0] * self.jacobian[1][1] - self.jacobian[1][0] * self.jacobian[0][1];
            return Ok(());
        }

        Err("geo_ndim and space_ndim combination is invalid for boundary normal")
    }

    /// Approximates the reference coordinates from given real coordinates (inverse mapping)
    ///
    /// **Note:** This function works with `geo_ndim == space_ndim` only
    ///
    /// # Input
    ///
    /// * `x` -- real coordinates (space_ndim = geo_ndim)
    /// * `nit_max` -- maximum number of iterations (e.g., 10)
    /// * `tol` -- tolerance for the norm of the difference x - x(ξ) (e.g., 1e-8)
    ///
    /// # Output
    ///
    /// * `ksi` -- ξ reference coordinates (geo_ndim = space_ndim)
    /// * returns the number of iterations
    ///
    /// # Warning
    ///
    /// You must set the coordinates matrix first, otherwise the computations
    /// will generate wrong results.
    pub fn approximate_ksi(
        &mut self,
        ksi: &mut Vector,
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
        if ksi.dim() != self.geo_ndim {
            return Err("ksi.dim() must equal geo_ndim");
        }

        // use linear scale to guess ksi
        for j in 0..self.geo_ndim {
            // ksi[j] = 0.0; // trial at center
            ksi[j] = 2.0 * (x[j] - self.min_coords[j]) / (self.max_coords[j] - self.min_coords[j]) - 1.0;
            if ksi[j] < -1.0 {
                ksi[j] = -1.0;
            }
            if ksi[j] > 1.0 {
                ksi[j] = 1.0;
            }
        }

        // perform iterations
        let mut residual = Vector::new(self.space_ndim);
        let mut x_at_ksi = Vector::new(self.space_ndim);
        let mut delta_ksi = Vector::new(self.geo_ndim);
        for it in 0..nit_max {
            self.calc_coords(&mut x_at_ksi, &ksi)?;
            for i in 0..self.space_ndim {
                residual[i] = x[i] - x_at_ksi[i];
            }
            if residual.norm(NormVec::Euc) <= tol {
                return Ok(it);
            }
            self.calc_jacobian(ksi)?;
            mat_vec_mul(&mut delta_ksi, 1.0, &self.inv_jacobian, &residual)?;
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
    /// * `ksi` -- ξ reference coordinate (geo_ndim = space_ndim)
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
        if ksi.dim() != self.geo_ndim {
            return Err("ksi.dim() must equal geo_ndim");
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
            GeoKind::Lin2 => 0,
            GeoKind::Lin3 => 0,
            GeoKind::Lin4 => 0,
            GeoKind::Lin5 => 0,
            GeoKind::Tri3 => Tri3::EDGE_POINT_IDS[e][i],
            GeoKind::Tri6 => Tri6::EDGE_POINT_IDS[e][i],
            GeoKind::Tri10 => Tri10::EDGE_POINT_IDS[e][i],
            GeoKind::Tri15 => Tri15::EDGE_POINT_IDS[e][i],
            GeoKind::Qua4 => Qua4::EDGE_POINT_IDS[e][i],
            GeoKind::Qua8 => Qua8::EDGE_POINT_IDS[e][i],
            GeoKind::Qua9 => Qua9::EDGE_POINT_IDS[e][i],
            GeoKind::Qua12 => Qua12::EDGE_POINT_IDS[e][i],
            GeoKind::Qua16 => Qua16::EDGE_POINT_IDS[e][i],
            GeoKind::Qua17 => Qua17::EDGE_POINT_IDS[e][i],
            GeoKind::Tet4 => Tet4::EDGE_POINT_IDS[e][i],
            GeoKind::Tet10 => Tet10::EDGE_POINT_IDS[e][i],
            GeoKind::Hex8 => Hex8::EDGE_POINT_IDS[e][i],
            GeoKind::Hex20 => Hex20::EDGE_POINT_IDS[e][i],
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
            GeoKind::Tet4 => Tet4::FACE_POINT_IDS[f][i],
            GeoKind::Tet10 => Tet10::FACE_POINT_IDS[f][i],
            GeoKind::Hex8 => Hex8::FACE_POINT_IDS[f][i],
            GeoKind::Hex20 => Hex20::FACE_POINT_IDS[f][i],
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
            GeoKind::Tet4 => Tet4::FACE_EDGE_POINT_IDS[f][k][i],
            GeoKind::Tet10 => Tet10::FACE_EDGE_POINT_IDS[f][k][i],
            GeoKind::Hex8 => Hex8::FACE_EDGE_POINT_IDS[f][k][i],
            GeoKind::Hex20 => Hex20::FACE_EDGE_POINT_IDS[f][k][i],
        }
    }

    /// Returns the reference coordinates at point m
    ///
    /// # Output
    ///
    /// * `ksi` -- (geo_ndim) reference coordinates `ξᵐ` at point m
    pub fn get_reference_coords(&self, ksi: &mut Vector, m: usize) {
        match self.kind {
            GeoKind::Lin2 => {
                ksi[0] = Lin2::POINT_REFERENCE_COORDS[m][0];
            }
            GeoKind::Lin3 => {
                ksi[0] = Lin3::POINT_REFERENCE_COORDS[m][0];
            }
            GeoKind::Lin4 => {
                ksi[0] = Lin4::POINT_REFERENCE_COORDS[m][0];
            }
            GeoKind::Lin5 => {
                ksi[0] = Lin5::POINT_REFERENCE_COORDS[m][0];
            }
            GeoKind::Tri3 => {
                ksi[0] = Tri3::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Tri3::POINT_REFERENCE_COORDS[m][1];
            }
            GeoKind::Tri6 => {
                ksi[0] = Tri6::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Tri6::POINT_REFERENCE_COORDS[m][1];
            }
            GeoKind::Tri10 => {
                ksi[0] = Tri10::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Tri10::POINT_REFERENCE_COORDS[m][1];
            }
            GeoKind::Tri15 => {
                ksi[0] = Tri15::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Tri15::POINT_REFERENCE_COORDS[m][1];
            }
            GeoKind::Qua4 => {
                ksi[0] = Qua4::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Qua4::POINT_REFERENCE_COORDS[m][1];
            }
            GeoKind::Qua8 => {
                ksi[0] = Qua8::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Qua8::POINT_REFERENCE_COORDS[m][1];
            }
            GeoKind::Qua9 => {
                ksi[0] = Qua9::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Qua9::POINT_REFERENCE_COORDS[m][1];
            }
            GeoKind::Qua12 => {
                ksi[0] = Qua12::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Qua12::POINT_REFERENCE_COORDS[m][1];
            }
            GeoKind::Qua16 => {
                ksi[0] = Qua16::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Qua16::POINT_REFERENCE_COORDS[m][1];
            }
            GeoKind::Qua17 => {
                ksi[0] = Qua17::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Qua17::POINT_REFERENCE_COORDS[m][1];
            }
            GeoKind::Tet4 => {
                ksi[0] = Tet4::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Tet4::POINT_REFERENCE_COORDS[m][1];
                ksi[2] = Tet4::POINT_REFERENCE_COORDS[m][2];
            }
            GeoKind::Tet10 => {
                ksi[0] = Tet10::POINT_REFERENCE_COORDS[m][0];
                ksi[1] = Tet10::POINT_REFERENCE_COORDS[m][1];
                ksi[2] = Tet10::POINT_REFERENCE_COORDS[m][2];
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
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::PI;
    use russell_chk::*;
    use russell_lab::copy_vector;
    use std::collections::HashMap;

    const RMIN: f64 = 1.0;
    const RMAX: f64 = 10.0;
    const AMIN: f64 = 30.0 * PI / 180.0;
    const AMAX: f64 = 60.0 * PI / 180.0;

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

    // Generates matrix of coordinates
    fn set_coords_matrix(shape: &mut Shape) {
        let mut x = Vector::new(shape.space_ndim);
        let mut ksi = Vector::new(shape.geo_ndim);
        let mut ksi_aux = Vector::new(shape.space_ndim);
        for m in 0..shape.npoint {
            shape.get_reference_coords(&mut ksi, m);
            if shape.geo_ndim == shape.space_ndim {
                gen_coords(&mut x, &ksi);
            } else if shape.geo_ndim == 1 && shape.space_ndim == 2 {
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = 1.0;
                gen_coords(&mut x, &ksi_aux);
            } else if shape.geo_ndim == 2 && shape.space_ndim == 3 {
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = ksi[1];
                ksi_aux[2] = 1.0;
                gen_coords(&mut x, &ksi_aux);
            } else {
                panic!("(geo_ndim,space_ndim) pair is invalid");
            }
            for j in 0..shape.space_ndim {
                shape.set_point_coord(m, j, x[j]).unwrap();
            }
        }
    }

    #[test]
    fn fn_interp_works() -> Result<(), StrError> {
        // define dims and number of points
        let pairs = vec![
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5),
            (2, 3),
            (2, 6),
            (2, 10),
            (2, 15),
            (2, 4),
            (2, 8),
            (2, 9),
            (2, 12),
            (2, 16),
            (2, 17),
            (3, 4),
            (3, 10),
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
        for (geo_ndim, npoint) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let shape = &mut Shape::new(space_ndim, geo_ndim, npoint)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();
            println!("{:?}: tol={:e}", shape.kind, tol);

            // loop over points of shape
            let mut ksi = Vector::new(shape.geo_ndim);
            for m in 0..shape.npoint {
                // get ξᵐ corresponding to point m
                shape.get_reference_coords(&mut ksi, m);

                // compute interpolation function Nⁿ(ξᵐ)
                (shape.fn_interp)(&mut shape.interp, &ksi);

                // check: Nⁿ(ξᵐ) = 1 if m==n; 0 otherwise
                for n in 0..shape.npoint {
                    if m == n {
                        assert_approx_eq!(shape.interp[n], 1.0, tol);
                    } else {
                        assert_approx_eq!(shape.interp[n], 0.0, tol);
                    }
                }
            }
        }
        Ok(())
    }

    // Holds arguments for numerical differentiation of N w.r.t ξ => L (deriv) matrix
    struct ArgsNumL {
        shape: Shape,   // auxiliary (copy) shape
        at_ksi: Vector, // at reference coord value
        ksi: Vector,    // temporary reference coord
        m: usize,       // point index from 0 to npoint
        j: usize,       // dimension index from 0 to geom_ndim
    }

    // Computes Nᵐ(ξ) with variable v := ξⱼ
    fn aux_deriv(v: f64, args: &mut ArgsNumL) -> f64 {
        copy_vector(&mut args.ksi, &args.at_ksi).unwrap();
        args.ksi[args.j] = v;
        (args.shape.fn_interp)(&mut args.shape.interp, &args.ksi);
        args.shape.interp[args.m]
    }

    #[test]
    fn fn_deriv_works() -> Result<(), StrError> {
        // define dims and number of points
        let pairs = vec![
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5),
            (2, 3),
            (2, 6),
            (2, 10),
            (2, 15),
            (2, 4),
            (2, 8),
            (2, 9),
            (2, 12),
            (2, 16),
            (2, 17),
            (3, 4),
            (3, 10),
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
        for (geo_ndim, npoint) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let shape = &mut Shape::new(space_ndim, geo_ndim, npoint)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();
            println!("{:?}: tol={:e}", shape.kind, tol);

            // set ξ within reference space
            let mut at_ksi = Vector::new(shape.geo_ndim);
            for j in 0..shape.geo_ndim {
                at_ksi[j] = 0.25;
            }

            // compute all derivatives of interpolation functions w.r.t ξ
            (shape.fn_deriv)(&mut shape.deriv, &at_ksi);

            // set arguments for numerical integration
            let args = &mut ArgsNumL {
                shape: Shape::new(shape.space_ndim, shape.geo_ndim, shape.npoint)?,
                at_ksi,
                ksi: Vector::new(shape.geo_ndim),
                m: 0,
                j: 0,
            };

            // check Lᵐ(ξ) = dNᵐ(ξ)/dξ
            for m in 0..shape.npoint {
                args.m = m;
                for j in 0..shape.geo_ndim {
                    args.j = j;
                    // Lᵐⱼ := dNᵐ/dξⱼ
                    assert_deriv_approx_eq!(shape.deriv[m][j], args.at_ksi[j], aux_deriv, args, tol);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn calc_coords_works() -> Result<(), StrError> {
        // define dims and number of points
        let pairs = vec![(2, 4), (2, 8), (2, 17), (3, 8), (3, 20)];

        // define tolerances
        let mut tols = HashMap::new();
        tols.insert(GeoKind::Qua4, 1e-15);
        tols.insert(GeoKind::Qua8, 1e-15);
        tols.insert(GeoKind::Qua17, 1e-15);
        tols.insert(GeoKind::Hex8, 1e-15);
        tols.insert(GeoKind::Hex20, 1e-15);

        // define tolerances for point at the middle of the reference domain
        let mut tols_mid = HashMap::new();
        tols_mid.insert(GeoKind::Qua4, 0.14); // linear maps are inaccurate for the circle wedge
        tols_mid.insert(GeoKind::Qua8, 1e-14);
        tols_mid.insert(GeoKind::Qua17, 1e-14);
        tols_mid.insert(GeoKind::Hex8, 0.14); // bi-linear maps are inaccurate for the circle wedge
        tols_mid.insert(GeoKind::Hex20, 1e-14);

        // loop over shapes
        for (geo_ndim, npoint) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let shape = &mut Shape::new(space_ndim, geo_ndim, npoint)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();
            let tol_mid = *tols_mid.get(&shape.kind).unwrap();
            println!("{:?}: tol={:e}, tol_mid={:e}", shape.kind, tol, tol_mid);

            // set coordinates matrix
            set_coords_matrix(shape);

            // loop over points of shape
            let mut ksi = Vector::new(shape.geo_ndim);
            let mut x = Vector::new(shape.space_ndim);
            let mut x_correct = Vector::new(shape.space_ndim);
            for m in 0..shape.npoint {
                // get ξᵐ corresponding to point m
                shape.get_reference_coords(&mut ksi, m);

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                shape.calc_coords(&mut x, &ksi)?;

                // compare xᵐ with generated coordinates
                gen_coords(&mut x_correct, &ksi);
                assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol);
            }

            // test again at middle of reference space with ξ := (0,0,0)
            let ksi_mid = Vector::new(shape.geo_ndim);
            shape.calc_coords(&mut x, &ksi_mid)?;
            gen_coords(&mut x_correct, &ksi_mid);
            assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol_mid);
        }
        Ok(())
    }

    // Holds arguments for numerical differentiation of x w.r.t ξ => Jacobian
    struct ArgsNumJ {
        shape: Shape,   // auxiliary (copy) shape
        at_ksi: Vector, // at reference coord value
        ksi: Vector,    // temporary reference coord
        x: Vector,      // (space_ndim) coordinates at ξ
        i: usize,       // dimension index from 0 to space_ndim
        j: usize,       // dimension index from 0 to geo_ndim
    }

    // Computes xᵢ(ξ) with variable v := ξⱼ
    fn aux_jacobian(v: f64, args: &mut ArgsNumJ) -> f64 {
        copy_vector(&mut args.ksi, &args.at_ksi).unwrap();
        args.ksi[args.j] = v;
        args.shape.calc_coords(&mut args.x, &args.ksi).unwrap();
        args.x[args.i]
    }

    #[test]
    fn calc_jacobian_works() -> Result<(), StrError> {
        // define dims and number of points
        let pairs = vec![(2, 4), (2, 8), (2, 17), (3, 8), (3, 20)];

        // define tolerances
        let mut tols = HashMap::new();
        tols.insert(GeoKind::Qua4, 1e-11);
        tols.insert(GeoKind::Qua8, 1e-11);
        tols.insert(GeoKind::Qua17, 1e-10);
        tols.insert(GeoKind::Hex8, 1e-11);
        tols.insert(GeoKind::Hex20, 1e-11);

        // loop over shapes
        for (geo_ndim, npoint) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let shape = &mut Shape::new(space_ndim, geo_ndim, npoint)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();
            println!("{:?}: tol={:e}", shape.kind, tol);

            // set coordinates matrix
            set_coords_matrix(shape);

            // set ξ within reference space
            let mut at_ksi = Vector::new(shape.geo_ndim);
            for j in 0..shape.geo_ndim {
                at_ksi[j] = 0.25;
            }

            // compute Jacobian, its inverse, and determinant
            let det_jac = shape.calc_jacobian(&at_ksi)?;
            assert!(det_jac > 0.0);

            // set arguments for numerical integration
            let args = &mut ArgsNumJ {
                shape: Shape::new(shape.space_ndim, shape.geo_ndim, shape.npoint)?,
                at_ksi,
                ksi: Vector::new(shape.geo_ndim),
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
                    assert_deriv_approx_eq!(shape.jacobian[i][j], args.at_ksi[j], aux_jacobian, args, tol);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn calc_boundary_normal_works() -> Result<(), StrError> {
        // allocate surface shape
        let mut surf = Shape::new(3, 2, 17)?;

        // set coordinates matrix
        set_coords_matrix(&mut surf);

        // compute boundary normal vector
        let at_ksi = Vector::new(surf.geo_ndim);
        let mut normal = Vector::new(surf.space_ndim);
        surf.calc_boundary_normal(&mut normal, &at_ksi)?;

        // check magnitude of normal vector
        let mag_normal = normal.norm(NormVec::Euc);
        let area = PI * (RMAX * RMAX - RMIN * RMIN) / 12.0;
        let ref_area = 4.0;
        let area_ratio = area / ref_area;
        println!("normal =\n{}", normal);
        println!("mag_normal = {}", mag_normal);
        println!("area_ratio = {}", area_ratio);
        assert_approx_eq!(mag_normal, area_ratio, 1e-4);

        // check direction of normal vector
        let mut unit_normal = Vector::from(normal.as_data());
        unit_normal.scale(1.0 / mag_normal);
        assert_vec_approx_eq!(unit_normal.as_data(), &[0.0, 0.0, 1.0], 1e-15);
        Ok(())
    }

    #[test]
    fn approximate_ksi_works() -> Result<(), StrError> {
        // define dims and number of points
        let pairs = vec![(2, 4), (2, 8), (3, 8), (3, 20)];

        // define tolerances
        let mut tols = HashMap::new();
        tols.insert(GeoKind::Qua4, 1e-15);
        tols.insert(GeoKind::Qua8, 1e-15);
        tols.insert(GeoKind::Hex8, 1e-15);
        tols.insert(GeoKind::Hex20, 1e-15);

        // loop over shapes
        for (geo_ndim, npoint) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let shape = &mut Shape::new(space_ndim, geo_ndim, npoint)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();
            println!("{:?}: tol={:e}", shape.kind, tol);

            // set coordinates matrix
            set_coords_matrix(shape);

            // loop over points of shape
            let mut x = Vector::new(shape.space_ndim);
            let mut ksi = Vector::new(shape.geo_ndim);
            let mut ksi_correct = Vector::new(shape.geo_ndim);
            for m in 0..shape.npoint {
                // get ξᵐ corresponding to point m
                shape.get_reference_coords(&mut ksi, m);
                copy_vector(&mut ksi_correct, &ksi)?;

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                shape.calc_coords(&mut x, &ksi)?;

                // compute approximation of the inverse mapping ξᵐ(xᵐ)
                let nit = shape.approximate_ksi(&mut ksi, &x, 10, 1e-14)?;

                // check (linear and bi-linear shapes converge with nit = 1)
                if shape.kind == GeoKind::Qua4 || shape.kind == GeoKind::Hex8 {
                    assert_eq!(nit, 1);
                }
                assert_vec_approx_eq!(ksi.as_data(), ksi_correct.as_data(), tol);
            }

            // test again at middle of reference space with ξ := (0,0,0)
            let ksi_mid = Vector::new(shape.geo_ndim);
            shape.calc_coords(&mut x, &ksi_mid)?;
            shape.approximate_ksi(&mut ksi, &x, 10, 1e-14)?;
            assert_vec_approx_eq!(ksi.as_data(), ksi_mid.as_data(), tol);
        }
        Ok(())
    }

    // Holds arguments for numerical differentiation of N w.r.t x => G (gradient) matrix
    struct ArgsNumG {
        shape: Shape, // auxiliary (copy) shape
        at_x: Vector, // at x coord value
        x: Vector,    // temporary x coord
        ksi: Vector,  // temporary reference coord
        m: usize,     // point index from 0 to npoint
        j: usize,     // dimension index from 0 to space_ndim
    }

    // Computes Nᵐ(ξ(x)) with variable v := xⱼ
    fn aux_grad(v: f64, args: &mut ArgsNumG) -> f64 {
        copy_vector(&mut args.x, &args.at_x).unwrap();
        args.x[args.j] = v;
        args.shape.approximate_ksi(&mut args.ksi, &args.x, 10, 1e-14).unwrap();
        (args.shape.fn_interp)(&mut args.shape.interp, &args.ksi);
        args.shape.interp[args.m]
    }

    #[test]
    fn fn_gradient_works() -> Result<(), StrError> {
        // define dims and number of points
        let pairs = vec![(2, 4), (2, 8), (3, 8), (3, 20)];

        // define tolerances
        let mut tols = HashMap::new();
        tols.insert(GeoKind::Qua4, 1e-11);
        tols.insert(GeoKind::Qua8, 1e-11);
        tols.insert(GeoKind::Hex8, 1e-11);
        tols.insert(GeoKind::Hex20, 1e-10);

        // loop over shapes
        for (geo_ndim, npoint) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let shape = &mut Shape::new(space_ndim, geo_ndim, npoint)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();
            println!("{:?}: tol={:e}", shape.kind, tol);

            // set coordinates matrix
            set_coords_matrix(shape);

            // set ξ at the middle of the reference domain
            let mut at_ksi = Vector::new(shape.geo_ndim);
            for j in 0..shape.geo_ndim {
                at_ksi[j] = 0.0;
            }

            // compute x corresponding to ξ using the isoparametric formula
            let mut at_x = Vector::new(shape.space_ndim);
            shape.calc_coords(&mut at_x, &at_ksi)?;

            // compute gradient
            let det_jac = shape.calc_gradient(&at_ksi)?;
            assert!(det_jac > 0.0);

            // set arguments for numerical integration
            let args = &mut ArgsNumG {
                shape: Shape::new(shape.space_ndim, shape.geo_ndim, shape.npoint)?,
                at_x,
                x: Vector::new(shape.space_ndim),
                ksi: Vector::new(shape.geo_ndim),
                m: 0,
                j: 0,
            };
            set_coords_matrix(&mut args.shape);

            // check Gᵐ(ξ(x)) = dNᵐ(ξ(x))/dx
            for m in 0..shape.npoint {
                args.m = m;
                for j in 0..shape.geo_ndim {
                    args.j = j;
                    // Gᵐⱼ := dNᵐ/dxⱼ
                    assert_deriv_approx_eq!(shape.gradient[m][j], args.at_x[j], aux_grad, args, tol);
                }
            }
        }
        Ok(())
    }
}
