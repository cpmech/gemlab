use super::*;
use crate::StrError;
use russell_lab::{inverse, mat_mat_mul, mat_vec_mul, Matrix, NormVec, Vector};

/// Defines the class of geometric shape
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum GeoClass {
    /// Lines (segments) class
    Lin,

    /// Triangles class
    Tri,

    /// Quadrilaterals class
    Qua,

    /// Tetrahedra class
    Tet,

    /// Hexahedra class
    Hex,
}

/// Defines the kind of geometric shape
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum GeoKind {
    /// Line (segment) with 2 points (linear functions)
    Lin2,

    /// Line (segment) with 3 points (quadratic functions)
    Lin3,

    /// Line (segment) with 4 points (cubic functions)
    Lin4,

    /// Line (segment) with 5 points (quartic functions)
    Lin5,

    /// Triangle with 3 points (linear edges)
    Tri3,

    /// Triangle with 6 points (quadratic edges)
    Tri6,

    /// Triangle with 10 points (cubic edges; interior points)
    Tri10,

    /// Triangle with 15 points (quartic edges; interior points)
    Tri15,

    /// Quadrilateral with 4 points (linear edges)
    Qua4,

    /// Quadrilateral with 8 points (quadratic edges)
    Qua8,

    /// Quadrilateral with 9 points (quadratic edges; interior point)
    Qua9,

    /// Quadrilateral with 12 points (cubic edges)
    Qua12,

    /// Quadrilateral with 16 points (cubic edges; interior points)
    Qua16,

    /// Quadrilateral with 17 points (quartic edges; interior point)
    Qua17,

    /// Tetrahedron with 4 points (linear faces)
    Tet4,

    /// Tetrahedron with 10 points (quadratic faces)
    Tet10,

    /// Hexahedron with 8 points (bilinear faces)
    Hex8,

    /// Hexahedron with 20 points (quadratic faces)
    Hex20,
}

/// Defines an alias for interpolation functions
type FnInterp = fn(&mut Vector, &[f64]);

/// Defines an alias for derivative of interpolation functions
type FnDeriv = fn(&mut Matrix, &[f64]);

// Defines an alias for integration points data (coordinates and weights)
pub type IpData = &'static [[f64; 4]];

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
/// →  →    dNᵐ(ξ)
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
/// # Line in multi-dimensions (geo_ndim == 1 and space_ndim > 1)
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
/// # Boundary line in 2D (geo_ndim == 1 and space_ndim == 2)
///
/// If the line defines a boundary in 2D, we compute a normal vector by means of
///
/// ```text
/// →   →    →
/// n = e3 × g1
/// ```
///
/// # Boundary surface (geo_ndim == 2 and space_ndim == 3)
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
/// # Note
///
/// All public members are **readonly** and should not be modified externally.
pub struct Shape {
    /// Geometry class
    pub class: GeoClass,

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
    fn_interp: FnInterp,

    /// Function to calculate local derivatives (w.r.t. ksi) of interpolation functions
    fn_deriv: FnDeriv,

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

    /// Matrix Xᵀ: (space_ndim,npoint) transposed coordinates matrix (real space)
    pub coords_transp: Matrix,

    /// Minimum (space_ndim) coordinates from the X matrix
    pub min_coords: Vec<f64>,

    /// Maximum (space_ndim) coordinates from the X matrix
    pub max_coords: Vec<f64>,

    /// Integration points data (coordinates and weights)
    pub ip_data: IpData,
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
    /// This can be accomplished by calling the `set_point` method.
    pub fn new(space_ndim: usize, geo_ndim: usize, npoint: usize) -> Result<Self, StrError> {
        // collect geometry data
        let (class, kind, nedge, nface, edge_npoint, face_npoint, face_nedge, fn_interp, fn_deriv, ip_data): (
            GeoClass,
            GeoKind,
            usize,
            usize,
            usize,
            usize,
            usize,
            FnInterp,
            FnDeriv,
            IpData,
        ) = match (geo_ndim, npoint) {
            // Lin
            (1, 2) => (
                GeoClass::Lin,
                GeoKind::Lin2,
                Lin2::NEDGE,
                Lin2::NFACE,
                Lin2::EDGE_NPOINT,
                Lin2::FACE_NPOINT,
                Lin2::FACE_NEDGE,
                Lin2::calc_interp,
                Lin2::calc_deriv,
                &IP_LIN_LEGENDRE_2,
            ),
            (1, 3) => (
                GeoClass::Lin,
                GeoKind::Lin3,
                Lin3::NEDGE,
                Lin3::NFACE,
                Lin3::EDGE_NPOINT,
                Lin3::FACE_NPOINT,
                Lin3::FACE_NEDGE,
                Lin3::calc_interp,
                Lin3::calc_deriv,
                &IP_LIN_LEGENDRE_3,
            ),
            (1, 4) => (
                GeoClass::Lin,
                GeoKind::Lin4,
                Lin4::NEDGE,
                Lin4::NFACE,
                Lin4::EDGE_NPOINT,
                Lin4::FACE_NPOINT,
                Lin4::FACE_NEDGE,
                Lin4::calc_interp,
                Lin4::calc_deriv,
                &IP_LIN_LEGENDRE_4,
            ),
            (1, 5) => (
                GeoClass::Lin,
                GeoKind::Lin5,
                Lin5::NEDGE,
                Lin5::NFACE,
                Lin5::EDGE_NPOINT,
                Lin5::FACE_NPOINT,
                Lin5::FACE_NEDGE,
                Lin5::calc_interp,
                Lin5::calc_deriv,
                &IP_LIN_LEGENDRE_5,
            ),

            // Tri
            (2, 3) => (
                GeoClass::Tri,
                GeoKind::Tri3,
                Tri3::NEDGE,
                Tri3::NFACE,
                Tri3::EDGE_NPOINT,
                Tri3::FACE_NPOINT,
                Tri3::FACE_NEDGE,
                Tri3::calc_interp,
                Tri3::calc_deriv,
                &IP_TRI_INTERNAL_1,
            ),
            (2, 6) => (
                GeoClass::Tri,
                GeoKind::Tri6,
                Tri6::NEDGE,
                Tri6::NFACE,
                Tri6::EDGE_NPOINT,
                Tri6::FACE_NPOINT,
                Tri6::FACE_NEDGE,
                Tri6::calc_interp,
                Tri6::calc_deriv,
                &IP_TRI_INTERNAL_3,
            ),
            (2, 10) => (
                GeoClass::Tri,
                GeoKind::Tri10,
                Tri10::NEDGE,
                Tri10::NFACE,
                Tri10::EDGE_NPOINT,
                Tri10::FACE_NPOINT,
                Tri10::FACE_NEDGE,
                Tri10::calc_interp,
                Tri10::calc_deriv,
                &IP_TRI_INTERNAL_12,
            ),
            (2, 15) => (
                GeoClass::Tri,
                GeoKind::Tri15,
                Tri15::NEDGE,
                Tri15::NFACE,
                Tri15::EDGE_NPOINT,
                Tri15::FACE_NPOINT,
                Tri15::FACE_NEDGE,
                Tri15::calc_interp,
                Tri15::calc_deriv,
                &IP_TRI_INTERNAL_12,
            ),

            // Qua
            (2, 4) => (
                GeoClass::Qua,
                GeoKind::Qua4,
                Qua4::NEDGE,
                Qua4::NFACE,
                Qua4::EDGE_NPOINT,
                Qua4::FACE_NPOINT,
                Qua4::FACE_NEDGE,
                Qua4::calc_interp,
                Qua4::calc_deriv,
                &IP_QUA_LEGENDRE_4,
            ),
            (2, 8) => (
                GeoClass::Qua,
                GeoKind::Qua8,
                Qua8::NEDGE,
                Qua8::NFACE,
                Qua8::EDGE_NPOINT,
                Qua8::FACE_NPOINT,
                Qua8::FACE_NEDGE,
                Qua8::calc_interp,
                Qua8::calc_deriv,
                &IP_QUA_LEGENDRE_9,
            ),
            (2, 9) => (
                GeoClass::Qua,
                GeoKind::Qua9,
                Qua9::NEDGE,
                Qua9::NFACE,
                Qua9::EDGE_NPOINT,
                Qua9::FACE_NPOINT,
                Qua9::FACE_NEDGE,
                Qua9::calc_interp,
                Qua9::calc_deriv,
                &IP_QUA_LEGENDRE_9,
            ),
            (2, 12) => (
                GeoClass::Qua,
                GeoKind::Qua12,
                Qua12::NEDGE,
                Qua12::NFACE,
                Qua12::EDGE_NPOINT,
                Qua12::FACE_NPOINT,
                Qua12::FACE_NEDGE,
                Qua12::calc_interp,
                Qua12::calc_deriv,
                &IP_QUA_LEGENDRE_9,
            ),
            (2, 16) => (
                GeoClass::Qua,
                GeoKind::Qua16,
                Qua16::NEDGE,
                Qua16::NFACE,
                Qua16::EDGE_NPOINT,
                Qua16::FACE_NPOINT,
                Qua16::FACE_NEDGE,
                Qua16::calc_interp,
                Qua16::calc_deriv,
                &IP_QUA_LEGENDRE_16,
            ),
            (2, 17) => (
                GeoClass::Qua,
                GeoKind::Qua17,
                Qua17::NEDGE,
                Qua17::NFACE,
                Qua17::EDGE_NPOINT,
                Qua17::FACE_NPOINT,
                Qua17::FACE_NEDGE,
                Qua17::calc_interp,
                Qua17::calc_deriv,
                &IP_QUA_LEGENDRE_16,
            ),

            // Tet
            (3, 4) => (
                GeoClass::Tet,
                GeoKind::Tet4,
                Tet4::NEDGE,
                Tet4::NFACE,
                Tet4::EDGE_NPOINT,
                Tet4::FACE_NPOINT,
                Tet4::FACE_NEDGE,
                Tet4::calc_interp,
                Tet4::calc_deriv,
                &IP_TET_INTERNAL_1,
            ),
            (3, 10) => (
                GeoClass::Tet,
                GeoKind::Tet10,
                Tet10::NEDGE,
                Tet10::NFACE,
                Tet10::EDGE_NPOINT,
                Tet10::FACE_NPOINT,
                Tet10::FACE_NEDGE,
                Tet10::calc_interp,
                Tet10::calc_deriv,
                &IP_TET_INTERNAL_4,
            ),

            // Hex
            (3, 8) => (
                GeoClass::Hex,
                GeoKind::Hex8,
                Hex8::NEDGE,
                Hex8::NFACE,
                Hex8::EDGE_NPOINT,
                Hex8::FACE_NPOINT,
                Hex8::FACE_NEDGE,
                Hex8::calc_interp,
                Hex8::calc_deriv,
                &IP_HEX_LEGENDRE_8,
            ),
            (3, 20) => (
                GeoClass::Hex,
                GeoKind::Hex20,
                Hex20::NEDGE,
                Hex20::NFACE,
                Hex20::EDGE_NPOINT,
                Hex20::FACE_NPOINT,
                Hex20::FACE_NEDGE,
                Hex20::calc_interp,
                Hex20::calc_deriv,
                &IP_HEX_LEGENDRE_27,
            ),
            _ => return Err("(geo_ndim,npoint) combination is invalid"),
        };

        // return new Shape
        Ok(Shape {
            class,
            kind,
            space_ndim,
            geo_ndim,
            npoint,
            nedge,
            nface,
            edge_npoint,
            face_npoint,
            face_nedge,
            fn_interp,
            fn_deriv,
            interp: Vector::new(npoint),
            deriv: Matrix::new(npoint, geo_ndim),
            jacobian: Matrix::new(space_ndim, geo_ndim),
            inv_jacobian: if geo_ndim == space_ndim {
                Matrix::new(space_ndim, space_ndim)
            } else {
                Matrix::new(0, 0)
            },
            gradient: if geo_ndim == space_ndim {
                Matrix::new(npoint, space_ndim)
            } else {
                Matrix::new(0, 0)
            },
            coords_transp: Matrix::new(space_ndim, npoint),
            min_coords: vec![f64::MAX; space_ndim],
            max_coords: vec![f64::MIN; space_ndim],
            ip_data,
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
    /// * `interp` -- interpolation functions (npoint)
    #[inline]
    pub fn calc_interp(&mut self, ksi: &[f64]) {
        (self.fn_interp)(&mut self.interp, ksi);
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
    /// * `deriv` -- interpolation functions (npoint)
    #[inline]
    pub fn calc_deriv(&mut self, ksi: &[f64]) {
        (self.fn_deriv)(&mut self.deriv, ksi);
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
    pub fn set_point(&mut self, m: usize, j: usize, value: f64) -> Result<(), StrError> {
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
    pub fn calc_coords(&mut self, x: &mut Vector, ksi: &[f64]) -> Result<(), StrError> {
        if x.dim() != self.space_ndim {
            return Err("x.dim() must equal space_ndim");
        }
        if ksi.len() != self.geo_ndim {
            return Err("ksi.len() must equal geo_ndim");
        }
        self.calc_interp(ksi);
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
    pub fn calc_jacobian(&mut self, ksi: &[f64]) -> Result<f64, StrError> {
        if ksi.len() != self.geo_ndim {
            return Err("ksi.len() must equal geo_ndim");
        }
        self.calc_deriv(ksi);
        mat_mat_mul(&mut self.jacobian, 1.0, &self.coords_transp, &self.deriv)?;
        if self.geo_ndim == self.space_ndim {
            inverse(&mut self.inv_jacobian, &self.jacobian)
        } else {
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
    pub fn calc_boundary_normal(&mut self, normal: &mut Vector, ksi: &[f64]) -> Result<(), StrError> {
        // check
        if self.geo_ndim == self.space_ndim {
            return Err("geo_ndim must be smaller than space_ndim");
        }
        if normal.dim() != self.space_ndim {
            return Err("normal.dim() must equal space_ndim");
        }
        if ksi.len() != self.geo_ndim {
            return Err("ksi.len() must equal geo_ndim");
        }

        // compute Jacobian
        self.calc_deriv(ksi);
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
    /// You must set the coordinates matrix first, otherwise the computations
    /// will generate wrong results.
    pub fn approximate_ksi(
        &mut self,
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
    /// →  →    dNᵐ(ξ)
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
    pub fn calc_gradient(&mut self, ksi: &[f64]) -> Result<f64, StrError> {
        if self.geo_ndim != self.space_ndim {
            return Err("geo_ndim must equal space_ndim");
        }
        if ksi.len() != self.geo_ndim {
            return Err("ksi.len() must equal geo_ndim");
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
    pub fn get_reference_coords(&self, m: usize) -> &'static [f64] {
        match self.kind {
            GeoKind::Lin2 => &Lin2::POINT_REFERENCE_COORDS[m],
            GeoKind::Lin3 => &Lin3::POINT_REFERENCE_COORDS[m],
            GeoKind::Lin4 => &Lin4::POINT_REFERENCE_COORDS[m],
            GeoKind::Lin5 => &Lin5::POINT_REFERENCE_COORDS[m],
            GeoKind::Tri3 => &Tri3::POINT_REFERENCE_COORDS[m],
            GeoKind::Tri6 => &Tri6::POINT_REFERENCE_COORDS[m],
            GeoKind::Tri10 => &Tri10::POINT_REFERENCE_COORDS[m],
            GeoKind::Tri15 => &Tri15::POINT_REFERENCE_COORDS[m],
            GeoKind::Qua4 => &Qua4::POINT_REFERENCE_COORDS[m],
            GeoKind::Qua8 => &Qua8::POINT_REFERENCE_COORDS[m],
            GeoKind::Qua9 => &Qua9::POINT_REFERENCE_COORDS[m],
            GeoKind::Qua12 => &Qua12::POINT_REFERENCE_COORDS[m],
            GeoKind::Qua16 => &Qua16::POINT_REFERENCE_COORDS[m],
            GeoKind::Qua17 => &Qua17::POINT_REFERENCE_COORDS[m],
            GeoKind::Tet4 => &Tet4::POINT_REFERENCE_COORDS[m],
            GeoKind::Tet10 => &Tet10::POINT_REFERENCE_COORDS[m],
            GeoKind::Hex8 => &Hex8::POINT_REFERENCE_COORDS[m],
            GeoKind::Hex20 => &Hex20::POINT_REFERENCE_COORDS[m],
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
    use super::*;
    use crate::util::PI;
    use russell_chk::*;
    use russell_lab::copy_vector;
    use std::collections::HashMap;

    const RMIN: f64 = 1.0;
    const RMAX: f64 = 10.0;
    const AMIN: f64 = 30.0 * PI / 180.0;
    const AMAX: f64 = 60.0 * PI / 180.0;

    /// Generate coordinates
    ///
    /// He shape is the area indicated with "?" or the edge with "%".
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
        for m in 0..shape.npoint {
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
                shape.set_point(m, j, x[j]).unwrap();
            }
        }
    }

    #[test]
    fn calc_interp_works() -> Result<(), StrError> {
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
            for m in 0..shape.npoint {
                // get ξᵐ corresponding to point m
                let ksi = shape.get_reference_coords(m);

                // compute interpolation function Nⁿ(ξᵐ)
                shape.calc_interp(ksi);

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
        shape: Shape,     // auxiliary (copy) shape
        at_ksi: Vec<f64>, // at reference coord value
        ksi: Vec<f64>,    // temporary reference coord
        m: usize,         // point index from 0 to npoint
        j: usize,         // dimension index from 0 to geom_ndim
    }

    // Computes Nᵐ(ξ) with variable v := ξⱼ
    fn aux_deriv(v: f64, args: &mut ArgsNumL) -> f64 {
        args.ksi.copy_from_slice(&args.at_ksi);
        args.ksi[args.j] = v;
        args.shape.calc_interp(&args.ksi);
        args.shape.interp[args.m]
    }

    #[test]
    fn calc_deriv_works() -> Result<(), StrError> {
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
            let at_ksi = vec![0.25; shape.geo_ndim];

            // compute all derivatives of interpolation functions w.r.t ξ
            shape.calc_deriv(&at_ksi);

            // set arguments for numerical integration
            let args = &mut ArgsNumL {
                shape: Shape::new(shape.space_ndim, shape.geo_ndim, shape.npoint)?,
                at_ksi,
                ksi: vec![0.0; shape.geo_ndim],
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

        // define tolerances for point in the reference domain
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
        for (geo_ndim, npoint) in pairs {
            // allocate shape
            let space_ndim = geo_ndim;
            let shape = &mut Shape::new(space_ndim, geo_ndim, npoint)?;

            // set tolerance
            let tol = *tols.get(&shape.kind).unwrap();
            let tol_in = *tols_in.get(&shape.kind).unwrap();
            println!("{:?}: tol={:e}, tol_in={:e}", shape.kind, tol, tol_in);

            // set coordinates matrix
            set_coords_matrix(shape);

            // loop over points of shape
            let mut x = Vector::new(shape.space_ndim);
            let mut x_correct = Vector::new(shape.space_ndim);
            for m in 0..shape.npoint {
                // get ξᵐ corresponding to point m
                let ksi = shape.get_reference_coords(m);

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                shape.calc_coords(&mut x, ksi)?;

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
            shape.calc_coords(&mut x, &ksi_in)?;
            gen_coords(&mut x_correct, &ksi_in, shape.class);
            assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol_in);
        }
        Ok(())
    }

    // Holds arguments for numerical differentiation of x w.r.t ξ => Jacobian
    struct ArgsNumJ {
        shape: Shape,     // auxiliary (copy) shape
        at_ksi: Vec<f64>, // at reference coord value
        ksi: Vec<f64>,    // temporary reference coord
        x: Vector,        // (space_ndim) coordinates at ξ
        i: usize,         // dimension index from 0 to space_ndim
        j: usize,         // dimension index from 0 to geo_ndim
    }

    // Computes xᵢ(ξ) with variable v := ξⱼ
    fn aux_jacobian(v: f64, args: &mut ArgsNumJ) -> f64 {
        args.ksi.copy_from_slice(&args.at_ksi);
        args.ksi[args.j] = v;
        args.shape.calc_coords(&mut args.x, &args.ksi).unwrap();
        args.x[args.i]
    }

    #[test]
    fn calc_jacobian_works() -> Result<(), StrError> {
        // define dims and number of points
        let pairs = vec![(2, 6), (2, 4), (2, 8), (2, 17), (3, 4), (3, 8), (3, 20)];

        // define tolerances
        let mut tols = HashMap::new();
        tols.insert(GeoKind::Tri6, 1e-11);
        tols.insert(GeoKind::Qua4, 1e-11);
        tols.insert(GeoKind::Qua8, 1e-11);
        tols.insert(GeoKind::Qua17, 1e-10);
        tols.insert(GeoKind::Hex8, 1e-11);
        tols.insert(GeoKind::Tet4, 1e-12);
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
            let at_ksi = vec![0.25; shape.geo_ndim];

            // compute Jacobian, its inverse, and determinant
            let det_jac = shape.calc_jacobian(&at_ksi)?;
            assert!(det_jac > 0.0);

            // set arguments for numerical integration
            let args = &mut ArgsNumJ {
                shape: Shape::new(shape.space_ndim, shape.geo_ndim, shape.npoint)?,
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
        let at_ksi = vec![0.0; surf.geo_ndim];
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
        let mut edge = Shape::new(2, 1, 5)?;

        // set coordinates matrix
        set_coords_matrix(&mut edge);

        // compute boundary normal vector
        let at_ksi = vec![0.0; edge.geo_ndim];
        let mut normal = Vector::new(edge.space_ndim);
        edge.calc_boundary_normal(&mut normal, &at_ksi)?;

        // check magnitude of normal vector
        let mag_normal = normal.norm(NormVec::Euc);
        let length = RMAX - RMIN;
        let ref_length = 2.0;
        let length_ratio = length / ref_length;
        println!("normal =\n{}", normal);
        println!("mag_normal = {}", mag_normal);
        println!("length_ratio = {}", length_ratio);
        assert_approx_eq!(mag_normal, length_ratio, 1e-15);

        // check direction of normal vector
        let mut unit_normal = Vector::from(normal.as_data());
        unit_normal.scale(1.0 / mag_normal);
        assert_vec_approx_eq!(unit_normal.as_data(), &[-f64::sin(AMAX), f64::cos(AMAX)], 1e-15);
        Ok(())
    }

    #[test]
    fn approximate_ksi_works() -> Result<(), StrError> {
        // define dims and number of points
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
            let mut ksi = vec![0.0; shape.geo_ndim];
            for m in 0..shape.npoint {
                // get ξᵐ corresponding to point m
                let ksi_ref = shape.get_reference_coords(m);

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                shape.calc_coords(&mut x, ksi_ref)?;

                // compute approximation of the inverse mapping ξᵐ(xᵐ)
                let nit = shape.approximate_ksi(&mut ksi, &x, 10, 1e-14)?;

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
            shape.calc_coords(&mut x, &ksi_in)?;
            shape.approximate_ksi(&mut ksi, &x, 10, 1e-14)?;
            assert_vec_approx_eq!(ksi, ksi_in, tol);
        }
        Ok(())
    }

    // Holds arguments for numerical differentiation of N w.r.t x => G (gradient) matrix
    struct ArgsNumG {
        shape: Shape,  // auxiliary (copy) shape
        at_x: Vector,  // at x coord value
        x: Vector,     // temporary x coord
        ksi: Vec<f64>, // temporary reference coord
        m: usize,      // point index from 0 to npoint
        j: usize,      // dimension index from 0 to space_ndim
    }

    // Computes Nᵐ(ξ(x)) with variable v := xⱼ
    fn aux_grad(v: f64, args: &mut ArgsNumG) -> f64 {
        copy_vector(&mut args.x, &args.at_x).unwrap();
        args.x[args.j] = v;
        args.shape.approximate_ksi(&mut args.ksi, &args.x, 10, 1e-14).unwrap();
        args.shape.calc_interp(&args.ksi);
        args.shape.interp[args.m]
    }

    #[test]
    fn fn_gradient_works() -> Result<(), StrError> {
        // define dims and number of points
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
            let at_ksi = vec![0.25; shape.geo_ndim];

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
                ksi: vec![0.0; shape.geo_ndim],
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
