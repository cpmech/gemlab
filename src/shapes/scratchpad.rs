use super::{FnDeriv, FnInterp, GeoKind};
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Holds the interpolation functions, gradients, and associated structures
///
/// This struct holds the variables calculated by the interpolation and derivative functions,
/// including the Jacobian, inverse Jacobian, and gradients. The **field** interpolation is given by
/// `kind`, while the **cell geometry** mapping is given by `geo_kind`. Both are equal for
/// single-order shapes (see [Scratchpad::new]); they may differ in order for mixed-order shapes
/// (see [Scratchpad::new_mixed]).
///
/// # Sequence of computations for the gradient
///
/// The derivatives `L` in step 1 are those of the field (`kind`), while the coordinates and the
/// Jacobian in step 2 belong to the cell geometry (`geo_kind`). For single-order shapes, the two
/// kinds coincide.
///
/// ## 1 Derivatives of interpolation functions @ ξ
///
/// ```text
///             →
/// →  →    dNᵐ(ξ)
/// Lᵐ(ξ) = ——————
///            →
///           dξ
/// ```
///
/// ## 2 Jacobian @ ξ
///
/// ```text
///         →
///   →    dx     →    →
/// J(ξ) = —— = Σ xᵐ ⊗ L_geoᵐ
///         →   m
///        dξ
/// ```
///
/// ```text
///            J         =        Xᵀ        ·       L_geo
/// (space_ndim,geo_ndim) (space_ndim,geo_nnode) (geo_nnode,geo_ndim)
/// ```
///
/// ## 3 Inverse Jacobian
///
/// Compute the inverse Jacobian matrix
///
/// ## 4 Gradient @ ξ
///
/// ```text
///             →
/// →  →    dNᵐ(ξ)
/// Bᵐ(ξ) = ——————
///            →
///           dx
/// ```
///
/// ```text
///        B          =       L        ·           J⁻¹
/// (nnode,space_ndim) (nnode,geo_ndim) (space_ndim,space_ndim)
/// ```
///
/// # Warning
///
/// All members are **readonly** and should not be modified externally. Use the member
/// functions to calculate the data members of this struct.
#[derive(Clone)]
pub struct Scratchpad {
    /// The kind of the field interpolation
    pub kind: GeoKind,

    /// The kind of the cell geometry mapping
    ///
    /// This member equals `kind` for single-order shapes. For mixed-order shapes, it defines the
    /// geometry that maps the reference coordinates to the real coordinates, while `kind` defines
    /// the field interpolation (see [Scratchpad::new_mixed]).
    pub geo_kind: GeoKind,

    /// Array N: (nnode) Nᵐ field interpolation functions at reference coordinate ξ (ksi)
    pub interp: Vector,

    /// Matrix L: (nnode,geo_ndim) dNᵐ/dξ field derivatives with respect to reference coordinate ξ (ksi)
    pub deriv: Matrix,

    /// Array N_geo: (geo_nnode) geometry interpolation functions at reference coordinate ξ (ksi)
    pub interp_geo: Vector,

    /// Matrix L_geo: (geo_nnode,geo_ndim) geometry derivatives with respect to reference coordinate ξ (ksi)
    pub deriv_geo: Matrix,

    /// Matrix J: (space_ndim,geo_ndim) dx/dξ Jacobian matrix
    ///
    /// ```text
    /// jacobian := J = Xᵀ · L
    /// jacobian := Jcable = Xᵀ · L
    /// jacobian := Jshell = Xᵀ · L
    /// ```
    pub jacobian: Matrix,

    /// Matrix inv(J): (space_ndim,space_ndim) Inverse Jacobian matrix (only if geo_ndim = space_ndim) at ξ (ksi)
    ///
    /// Only available if `geo_ndim = space_ndim` (otherwise, the matrix is set to empty; 0 x 0 matrix)
    pub inv_jacobian: Matrix,

    /// Matrix B: (nnode,space_ndim) dNᵐ/dx Gradient of shape functions (only if geo_ndim = space_ndim) at ξ (ksi)
    ///
    /// Only available if `geo_ndim = space_ndim` (otherwise, the matrix is set to empty; 0 x 0 matrix)
    ///
    /// ```text
    /// B = L · J⁻¹
    /// ```
    pub gradient: Matrix,

    /// Matrix Xᵀ: (space_ndim,geo_nnode) transposed coordinates matrix (real space)
    ///
    /// ```text
    ///      ┌                              ┐  superscript = node
    ///      | x⁰₀  x¹₀  x²₀  x³₀       xᴹ₀ |  subscript = space dimension
    /// Xᵀ = | x⁰₁  x¹₁  x²₁  x³₁  ...  xᴹ₁ |
    ///      | x⁰₂  x¹₂  x²₂  x³₂       xᴹ₂ |
    ///      └                              ┘_(space_ndim,geo_nnode)
    /// where `M = geo_nnode - 1`
    /// ```
    ///
    /// You must call `set_xx` to set the components, otherwise calculations will be **incorrect.**
    pub xxt: Matrix,

    /// Indicates that all coordinates in the X matrix have been set and the min,max values computed
    pub ok_xxt: bool,

    /// Function to calculate the field interpolation functions
    pub fn_interp: FnInterp,

    /// Function to calculate derivatives of the field interpolation functions
    pub fn_deriv: FnDeriv,

    /// Function to calculate the geometry interpolation functions
    pub fn_interp_geo: FnInterp,

    /// Function to calculate derivatives of the geometry interpolation functions
    pub fn_deriv_geo: FnDeriv,
}

impl Scratchpad {
    /// Allocates a new instance with a single-order shape (`geo_kind = kind`)
    ///
    /// # Notes
    ///
    /// 1. `space_ndim` must be 2 or 3
    /// 2. `space_ndim` must be greater than or equal to `geo_ndim`;
    ///    e.g., you cannot have a 3D shape in a 2D space.
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Scratchpad};
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
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
    ///     let space_ndim = 3;
    ///     let pad = Scratchpad::new(space_ndim, GeoKind::Hex8)?;
    ///     assert_eq!(pad.kind, GeoKind::Hex8);
    ///     assert_eq!(pad.geo_kind, GeoKind::Hex8);
    ///     assert_eq!(pad.ok_xxt, false);
    ///     Ok(())
    /// }
    /// ```
    pub fn new(space_ndim: usize, kind: GeoKind) -> Result<Self, StrError> {
        if space_ndim < 2 || space_ndim > 3 {
            return Err("space_ndim must be 2 or 3");
        }
        let geo_ndim = kind.ndim();
        if geo_ndim > space_ndim {
            return Err("geo_ndim cannot be greater than space_ndim");
        }
        Ok(Self::alloc(space_ndim, kind, kind))
    }

    /// Allocates a new instance with a field interpolation over a cell geometry mapping
    ///
    /// This function supports **mixed-order** shapes: the field functions `kind` are evaluated over
    /// the geometry of `geo_kind`, which may be of lower, equal, or higher order than `kind`. For
    /// example, a `Qua4` field over a `Qua8` cell, or a `Qua8` field over a `Qua4` cell. The
    /// geometric Jacobian, inverse Jacobian, and integration measure always come from `geo_kind`,
    /// while `interp`, `deriv`, and `gradient` belong to `kind`.
    ///
    /// # Notes
    ///
    /// 1. `geo_kind` must match `space_ndim` (e.g., a `Qua4` or `Qua8` cell in 2D)
    /// 2. `kind` and `geo_kind` must be of the same class and geometric dimension
    /// 3. `space_ndim` must be 2 or 3
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Scratchpad};
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     // A Qua4 field over a Qua8 cell geometry
    ///     let pad = Scratchpad::new_mixed(2, GeoKind::Qua4, GeoKind::Qua8)?;
    ///     assert_eq!(pad.kind, GeoKind::Qua4);
    ///     assert_eq!(pad.geo_kind, GeoKind::Qua8);
    ///     assert_eq!(pad.nnode(), 4);
    ///     assert_eq!(pad.geo_nnode(), 8);
    ///     Ok(())
    /// }
    /// ```
    pub fn new_mixed(space_ndim: usize, kind: GeoKind, geo_kind: GeoKind) -> Result<Self, StrError> {
        if geo_kind.ndim() != space_ndim {
            return Err("new_mixed requires geo_kind with geo_ndim = space_ndim");
        }
        if kind.class() != geo_kind.class() || kind.ndim() != geo_kind.ndim() {
            return Err("kind and geo_kind must be of the same class and geo_ndim (e.g. Qua4 with Qua8)");
        }
        if space_ndim < 2 || space_ndim > 3 {
            return Err("space_ndim must be 2 or 3");
        }
        Ok(Self::alloc(space_ndim, kind, geo_kind))
    }

    /// Allocates the arrays for a given field kind and geometry kind
    fn alloc(space_ndim: usize, kind: GeoKind, geo_kind: GeoKind) -> Self {
        let geo_ndim = geo_kind.ndim();
        let nnode = kind.nnode();
        let geo_nnode = geo_kind.nnode();
        let (fn_interp, fn_deriv) = kind.functions();
        let (fn_interp_geo, fn_deriv_geo) = geo_kind.functions();
        Scratchpad {
            kind,
            geo_kind,
            interp: Vector::new(nnode),
            deriv: Matrix::new(nnode, geo_ndim),
            interp_geo: Vector::new(geo_nnode),
            deriv_geo: Matrix::new(geo_nnode, geo_ndim),
            jacobian: Matrix::new(space_ndim, geo_ndim),
            inv_jacobian: if geo_ndim == space_ndim {
                Matrix::new(space_ndim, space_ndim)
            } else {
                Matrix::new(0, 0)
            },
            gradient: if geo_ndim == space_ndim {
                Matrix::new(nnode, space_ndim)
            } else {
                Matrix::new(0, 0)
            },
            xxt: Matrix::new(space_ndim, geo_nnode),
            ok_xxt: false,
            fn_interp,
            fn_deriv,
            fn_interp_geo,
            fn_deriv_geo,
        }
    }

    /// Returns space_ndim, geo_ndim, and the field's nnode
    #[inline]
    pub fn dims(&self) -> (usize, usize, usize) {
        (self.xxt.nrow(), self.kind.ndim(), self.interp.dim())
    }

    /// Returns the space ndim constant (2 or 3)
    #[inline]
    pub fn space_ndim(&self) -> usize {
        self.xxt.nrow()
    }

    /// Returns the number of nodes of the field interpolation (`kind`)
    #[inline]
    pub fn nnode(&self) -> usize {
        self.interp.dim()
    }

    /// Returns the number of nodes of the cell geometry (`geo_kind`)
    #[inline]
    pub fn geo_nnode(&self) -> usize {
        self.xxt.ncol()
    }

    /// Sets the component of the coordinates matrix corresponding to node-m, dimension-j
    ///
    /// # Input
    ///
    /// * `m` -- (local) node index
    /// * `j` -- index of space dimension
    ///
    /// # Panics
    ///
    /// This function will panic if either m or j is out of bounds
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// use gemlab::shapes::{GeoKind, Scratchpad};
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
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
    ///     let space_ndim = 3;
    ///     let mut pad = Scratchpad::new(space_ndim, GeoKind::Hex8)?;
    ///     assert_eq!(pad.ok_xxt, false);
    ///     // set coordinates of node 0
    ///     pad.set_xx(0, 0, 0.0);
    ///     pad.set_xx(0, 1, 0.0);
    ///     pad.set_xx(0, 2, 0.0);
    ///     // set coordinates of node 1
    ///     pad.set_xx(1, 0, 5.0);
    ///     pad.set_xx(1, 1, 0.0);
    ///     pad.set_xx(1, 2, 0.0);
    ///     // set coordinates of node 2
    ///     pad.set_xx(2, 0, 5.0);
    ///     pad.set_xx(2, 1, 5.0);
    ///     pad.set_xx(2, 2, 0.0);
    ///     // set coordinates of node 3
    ///     pad.set_xx(3, 0, 0.0);
    ///     pad.set_xx(3, 1, 5.0);
    ///     pad.set_xx(3, 2, 0.0);
    ///     // set coordinates of node 4
    ///     pad.set_xx(4, 0, 0.0);
    ///     pad.set_xx(4, 1, 0.0);
    ///     pad.set_xx(4, 2, 0.0);
    ///     // set coordinates of node 5
    ///     pad.set_xx(5, 0, 5.0);
    ///     pad.set_xx(5, 1, 0.0);
    ///     pad.set_xx(5, 2, 0.0);
    ///     // set coordinates of node 6
    ///     pad.set_xx(6, 0, 5.0);
    ///     pad.set_xx(6, 1, 5.0);
    ///     pad.set_xx(6, 2, 0.0);
    ///     // set coordinates of node 7
    ///     pad.set_xx(7, 0, 0.0);
    ///     pad.set_xx(7, 1, 5.0);
    ///     assert_eq!(pad.ok_xxt, false);
    ///     pad.set_xx(7, 2, 0.0);
    ///     assert_eq!(pad.ok_xxt, true);
    ///     pad.set_xx(0, 0, 0.0);
    ///     assert_eq!(pad.ok_xxt, false);
    ///     Ok(())
    /// }
    /// ```
    pub fn set_xx(&mut self, m: usize, j: usize, value: f64) {
        self.ok_xxt = false;
        self.xxt.set(j, m, value);
        let space_ndim = self.space_ndim();
        let geo_nnode = self.geo_nnode();
        if m == geo_nnode - 1 && j == space_ndim - 1 {
            self.ok_xxt = true;
        }
    }

    /// Calculates the interpolation functions of the field at a reference coordinate
    ///
    /// The interpolation functions are such that:
    ///
    /// ```text
    /// → →         →  →
    /// u(ξ) = Σ Nᵐ(ξ) uᵐ
    ///        m
    /// ```
    ///
    /// for any quantity `uᵐ` specified at the nodes of an element/shape. Above, `ξ` is the (geo_ndim)
    /// vector of reference coordinates, and `Nᵐ` are the (nnode) interpolation functions of `kind`.
    ///
    /// # Input
    ///
    /// * `ksi` -- set of coordinates in the reference space; its length must be ≥ `geo_ndim`
    ///
    /// # Output
    ///
    /// * `interp` -- the results are available in the `interp` data member
    ///
    /// # Panics
    ///
    /// None of the dimensions are checked, thus a panic my occur if the
    /// arguments have dimensions incompatible with the related [GeoKind].
    #[inline]
    pub fn calc_interp(&mut self, ksi: &[f64]) {
        (self.fn_interp)(&mut self.interp, ksi);
    }

    /// Calculates the axisymmetric radius at a reference coordinate
    ///
    /// The radius is the first real coordinate `x₀(ξ)` of the cell geometry. When `geo_kind` and
    /// `kind` are equal, this function uses `interp`; therefore, call [Scratchpad::calc_gradient] or
    /// [Scratchpad::calc_interp] at the same `ksi` first. When the kinds differ, this function
    /// evaluates the geometry interpolation functions and uses `interp_geo`.
    ///
    /// # Input
    ///
    /// * `ksi` -- reference coordinates ξ with len ≥ `geo_ndim`
    ///
    /// # Returns
    ///
    /// Returns the radius or an error if the coordinates have not been set.
    pub fn calc_radius(&mut self, ksi: &[f64]) -> Result<f64, StrError> {
        if !self.ok_xxt {
            return Err("all components of the coordinates matrix must be set first");
        }
        let mut radius = 0.0;
        if self.geo_kind == self.kind {
            for k in 0..self.nnode() {
                radius += self.interp[k] * self.xxt.get(0, k);
            }
        } else {
            (self.fn_interp_geo)(&mut self.interp_geo, ksi);
            for k in 0..self.geo_nnode() {
                radius += self.interp_geo[k] * self.xxt.get(0, k);
            }
        }
        Ok(radius)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Scratchpad;
    use crate::integ::Gauss;
    use crate::mesh::Samples;
    use crate::recovery::get_points_coords;
    use crate::shapes::{GeoClass, GeoKind};
    use russell_lab::{approx_eq, mat_approx_eq, mat_mat_mul, Matrix};
    use russell_tensor::{t1_approx_eq, Tensor1};

    #[test]
    fn derive_works() {
        let pad = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let pad_clone = pad.clone();
        assert_eq!(pad_clone.kind, pad.kind);
        assert_eq!(pad_clone.geo_kind, pad.geo_kind);
    }

    #[test]
    fn new_fails_on_wrong_input() {
        assert_eq!(
            Scratchpad::new(1, GeoKind::Lin2).err(),
            Some("space_ndim must be 2 or 3")
        );
        assert_eq!(
            Scratchpad::new(2, GeoKind::Hex8).err(),
            Some("geo_ndim cannot be greater than space_ndim")
        );
    }

    #[test]
    fn new_works() {
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
            let kind = GeoKind::VALUES[i];
            let space_ndim = usize::max(2, kind.ndim());
            let pad = Scratchpad::new(space_ndim, kind).unwrap();
            assert_eq!(pad.dims(), (space_ndim, geo_ndim[i], nnode[i]));
            assert_eq!(pad.space_ndim(), space_ndim);
            assert_eq!(pad.nnode(), nnode[i]);
            assert_eq!(pad.geo_nnode(), nnode[i]);
            assert_eq!(pad.kind, kind);
            assert_eq!(pad.geo_kind, kind);
            assert_eq!(pad.interp.dim(), nnode[i]);
            assert_eq!(pad.interp_geo.dim(), nnode[i]);
            assert_eq!(pad.deriv.dims(), (nnode[i], geo_ndim[i]));
            assert_eq!(pad.deriv_geo.dims(), (nnode[i], geo_ndim[i]));
            assert_eq!(pad.jacobian.dims(), (space_ndim, geo_ndim[i]));
            if kind.class() == GeoClass::Lin {
                assert_eq!(pad.inv_jacobian.dims(), (0, 0));
                assert_eq!(pad.gradient.dims(), (0, 0));
            } else {
                assert_eq!(pad.inv_jacobian.dims(), (space_ndim, space_ndim));
                assert_eq!(pad.gradient.dims(), (nnode[i], space_ndim));
            }
            assert_eq!(pad.xxt.dims(), (space_ndim, nnode[i]));
            assert_eq!(pad.ok_xxt, false);
        }
    }

    #[test]
    fn new_mixed_works() {
        // space_ndim, field kind, geometry kind
        let problem = [
            (2, GeoKind::Qua4, GeoKind::Qua8),
            (2, GeoKind::Qua8, GeoKind::Qua4),
            (2, GeoKind::Qua4, GeoKind::Qua4),
            (2, GeoKind::Tri3, GeoKind::Tri6),
            (2, GeoKind::Tri6, GeoKind::Tri3),
            (3, GeoKind::Tet4, GeoKind::Tet10),
            (3, GeoKind::Tet10, GeoKind::Tet4),
            (3, GeoKind::Hex8, GeoKind::Hex20),
            (3, GeoKind::Hex20, GeoKind::Hex8),
        ];
        for (space_ndim, kind, geo_kind) in problem {
            let pad = Scratchpad::new_mixed(space_ndim, kind, geo_kind).unwrap();
            let nnode = kind.nnode();
            let geo_nnode = geo_kind.nnode();
            assert_eq!(pad.kind, kind);
            assert_eq!(pad.geo_kind, geo_kind);
            assert_eq!(pad.dims(), (space_ndim, geo_kind.ndim(), nnode));
            assert_eq!(pad.space_ndim(), space_ndim);
            assert_eq!(pad.nnode(), nnode);
            assert_eq!(pad.geo_nnode(), geo_nnode);
            assert_eq!(pad.interp.dim(), nnode);
            assert_eq!(pad.deriv.dims(), (nnode, geo_kind.ndim()));
            assert_eq!(pad.interp_geo.dim(), geo_nnode);
            assert_eq!(pad.deriv_geo.dims(), (geo_nnode, geo_kind.ndim()));
            assert_eq!(pad.jacobian.dims(), (space_ndim, geo_kind.ndim()));
            assert_eq!(pad.inv_jacobian.dims(), (space_ndim, space_ndim));
            assert_eq!(pad.gradient.dims(), (nnode, space_ndim));
            assert_eq!(pad.xxt.dims(), (space_ndim, geo_nnode));
            assert_eq!(pad.ok_xxt, false);
        }
    }

    #[test]
    fn new_mixed_handles_errors() {
        // the geometry kind must match the space dimension
        assert_eq!(
            Scratchpad::new_mixed(2, GeoKind::Lin2, GeoKind::Lin2).err(),
            Some("new_mixed requires geo_kind with geo_ndim = space_ndim")
        );
        assert_eq!(
            Scratchpad::new_mixed(3, GeoKind::Tri3, GeoKind::Tri3).err(),
            Some("new_mixed requires geo_kind with geo_ndim = space_ndim")
        );
        assert_eq!(
            Scratchpad::new_mixed(4, GeoKind::Hex8, GeoKind::Hex8).err(),
            Some("new_mixed requires geo_kind with geo_ndim = space_ndim")
        );

        // the field must be of the same class and geometric dimension
        assert_eq!(
            Scratchpad::new_mixed(2, GeoKind::Tri3, GeoKind::Qua4).err(),
            Some("kind and geo_kind must be of the same class and geo_ndim (e.g. Qua4 with Qua8)")
        );

        // error from the space dimension check
        assert_eq!(
            Scratchpad::new_mixed(1, GeoKind::Lin2, GeoKind::Lin2).err(),
            Some("space_ndim must be 2 or 3")
        );

        // coordinates not set
        let mut pad = Scratchpad::new_mixed(2, GeoKind::Qua4, GeoKind::Qua8).unwrap();
        assert_eq!(
            pad.calc_gradient(&[0.0, 0.0]).err(),
            Some("all components of the coordinates matrix must be set first")
        );
    }

    #[test]
    fn set_xx_works() {
        let mut pad = Scratchpad::new(2, GeoKind::Lin2).unwrap();
        assert_eq!(pad.ok_xxt, false);
        pad.set_xx(0, 0, 1.0);
        assert_eq!(pad.ok_xxt, false);
        pad.set_xx(0, 1, 2.0);
        assert_eq!(pad.ok_xxt, false);
        pad.set_xx(1, 0, 3.0);
        assert_eq!(pad.ok_xxt, false);
        pad.set_xx(1, 1, 4.0);
        assert_eq!(pad.ok_xxt, true);
        pad.set_xx(0, 0, 1.0); // reset
        assert_eq!(pad.ok_xxt, false);
    }

    #[test]
    fn set_xx_uses_geo_nnode() {
        // the last field node is not the last geometry node
        let mut pad = Scratchpad::new_mixed(2, GeoKind::Qua4, GeoKind::Qua8).unwrap();
        pad.set_xx(3, 0, 0.0);
        pad.set_xx(3, 1, 0.0);
        assert_eq!(pad.ok_xxt, false);
        pad.set_xx(7, 0, 1.0);
        assert_eq!(pad.ok_xxt, false);
        pad.set_xx(7, 1, 1.0);
        assert_eq!(pad.ok_xxt, true);
    }

    /// Coefficients of the linear field used in the patch tests: a₀ + a · x
    const A: [f64; 3] = [2.0, -3.0, 0.5];

    /// Returns the linear field a₀ + a · x at the point with coordinates x
    fn linear_field(x: &[f64]) -> f64 {
        let mut value = 1.0;
        for j in 0..x.len() {
            value += A[j] * x[j];
        }
        value
    }

    /// Computes the area of a straight-edged 2D element with the shoelace formula on its `ncorners`
    fn shoelace(pad: &Scratchpad, ncorners: usize) -> f64 {
        let mut sum = 0.0;
        for m in 0..ncorners {
            let next = if m + 1 == ncorners { 0 } else { m + 1 };
            sum += pad.xxt.get(0, m) * pad.xxt.get(1, next) - pad.xxt.get(0, next) * pad.xxt.get(1, m);
        }
        0.5 * sum
    }

    /// Checks the invariants of a mixed-order pad at every Gauss point of the geometry kind
    ///
    /// Requires the field nodes to be the first nodes of the geometry (the usual mixed element).
    ///
    /// 1. `Σ N = 1` for the field and geometry shapes (machine precision)
    /// 2. `Σ B = 0` per column (rigid-body translation)
    /// 3. a linear field is reproduced by the field shape and its gradient
    /// 4. `Σ_i Σ_j ∫ N_i N_j dΩ = measure` and `Σ_i Σ_j ∫ N_i B_geo_j dΩ = 0`
    fn check_invariants(pad: &mut Scratchpad, ngauss: usize, measure: f64) {
        let space_ndim = pad.space_ndim();
        let n_geo = pad.geo_nnode();
        let k = pad.nnode();
        assert!(k <= n_geo);
        let gauss = Gauss::new_or_sized(pad.geo_kind, Some(ngauss)).unwrap();

        // nodal coordinates of the field shape (its nodes are the first k nodes of the geometry)
        let mut xr = vec![vec![0.0; space_ndim]; k];
        for m in 0..k {
            for j in 0..space_ndim {
                xr[m][j] = pad.xxt.get(j, m);
            }
        }

        // values of the linear field at the field nodes
        let phi: Vec<f64> = (0..k).map(|m| linear_field(&xr[m])).collect();

        // integral accumulators: Σ_ij ∫ N_i N_j dΩ and Σ_ij ∫ N_i B_geo_j dΩ
        let mut norm_int = 0.0;
        let mut mixed_int = vec![0.0; space_ndim];

        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad.calc_gradient(iota).unwrap();
            assert!(det_jac > 0.0);
            let w = det_jac * gauss.weight(p);

            // coordinates and geometry interpolation at the current point
            let mut x = Tensor1::new();
            pad.calc_coords(&mut x, iota).unwrap();
            let xyz = [x.get(0), x.get(1), x.get(2)];
            let interp_geo = if pad.geo_kind == pad.kind {
                &pad.interp
            } else {
                &pad.interp_geo
            };

            // 1. partition of unity of the geometry and field shapes
            let mut sum_geo = 0.0;
            for m in 0..n_geo {
                sum_geo += interp_geo[m];
            }
            approx_eq(sum_geo, 1.0, 1e-14);
            let mut sum_field = 0.0;
            for m in 0..k {
                sum_field += pad.interp[m];
            }
            approx_eq(sum_field, 1.0, 1e-14);

            // 2. the field gradient sums to zero over the nodes
            for j in 0..space_ndim {
                let mut sum = 0.0;
                for m in 0..k {
                    sum += pad.gradient.get(m, j);
                }
                approx_eq(sum, 0.0, 1e-14);
            }

            // 3. the linear field is reproduced at x(ξ), and its gradient is exact
            let mut value = 0.0;
            for m in 0..k {
                value += pad.interp[m] * phi[m];
            }
            approx_eq(value, linear_field(&xyz[..space_ndim]), 1e-13);
            for j in 0..space_ndim {
                let mut g = 0.0;
                for m in 0..k {
                    g += pad.gradient.get(m, j) * phi[m];
                }
                approx_eq(g, A[j], 1e-13);
            }

            // geometry gradient B_geo = dN_geo/dξ · J_geo⁻¹
            let b_geo = if pad.geo_kind == pad.kind {
                pad.gradient.clone()
            } else {
                (pad.fn_deriv_geo)(&mut pad.deriv_geo, iota);
                let mut b = Matrix::new(n_geo, space_ndim);
                mat_mat_mul(&mut b, 1.0, &pad.deriv_geo, &pad.inv_jacobian, 0.0).unwrap();
                b
            };

            // 4. accumulate the two integrals
            for i in 0..k {
                for j in 0..k {
                    norm_int += pad.interp[i] * pad.interp[j] * w;
                }
                for j in 0..n_geo {
                    for d in 0..space_ndim {
                        mixed_int[d] += pad.interp[i] * b_geo.get(j, d) * w;
                    }
                }
            }
        }

        // Σ_ij ∫ N_i N_j dΩ = ∫ (Σ N)² dΩ = measure
        approx_eq(norm_int, measure, 1e-12);

        // Σ_ij ∫ N_i B_geo_j dΩ = ∫ (Σ N) (Σ B_geo) dΩ = 0
        for d in 0..space_ndim {
            approx_eq(mixed_int[d], 0.0, 1e-12);
        }
    }

    #[test]
    fn get_pad_alt_qua8() {
        // unit square: Qua8 geometry with a Qua4 field
        let mesh = Samples::one_qua8();
        let mut pad = mesh.get_pad_alt(0, GeoKind::Qua4);
        assert_eq!(pad.kind, GeoKind::Qua4);
        assert_eq!(pad.geo_kind, GeoKind::Qua8);
        assert_eq!(pad.geo_nnode(), 8);
        assert_eq!(pad.nnode(), 4);
        check_invariants(&mut pad, 9, 1.0);
    }

    #[test]
    fn get_pad_alt_tri6() {
        // triangle with base 1 and height 0.85: Tri6 geometry with a Tri3 field, area = 0.425
        let mesh = Samples::one_tri6();
        let mut pad = mesh.get_pad_alt(0, GeoKind::Tri3);
        assert_eq!(pad.kind, GeoKind::Tri3);
        assert_eq!(pad.geo_kind, GeoKind::Tri6);
        assert_eq!(pad.geo_nnode(), 6);
        assert_eq!(pad.nnode(), 3);
        check_invariants(&mut pad, 7, 0.425);
    }

    #[test]
    fn get_pad_alt_hex20() {
        // unit cube: Hex20 geometry with a Hex8 field, volume = 1
        let mesh = Samples::one_hex20(0.0, 0.0, 0.0);
        let mut pad = mesh.get_pad_alt(0, GeoKind::Hex8);
        assert_eq!(pad.kind, GeoKind::Hex8);
        assert_eq!(pad.geo_kind, GeoKind::Hex20);
        assert_eq!(pad.geo_nnode(), 20);
        assert_eq!(pad.nnode(), 8);
        check_invariants(&mut pad, 27, 1.0);
    }

    #[test]
    fn distorted_qua8_matches_the_reference_values() {
        // Distorted quadrilateral with straight edges: the midside nodes are the edge midpoints, so
        // the area is exactly the shoelace area of the corners (2.195), as measured in the
        // investigation recorded in the plan note.
        let corners = [(0.0, 0.0), (2.0, 0.1), (1.7, 1.3), (-0.2, 1.0)];
        let mut pad = Scratchpad::new_mixed(2, GeoKind::Qua4, GeoKind::Qua8).unwrap();
        for m in 0..4 {
            pad.set_xx(m, 0, corners[m].0);
            pad.set_xx(m, 1, corners[m].1);
        }
        // midside nodes: 4 = middle(0,1), 5 = middle(1,2), 6 = middle(2,3), 7 = middle(3,0)
        let edges = [(0, 1), (1, 2), (2, 3), (3, 0)];
        for (i, (a, b)) in edges.iter().enumerate() {
            let m = 4 + i;
            pad.set_xx(m, 0, 0.5 * (corners[*a].0 + corners[*b].0));
            pad.set_xx(m, 1, 0.5 * (corners[*a].1 + corners[*b].1));
        }

        // the area comes from the corners only (straight edges)
        approx_eq(shoelace(&pad, 4), 2.195, 1e-14);

        // 9 Gauss points, as in the investigation
        check_invariants(&mut pad, 9, 2.195);
    }

    #[test]
    fn equal_order_works() {
        // Qua4 with Qua4: the field and geometry shapes coincide
        let mut pad = Scratchpad::new_mixed(2, GeoKind::Qua4, GeoKind::Qua4).unwrap();
        assert_eq!(pad.kind, GeoKind::Qua4);
        assert_eq!(pad.geo_kind, GeoKind::Qua4);
        assert_eq!(pad.geo_nnode(), 4);
        assert_eq!(pad.nnode(), 4);

        // skewed quadrilateral
        let xx = [(0.0, 0.0), (2.0, 0.1), (1.7, 1.3), (-0.2, 1.0)];
        for m in 0..4 {
            pad.set_xx(m, 0, xx[m].0);
            pad.set_xx(m, 1, xx[m].1);
        }

        // compare with a plain scratchpad at the same ξ
        let at_ksi = [-0.5, 0.3];
        let det_jac = pad.calc_gradient(&at_ksi).unwrap();
        let mut pad_ref = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad_ref.set_xx(m, 0, xx[m].0);
            pad_ref.set_xx(m, 1, xx[m].1);
        }
        let det_jac_ref = pad_ref.calc_gradient(&at_ksi).unwrap();
        approx_eq(det_jac, det_jac_ref, 1e-15);
        mat_approx_eq(&pad.gradient, &pad_ref.gradient, 1e-15);

        // the geometry gradient also matches the reference
        (pad.fn_deriv_geo)(&mut pad.deriv_geo, &at_ksi);
        let mut b_geo = Matrix::new(4, 2);
        mat_mat_mul(&mut b_geo, 1.0, &pad.deriv_geo, &pad.inv_jacobian, 0.0).unwrap();
        mat_approx_eq(&b_geo, &pad_ref.gradient, 1e-15);

        for m in 0..4 {
            approx_eq(pad.interp[m], pad_ref.interp[m], 1e-15);
        }

        // the invariants also hold for the equal-order case
        check_invariants(&mut pad, 9, 2.195);
    }

    #[test]
    fn reversed_order_field_qua8_over_qua4_geometry() {
        // The geometry is Qua4; the field is quadratic. Its gradient must be
        // dN_field/dξ · J_geo⁻¹, and the field nodes beyond the geometry nodes have no coordinates.
        let corners = [(0.0, 0.0), (2.0, 0.1), (1.7, 1.3), (-0.2, 1.0)];
        let mut pad = Scratchpad::new_mixed(2, GeoKind::Qua8, GeoKind::Qua4).unwrap();
        assert_eq!(pad.kind, GeoKind::Qua8);
        assert_eq!(pad.geo_kind, GeoKind::Qua4);
        assert_eq!(pad.geo_nnode(), 4);
        assert_eq!(pad.nnode(), 8);
        for m in 0..4 {
            pad.set_xx(m, 0, corners[m].0);
            pad.set_xx(m, 1, corners[m].1);
        }

        // the field derivatives are dN/dξ; the reference multiplies them by the geometry's J⁻¹
        let mut pad_field = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        let gauss = Gauss::new_or_sized(GeoKind::Qua4, Some(9)).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            pad.calc_gradient(iota).unwrap();
            (pad_field.fn_deriv)(&mut pad_field.deriv, iota);

            let mut sum_field = 0.0;
            for m in 0..8 {
                sum_field += pad.interp[m];
            }
            approx_eq(sum_field, 1.0, 1e-14);

            for m in 0..8 {
                for d in 0..2 {
                    let mut b = 0.0;
                    for j in 0..2 {
                        b += pad_field.deriv.get(m, j) * pad.inv_jacobian.get(j, d);
                    }
                    approx_eq(pad.gradient.get(m, d), b, 1e-15);
                }
            }
            for d in 0..2 {
                let mut sum = 0.0;
                for m in 0..8 {
                    sum += pad.gradient.get(m, d);
                }
                approx_eq(sum, 0.0, 1e-13);
            }
        }
    }

    #[test]
    fn calc_radius_works() {
        // single-order shape
        let mut pad = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let xx = [(2.0, 1.0), (6.0, 1.0), (6.0, 3.0), (2.0, 3.0)];
        for m in 0..4 {
            pad.set_xx(m, 0, xx[m].0);
            pad.set_xx(m, 1, xx[m].1);
        }
        let at_ksi = [0.3, -0.2];
        pad.calc_gradient(&at_ksi).unwrap();
        let radius = pad.calc_radius(&at_ksi).unwrap();
        let mut x = Tensor1::new();
        pad.calc_coords(&mut x, &at_ksi).unwrap();
        approx_eq(radius, x.get(0), 1e-15);

        // mixed-order shape with a curved Qua8 geometry and a Qua4 field
        let mut pad = Scratchpad::new_mixed(2, GeoKind::Qua4, GeoKind::Qua8).unwrap();
        let coords = [
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.5, -0.25), // curved bottom edge
            (1.0, 0.5),
            (0.5, 1.0),
            (0.0, 0.5),
        ];
        for m in 0..8 {
            pad.set_xx(m, 0, coords[m].0);
            pad.set_xx(m, 1, coords[m].1);
        }
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            pad.calc_gradient(iota).unwrap();
            let radius = pad.calc_radius(iota).unwrap();
            pad.calc_coords(&mut x, iota).unwrap();
            approx_eq(radius, x.get(0), 1e-15);
        }

        // coordinates not set
        let mut pad = Scratchpad::new_mixed(2, GeoKind::Qua4, GeoKind::Qua8).unwrap();
        assert_eq!(
            pad.calc_radius(&[0.0, 0.0]).err(),
            Some("all components of the coordinates matrix must be set first")
        );
    }

    #[test]
    fn mixed_order_uses_the_cell_geometry() {
        // curved Qua8 geometry with a Qua4 field
        let coords = [
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.5, -0.25), // curved bottom edge
            (1.0, 0.5),
            (0.5, 1.0),
            (0.0, 0.5),
        ];
        let mut pad = Scratchpad::new_mixed(2, GeoKind::Qua4, GeoKind::Qua8).unwrap();
        let mut pad_geo = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for m in 0..8 {
            pad.set_xx(m, 0, coords[m].0);
            pad.set_xx(m, 1, coords[m].1);
            pad_geo.set_xx(m, 0, coords[m].0);
            pad_geo.set_xx(m, 1, coords[m].1);
        }

        // the Jacobian must match the single-order geometry
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad.calc_jacobian(iota).unwrap();
            let det_jac_geo = pad_geo.calc_jacobian(iota).unwrap();
            approx_eq(det_jac, det_jac_geo, 1e-15);
            mat_approx_eq(&pad.jacobian, &pad_geo.jacobian, 1e-15);
        }

        // get_points_coords must return the exact Gauss coordinates of the curved cell
        let x_mixed = get_points_coords(&mut pad, &gauss).unwrap();
        let x_geo = get_points_coords(&mut pad_geo, &gauss).unwrap();
        for p in 0..gauss.npoint() {
            t1_approx_eq(&x_mixed[p], &x_geo[p], 1e-15);
        }
    }
}
