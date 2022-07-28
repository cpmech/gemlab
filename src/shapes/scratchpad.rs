use super::{FnDeriv, FnInterp, GeoKind};
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Holds the variables used by the functions in the sub-module named op
///
/// Basically, this struct holds the variables calculated by the interpolation
/// and derivative functions, including the Jacobian, inverse Jacobian, and gradients.
///
/// # Sequence of computations for the gradient
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
/// J(ξ) = —— = Σ xᵐ ⊗ Lᵐ
///         →   m
///        dξ
/// ```
///
/// ```text
///            J         =        Xᵀ        ·       L
/// (space_ndim,geo_ndim) (space_ndim,nnode) (nnode,geo_ndim)
/// ```
///
/// ## 3 Inverse Jacobian
///
/// ## 4 Gradient @ ξ
///
/// ```text
///             →
/// →  →    dNᵐ(ξ)
/// Gᵐ(ξ) = ——————
///            →
///           dx
/// ```
///
/// ```text
///        G          =       L        ·           J⁻¹
/// (nnode,space_ndim) (nnode,geo_ndim) (space_ndim,space_ndim)
/// ```
///
/// # Warning
///
/// All members are **readonly** and should not be modified externally. The functions
/// in the op (operators) sub-module are the only ones allowed to modify the scratchpad.
#[derive(Clone)]
pub struct Scratchpad {
    /// The kind of the shape
    pub kind: GeoKind,

    /// Array N: (nnode) Nᵐ interpolation functions at reference coordinate ξ (ksi)
    pub interp: Vector,

    /// Matrix L: (nnode,geo_ndim) dNᵐ/dξ derivatives of interpolation functions with respect to reference coordinate ξ (ksi)
    pub deriv: Matrix,

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

    /// Matrix G: (nnode,space_ndim) dNᵐ/dx Gradient of shape functions (only if geo_ndim = space_ndim) at ξ (ksi)
    ///
    /// Only available if `geo_ndim = space_ndim` (otherwise, the matrix is set to empty; 0 x 0 matrix)
    ///
    /// ```text
    /// G = L · J⁻¹
    /// ```
    pub gradient: Matrix,

    /// Matrix Xᵀ: (space_ndim,nnode) transposed coordinates matrix (real space)
    ///
    /// ```text
    ///      ┌                              ┐  superscript = node
    ///      | x⁰₀  x¹₀  x²₀  x³₀       xᴹ₀ |  subscript = space dimension
    /// Xᵀ = | x⁰₁  x¹₁  x²₁  x³₁  ...  xᴹ₁ |
    ///      | x⁰₂  x¹₂  x²₂  x³₂       xᴹ₂ |
    ///      └                              ┘_(space_ndim,nnode)
    /// where `M = nnode - 1`
    /// ```
    ///
    /// You must call `set_xx` to set the components, otherwise calculations will be **incorrect.**
    pub xxt: Matrix,

    /// Minimum coordinates (space_ndim)
    pub xmin: Vec<f64>,

    /// Maximum coordinates (space_ndim)
    pub xmax: Vec<f64>,

    /// Indicates that all coordinates in the X matrix have been set and the min,max values computed
    pub ok_xxt: bool,

    /// Function to calculate interpolation functions
    pub fn_interp: FnInterp,

    /// Function to calculate derivatives of interpolation functions
    pub fn_deriv: FnDeriv,
}

impl Scratchpad {
    /// Allocates a new instance
    ///
    /// # Notes
    ///
    /// 1. `space_ndim` must be 2 or 3
    /// 2. `space_ndim` must be greater than or equal to `geo_ndim`;
    ///     e.g., you cannot have a 3D shape in a 2D space.
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
        let nnode = kind.nnode();
        let (fn_interp, fn_deriv) = kind.functions();
        Ok(Scratchpad {
            kind,
            interp: Vector::new(nnode),
            deriv: Matrix::new(nnode, geo_ndim),
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
            xxt: Matrix::new(space_ndim, nnode),
            xmin: vec![f64::MAX; space_ndim],
            xmax: vec![f64::MIN; space_ndim],
            ok_xxt: false,
            fn_interp,
            fn_deriv,
        })
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
        self.xxt[j][m] = value;
        self.xmin[j] = f64::min(self.xmin[j], self.xxt[j][m]);
        self.xmax[j] = f64::max(self.xmax[j], self.xxt[j][m]);
        let (space_ndim, nnode) = self.xxt.dims();
        if m == nnode - 1 && j == space_ndim - 1 {
            self.ok_xxt = true;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Scratchpad;
    use crate::shapes::{GeoClass, GeoKind};

    #[test]
    fn derive_works() {
        let pad = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let pad_clone = pad.clone();
        assert_eq!(pad_clone.kind, pad.kind);
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
            assert_eq!(pad.kind, kind);
            assert_eq!(pad.interp.dim(), nnode[i]);
            assert_eq!(pad.deriv.dims(), (nnode[i], geo_ndim[i]));
            assert_eq!(pad.jacobian.dims(), (space_ndim, geo_ndim[i]));
            if kind.class() == GeoClass::Lin {
                assert_eq!(pad.inv_jacobian.dims(), (0, 0));
                assert_eq!(pad.gradient.dims(), (0, 0));
            } else {
                assert_eq!(pad.inv_jacobian.dims(), (space_ndim, space_ndim));
                assert_eq!(pad.gradient.dims(), (nnode[i], space_ndim));
            }
            assert_eq!(pad.xxt.dims(), (space_ndim, nnode[i]));
            assert_eq!(pad.xmin.len(), space_ndim);
            assert_eq!(pad.xmax.len(), space_ndim);
            assert_eq!(pad.ok_xxt, false);
        }
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
}
