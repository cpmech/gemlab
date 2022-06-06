use super::GeoKind;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// The Matrix of coordinates is:
///
/// ```text
///     ┌               ┐
///     | x⁰₀  x⁰₁  x⁰₂ |  
///     | x¹₀  x¹₁  x¹₂ |  superscript = node
/// X = | x²₀  x²₁  x²₂ |  subscript = dimension
///     | x³₀  x³₁  x³₂ |
///     |      ...      |
///     | xᴹ₀  xᴹ₁  xᴹ₂ |
///     └               ┘_(nnode,space_ndim)
/// ```
///
/// where `M = nnode - 1`
///
/// # Warnings
///
/// 1. All public members here are supposed to be **readonly**.
/// 2. The Scratchpad must be initialized by setting the coordinates matrix `X` first.
#[derive(Clone, Debug)]
pub struct Scratchpad {
    /// Array N: (nnode) interpolation functions at reference coordinate ksi
    pub interp: Vector,

    /// Matrix L: (nnode,geo_ndim) derivatives of interpolation functions w.r.t reference coordinate ksi
    pub deriv: Matrix,

    /// Matrix J: (space_ndim,geo_ndim) Jacobian matrix
    pub jacobian: Matrix,

    /// Matrix inv(J): (space_ndim,space_ndim) Inverse Jacobian matrix (only if geo_ndim == space_ndim) at ksi
    ///
    /// Only available if `geo_ndim == space_ndim` (otherwise, the matrix is set to empty; 0 x 0 matrix)
    pub inv_jacobian: Matrix,

    /// Matrix G: (nnode,space_ndim) Gradient of shape functions (only if geo_ndim == space_ndim) at ksi
    ///
    /// Only available if `geo_ndim == space_ndim` (otherwise, the matrix is set to empty; 0 x 0 matrix)
    pub gradient: Matrix,

    /// Matrix Xᵀ: (space_ndim,nnode) transposed coordinates matrix (real space)
    ///
    /// ```text
    ///      ┌                              ┐  superscript = node
    ///      | x⁰₀  x¹₀  x²₀  x³₀       xᴹ₀ |  subscript = dimension
    /// Xᵀ = | x⁰₁  x¹₁  x²₁  x³₁  ...  xᴹ₁ |
    ///      | x⁰₂  x¹₂  x²₂  x³₂       xᴹ₂ |
    ///      └                              ┘_(space_ndim,nnode)
    /// ```
    ///
    /// where `M = nnode - 1`
    xx_transp: Matrix,

    /// Indicates that all coordinates in the X matrix have been set (at least once)
    initialized: bool,
}

impl Scratchpad {
    /// Allocates a new instance
    ///
    /// # Notes
    ///
    /// 1. `space_ndim` must be 2 or 3
    /// 2. `space_ndim` must be greater than or equal to `geo_ndim`;
    ///     e.g., you cannot have a 3D shape in a 2D space.
    pub fn new(space_ndim: usize, kind: GeoKind) -> Result<Self, StrError> {
        if space_ndim < 2 || space_ndim > 3 {
            return Err("space_ndim must be 2 or 3");
        }
        let geo_ndim = kind.ndim();
        if geo_ndim > space_ndim {
            return Err("space_ndim must be ≥ geo_ndim");
        }
        let nnode = kind.nnode();
        Ok(Scratchpad {
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
            xx_transp: Matrix::new(space_ndim, nnode),
            initialized: false,
        })
    }

    /// Returns whether the matrix of coordinates has been fully initialized or not (at least once)
    pub fn initialized(&self) -> bool {
        self.initialized
    }

    /// Sets the m-node's coordinate j
    ///
    /// # Input
    ///
    /// * `m` -- (local) node index
    /// * `j` -- index of space dimension
    ///
    /// # Panics
    ///
    /// This function will panic if the values m or j are out of bounds
    pub fn set_xx(&mut self, m: usize, j: usize, value: f64) {
        self.initialized = false;
        self.xx_transp[j][m] = value;
        let (space_ndim, nnode) = self.xx_transp.dims();
        if m == nnode - 1 && j == space_ndim - 1 {
            self.initialized = true;
        }
    }

    /// Returns the m-node's coordinate j
    ///
    /// # Input
    ///
    /// * `m` -- (local) node index
    /// * `j` -- index of space dimension
    ///
    /// # Panics
    ///
    /// This function will panic if the values m or j are out of bounds
    pub fn xx(&self, m: usize, j: usize) -> f64 {
        self.xx_transp[j][m]
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Scratchpad;
    use crate::shapes::GeoKind;

    #[test]
    fn new_fails_on_wrong_input() {
        assert_eq!(
            Scratchpad::new(1, GeoKind::Lin2).err(),
            Some("space_ndim must be 2 or 3")
        );
        assert_eq!(
            Scratchpad::new(2, GeoKind::Hex8).err(),
            Some("space_ndim must be ≥ geo_ndim")
        );
    }
}
