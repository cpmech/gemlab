use super::GeoKind;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Holds mutable data used by the functions in the sub-module named op
///
/// Basically, this struct holds the variables calculated by the interpolation
/// and derivative functions, including the Jacobian, inverse Jacobian, and gradients.
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
        })
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
