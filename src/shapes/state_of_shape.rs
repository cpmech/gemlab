use super::geo_class_and_kind;
use crate::util::AsArray2D;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Holds the state of a shape entity
#[derive(Clone, Debug)]
pub struct StateOfShape {
    // -- steady ----------------------------------------------------------------------------------------------------
    //
    /// Matrix Xᵀ: (space_ndim,nnode) transposed coordinates matrix (real space)
    ///
    /// ```text
    /// The (non-transposed) matrix X is such that:
    ///
    ///     ┌               ┐
    ///     | x⁰₀  x⁰₁  x⁰₂ |  subscript = dimension
    ///     | x¹₀  x¹₁  x¹₂ |  superscript = node
    /// X = | x²₀  x²₁  x²₂ |
    ///     |      ...      |
    ///     | xᴹ₀  xᴹ₁  xᴹ₂ |
    ///     └               ┘_(nnode,space_ndim)
    ///
    /// where `M = nnode - 1`
    /// ```
    pub coords_transp: Matrix,

    /// Minimum (space_ndim) coordinates from the X matrix
    pub coords_min: Vec<f64>,

    /// Maximum (space_ndim) coordinates from the X matrix
    pub coords_max: Vec<f64>,

    // -- transient -------------------------------------------------------------------------------------------------
    //
    /// Array N: (nnode) interpolation functions at reference coordinate ksi
    pub interp: Vector,

    /// Matrix L: (nnode,geo_ndim) derivatives of interpolation functions w.r.t reference coordinate ksi
    pub deriv: Matrix,

    /// Matrix J: (space_ndim,geo_ndim) Jacobian matrix
    pub jacobian: Matrix,

    /// Matrix inv(J): (space_ndim,space_ndim) Inverse Jacobian matrix (only if geo_ndim == space_ndim) at ksi
    pub inv_jacobian: Matrix,

    /// Matrix G: (nnode,space_ndim) Gradient of shape functions (only if geo_ndim == space_ndim) at ksi
    pub gradient: Matrix,
}

impl StateOfShape {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `geo_ndim` -- the dimension of the shape; e.g. a 1D line in a 3D space
    /// * `coords` -- the (non-transposed) matrix of coordinates (`X` matrix); see below
    ///
    /// ```text
    ///          ┌               ┐
    ///          | x⁰₀  x⁰₁  x⁰₂ |  subscript = dimension
    ///          | x¹₀  x¹₁  x¹₂ |  superscript = node
    /// coords = | x²₀  x²₁  x²₂ |
    ///          |      ...      |
    ///          | xᴹ₀  xᴹ₁  xᴹ₂ |
    ///          └               ┘_(nnode,space_ndim)
    /// ```
    pub fn new<'a, T>(geo_ndim: usize, coords: &'a T) -> Result<Self, StrError>
    where
        T: AsArray2D<'a, f64>,
    {
        let (nnode, space_ndim) = coords.size();
        if space_ndim < 2 || space_ndim > 3 {
            return Err("space_ndim must be 2 or 3");
        }
        geo_class_and_kind(geo_ndim, nnode)?; // just to check combination
        let mut state = StateOfShape {
            coords_transp: Matrix::new(space_ndim, nnode),
            coords_min: vec![f64::MAX; space_ndim],
            coords_max: vec![f64::MIN; space_ndim],
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
        };
        for m in 0..nnode {
            for j in 0..space_ndim {
                let value = coords.at(m, j);
                state.coords_transp[j][m] = value;
                if value < state.coords_min[j] {
                    state.coords_min[j] = value;
                }
                if value > state.coords_max[j] {
                    state.coords_max[j] = value;
                }
            }
        }
        Ok(state)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::shapes::StateOfShape;
    use crate::StrError;

    #[test]
    fn new_fails_on_wrong_input() {
        assert_eq!(
            StateOfShape::new(100, &[[1.0, 1.0], [2.0, 2.0]]).err(),
            Some("(geo_ndim,nnode) combination is invalid")
        );
        assert_eq!(StateOfShape::new(1, &[[], []]).err(), Some("space_ndim must be 2 or 3"));
        assert_eq!(
            StateOfShape::new(1, &[[0.0], [0.0]]).err(),
            Some("space_ndim must be 2 or 3")
        );
        assert_eq!(
            StateOfShape::new(1, &[[1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]]).err(),
            Some("space_ndim must be 2 or 3")
        );
    }

    #[test]
    fn new_works() -> Result<(), StrError> {
        let space_ndim = 2;
        let geo_ndim = 1;
        let nnode = 2;
        let state = StateOfShape::new(geo_ndim, &[[-1.23, 1.23], [-4.56, 4.56]])?;
        assert_eq!(
            format!("{}", state.coords_transp),
            "┌             ┐\n\
             │ -1.23 -4.56 │\n\
             │  1.23  4.56 │\n\
             └             ┘"
        );
        assert_eq!(state.coords_min, &[-4.56, 1.23]);
        assert_eq!(state.coords_max, &[-1.23, 4.56]);
        assert_eq!(state.interp.dim(), nnode);
        assert_eq!(state.deriv.dims(), (nnode, geo_ndim));
        assert_eq!(state.jacobian.dims(), (space_ndim, geo_ndim));
        assert_eq!(state.inv_jacobian.dims(), (0, 0));
        assert_eq!(state.gradient.dims(), (0, 0));

        let space_ndim = 2;
        let geo_ndim = 2;
        let nnode = 3;
        let state = StateOfShape::new(geo_ndim, &[[-1.23, 1.23], [-4.56, 4.56], [-7.89, 7.89]])?;
        assert_eq!(
            format!("{}", state.coords_transp),
            "┌                   ┐\n\
             │ -1.23 -4.56 -7.89 │\n\
             │  1.23  4.56  7.89 │\n\
             └                   ┘"
        );
        assert_eq!(state.coords_min, &[-7.89, 1.23]);
        assert_eq!(state.coords_max, &[-1.23, 7.89]);
        assert_eq!(state.interp.dim(), nnode);
        assert_eq!(state.deriv.dims(), (nnode, geo_ndim));
        assert_eq!(state.jacobian.dims(), (space_ndim, geo_ndim));
        assert_eq!(state.inv_jacobian.dims(), (space_ndim, space_ndim));
        assert_eq!(state.gradient.dims(), (nnode, space_ndim));
        Ok(())
    }

    #[test]
    fn derive_works() -> Result<(), StrError> {
        let state = StateOfShape::new(2, &[[-1.23, 1.23], [-4.56, 4.56], [-7.89, 7.89]])?;
        let state_clone = state.clone();
        assert_eq!(format!("{:?}", state), "StateOfShape { coords_transp: NumMatrix { nrow: 2, ncol: 3, data: [-1.23, -4.56, -7.89, 1.23, 4.56, 7.89] }, coords_min: [-7.89, 1.23], coords_max: [-1.23, 7.89], interp: NumVector { data: [0.0, 0.0, 0.0] }, deriv: NumMatrix { nrow: 3, ncol: 2, data: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] }, jacobian: NumMatrix { nrow: 2, ncol: 2, data: [0.0, 0.0, 0.0, 0.0] }, inv_jacobian: NumMatrix { nrow: 2, ncol: 2, data: [0.0, 0.0, 0.0, 0.0] }, gradient: NumMatrix { nrow: 3, ncol: 2, data: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] } }");
        assert_eq!(state_clone.coords_min, &[-7.89, 1.23]);
        assert_eq!(state_clone.coords_max, &[-1.23, 7.89]);
        Ok(())
    }
}
