use super::geo_class_and_kind;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Holds the state of a shape entity
#[derive(Clone, Debug)]
pub struct StateOfShape {
    // -- steady ----------------------------------------------------------------------------------------------------
    //
    /// Maps local node index to global point id (aka "connectivity"); corresponds to coords_transp (nnode)
    pub node_to_point: Vec<usize>,

    /// Matrix Xᵀ: (space_ndim,nnode) transposed coordinates matrix (real space)
    ///
    /// ```text
    /// The (non-transposed) matrix X is such that:
    ///
    ///     ┌               ┐
    ///     | x⁰₀  x⁰₁  x⁰₂ |
    ///     | x¹₀  x¹₁  x¹₂ |
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

    // -- internal --------------------------------------------------------------------------------------------------
    //
    /// Tells whether the last entry of the coordinates matrix has been given or not
    pub(super) last_coord_given: bool,
}

impl StateOfShape {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `space_ndim` -- the space dimension (1, 2, or 3)
    /// * `geo_ndim` -- the dimension of the shape; e.g. a 1D line in a 3D space
    /// * `nnode` -- the number of points defining the shape (number of nodes)
    pub fn new(space_ndim: usize, geo_ndim: usize, nnode: usize) -> Result<Self, StrError> {
        geo_class_and_kind(geo_ndim, nnode)?; // just to check combination
        Ok(StateOfShape {
            node_to_point: vec![0; nnode],
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
            last_coord_given: false,
        })
    }

    /// Sets a component of the coordinates matrix
    ///
    /// ```text
    /// The (non-transposed) matrix X is such that:
    ///
    ///     ┌               ┐
    ///     | x⁰₀  x⁰₁  x⁰₂ |
    ///     | x¹₀  x¹₁  x¹₂ |
    /// X = | x²₀  x²₁  x²₂ |
    ///     |      ...      |
    ///     | xᴹ₀  xᴹ₁  xᴹ₂ |
    ///     └               ┘_(nnode,space_ndim)
    ///
    /// where `M = nnode - 1`
    /// ```
    ///
    /// # Input
    ///
    /// * `point_id` -- the global point id to configure the `node_to_point` map; aka "connectivity"
    /// * `m` -- node index in `0 ≤ m ≤ nnode - 1`
    /// * `j` -- dimension index in `0 ≤ j ≤ space_ndim - 1`
    /// * `value` -- the X(m,j) component
    ///
    /// # Updated variables
    ///
    /// * `node_to_point` -- (nnode) array of point ids that works as a map from "local node index" to "global point id"
    /// * `coords_transp` -- Matrix Xᵀ: (space_ndim,nnode) transposed coordinates matrix (real space)
    /// * `coords_min` -- minimum (space_ndim) coordinates from the X matrix
    /// * `coords_max` -- maximum (space_ndim) coordinates from the X matrix
    pub fn set_node(&mut self, point_id: usize, m: usize, j: usize, value: f64) -> Result<(), StrError> {
        let (space_ndim, nnode) = self.coords_transp.dims();
        if m >= nnode {
            return Err("index of node must be in 0 ≤ m ≤ nnode - 1");
        }
        if j >= space_ndim {
            return Err("index of space dimension must be in 0 ≤ j ≤ space_ndim - 1");
        }
        self.node_to_point[m] = point_id;
        self.coords_transp[j][m] = value;
        if value < self.coords_min[j] {
            self.coords_min[j] = value;
        }
        if value > self.coords_max[j] {
            self.coords_max[j] = value;
        }
        if m == nnode - 1 && j == space_ndim - 1 {
            self.last_coord_given = true;
        } else {
            self.last_coord_given = false;
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::shapes::StateOfShape;
    use crate::StrError;

    #[test]
    fn set_node_captures_wrong_input() -> Result<(), StrError> {
        let mut state = StateOfShape::new(1, 1, 2)?;
        let point_id = 123;
        assert_eq!(
            state.set_node(point_id, 3, 0, 0.0).err(),
            Some("index of node must be in 0 ≤ m ≤ nnode - 1")
        );
        assert_eq!(
            state.set_node(point_id, 0, 2, 0.0).err(),
            Some("index of space dimension must be in 0 ≤ j ≤ space_ndim - 1")
        );
        Ok(())
    }

    #[test]
    fn set_node_resets_last_coord_given() -> Result<(), StrError> {
        let point_id = 123;
        let mut state = StateOfShape::new(1, 1, 2)?;
        assert_eq!(state.last_coord_given, false);
        state.set_node(point_id, 0, 0, 0.0)?;
        assert_eq!(state.last_coord_given, false);
        state.set_node(point_id, 1, 0, 0.0)?;
        assert_eq!(state.last_coord_given, true);
        state.set_node(point_id, 0, 0, 0.0)?; // must reset
        assert_eq!(state.last_coord_given, false);
        Ok(())
    }

    #[test]
    fn set_node_works() -> Result<(), StrError> {
        let mut state = StateOfShape::new(2, 2, 3)?;
        state.set_node(123, 0, 0, -1.23)?;
        state.set_node(123, 0, 1, 1.23)?;
        state.set_node(456, 1, 0, -4.56)?;
        state.set_node(456, 1, 1, 4.56)?;
        state.set_node(789, 2, 0, -7.89)?;
        state.set_node(789, 2, 1, 7.89)?;
        assert_eq!(state.node_to_point, &[123, 456, 789]);
        assert_eq!(
            format!("{}", state.coords_transp),
            "┌                   ┐\n\
             │ -1.23 -4.56 -7.89 │\n\
             │  1.23  4.56  7.89 │\n\
             └                   ┘"
        );
        assert_eq!(state.last_coord_given, true);
        Ok(())
    }
}
