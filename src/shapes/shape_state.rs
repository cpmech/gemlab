use super::*;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Defines an alias for integration points' constants (coordinates and weights)
pub type IntegPointConstants = &'static [[f64; 4]];

// Stores state variables for Shape
pub struct ShapeState {
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

    /// Integration points data (coordinates and weights)
    pub integ_point_constants: IntegPointConstants,

    /// Store a copy of GeoClass to aid in setting integration points
    class: GeoClass,
}

impl ShapeState {
    /// Allocates state variables
    pub fn new(shape: &Shape) -> Self {
        let integ_point_constants: IntegPointConstants = match shape.kind {
            // Lin
            GeoKind::Lin2 => &IP_LIN_LEGENDRE_2,
            GeoKind::Lin3 => &IP_LIN_LEGENDRE_3,
            GeoKind::Lin4 => &IP_LIN_LEGENDRE_4,
            GeoKind::Lin5 => &IP_LIN_LEGENDRE_5,
            // Tri
            GeoKind::Tri3 => &IP_TRI_INTERNAL_1,
            GeoKind::Tri6 => &IP_TRI_INTERNAL_3,
            GeoKind::Tri10 => &IP_TRI_INTERNAL_12,
            GeoKind::Tri15 => &IP_TRI_INTERNAL_12,
            // Qua
            GeoKind::Qua4 => &IP_QUA_LEGENDRE_4,
            GeoKind::Qua8 => &IP_QUA_LEGENDRE_9,
            GeoKind::Qua9 => &IP_QUA_LEGENDRE_9,
            GeoKind::Qua12 => &IP_QUA_LEGENDRE_9,
            GeoKind::Qua16 => &IP_QUA_LEGENDRE_16,
            GeoKind::Qua17 => &IP_QUA_LEGENDRE_16,
            // Tet
            GeoKind::Tet4 => &IP_TET_INTERNAL_1,
            GeoKind::Tet10 => &IP_TET_INTERNAL_4,
            // Hex
            GeoKind::Hex8 => &IP_HEX_LEGENDRE_8,
            GeoKind::Hex20 => &IP_HEX_LEGENDRE_27,
        };
        ShapeState {
            interp: Vector::new(shape.nnode),
            deriv: Matrix::new(shape.nnode, shape.geo_ndim),
            jacobian: Matrix::new(shape.space_ndim, shape.geo_ndim),
            inv_jacobian: if shape.geo_ndim == shape.space_ndim {
                Matrix::new(shape.space_ndim, shape.space_ndim)
            } else {
                Matrix::new(0, 0)
            },
            gradient: if shape.geo_ndim == shape.space_ndim {
                Matrix::new(shape.nnode, shape.space_ndim)
            } else {
                Matrix::new(0, 0)
            },
            integ_point_constants,
            class: shape.class,
        }
    }

    /// Selects integrations points and weights
    ///
    /// # Options
    ///
    /// ## n_integ_point for Lin class
    ///
    /// * `1` -- Conventional Legendre integration points and weights
    /// * `2` -- Conventional Legendre integration points and weights
    /// * `3` -- Conventional Legendre integration points and weights
    /// * `4` -- Conventional Legendre integration points and weights
    /// * `5` -- Conventional Legendre integration points and weights
    ///
    /// ## n_integ_point for Tri class
    ///
    /// * `1` -- Internal integration points and weights
    /// * `3` -- Internal integration points and weights
    /// * `1_003` -- Edge integration points and weights
    /// * `4` -- Internal integration points and weights
    /// * `12` -- Internal integration points and weights
    /// * `16` -- Internal integration points and weights
    ///
    /// ## n_integ_point for Qua class
    ///
    /// * `1` -- Conventional Legendre integration points and weights
    /// * `4` -- Conventional Legendre integration points and weights
    /// * `5` -- Wilson's integration points and weights. "Corner" version
    /// * `1_005` -- 5 points. Wilson's integration points and weights. "Stable" version version with w0=0.004 and wa=0.999 to mimic 4-point rule
    /// * `8` -- Wilson's integration points and weights.
    /// * `9` -- Conventional Legendre integration points and weights
    /// * `16` -- Conventional Legendre integration points and weights
    ///
    /// ## n_integ_point for Tet class
    ///
    /// * `1` -- Internal integration points and weights
    /// * `4` -- Internal integration points and weights
    /// * `5` -- Internal integration points and weights
    /// * `6` -- Internal integration points and weights
    ///
    /// ## n_integ_point for Hex class
    ///
    /// * `6` -- Iron's integration points and weights
    /// * `8` -- Conventional Legendre integration points and weights
    /// * `9` -- Wilson's integration points and weights. "Corner" version
    /// * `1_009` -- Wilson's integration points and weights. "Stable" version
    /// * `14` -- Iron's integration points and weights
    /// * `27` -- Conventional Legendre integration points and weights
    pub fn select_int_points(&mut self, n_integ_point: usize) -> Result<(), StrError> {
        self.integ_point_constants = match self.class {
            // Lin
            GeoClass::Lin => match n_integ_point {
                1 => &IP_LIN_LEGENDRE_1,
                2 => &IP_LIN_LEGENDRE_2,
                3 => &IP_LIN_LEGENDRE_3,
                4 => &IP_LIN_LEGENDRE_4,
                5 => &IP_LIN_LEGENDRE_5,
                _ => return Err("number of integration points is not available for Lin class"),
            },
            // Tri
            GeoClass::Tri => match n_integ_point {
                1 => &IP_TRI_INTERNAL_1,
                3 => &IP_TRI_INTERNAL_3,
                1_003 => &IP_TRI_EDGE_3,
                4 => &IP_TRI_INTERNAL_4,
                12 => &IP_TRI_INTERNAL_12,
                16 => &IP_TRI_INTERNAL_16,
                _ => return Err("number of integration points is not available for Tri class"),
            },
            // Qua
            GeoClass::Qua => match n_integ_point {
                1 => &IP_QUA_LEGENDRE_1,
                4 => &IP_QUA_LEGENDRE_4,
                5 => &IP_QUA_WILSON_CORNER_5,
                1_005 => &IP_QUA_WILSON_STABLE_5,
                8 => &IP_QUA_WILSON_8,
                9 => &IP_QUA_LEGENDRE_9,
                16 => &IP_QUA_LEGENDRE_16,
                _ => return Err("number of integration points is not available for Qua class"),
            },
            // Tet
            GeoClass::Tet => match n_integ_point {
                1 => &IP_TET_INTERNAL_1,
                4 => &IP_TET_INTERNAL_4,
                5 => &IP_TET_INTERNAL_5,
                6 => &IP_TET_INTERNAL_6,
                _ => return Err("number of integration points is not available for Tet class"),
            },
            // Hex
            GeoClass::Hex => match n_integ_point {
                6 => &IP_HEX_IRONS_6,
                8 => &IP_HEX_LEGENDRE_8,
                9 => &IP_HEX_WILSON_CORNER_9,
                1_009 => &IP_HEX_WILSON_STABLE_9,
                14 => &IP_HEX_IRONS_14,
                27 => &IP_HEX_LEGENDRE_27,
                _ => return Err("number of integration points is not available for Hex class"),
            },
        };
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::shapes::{Shape, ShapeState};
    use crate::StrError;

    #[test]
    fn new_works() -> Result<(), StrError> {
        let (space_ndim, geo_ndim, nnode) = (1, 1, 2);
        let shape = Shape::new(space_ndim, geo_ndim, nnode)?;
        let state = ShapeState::new(&shape);
        assert_eq!(state.interp.dim(), nnode);
        assert_eq!(state.deriv.dims(), (nnode, space_ndim));
        assert_eq!(state.jacobian.dims(), (space_ndim, geo_ndim));
        assert_eq!(state.inv_jacobian.dims(), (space_ndim, space_ndim));
        assert_eq!(state.gradient.dims(), (nnode, space_ndim));
        assert_eq!(state.integ_point_constants.len(), 2);
        Ok(())
    }

    #[test]
    fn select_int_points_works() -> Result<(), StrError> {
        // Lin
        let shape = Shape::new(1, 1, 2)?;
        let mut state = ShapeState::new(&shape);
        for n_integ_point in [1, 2, 3, 4, 5] {
            state.select_int_points(n_integ_point)?;
            assert_eq!(state.integ_point_constants.len(), n_integ_point);
        }
        assert_eq!(
            state.select_int_points(100).err(),
            Some("number of integration points is not available for Lin class")
        );

        // Tri
        let shape = Shape::new(2, 2, 3)?;
        let mut state = ShapeState::new(&shape);
        for n_integ_point in [1, 3, 4, 12, 16] {
            state.select_int_points(n_integ_point)?;
            assert_eq!(state.integ_point_constants.len(), n_integ_point);
        }
        state.select_int_points(1_003)?;
        assert_eq!(state.integ_point_constants.len(), 3);
        assert_eq!(
            state.select_int_points(100).err(),
            Some("number of integration points is not available for Tri class")
        );

        // Qua
        let shape = Shape::new(2, 2, 4)?;
        let mut state = ShapeState::new(&shape);
        for n_integ_point in [1, 4, 5, 8, 9, 16] {
            state.select_int_points(n_integ_point)?;
            assert_eq!(state.integ_point_constants.len(), n_integ_point);
        }
        state.select_int_points(1_005)?;
        assert_eq!(state.integ_point_constants.len(), 5);
        assert_eq!(
            state.select_int_points(100).err(),
            Some("number of integration points is not available for Qua class")
        );

        // Tet
        let shape = Shape::new(3, 3, 4)?;
        let mut state = ShapeState::new(&shape);
        for n_integ_point in [1, 4, 5, 6] {
            state.select_int_points(n_integ_point)?;
            assert_eq!(state.integ_point_constants.len(), n_integ_point);
        }
        assert_eq!(
            state.select_int_points(100).err(),
            Some("number of integration points is not available for Tet class")
        );

        // Hex
        let shape = Shape::new(3, 3, 8)?;
        let mut state = ShapeState::new(&shape);
        for n_integ_point in [6, 8, 9, 14, 27] {
            state.select_int_points(n_integ_point)?;
            assert_eq!(state.integ_point_constants.len(), n_integ_point);
        }
        state.select_int_points(1_009)?;
        assert_eq!(state.integ_point_constants.len(), 9);
        assert_eq!(
            state.select_int_points(100).err(),
            Some("number of integration points is not available for Hex class")
        );
        Ok(())
    }
}
