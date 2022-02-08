use super::*;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Defines an alias for integration points data (coordinates and weights)
pub type IpData = &'static [[f64; 4]];

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
    pub ip_data: IpData,

    /// Store a copy of GeoClass to aid in setting integration points
    class: GeoClass,
}

impl ShapeState {
    /// Allocates state variables
    pub fn new(shape: &Shape) -> Self {
        let ip_data: IpData = match shape.kind {
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
            ip_data,
            class: shape.class,
        }
    }

    /// Selects integrations points and weights
    ///
    /// # Options
    ///
    /// ## nip for Lin class
    ///
    /// * `1` -- Legendre points
    /// * `2` -- Legendre points
    /// * `3` -- Legendre points
    /// * `4` -- Legendre points
    /// * `5` -- Legendre points
    ///
    /// ## nip for Tri class
    ///
    /// * `1` -- Internal points
    /// * `3`:
    ///     - `edge` -- Points on edge
    ///     - otherwise, internal points
    /// * `4` -- Internal points
    /// * `12` -- Internal points
    /// * `16` -- Internal points
    ///
    /// ## nip for Qua class
    ///
    /// * `1` -- Legendre points
    /// * `4` -- Legendre points
    /// * `5`:
    ///     - `ws` -- Wilson's "Stable" version with w0=0.004 and wa=0.999 to mimic 4-point rule
    ///     - otherwise, standard Wilson's formula
    /// * `9` -- Legendre points
    /// * `16` -- Legendre points
    ///
    /// ## nip for Tet class
    ///
    /// * `1` -- Internal points
    /// * `4` -- Internal points
    /// * `5` -- Internal points
    /// * `6` -- Internal points
    ///
    /// ## nip for Hex class
    ///
    /// * `6` -- Iron's formula
    /// * `8` -- Legendre points
    /// * `9`:
    ///     - `ws` -- Wilson's "Stable" version
    ///     - otherwise, Wilson's standard formula
    /// * `14` -- Iron's formula
    /// * `27` -- Legendre points
    pub fn select_int_points(&mut self, nip: usize, edge: bool, ws: bool) -> Result<(), StrError> {
        match self.class {
            // Lin
            GeoClass::Lin => match nip {
                1 => self.ip_data = &IP_LIN_LEGENDRE_1,
                2 => self.ip_data = &IP_LIN_LEGENDRE_2,
                3 => self.ip_data = &IP_LIN_LEGENDRE_3,
                4 => self.ip_data = &IP_LIN_LEGENDRE_4,
                5 => self.ip_data = &IP_LIN_LEGENDRE_5,
                _ => return Err("number of integration points is not available for Lin class"),
            },
            // Tri
            GeoClass::Tri => match nip {
                1 => self.ip_data = &IP_TRI_INTERNAL_1,
                3 => {
                    if edge {
                        self.ip_data = &IP_TRI_EDGE_3
                    } else {
                        self.ip_data = &IP_TRI_INTERNAL_3
                    }
                }
                4 => self.ip_data = &IP_TRI_INTERNAL_4,
                12 => self.ip_data = &IP_TRI_INTERNAL_12,
                16 => self.ip_data = &IP_TRI_INTERNAL_16,
                _ => return Err("number of integration points is not available for Tri class"),
            },
            // Qua
            GeoClass::Qua => match nip {
                1 => self.ip_data = &IP_QUA_LEGENDRE_1,
                4 => self.ip_data = &IP_QUA_LEGENDRE_4,
                5 => {
                    if ws {
                        self.ip_data = &IP_QUA_WILSON_STABLE_5
                    } else {
                        self.ip_data = &IP_QUA_WILSON_CORNER_5
                    }
                }
                8 => self.ip_data = &IP_QUA_WILSON_8,
                9 => self.ip_data = &IP_QUA_LEGENDRE_9,
                16 => self.ip_data = &IP_QUA_LEGENDRE_16,
                _ => return Err("number of integration points is not available for Qua class"),
            },
            // Tet
            GeoClass::Tet => match nip {
                1 => self.ip_data = &IP_TET_INTERNAL_1,
                4 => self.ip_data = &IP_TET_INTERNAL_4,
                5 => self.ip_data = &IP_TET_INTERNAL_5,
                6 => self.ip_data = &IP_TET_INTERNAL_6,
                _ => return Err("number of integration points is not available for Tet class"),
            },
            // Hex
            GeoClass::Hex => match nip {
                6 => self.ip_data = &IP_HEX_IRONS_6,
                8 => self.ip_data = &IP_HEX_LEGENDRE_8,
                9 => {
                    if ws {
                        self.ip_data = &IP_HEX_WILSON_STABLE_9
                    } else {
                        self.ip_data = &IP_HEX_WILSON_CORNER_9
                    }
                }
                14 => self.ip_data = &IP_HEX_IRONS_14,
                27 => self.ip_data = &IP_HEX_LEGENDRE_27,
                _ => return Err("number of integration points is not available for Hex class"),
            },
        }
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
        assert_eq!(state.ip_data.len(), 2);
        Ok(())
    }

    #[test]
    fn select_int_points_works() -> Result<(), StrError> {
        // Lin
        let shape = Shape::new(1, 1, 2)?;
        let mut state = ShapeState::new(&shape);
        for nip in 1..6 {
            state.select_int_points(nip, false, false)?;
            assert_eq!(state.ip_data.len(), nip);
        }
        assert_eq!(
            state.select_int_points(100, false, false).err(),
            Some("number of integration points is not available for Lin class")
        );

        // Tri
        let shape = Shape::new(2, 2, 3)?;
        let mut state = ShapeState::new(&shape);
        for nip in [1, 3, 4, 12, 16] {
            state.select_int_points(nip, false, false)?;
            assert_eq!(state.ip_data.len(), nip);
        }
        state.select_int_points(3, true, false)?;
        assert_eq!(state.ip_data.len(), 3);
        assert_eq!(
            state.select_int_points(100, false, false).err(),
            Some("number of integration points is not available for Tri class")
        );

        // Qua
        let shape = Shape::new(2, 2, 4)?;
        let mut state = ShapeState::new(&shape);
        for nip in [1, 4, 5, 8, 9, 16] {
            state.select_int_points(nip, false, false)?;
            assert_eq!(state.ip_data.len(), nip);
        }
        state.select_int_points(5, false, true)?;
        assert_eq!(state.ip_data.len(), 5);
        assert_eq!(
            state.select_int_points(100, false, false).err(),
            Some("number of integration points is not available for Qua class")
        );

        // Tet
        let shape = Shape::new(3, 3, 4)?;
        let mut state = ShapeState::new(&shape);
        for nip in [1, 4, 5, 6] {
            state.select_int_points(nip, false, false)?;
            assert_eq!(state.ip_data.len(), nip);
        }
        assert_eq!(
            state.select_int_points(100, false, false).err(),
            Some("number of integration points is not available for Tet class")
        );

        // Hex
        let shape = Shape::new(3, 3, 8)?;
        let mut state = ShapeState::new(&shape);
        for nip in [6, 8, 9, 14, 27] {
            state.select_int_points(nip, false, false)?;
            assert_eq!(state.ip_data.len(), nip);
        }
        state.select_int_points(9, false, true)?;
        assert_eq!(state.ip_data.len(), 9);
        assert_eq!(
            state.select_int_points(100, false, false).err(),
            Some("number of integration points is not available for Hex class")
        );
        Ok(())
    }
}
