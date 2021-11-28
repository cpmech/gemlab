use super::*;
use crate::shapes::{GeoClass, GeoKind, Shape};
use crate::StrError;

/// Defines options to select integration points
#[derive(Eq, PartialEq)]
pub enum OptionIntPoint {
    /// Option for Tri (points on edge)
    Edge,

    /// Option for Qua and Hex (Wilson formula)
    WilsonStable,
}

// Defines an alias for integration points data
pub type IpData = &'static [[f64; 4]];

pub struct Integrator {
    /// Number of integration points
    pub nip: usize,

    /// Integration points and weights
    pub points: IpData,
    // pub int_points: Vec<Vector>,
}

/// Defines scalar function (index,ksi)
pub type FnScalar = fn(usize, &[f64]) -> f64;

impl Integrator {
    /// Creates new Integrator
    pub fn new(kind: GeoKind) -> Self {
        match kind {
            // Lin
            GeoKind::Lin2 => Integrator {
                nip: 2,
                points: &IP_LIN_LEGENDRE_2,
            },
            GeoKind::Lin3 => Integrator {
                nip: 3,
                points: &IP_LIN_LEGENDRE_3,
            },
            GeoKind::Lin4 => Integrator {
                nip: 4,
                points: &IP_LIN_LEGENDRE_4,
            },
            GeoKind::Lin5 => Integrator {
                nip: 5,
                points: &IP_LIN_LEGENDRE_5,
            },

            // Tri
            GeoKind::Tri3 => Integrator {
                nip: 1,
                points: &IP_TRI_INTERNAL_1,
            },
            GeoKind::Tri6 => Integrator {
                nip: 3,
                points: &IP_TRI_INTERNAL_3,
            },
            GeoKind::Tri10 => Integrator {
                nip: 12,
                points: &IP_TRI_INTERNAL_12,
            },
            GeoKind::Tri15 => Integrator {
                nip: 12,
                points: &IP_TRI_INTERNAL_12,
            },

            // Qua
            GeoKind::Qua4 => Integrator {
                nip: 4,
                points: &IP_QUA_LEGENDRE_4,
            },
            GeoKind::Qua8 => Integrator {
                nip: 9,
                points: &IP_QUA_LEGENDRE_9,
            },
            GeoKind::Qua9 => Integrator {
                nip: 9,
                points: &IP_QUA_LEGENDRE_9,
            },
            GeoKind::Qua12 => Integrator {
                nip: 9,
                points: &IP_QUA_LEGENDRE_9,
            },
            GeoKind::Qua16 => Integrator {
                nip: 16,
                points: &IP_QUA_LEGENDRE_16,
            },
            GeoKind::Qua17 => Integrator {
                nip: 16,
                points: &IP_QUA_LEGENDRE_16,
            },

            // Tet
            GeoKind::Tet4 => Integrator {
                nip: 1,
                points: &IP_TET_INTERNAL_1,
            },
            GeoKind::Tet10 => Integrator {
                nip: 4,
                points: &IP_TET_INTERNAL_4,
            },

            // Hex
            GeoKind::Hex8 => Integrator {
                nip: 8,
                points: &IP_HEX_LEGENDRE_8,
            },
            GeoKind::Hex20 => Integrator {
                nip: 27,
                points: &IP_HEX_LEGENDRE_27,
            },
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
    ///     - `Edge` -- Points on edge
    ///     - Otherwise, internal points
    /// * `4` -- Internal points
    /// * `12` -- Internal points
    /// * `16` -- Internal points
    ///
    /// ## nip for Qua class
    ///
    /// * `1` -- Legendre points
    /// * `4` -- Legendre points
    /// * `5`:
    ///     - `WilsonStable` -- Wilson's "Stable" version with w0=0.004 and wa=0.999 to mimic 4-point rule
    ///     - Otherwise, Wilson's formula
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
    ///     - `WilsonStable` -- Wilson's "Stable" version
    ///     - Otherwise, Wilson's formula
    /// * `14` -- Iron's formula
    /// * `27` -- Legendre points
    pub fn select_int_points(&mut self, class: GeoClass, nip: usize, option: OptionIntPoint) -> Result<(), StrError> {
        match class {
            // Lin
            GeoClass::Lin => match nip {
                1 => self.points = &IP_LIN_LEGENDRE_1,
                2 => self.points = &IP_LIN_LEGENDRE_2,
                3 => self.points = &IP_LIN_LEGENDRE_3,
                4 => self.points = &IP_LIN_LEGENDRE_4,
                5 => self.points = &IP_LIN_LEGENDRE_5,
                _ => return Err("number of integration points is not available for Lin class"),
            },

            // Tri
            GeoClass::Tri => match nip {
                1 => self.points = &IP_TRI_INTERNAL_1,
                3 => {
                    if option == OptionIntPoint::Edge {
                        self.points = &IP_TRI_EDGE_3
                    } else {
                        self.points = &IP_TRI_INTERNAL_3
                    }
                }
                4 => self.points = &IP_TRI_INTERNAL_4,
                12 => self.points = &IP_TRI_INTERNAL_12,
                16 => self.points = &IP_TRI_INTERNAL_16,
                _ => return Err("number of integration points is not available for Tri class"),
            },

            // Qua
            GeoClass::Qua => match nip {
                1 => self.points = &IP_QUA_LEGENDRE_1,
                4 => self.points = &IP_QUA_LEGENDRE_4,
                5 => {
                    if option == OptionIntPoint::WilsonStable {
                        self.points = &IP_QUA_WILSON_STABLE_5
                    } else {
                        self.points = &IP_QUA_WILSON_CORNER_5
                    }
                }
                8 => self.points = &IP_QUA_WILSON_8,
                9 => self.points = &IP_QUA_LEGENDRE_9,
                16 => self.points = &IP_QUA_LEGENDRE_16,
                _ => return Err("number of integration points is not available for Qua class"),
            },

            // Tet
            GeoClass::Tet => match nip {
                1 => self.points = &IP_TET_INTERNAL_1,
                4 => self.points = &IP_TET_INTERNAL_4,
                5 => self.points = &IP_TET_INTERNAL_5,
                6 => self.points = &IP_TET_INTERNAL_6,
                _ => return Err("number of integration points is not available for Tet class"),
            },

            // Hex
            GeoClass::Hex => match nip {
                6 => self.points = &IP_HEX_IRONS_6,
                8 => self.points = &IP_HEX_LEGENDRE_8,
                9 => {
                    if option == OptionIntPoint::WilsonStable {
                        self.points = &IP_HEX_WILSON_STABLE_9
                    } else {
                        self.points = &IP_HEX_WILSON_CORNER_9
                    }
                }
                14 => self.points = &IP_HEX_IRONS_14,
                27 => self.points = &IP_HEX_LEGENDRE_27,
                _ => return Err("number of integration points is not available for Hex class"),
            },
        }
        self.nip = nip;
        Ok(())
    }

    pub fn case_a(&self, a: &mut [f64], shape: &mut Shape, s: FnScalar) -> Result<(), StrError> {
        a.fill(0.0);
        for (index, iota) in self.points.iter().enumerate() {
            let ksi = &iota[0..shape.geo_ndim];
            let weight = iota[3];
            shape.calc_interp(ksi);
            let det_jac = shape.calc_jacobian(ksi)?;
            for m in 0..shape.npoint {
                a[m] += shape.interp[m] * s(index, ksi) * det_jac * weight;
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    const NOISE: f64 = 1234.56;

    #[test]
    fn new_works() -> Result<(), StrError> {
        let shape = Shape::new(2, 2, 4)?;
        let integ = Integrator::new(shape.kind);
        assert_eq!(integ.nip, 4);
        Ok(())
    }

    #[test]
    fn case_a_works() -> Result<(), StrError> {
        let mut shape = Shape::new(2, 2, 4)?;
        shape.set_point(0, 0, 0.0)?;
        shape.set_point(0, 1, 0.0)?;
        shape.set_point(1, 0, 1.0)?;
        shape.set_point(1, 1, 0.0)?;
        shape.set_point(2, 0, 1.0)?;
        shape.set_point(2, 1, 1.0)?;
        shape.set_point(3, 0, 0.0)?;
        shape.set_point(3, 1, 1.0)?;

        let integ = Integrator::new(shape.kind);
        let mut a = vec![NOISE; shape.npoint];
        integ.case_a(&mut a, &mut shape, |_, _| 0.0)?;
        println!("a = {:?}", a);

        Ok(())
    }
}
