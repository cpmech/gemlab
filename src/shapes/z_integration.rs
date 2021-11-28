use super::*;
use crate::shapes::{GeoClass, Shape};
use crate::StrError;
use russell_lab::Vector;

/// Defines scalar function (the argument is index of the integration point)
pub type FnScalar = fn(usize) -> f64;

/// Defines vector function (the second argument is index of the integration point)
pub type FnVector = fn(&mut Vector, usize);

impl Shape {
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

    pub fn integ_case_a(&mut self, a: &mut [f64], fn_s: FnScalar) -> Result<(), StrError> {
        // clear results
        a.fill(0.0);

        // loop over integration points
        for index in 0..self.ip_data.len() {
            // ksi coordinates and weight
            let iota = &self.ip_data[index][0..self.geo_ndim];
            let weight = self.ip_data[index][3];

            // calculate interpolation functions and Jacobian
            self.calc_interp(iota);
            let det_jac = self.calc_jacobian(iota)?;

            // loop over points and perform summation
            for m in 0..self.npoint {
                a[m] += self.interp[m] * fn_s(index) * det_jac * weight;
            }
        }
        Ok(())
    }

    pub fn integ_case_b(&mut self, b: &mut Vec<Vector>, fn_v: FnVector) -> Result<(), StrError> {
        // check
        if b.len() != self.npoint {
            return Err("b.len() must equal npoint");
        }

        // check dims and clear results
        for m in 0..self.npoint {
            if b[m].dim() != self.space_ndim {
                return Err("b[m].dim() must equal space_ndim");
            }
            b[m].fill(0.0);
        }

        // auxiliary vector holding the output of fn_v
        let mut v = Vector::new(self.space_ndim);

        // loop over integration points
        for index in 0..self.ip_data.len() {
            // ksi coordinates and weight
            let iota = &self.ip_data[index][0..self.geo_ndim];
            let weight = self.ip_data[index][3];

            // calculate interpolation functions and Jacobian
            self.calc_interp(iota);
            let det_jac = self.calc_jacobian(iota)?;

            // loop over points and perform summation
            for m in 0..self.npoint {
                fn_v(&mut v, index);
                for i in 0..b[m].dim() {
                    b[m][i] += self.interp[m] * v[i] * det_jac * weight;
                }
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    // to test if variables are cleared before summation
    const NOISE: f64 = 1234.56;

    #[test]
    fn integ_case_a_works() -> Result<(), StrError> {
        let mut shape = Shape::new(2, 2, 4)?;
        shape.set_point(0, 0, 0.0)?;
        shape.set_point(0, 1, 0.0)?;
        shape.set_point(1, 0, 1.0)?;
        shape.set_point(1, 1, 0.0)?;
        shape.set_point(2, 0, 1.0)?;
        shape.set_point(2, 1, 1.0)?;
        shape.set_point(3, 0, 0.0)?;
        shape.set_point(3, 1, 1.0)?;

        let mut a = vec![NOISE; shape.npoint];

        shape.integ_case_a(&mut a, |_| 0.0)?;

        println!("a = {:?}", a);

        Ok(())
    }
}
