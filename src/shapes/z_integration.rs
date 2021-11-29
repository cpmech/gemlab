use super::*;
use crate::shapes::{GeoClass, Shape};
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::Tensor2;

/// Defines scalar function (the argument is index of the integration point)
pub type FnScalar = fn(usize) -> f64;

/// Defines vector function (the second argument is index of the integration point)
pub type FnVector = fn(&mut Vector, usize);

/// Defines tensor function (the second argument is index of the integration point)
pub type FnTensor = fn(&mut Tensor2, usize);

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

    /// Performs the Case-A integration
    ///
    /// Case-A, interpolation functions times scalar field:
    ///
    /// ```text
    ///      ⌠    → →     →
    /// aᵐ = │ Nᵐ(x(ξ)) s(x) dΩ
    ///      ⌡
    ///      Ω
    /// ```
    ///
    /// or, for lines in multi-dimensions,
    ///
    /// ```text
    ///      ⌠
    /// aᵐ = │ Nᵐ(ℓ(ξ)) s(ℓ) dℓ
    ///      ⌡
    ///      Ω
    /// ```
    ///
    /// The numerical integration is:
    ///
    /// ```text
    ///     nip-1    →     →       →
    /// aᵐ ≈  Σ   Nᵐ(ιp) s(ιp) |J|(ιp) wp
    ///      p=0
    /// ```
    ///
    /// # Input
    ///
    /// * `fn_s` -- s(x(ξ)) or s(ℓ) scalar function, however written as the
    ///             function of the index of the integration point.
    ///
    /// # Output
    ///
    /// * `a` -- (npoint) aᵐ scalars; result from the integration
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

    /// Performs the Case-B integration
    ///
    /// Case-B, interpolation functions times vector field:
    ///
    /// ```text
    /// →    ⌠    → →   → →
    /// bᵐ = │ Nᵐ(x(ξ)) v(x) dΩ
    ///      ⌡
    ///      Ω
    /// ```
    ///
    /// The numerical integration is:
    ///
    /// ```text
    /// →   nip-1    →   → →       →
    /// bᵐ ≈  Σ   Nᵐ(ιp) v(ιp) |J|(ιp) wp
    ///      p=0
    /// ```
    ///
    /// # Input
    ///
    /// * `fn_v` -- v(x(ξ)) vector function, however written as the
    ///             function of the index of the integration point.
    ///
    /// # Output
    ///
    /// * `b` -- (npoint) bᵐ vectors; result from the integration
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
                for i in 0..self.space_ndim {
                    b[m][i] += self.interp[m] * v[i] * det_jac * weight;
                }
            }
        }
        Ok(())
    }

    /// Performs the Case-C integration
    ///
    /// Case-C, vector dot gradient:
    ///
    /// ```text
    ///      ⌠ → →    →  → →
    /// cᵐ = │ w(x) · Gᵐ(x(ξ)) dΩ
    ///      ⌡
    ///      Ω
    /// ```
    ///
    /// The numerical integration is:
    ///
    /// ```text
    ///     nip-1 → →     →  →       →
    /// cᵐ ≈  Σ   w(ιp) · Gᵐ(ιp) |J|(ιp) wp
    ///      p=0
    /// ```
    ///
    /// # Input
    ///
    /// * `fn_w` -- w(x(ξ)) vector function, however written as the
    ///             function of the index of the integration point.
    ///
    /// # Output
    ///
    /// * `c` -- (npoint) cᵐ scalars; result from the integration
    pub fn integ_case_c(&mut self, c: &mut [f64], fn_w: FnVector) -> Result<(), StrError> {
        // clear results
        c.fill(0.0);

        // auxiliary vector holding the output of fn_w
        let mut w = Vector::new(self.space_ndim);

        // loop over integration points
        for index in 0..self.ip_data.len() {
            // ksi coordinates and weight
            let iota = &self.ip_data[index][0..self.geo_ndim];
            let weight = self.ip_data[index][3];

            // calculate Jacobian and Gradient
            let det_jac = self.calc_gradient(iota)?;

            // loop over points and perform summation
            for m in 0..self.npoint {
                fn_w(&mut w, index);
                let w_dot_grad = self.vec_dot_grad(m, &w);
                c[m] += w_dot_grad * det_jac * weight;
            }
        }
        Ok(())
    }

    /// Performs the Case-D integration
    ///
    /// Case-D, tensor dot gradient:
    ///
    /// ```text
    /// →    ⌠   →    →  → →   →
    /// dᵐ = │ σ(x) · Gᵐ(x(ξ)) dΩ
    ///      ⌡ ▔
    ///      Ω
    /// ```
    ///
    /// The numerical integration is:
    ///
    /// ```text
    /// →   nip-1   →     →  →       →
    /// dᵐ ≈  Σ   σ(ιp) · Gᵐ(ιp) |J|(ιp) wp
    ///      p=0  ▔
    /// ```
    ///
    /// # Input
    ///
    /// * `fn_sig` -- sig(x(ξ)) tensor function, however written as the
    ///               function of the index of the integration point.
    ///
    /// # Output
    ///
    /// * `d` -- (npoint) dᵐ vectors; result from the integration
    pub fn integ_case_d(&mut self, d: &mut Vec<Vector>, fn_sig: FnTensor) -> Result<(), StrError> {
        // check
        if d.len() != self.npoint {
            return Err("b.len() must equal npoint");
        }

        // check dims and clear results
        for m in 0..self.npoint {
            if d[m].dim() != self.space_ndim {
                return Err("d[m].dim() must equal space_ndim");
            }
            d[m].fill(0.0);
        }

        // auxiliary tensor holding the output of fn_sig
        let mut sig = Tensor2::new(true);

        // auxiliary vector equal to σ · G
        let mut sig_dot_grad = vec![0.0; self.space_ndim];

        // loop over integration points
        for index in 0..self.ip_data.len() {
            // ksi coordinates and weight
            let iota = &self.ip_data[index][0..self.geo_ndim];
            let weight = self.ip_data[index][3];

            // calculate Jacobian and Gradient
            let det_jac = self.calc_gradient(iota)?;

            // loop over points and perform summation
            for m in 0..self.npoint {
                fn_sig(&mut sig, index);
                self.tensor_dot_grad(&mut sig_dot_grad, m, &sig);
                for i in 0..self.space_ndim {
                    d[m][i] += sig_dot_grad[i] * det_jac * weight;
                }
            }
        }
        Ok(())
    }

    /// Computes vector dot the gradient at point m
    fn vec_dot_grad(&self, m: usize, w: &Vector) -> f64 {
        let mut res = 0.0;
        for i in 0..self.space_ndim {
            res += w[i] * self.gradient[m][i];
        }
        res
    }

    /// Computes tensor dot the gradient at point m
    fn tensor_dot_grad(&self, res: &mut [f64], m: usize, sig: &Tensor2) {
        for i in 0..self.space_ndim {
            res[i] = 0.0;
            for j in 0..self.space_ndim {
                res[i] += sig.get(i, j) * self.gradient[m][j];
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::SQRT_3;
    use russell_chk::assert_vec_approx_eq;

    // to test if variables are cleared before summation
    const NOISE: f64 = 1234.56;

    fn gen_eq_triangle() -> (Shape, f64) {
        // equilateral triangle with sides equal to l
        //       /\
        //      /  \
        //   l /    \ l
        //    /      \
        //   /________\
        //        l
        let l = 5.0;
        let h = l * SQRT_3 / 2.0;
        let area = l * h / 2.0;
        let mut shape = Shape::new(2, 2, 3).unwrap();
        let (xmin, ymin) = (3.0, 4.0);
        shape.set_point(0, 0, xmin).unwrap();
        shape.set_point(0, 1, ymin).unwrap();
        shape.set_point(1, 0, xmin + l).unwrap();
        shape.set_point(1, 1, ymin).unwrap();
        shape.set_point(2, 0, xmin + l / 2.0).unwrap();
        shape.set_point(2, 1, ymin + h).unwrap();
        (shape, area)
    }

    #[test]
    fn integ_case_a_works() -> Result<(), StrError> {
        // with a constant source term:
        //
        // s(x) = cₛ
        //
        // we get:
        //           ┌   ┐
        //      cₛ A │ 1 │
        // Fₛ = ———— │ 1 │
        //        3  │ 1 │
        //           └   ┘
        let (mut shape, area) = gen_eq_triangle();
        let mut a = vec![NOISE; shape.npoint];
        const CS: f64 = 3.0;
        let fn_s = |_| CS;
        shape.integ_case_a(&mut a, fn_s)?;
        let cf = CS * area / 3.0;
        let a_correct = &[cf, cf, cf];
        assert_vec_approx_eq!(a, a_correct, 1e-14);
        Ok(())
    }
}
