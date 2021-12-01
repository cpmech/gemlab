use super::*;
use crate::shapes::{GeoClass, Shape};
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::Tensor2;

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
    /// aᵐ ≈  Σ   Nᵐ(ιᵖ) s(ιᵖ) |J|(ιᵖ) wᵖ
    ///      p=0
    /// ```
    ///
    /// # Output
    ///
    /// ```text
    ///     ┌     ┐
    ///     |  a⁰ |
    ///     |  a¹ |
    /// a = |  a² |
    ///     | ··· |
    ///     |  aᵐ |
    ///     └     ┘
    /// ```
    ///
    /// * `a` -- A vector containing all `aᵐ` values, one after another, and
    ///          sequentially placed as shown above. `m` is the index of the point.
    ///          The length of `a` must be equal to `npoint`.
    ///
    /// # Input
    ///
    /// * `fn_s(index: usize) -> f64` -- s(x(ξ)) or s(ℓ) scalar function, however written as
    ///                                  a function of the index of the integration point.
    pub fn integ_case_a<F>(&mut self, a: &mut Vector, fn_s: F) -> Result<(), StrError>
    where
        F: Fn(usize) -> f64,
    {
        // check
        if a.dim() != self.npoint {
            return Err("the length of vector 'a' must be equal to npoint");
        }

        // clear output vector
        a.fill(0.0);

        // loop over integration points
        for index in 0..self.ip_data.len() {
            // ksi coordinates and weight
            let iota = &self.ip_data[index];
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
    /// bᵐ ≈  Σ   Nᵐ(ιᵖ) v(ιᵖ) |J|(ιᵖ) wᵖ
    ///      p=0
    /// ```
    ///
    /// # Output
    ///
    /// ```text
    ///     ┌     ┐
    ///     | b⁰₀ |
    ///     | b⁰₁ |
    ///     | b¹₀ |
    /// b = | b¹₁ |
    ///     | b²₀ |
    ///     | b²₁ |
    ///     | ··· |
    ///     | bᵐᵢ |
    ///     └     ┘
    /// ```
    ///
    /// * `b` -- A vector containing all `bᵐᵢ` values, one after another, and sequentially placed
    ///          as shown above (in 2D). `m` is the index of the point and `i` corresponds to `space_ndim`.
    ///          The length of `b` must equal to `npoint * space_ndim`.
    ///
    /// # Input
    ///
    /// * `fn_v(v: &mut Vector, index: usize)` -- v(x(ξ)) vector function with `v.dim() == space_ndim`, however written as
    ///                                           a function of the index of the integration point.
    /// * `aux_v` -- is an auxiliary vector with size equal to `space_ndim`.
    pub fn integ_case_b<F>(&mut self, b: &mut Vector, fn_v: F, aux_v: &mut Vector) -> Result<(), StrError>
    where
        F: Fn(&mut Vector, usize),
    {
        // check
        if b.dim() != self.npoint * self.space_ndim {
            return Err("the length of vector 'b' must be equal to npoint * space_ndim");
        }
        if aux_v.dim() != self.space_ndim {
            return Err("the length of vector 'aux_v' must be equal to space_ndim");
        }

        // clear output vector
        b.fill(0.0);

        // loop over integration points
        for index in 0..self.ip_data.len() {
            // ksi coordinates and weight
            let iota = &self.ip_data[index];
            let weight = self.ip_data[index][3];

            // calculate interpolation functions and Jacobian
            self.calc_interp(iota);
            let det_jac = self.calc_jacobian(iota)?;

            // loop over points and perform summation
            for m in 0..self.npoint {
                fn_v(aux_v, index);
                for i in 0..self.space_ndim {
                    let ii = i + m * self.space_ndim;
                    b[ii] += self.interp[m] * aux_v[i] * det_jac * weight;
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
    /// cᵐ ≈  Σ   w(ιᵖ) · Gᵐ(ιᵖ) |J|(ιᵖ) wᵖ
    ///      p=0
    /// ```
    ///
    /// # Input
    ///
    /// * `fn_w(w: &mut Vector, index: usize)` -- w(x(ξ)) vector function, however written as
    ///                                           a function of the index of the integration point.
    ///
    /// # Output
    ///
    /// * `c` -- (npoint) cᵐ scalars; result from the integration
    pub fn integ_case_c<F>(&mut self, c: &mut [f64], fn_w: F) -> Result<(), StrError>
    where
        F: Fn(&mut Vector, usize),
    {
        // clear results
        c.fill(0.0);

        // auxiliary vector holding the output of fn_w
        let mut w = Vector::new(self.space_ndim);

        // loop over integration points
        for index in 0..self.ip_data.len() {
            // ksi coordinates and weight
            let iota = &self.ip_data[index];
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
    /// dᵐ ≈  Σ   σ(ιᵖ) · Gᵐ(ιᵖ) |J|(ιᵖ) wᵖ
    ///      p=0  ▔
    /// ```
    ///
    /// # Input
    ///
    /// * `fn_sig(sig: &mut Tensor2, index: usize)` -- σ(x(ξ)) tensor function, however written as
    ///                                                a function of the index of the integration point.
    ///
    /// # Output
    ///
    /// * `d` -- (npoint) dᵐ vectors; result from the integration
    pub fn integ_case_d<F>(&mut self, d: &mut Vec<Vector>, fn_sig: F) -> Result<(), StrError>
    where
        F: Fn(&mut Tensor2, usize),
    {
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
            let iota = &self.ip_data[index];
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

    // equilateral triangle with sides equal to l
    //       /\
    //      /  \
    //   l /    \ l
    //    /      \
    //   /________\
    //        l
    fn gen_tri3() -> (Shape, f64) {
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

    // line segment from xa to xb (geo_ndim == space_ndim)
    fn gen_lin2() -> (Shape, f64, f64) {
        let mut shape = Shape::new(1, 1, 2).unwrap();
        let (xa, xb) = (3.0, 9.0);
        shape.set_point(0, 0, xa).unwrap();
        shape.set_point(1, 0, xb).unwrap();
        (shape, xa, xb)
    }

    #[test]
    fn integ_case_a_works() -> Result<(), StrError> {
        // tri3 with a constant source term:
        //
        // s(x) = cₛ
        //
        // we get:
        //           ┌   ┐
        //      cₛ A │ 1 │
        // Fₛ = ———— │ 1 │
        //        3  │ 1 │
        //           └   ┘
        let (mut tri3, area) = gen_tri3();
        const CS: f64 = 3.0;
        let fn_s = |_| CS;
        let mut a = Vector::filled(tri3.npoint, NOISE);
        tri3.integ_case_a(&mut a, fn_s)?;
        let cf = CS * area / 3.0;
        let a_correct = &[cf, cf, cf];
        assert_vec_approx_eq!(a.as_data(), a_correct, 1e-14);

        // lin2 with linear source term:
        //
        // s(x) = x
        //
        //        ┌           ┐
        //      l │ 2 xa + xb │
        // Fₛ = — │           │
        //      6 │ xa + 2 xb │
        //        └           ┘
        let (mut lin2, xa, xb) = gen_lin2();
        let all_int_points = lin2.calc_int_points_coords()?;
        let fn_s = |index: usize| all_int_points[index][0];
        let mut a = Vector::new(lin2.npoint);
        lin2.integ_case_a(&mut a, fn_s)?;
        let cf = (xb - xa) / 6.0;
        let a_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];
        assert_vec_approx_eq!(a.as_data(), a_correct, 1e-15);
        Ok(())
    }

    #[test]
    fn integ_case_b_works() -> Result<(), StrError> {
        // This test is similar to the case_a with tri3, however using a vector
        // So, each component of `b` equals `Fₛ`
        let (mut tri3, area) = gen_tri3();
        const CS: f64 = 3.0;
        let fn_v = |v: &mut Vector, _: usize| v.fill(CS);
        let mut b = Vector::filled(tri3.npoint * tri3.space_ndim, NOISE);
        let mut aux_v = Vector::new(tri3.space_ndim);
        tri3.integ_case_b(&mut b, fn_v, &mut aux_v)?;
        let cf = CS * area / 3.0;
        let b_correct = &[cf, cf, cf, cf, cf, cf];
        assert_vec_approx_eq!(b.as_data(), b_correct, 1e-14);

        // Likewise, this test is similar to case_a with lin2, however using a vector
        // with a single component. So, each component of `b` equals `Fₛ`
        let (mut lin2, xa, xb) = gen_lin2();
        let all_int_points = lin2.calc_int_points_coords()?;
        let fn_v = |v: &mut Vector, index: usize| {
            v.fill(all_int_points[index][0]);
        };
        let mut b = Vector::filled(lin2.npoint * lin2.space_ndim, NOISE);
        let mut aux_v = Vector::new(lin2.space_ndim);
        lin2.integ_case_b(&mut b, fn_v, &mut aux_v)?;
        let cf = (xb - xa) / 6.0;
        let b_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];
        assert_vec_approx_eq!(b.as_data(), b_correct, 1e-15);
        Ok(())
    }
}
