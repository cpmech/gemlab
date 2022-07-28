use crate::shapes::{GeoKind, Scratchpad};
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::{mat_mat_mul, mat_t_mat_mul, Matrix};
use russell_tensor::LinElasticity;

/// Performs analytical integrations on a Tri3
pub struct AnalyticalTri3 {
    /// Holds the area of the triangle
    pub area: f64,

    /// Holds the gradients (G-matrix)
    ///
    /// ```text
    ///             →
    /// →  →    dNᵐ(ξ)
    /// Gᵐ(ξ) = ——————
    ///            →
    ///           dx
    /// ```
    ///
    /// Organized as the G matrix (nnode=3, space_ndim=2)
    pub gg: Matrix,

    /// Holds the B-matrix (4, 6)
    pub bb: Matrix,
}

impl AnalyticalTri3 {
    pub fn new(pad: &Scratchpad) -> Self {
        assert_eq!(pad.kind, GeoKind::Tri3);

        // coefficients
        let (x0, y0) = (pad.xxt[0][0], pad.xxt[1][0]);
        let (x1, y1) = (pad.xxt[0][1], pad.xxt[1][1]);
        let (x2, y2) = (pad.xxt[0][2], pad.xxt[1][2]);
        let (a0, a1, a2) = (y1 - y2, y2 - y0, y0 - y1);
        let (b0, b1, b2) = (x2 - x1, x0 - x2, x1 - x0);
        let (f0, f1, f2) = (x1 * y2 - x2 * y1, x2 * y0 - x0 * y2, x0 * y1 - x1 * y0);

        // area
        let area = (f0 + f1 + f2) / 2.0;

        // auxiliary
        let r = 2.0 * area;
        let s = r * SQRT_2;

        // gradients
        #[rustfmt::skip]
        let gg = Matrix::from(&[
            [a0/r, b0/r],
            [a1/r, b1/r],
            [a2/r, b2/r],
        ]);

        // B-matrix
        #[rustfmt::skip]
        let bb = Matrix::from(&[
            [a0/r,  0.0, a1/r,  0.0, a2/r,  0.0],
            [ 0.0, b0/r,  0.0, b1/r,  0.0, b2/r],
            [ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
            [b0/s, a0/s, b1/s, a1/s, b2/s, a2/s],
        ]);

        // results
        AnalyticalTri3 { area, gg, bb }
    }

    /// Integrates shape times scalar with constant function s(x) = cₛ
    ///
    /// solution:
    ///
    /// ```text
    /// aᵐ = cₛ A / 3
    /// ```
    pub fn integ_vec_a_constant(&self, cs: f64) -> Vec<f64> {
        vec![cs * self.area / 3.0, cs * self.area / 3.0, cs * self.area / 3.0]
    }

    /// Integrates shape times vector with constant vector v(x) = {v₀, v₁}
    ///
    /// solution:
    ///
    /// ```text
    /// bᵐ₀ = v₀ A / 3
    /// bᵐ₁ = v₁ A / 3
    ///       ┌    ┐
    ///       │ v0 │
    ///       │ v1 │
    ///     A │ v0 │
    /// b = — │ v1 │
    ///     3 │ v0 │
    ///       │ v1 │
    ///       └    ┘
    /// ```
    pub fn integ_vec_b_constant(&self, v0: f64, v1: f64) -> Vec<f64> {
        vec![
            v0 * self.area / 3.0,
            v1 * self.area / 3.0,
            v0 * self.area / 3.0,
            v1 * self.area / 3.0,
            v0 * self.area / 3.0,
            v1 * self.area / 3.0,
        ]
    }

    /// Integrates vector dot gradient with constant vector function w(x) = {w₀, w₁}
    ///
    /// solution:
    ///
    /// ```text
    /// cᵐ = (w₀ Gᵐ₀ + w₁ Gᵐ₁) A
    /// ```
    pub fn integ_vec_c_constant(&self, w0: f64, w1: f64) -> Vec<f64> {
        vec![
            (w0 * self.gg[0][0] + w1 * self.gg[0][1]) * self.area,
            (w0 * self.gg[1][0] + w1 * self.gg[1][1]) * self.area,
            (w0 * self.gg[2][0] + w1 * self.gg[2][1]) * self.area,
        ]
    }

    /// Integrates vector dot gradient with bilinear vector function w(x) = {x, y}
    ///
    /// solution:
    ///
    /// ```text
    /// cᵐ = ((x₀+x₁+x₂) Gᵐ₀ + (y₀+y₁+y₂) Gᵐ₁) A / 3
    /// ```
    ///
    /// # Input
    ///
    /// * `pad` -- The same shape used in `new` because we need the nodal coordinates here
    ///            Do not change the coordinates, otherwise the values will be wrong.
    pub fn integ_vec_c_bilinear(&self, pad: &Scratchpad) -> Vec<f64> {
        let (x0, x1, x2) = (pad.xxt[0][0], pad.xxt[0][1], pad.xxt[0][2]);
        let (y0, y1, y2) = (pad.xxt[1][0], pad.xxt[1][1], pad.xxt[1][2]);
        vec![
            ((x0 + x1 + x2) * self.gg[0][0] + (y0 + y1 + y2) * self.gg[0][1]) * self.area / 3.0,
            ((x0 + x1 + x2) * self.gg[1][0] + (y0 + y1 + y2) * self.gg[1][1]) * self.area / 3.0,
            ((x0 + x1 + x2) * self.gg[2][0] + (y0 + y1 + y2) * self.gg[2][1]) * self.area / 3.0,
        ]
    }

    /// Integrates tensor dot gradient with constant tensor function σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2}
    ///
    /// solution:
    ///
    /// ```text
    /// dᵐ₀ = (σ₀₀ Gᵐ₀ + σ₀₁ Gᵐ₁) A
    /// dᵐ₁ = (σ₁₀ Gᵐ₀ + σ₁₁ Gᵐ₁) A
    /// ```
    ///
    /// σ₂₂ is ignored.
    ///
    /// # Input
    ///
    /// * `s₀₀, s₁₁, s₀₁` -- components of the constant tensor function: σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2}
    pub fn integ_vec_d_constant(&self, s00: f64, s11: f64, s01: f64) -> Vec<f64> {
        vec![
            (s00 * self.gg[0][0] + s01 * self.gg[0][1]) * self.area,
            (s01 * self.gg[0][0] + s11 * self.gg[0][1]) * self.area,
            (s00 * self.gg[1][0] + s01 * self.gg[1][1]) * self.area,
            (s01 * self.gg[1][0] + s11 * self.gg[1][1]) * self.area,
            (s00 * self.gg[2][0] + s01 * self.gg[2][1]) * self.area,
            (s01 * self.gg[2][0] + s11 * self.gg[2][1]) * self.area,
        ]
    }

    /// Performs the nsn integration
    pub fn integ_nsn(&self, s: f64, th: f64) -> Matrix {
        let c = th * s * self.area / 12.0;
        Matrix::from(&[[2.0 * c, c, c], [c, 2.0 * c, c], [c, c, 2.0 * c]])
    }

    /// Performs the gtg integration
    pub fn integ_gtg(&self, kx: f64, ky: f64) -> Matrix {
        let g = &self.gg;
        let k00 = self.area * (g[0][0] * g[0][0] * kx + g[0][1] * g[0][1] * ky);
        let k11 = self.area * (g[1][0] * g[1][0] * kx + g[1][1] * g[1][1] * ky);
        let k01 = self.area * (g[0][0] * g[1][0] * kx + g[0][1] * g[1][1] * ky);
        let k12 = self.area * (g[1][0] * g[2][0] * kx + g[1][1] * g[2][1] * ky);
        let k02 = self.area * (g[0][0] * g[2][0] * kx + g[0][1] * g[2][1] * ky);
        let k22 = self.area * (g[2][0] * g[2][0] * kx + g[2][1] * g[2][1] * ky);
        Matrix::from(&[[k00, k01, k02], [k01, k11, k12], [k02, k12, k22]])
    }

    /// Performs the ntn integration
    pub fn integ_ntn(&self, rho: f64, th: f64) -> Matrix {
        let c = th * rho * self.area / 12.0;
        Matrix::from(&[
            [c * 2.0, c * 0.0, c * 1.0, c * 0.0, c * 1.0, c * 0.0],
            [c * 0.0, c * 2.0, c * 0.0, c * 1.0, c * 0.0, c * 1.0],
            [c * 1.0, c * 0.0, c * 2.0, c * 0.0, c * 1.0, c * 0.0],
            [c * 0.0, c * 1.0, c * 0.0, c * 2.0, c * 0.0, c * 1.0],
            [c * 1.0, c * 0.0, c * 1.0, c * 0.0, c * 2.0, c * 0.0],
            [c * 0.0, c * 1.0, c * 0.0, c * 1.0, c * 0.0, c * 2.0],
        ])
    }

    /// Performs the gvn integration with constant vector
    #[rustfmt::skip]
    pub fn integ_gvn_constant(&self, v0: f64, v1: f64) -> Matrix {
        let aa = self.area;
        let g = &self.gg;
        Matrix::from(&[
            [aa*(g[0][0]*v0 + g[0][1]*v1)/3.0, aa*(g[0][0]*v0 + g[0][1]*v1)/3.0, aa*(g[0][0]*v0 + g[0][1]*v1)/3.0],
            [aa*(g[1][0]*v0 + g[1][1]*v1)/3.0, aa*(g[1][0]*v0 + g[1][1]*v1)/3.0, aa*(g[1][0]*v0 + g[1][1]*v1)/3.0],
            [aa*(g[2][0]*v0 + g[2][1]*v1)/3.0, aa*(g[2][0]*v0 + g[2][1]*v1)/3.0, aa*(g[2][0]*v0 + g[2][1]*v1)/3.0],
        ])
    }

    /// Performs the gvn integration with v = {x, y}
    #[rustfmt::skip]
    pub fn integ_gvn_bilinear(&self, pad: &Scratchpad) -> Matrix {
        let aa = self.area;
        let (g00, g01) = (self.gg[0][0], self.gg[0][1]);
        let (g10, g11) = (self.gg[1][0], self.gg[1][1]);
        let (g20, g21) = (self.gg[2][0], self.gg[2][1]);
        let (x0, y0) = (pad.xxt[0][0], pad.xxt[1][0]);
        let (x1, y1) = (pad.xxt[0][1], pad.xxt[1][1]);
        let (x2, y2) = (pad.xxt[0][2], pad.xxt[1][2]);
        Matrix::from(&[
            [(aa*(2.0*g00*x0 + g00*x1 + g00*x2 + 2.0*g01*y0 + g01*y1 + g01*y2))/12.0,(aa*(g00*x0 + 2.0*g00*x1 + g00*x2 + g01*y0 + 2.0*g01*y1 + g01*y2))/12.0, (aa*(g00*x0 + g00*x1 + 2.0*g00*x2 + g01*y0 + g01*y1 + 2.0*g01*y2))/12.0],
            [(aa*(2.0*g10*x0 + g10*x1 + g10*x2 + 2.0*g11*y0 + g11*y1 + g11*y2))/12.0,(aa*(g10*x0 + 2.0*g10*x1 + g10*x2 + g11*y0 + 2.0*g11*y1 + g11*y2))/12.0, (aa*(g10*x0 + g10*x1 + 2.0*g10*x2 + g11*y0 + g11*y1 + 2.0*g11*y2))/12.0],
            [(aa*(2.0*g20*x0 + g20*x1 + g20*x2 + 2.0*g21*y0 + g21*y1 + g21*y2))/12.0,(aa*(g20*x0 + 2.0*g20*x1 + g20*x2 + g21*y0 + 2.0*g21*y1 + g21*y2))/12.0, (aa*(g20*x0 + g20*x1 + 2.0*g20*x2 + g21*y0 + g21*y1 + 2.0*g21*y2))/12.0],
        ])
    }

    /// Performs the nvg integration with constant vector
    #[rustfmt::skip]
    pub fn integ_nvg_constant(&self, v0: f64, v1: f64) -> Matrix {
        let (g00, g01) = (self.gg[0][0], self.gg[0][1]);
        let (g10, g11) = (self.gg[1][0], self.gg[1][1]);
        let (g20, g21) = (self.gg[2][0], self.gg[2][1]);
        let c = self.area / 3.0;
        Matrix::from(&[
            [c*g00*v0, c*g01*v0, c*g10*v0, c*g11*v0, c*g20*v0, c*g21*v0],
            [c*g00*v1, c*g01*v1, c*g10*v1, c*g11*v1, c*g20*v1, c*g21*v1],
            [c*g00*v0, c*g01*v0, c*g10*v0, c*g11*v0, c*g20*v0, c*g21*v0],
            [c*g00*v1, c*g01*v1, c*g10*v1, c*g11*v1, c*g20*v1, c*g21*v1],
            [c*g00*v0, c*g01*v0, c*g10*v0, c*g11*v0, c*g20*v0, c*g21*v0],
            [c*g00*v1, c*g01*v1, c*g10*v1, c*g11*v1, c*g20*v1, c*g21*v1],
        ])
    }

    /// Calculates the stiffness matrix
    ///
    /// solution:
    ///
    /// ```text
    /// K = Bᵀ ⋅ D ⋅ B ⋅ th ⋅ area
    /// ```
    pub fn integ_stiffness(&self, young: f64, poisson: f64, plane_stress: bool, th: f64) -> Result<Matrix, StrError> {
        let ela = LinElasticity::new(young, poisson, true, plane_stress);
        let dd = ela.get_modulus();
        let dim_dd = 4;
        let dim_kk = 6;
        let mut bb_t_dd = Matrix::new(dim_kk, dim_dd);
        let mut kk = Matrix::new(dim_kk, dim_kk);
        mat_t_mat_mul(&mut bb_t_dd, 1.0, &self.bb, &dd.mat)?;
        mat_mat_mul(&mut kk, th * self.area, &bb_t_dd, &self.bb)?;
        Ok(kk)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::AnalyticalTri3;
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::StrError;
    use russell_chk::assert_approx_eq;

    #[test]
    fn new_works() -> Result<(), StrError> {
        // Tri3 # 1 from Figure 1.18, page 29 of [@bhatti]
        //
        // [@bhatti] Bhatti, M.A. (2005) Fundamental Finite Element Analysis
        //          and Applications, Wiley, 700p.
        //
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Tri3)?;
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, 0.2);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(2, 0, 0.1);
        pad.set_xx(2, 1, 0.1);
        let ana = AnalyticalTri3::new(&pad);
        assert_approx_eq!(ana.area, 0.01, 1e-15);
        // a = [-0.1, 0.1, 0.0]
        // b = [-0.1, -0.1, 0.2]
        // A = 0.01, Gmi = ai/(2 A) = -0.1 / 0.02 = -5
        assert_eq!(
            format!("{:.2}", ana.gg),
            "┌             ┐\n\
             │ -5.00 -5.00 │\n\
             │  5.00 -5.00 │\n\
             │  0.00 10.00 │\n\
             └             ┘"
        );
        Ok(())
    }
}
