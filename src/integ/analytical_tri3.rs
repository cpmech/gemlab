use crate::shapes::{GeoKind, Scratchpad};
use crate::StrError;
use russell_lab::math::SQRT_2;
use russell_lab::{mat_mat_mul, mat_t_mat_mul, Matrix, Vector};
use russell_tensor::{LinElasticity, Tensor2};

/// Performs analytical integrations on a Tri3
pub struct AnalyticalTri3 {
    /// Holds the area of the triangle
    pub area: f64,

    /// Holds the gradients (B-matrix)
    ///
    /// ```text
    ///             →
    /// →  →    dNᵐ(ξ)
    /// Bᵐ(ξ) = ——————
    ///            →
    ///           dx
    /// ```
    ///
    /// Organized as the B matrix (nnode=3, space_ndim=2)
    pub bb: Matrix,

    /// Holds the (element) Be-matrix (4, 6)
    pub bbe: Matrix,

    // Holds the lengths of each edge
    pub ll: Vec<f64>,

    // Coordinates
    x0: f64,
    x1: f64,
    x2: f64,
}

impl AnalyticalTri3 {
    /// Allocates a new instance
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
        let bb = Matrix::from(&[
            [a0/r, b0/r],
            [a1/r, b1/r],
            [a2/r, b2/r],
        ]);

        // element Be-matrix (plane-strain or plane-stress)
        #[rustfmt::skip]
        let bbe = Matrix::from(&[
            [a0/r,  0.0, a1/r,  0.0, a2/r,  0.0],
            [ 0.0, b0/r,  0.0, b1/r,  0.0, b2/r],
            [ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
            [b0/s, a0/s, b1/s, a1/s, b2/s, a2/s],
        ]);

        // edges lengths
        let ll = vec![
            f64::sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1)), // L01
            f64::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)), // L12
            f64::sqrt((x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0)), // L20
        ];

        // results
        AnalyticalTri3 {
            area,
            bb,
            bbe,
            ll,
            x0,
            x1,
            x2,
        }
    }

    /// Integrates shape times scalar with constant function s(x)
    pub fn vec_01_ns(&self, s: f64, axisymmetric: bool) -> Vector {
        if axisymmetric {
            let c = s * self.area / 12.0;
            Vector::from(&[
                c * (2.0 * self.x0 + self.x1 + self.x2),
                c * (self.x0 + 2.0 * self.x1 + self.x2),
                c * (self.x0 + self.x1 + 2.0 * self.x2),
            ])
        } else {
            let c = s * self.area / 3.0;
            Vector::from(&[c, c, c])
        }
    }

    /// Integrates shape times scalar with constant function s(x)
    ///
    /// **Important:** `side` must be 0, 1, or 2
    pub fn vec_01_ns_bry(&self, side: usize, s: f64, axisymmetric: bool) -> Vector {
        if axisymmetric {
            let (a0, a1) = match side {
                0 => (self.x1, self.x0),
                1 => (self.x2, self.x1),
                2 => (self.x0, self.x2),
                _ => panic!("invalid side"),
            };
            let c = s * self.ll[side] / 6.0;
            Vector::from(&[c * (2.0 * a0 + a1), c * (a0 + 2.0 * a1)])
        } else {
            let c = s * self.ll[side] / 2.0;
            Vector::from(&[c, c])
        }
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
    pub fn vec_02_nv(&self, v0: f64, v1: f64) -> Vector {
        Vector::from(&[
            v0 * self.area / 3.0,
            v1 * self.area / 3.0,
            v0 * self.area / 3.0,
            v1 * self.area / 3.0,
            v0 * self.area / 3.0,
            v1 * self.area / 3.0,
        ])
    }

    /// Integrates vector dot gradient with constant vector function w(x) = {w₀, w₁}
    ///
    /// solution:
    ///
    /// ```text
    /// cᵐ = (w₀ Bᵐ₀ + w₁ Bᵐ₁) A
    /// ```
    pub fn vec_03_vb(&self, w0: f64, w1: f64) -> Vector {
        Vector::from(&[
            (w0 * self.bb[0][0] + w1 * self.bb[0][1]) * self.area,
            (w0 * self.bb[1][0] + w1 * self.bb[1][1]) * self.area,
            (w0 * self.bb[2][0] + w1 * self.bb[2][1]) * self.area,
        ])
    }

    /// Integrates vector dot gradient with bilinear vector function w(x) = {x, y}
    ///
    /// solution:
    ///
    /// ```text
    /// cᵐ = ((x₀+x₁+x₂) Bᵐ₀ + (y₀+y₁+y₂) Bᵐ₁) A / 3
    /// ```
    ///
    /// # Input
    ///
    /// * `pad` -- The same shape used in `new` because we need the nodal coordinates here
    ///            Do not change the coordinates, otherwise the values will be wrong.
    pub fn vec_03_vb_bilinear(&self, pad: &Scratchpad) -> Vector {
        let (x0, x1, x2) = (pad.xxt[0][0], pad.xxt[0][1], pad.xxt[0][2]);
        let (y0, y1, y2) = (pad.xxt[1][0], pad.xxt[1][1], pad.xxt[1][2]);
        Vector::from(&[
            ((x0 + x1 + x2) * self.bb[0][0] + (y0 + y1 + y2) * self.bb[0][1]) * self.area / 3.0,
            ((x0 + x1 + x2) * self.bb[1][0] + (y0 + y1 + y2) * self.bb[1][1]) * self.area / 3.0,
            ((x0 + x1 + x2) * self.bb[2][0] + (y0 + y1 + y2) * self.bb[2][1]) * self.area / 3.0,
        ])
    }

    /// Integrates tensor dot gradient with constant tensor function σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2}
    #[rustfmt::skip]
    pub fn vec_04_tb(&self, tt: &Tensor2, axisymmetric: bool) -> Vector {
        let (x0, x1, x2) = (self.x0, self.x1, self.x2);
        let (b00, b01) = (self.bb[0][0], self.bb[0][1]);
        let (b10, b11) = (self.bb[1][0], self.bb[1][1]);
        let (b20, b21) = (self.bb[2][0], self.bb[2][1]);
        let (t0, t1, t2, t3) = (tt.vec[0], tt.vec[1], tt.vec[2], tt.vec[3]);
        if axisymmetric {
            let c = self.area / 3.0;
            Vector::from(&[
                c * (b00 * t0 + b01 * t3 / SQRT_2) * (x0 + x1 + x2) + c * t2,
                c * (b01 * t1 + b00 * t3 / SQRT_2) * (x0 + x1 + x2),
                c * (b10 * t0 + b11 * t3 / SQRT_2) * (x0 + x1 + x2) + c * t2,
                c * (b11 * t1 + b10 * t3 / SQRT_2) * (x0 + x1 + x2),
                c * (b20 * t0 + b21 * t3 / SQRT_2) * (x0 + x1 + x2) + c * t2,
                c * (b21 * t1 + b20 * t3 / SQRT_2) * (x0 + x1 + x2),
            ])
        } else {
            let c = self.area;
            Vector::from(&[
                c * (b00 * t0 + b01 * t3 / SQRT_2),
                c * (b01 * t1 + b00 * t3 / SQRT_2),
                c * (b10 * t0 + b11 * t3 / SQRT_2),
                c * (b11 * t1 + b10 * t3 / SQRT_2),
                c * (b20 * t0 + b21 * t3 / SQRT_2),
                c * (b21 * t1 + b20 * t3 / SQRT_2),
            ])
        }
    }

    /// Performs the n-s-n integration with constant s(x) field
    pub fn mat_01_nsn(&self, s: f64, th: f64) -> Matrix {
        let c = th * s * self.area / 12.0;
        Matrix::from(&[[2.0 * c, c, c], [c, 2.0 * c, c], [c, c, 2.0 * c]])
    }

    /// Performs the n-s-n integration with constant s(x) field (boundary integral version)
    ///
    /// **Important:** `side` must be 0, 1, or 2
    pub fn mat_01_nsn_bry(&self, side: usize, s: f64, axisymmetric: bool) -> Matrix {
        if axisymmetric {
            let (a0, a1) = match side {
                0 => (self.x1, self.x0),
                1 => (self.x2, self.x1),
                2 => (self.x0, self.x2),
                _ => panic!("invalid side"),
            };
            let c = s * self.ll[side] / 12.0;
            Matrix::from(&[
                [c * (3.0 * a0 + a1), c * (a0 + a1)],
                [c * (a0 + a1), c * (a0 + 3.0 * a1)],
            ])
        } else {
            let c = s * self.ll[side] / 6.0;
            Matrix::from(&[[2.0 * c, c], [c, 2.0 * c]])
        }
    }

    /// Performs the b-v-n integration with constant vector
    #[rustfmt::skip]
    pub fn mat_02_bvn(&self, v0: f64, v1: f64) -> Matrix {
        let aa = self.area;
        let g = &self.bb;
        Matrix::from(&[
            [aa*(g[0][0]*v0 + g[0][1]*v1)/3.0, aa*(g[0][0]*v0 + g[0][1]*v1)/3.0, aa*(g[0][0]*v0 + g[0][1]*v1)/3.0],
            [aa*(g[1][0]*v0 + g[1][1]*v1)/3.0, aa*(g[1][0]*v0 + g[1][1]*v1)/3.0, aa*(g[1][0]*v0 + g[1][1]*v1)/3.0],
            [aa*(g[2][0]*v0 + g[2][1]*v1)/3.0, aa*(g[2][0]*v0 + g[2][1]*v1)/3.0, aa*(g[2][0]*v0 + g[2][1]*v1)/3.0],
        ])
    }

    /// Performs the b-v-n integration with v = {x, y}
    #[rustfmt::skip]
    pub fn mat_02_bvn_bilinear(&self, pad: &Scratchpad) -> Matrix {
        let aa = self.area;
        let (b00, b01) = (self.bb[0][0], self.bb[0][1]);
        let (b10, b11) = (self.bb[1][0], self.bb[1][1]);
        let (b20, b21) = (self.bb[2][0], self.bb[2][1]);
        let (x0, y0) = (pad.xxt[0][0], pad.xxt[1][0]);
        let (x1, y1) = (pad.xxt[0][1], pad.xxt[1][1]);
        let (x2, y2) = (pad.xxt[0][2], pad.xxt[1][2]);
        Matrix::from(&[
            [(aa*(2.0*b00*x0 + b00*x1 + b00*x2 + 2.0*b01*y0 + b01*y1 + b01*y2))/12.0,(aa*(b00*x0 + 2.0*b00*x1 + b00*x2 + b01*y0 + 2.0*b01*y1 + b01*y2))/12.0, (aa*(b00*x0 + b00*x1 + 2.0*b00*x2 + b01*y0 + b01*y1 + 2.0*b01*y2))/12.0],
            [(aa*(2.0*b10*x0 + b10*x1 + b10*x2 + 2.0*b11*y0 + b11*y1 + b11*y2))/12.0,(aa*(b10*x0 + 2.0*b10*x1 + b10*x2 + b11*y0 + 2.0*b11*y1 + b11*y2))/12.0, (aa*(b10*x0 + b10*x1 + 2.0*b10*x2 + b11*y0 + b11*y1 + 2.0*b11*y2))/12.0],
            [(aa*(2.0*b20*x0 + b20*x1 + b20*x2 + 2.0*b21*y0 + b21*y1 + b21*y2))/12.0,(aa*(b20*x0 + 2.0*b20*x1 + b20*x2 + b21*y0 + 2.0*b21*y1 + b21*y2))/12.0, (aa*(b20*x0 + b20*x1 + 2.0*b20*x2 + b21*y0 + b21*y1 + 2.0*b21*y2))/12.0],
        ])
    }

    /// Performs the b-t-b integration with constant tensor
    pub fn mat_03_btb(&self, kx: f64, ky: f64, axisymmetric: bool) -> Matrix {
        let c = if axisymmetric {
            self.area * (self.x0 + self.x1 + self.x2) / 3.0
        } else {
            self.area
        };
        let (b00, b01) = (self.bb[0][0], self.bb[0][1]);
        let (b10, b11) = (self.bb[1][0], self.bb[1][1]);
        let (b20, b21) = (self.bb[2][0], self.bb[2][1]);
        let k00 = c * (b00 * b00 * kx + b01 * b01 * ky);
        let k11 = c * (b10 * b10 * kx + b11 * b11 * ky);
        let k01 = c * (b00 * b10 * kx + b01 * b11 * ky);
        let k12 = c * (b10 * b20 * kx + b11 * b21 * ky);
        let k02 = c * (b00 * b20 * kx + b01 * b21 * ky);
        let k22 = c * (b20 * b20 * kx + b21 * b21 * ky);
        Matrix::from(&[[k00, k01, k02], [k01, k11, k12], [k02, k12, k22]])
    }

    /// Performs the n-t-n integration with constant tensor field rho*delta
    pub fn mat_08_ntn(&self, rho: f64, th: f64) -> Matrix {
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

    /// Performs the n-v-b integration with constant vector
    #[rustfmt::skip]
    pub fn mat_09_nvb(&self, v0: f64, v1: f64) -> Matrix {
        let (b00, b01) = (self.bb[0][0], self.bb[0][1]);
        let (b10, b11) = (self.bb[1][0], self.bb[1][1]);
        let (b20, b21) = (self.bb[2][0], self.bb[2][1]);
        let c = self.area / 3.0;
        Matrix::from(&[
            [c*b00*v0, c*b01*v0, c*b10*v0, c*b11*v0, c*b20*v0, c*b21*v0],
            [c*b00*v1, c*b01*v1, c*b10*v1, c*b11*v1, c*b20*v1, c*b21*v1],
            [c*b00*v0, c*b01*v0, c*b10*v0, c*b11*v0, c*b20*v0, c*b21*v0],
            [c*b00*v1, c*b01*v1, c*b10*v1, c*b11*v1, c*b20*v1, c*b21*v1],
            [c*b00*v0, c*b01*v0, c*b10*v0, c*b11*v0, c*b20*v0, c*b21*v0],
            [c*b00*v1, c*b01*v1, c*b10*v1, c*b11*v1, c*b20*v1, c*b21*v1],
        ])
    }

    /// Performs the b-d-b integration with constant tensor field (calculates the stiffness matrix)
    ///
    /// solution:
    ///
    /// ```text
    /// K = Bᵀ ⋅ D ⋅ B ⋅ th ⋅ area
    /// ```
    pub fn mat_10_bdb(&self, young: f64, poisson: f64, plane_stress: bool, th: f64) -> Result<Matrix, StrError> {
        let ela = LinElasticity::new(young, poisson, true, plane_stress);
        let dd = ela.get_modulus();
        let dim_dd = 4;
        let dim_kk = 6;
        let mut bb_t_dd = Matrix::new(dim_kk, dim_dd);
        let mut kk = Matrix::new(dim_kk, dim_kk);
        mat_t_mat_mul(&mut bb_t_dd, 1.0, &self.bbe, &dd.mat).unwrap();
        mat_mat_mul(&mut kk, th * self.area, &bb_t_dd, &self.bbe).unwrap();
        Ok(kk)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::AnalyticalTri3;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_chk::approx_eq;

    #[test]
    fn new_works() {
        // Tri3 # 1 from Figure 1.18, page 29 of @bhatti
        // @bhatti Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Tri3).unwrap();
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, 0.2);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(2, 0, 0.1);
        pad.set_xx(2, 1, 0.1);
        let ana = AnalyticalTri3::new(&pad);
        approx_eq(ana.area, 0.01, 1e-15);
        // a = [-0.1, 0.1, 0.0]
        // b = [-0.1, -0.1, 0.2]
        // A = 0.01, Gmi = ai/(2 A) = -0.1 / 0.02 = -5
        assert_eq!(
            format!("{:.2}", ana.bb),
            "┌             ┐\n\
             │ -5.00 -5.00 │\n\
             │  5.00 -5.00 │\n\
             │  0.00 10.00 │\n\
             └             ┘"
        );
    }
}
