use crate::shapes::{GeoKind, Scratchpad};
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::{mat_mat_mul, mat_t_mat_mul, Matrix};
use russell_tensor::{LinElasticity, Tensor2};

/// Performs analytical integrations on a Tet4
pub struct AnalyticalTet4 {
    /// Holds the volume of the tetrahedron
    pub volume: f64,

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
    /// Organized as the G matrix (nnode=4, space_ndim=3)
    pub gg: Matrix,

    /// Holds the B-matrix (6, 12)
    pub bb: Matrix,
}

impl AnalyticalTet4 {
    /// Allocates a new instance
    pub fn new(pad: &Scratchpad) -> Self {
        assert_eq!(pad.kind, GeoKind::Tet4);

        let x1 = pad.xxt[0][0];
        let x2 = pad.xxt[0][1];
        let x3 = pad.xxt[0][2];
        let x4 = pad.xxt[0][3];

        let y1 = pad.xxt[1][0];
        let y2 = pad.xxt[1][1];
        let y3 = pad.xxt[1][2];
        let y4 = pad.xxt[1][3];

        let z1 = pad.xxt[2][0];
        let z2 = pad.xxt[2][1];
        let z3 = pad.xxt[2][2];
        let z4 = pad.xxt[2][3];

        let x12 = x1 - x2;
        let x13 = x1 - x3;
        let x14 = x1 - x4;
        let x23 = x2 - x3;
        let x24 = x2 - x4;
        let x34 = x3 - x4;
        let x21 = -x12;
        let x31 = -x13;
        // let x41 = -x14; // no needed
        let x32 = -x23;
        let x42 = -x24;
        let x43 = -x34;

        let y12 = y1 - y2;
        let y13 = y1 - y3;
        let y14 = y1 - y4;
        let y23 = y2 - y3;
        let y24 = y2 - y4;
        let y34 = y3 - y4;
        let y21 = -y12;
        let y31 = -y13;
        // let y41 = -y14; // no needed
        let y32 = -y23;
        let y42 = -y24;
        let y43 = -y34;

        let z12 = z1 - z2;
        let z13 = z1 - z3;
        let z14 = z1 - z4;
        let z23 = z2 - z3;
        let z24 = z2 - z4;
        let z34 = z3 - z4;
        let z21 = -z12;
        let z31 = -z13;
        // let z41 = -z14; // no needed
        let z32 = -z23;
        let z42 = -z24;
        let z43 = -z34;

        let jj_det = x21 * (y23 * z34 - y34 * z23) + x32 * (y34 * z12 - y12 * z34) + x43 * (y12 * z23 - y23 * z12);

        let a1 = y42 * z32 - y32 * z42;
        let b1 = x32 * z42 - x42 * z32;
        let c1 = x42 * y32 - x32 * y42;
        let a2 = y31 * z43 - y34 * z13;
        let b2 = x43 * z31 - x13 * z34;
        let c2 = x31 * y43 - x34 * y13;
        let a3 = y24 * z14 - y14 * z24;
        let b3 = x14 * z24 - x24 * z14;
        let c3 = x24 * y14 - x14 * y24;
        let a4 = y13 * z21 - y12 * z31;
        let b4 = x21 * z13 - x31 * z12;
        let c4 = x13 * y21 - x12 * y31;

        // auxiliary
        let r = jj_det; // == 6 * V
        let s = r * SQRT_2;

        // gradients
        #[rustfmt::skip]
        let gg = Matrix::from(&[
            [a1/r, b1/r, c1/r],
            [a2/r, b2/r, c2/r],
            [a3/r, b3/r, c3/r],
            [a4/r, b4/r, c4/r],
        ]);

        // B-matrix
        #[rustfmt::skip]
        let bb = Matrix::from(&[
            [a1/r,  0.0,  0.0, a2/r,  0.0,  0.0, a3/r,  0.0,  0.0, a4/r,  0.0,  0.0],
            [ 0.0, b1/r,  0.0,  0.0, b2/r,  0.0,  0.0, b3/r,  0.0,  0.0, b4/r,  0.0],
            [ 0.0,  0.0, c1/r,  0.0,  0.0, c2/r,  0.0,  0.0, c3/r,  0.0,  0.0, c4/r],
            [b1/s, a1/s,  0.0, b2/s, a2/s,  0.0, b3/s, a3/s,  0.0, b4/s, a4/s,  0.0],
            [ 0.0, c1/s, b1/s,  0.0, c2/s, b2/s,  0.0, c3/s, b3/s,  0.0, c4/s, b4/s],
            [c1/s,  0.0, a1/s, c2/s,  0.0, a2/s, c3/s,  0.0, a3/s, c4/s,  0.0, a4/s],
        ]);

        AnalyticalTet4 {
            volume: jj_det / 6.0,
            gg,
            bb,
        }
    }

    /// Integrates shape times scalar with linear scalar function s(x) = x₂ = z
    ///
    /// # Input
    ///
    /// * `pad` -- The same pad used in `new` because we need the nodal coordinates here.
    ///            Do not change the coordinates, otherwise the values will be wrong.
    pub fn integ_vec_a_linear_along_z(&self, pad: &Scratchpad) -> Vec<f64> {
        let (z1, z2, z3, z4) = (pad.xxt[2][0], pad.xxt[2][1], pad.xxt[2][2], pad.xxt[2][3]);
        vec![
            (self.volume * (2.0 * z1 + z2 + z3 + z4)) / 20.0,
            (self.volume * (z1 + 2.0 * z2 + z3 + z4)) / 20.0,
            (self.volume * (z1 + z2 + 2.0 * z3 + z4)) / 20.0,
            (self.volume * (z1 + z2 + z3 + 2.0 * z4)) / 20.0,
        ]
    }

    /// Integrates shape times vector with constant vector v(x) = {v0, v1, v2}
    ///
    /// solution:
    ///
    /// ```text
    /// bᵐ₀ = v₀ V / 4
    /// bᵐ₁ = v₁ V / 4
    /// bᵐ₂ = v₂ V / 4
    /// ```
    pub fn integ_vec_b(&self, v0: f64, v1: f64, v2: f64) -> Vec<f64> {
        vec![
            v0 * self.volume / 4.0,
            v1 * self.volume / 4.0,
            v2 * self.volume / 4.0,
            v0 * self.volume / 4.0,
            v1 * self.volume / 4.0,
            v2 * self.volume / 4.0,
            v0 * self.volume / 4.0,
            v1 * self.volume / 4.0,
            v2 * self.volume / 4.0,
            v0 * self.volume / 4.0,
            v1 * self.volume / 4.0,
            v2 * self.volume / 4.0,
        ]
    }

    /// Integrates vector dot gradient with constant vector w(x) = {w0, w1, w2}
    ///
    /// Solution:
    ///
    /// ```text
    /// cᵐ = (w₀ Gᵐ₀ + w₁ Gᵐ₁ + w₂ Gᵐ₂) V
    /// ```
    pub fn integ_vec_c(&self, w0: f64, w1: f64, w2: f64) -> Vec<f64> {
        vec![
            (w0 * self.gg[0][0] + w1 * self.gg[0][1] + w2 * self.gg[0][2]) * self.volume,
            (w0 * self.gg[1][0] + w1 * self.gg[1][1] + w2 * self.gg[1][2]) * self.volume,
            (w0 * self.gg[2][0] + w1 * self.gg[2][1] + w2 * self.gg[2][2]) * self.volume,
            (w0 * self.gg[3][0] + w1 * self.gg[3][1] + w2 * self.gg[3][2]) * self.volume,
        ]
    }

    /// Integrates tensor dot gradient with constant tensor function σ(x)
    ///
    /// Solution:
    ///
    /// ```text
    /// σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2, σ₁₂√2, σ₀₂√2}
    ///
    /// dᵐ₀ = (σ₀₀ Gᵐ₀ + σ₀₁ Gᵐ₁ + σ₀₂ Gᵐ₂) V
    /// dᵐ₁ = (σ₁₀ Gᵐ₀ + σ₁₁ Gᵐ₁ + σ₁₂ Gᵐ₂) V
    /// dᵐ₂ = (σ₂₀ Gᵐ₀ + σ₂₁ Gᵐ₁ + σ₂₂ Gᵐ₂) V
    /// ```
    pub fn integ_vec_d(&self, tt: &Tensor2) -> Vec<f64> {
        let c = self.volume;
        let mat = tt.to_matrix();
        let (a00, a01, a02) = (mat[0][0], mat[0][1], mat[0][2]);
        let (a11, a12) = (mat[1][1], mat[1][2]);
        let a22 = mat[2][2];
        vec![
            c * (a00 * self.gg[0][0] + a01 * self.gg[0][1] + a02 * self.gg[0][2]),
            c * (a01 * self.gg[0][0] + a11 * self.gg[0][1] + a12 * self.gg[0][2]),
            c * (a02 * self.gg[0][0] + a12 * self.gg[0][1] + a22 * self.gg[0][2]),
            c * (a00 * self.gg[1][0] + a01 * self.gg[1][1] + a02 * self.gg[1][2]),
            c * (a01 * self.gg[1][0] + a11 * self.gg[1][1] + a12 * self.gg[1][2]),
            c * (a02 * self.gg[1][0] + a12 * self.gg[1][1] + a22 * self.gg[1][2]),
            c * (a00 * self.gg[2][0] + a01 * self.gg[2][1] + a02 * self.gg[2][2]),
            c * (a01 * self.gg[2][0] + a11 * self.gg[2][1] + a12 * self.gg[2][2]),
            c * (a02 * self.gg[2][0] + a12 * self.gg[2][1] + a22 * self.gg[2][2]),
            c * (a00 * self.gg[3][0] + a01 * self.gg[3][1] + a02 * self.gg[3][2]),
            c * (a01 * self.gg[3][0] + a11 * self.gg[3][1] + a12 * self.gg[3][2]),
            c * (a02 * self.gg[3][0] + a12 * self.gg[3][1] + a22 * self.gg[3][2]),
        ]
    }

    /// Performs the n-s-n integration with constant scalar function
    #[rustfmt::skip]
    pub fn integ_nsn(&self, s: f64) -> Matrix {
        let c = self.volume / 20.0;
        Matrix::from(&[
            [2.0 * s * c, s * c, s * c, s * c],
            [s * c, 2.0 * s * c, s * c, s * c],
            [s * c, s * c, 2.0 * s * c, s * c],
            [s * c, s * c, s * c, 2.0 * s * c],
        ])
    }

    /// Performs the g-v-n integration with constant vector field
    #[rustfmt::skip]
    pub fn integ_gvn(&self, v0: f64, v1: f64, v2: f64) -> Matrix {
        let c = self.volume / 4.0;
        let (g00, g01, g02) = (self.gg[0][0], self.gg[0][1], self.gg[0][2]);
        let (g10, g11, g12) = (self.gg[1][0], self.gg[1][1], self.gg[1][2]);
        let (g20, g21, g22) = (self.gg[2][0], self.gg[2][1], self.gg[2][2]);
        let (g30, g31, g32) = (self.gg[3][0], self.gg[3][1], self.gg[3][2]);
        Matrix::from(&[
            [c*(g00*v0 + g01*v1 + g02*v2), c*(g00*v0 + g01*v1 + g02*v2), c*(g00*v0 + g01*v1 + g02*v2), c*(g00*v0 + g01*v1 + g02*v2)],
            [c*(g10*v0 + g11*v1 + g12*v2), c*(g10*v0 + g11*v1 + g12*v2), c*(g10*v0 + g11*v1 + g12*v2), c*(g10*v0 + g11*v1 + g12*v2)],
            [c*(g20*v0 + g21*v1 + g22*v2), c*(g20*v0 + g21*v1 + g22*v2), c*(g20*v0 + g21*v1 + g22*v2), c*(g20*v0 + g21*v1 + g22*v2)],
            [c*(g30*v0 + g31*v1 + g32*v2), c*(g30*v0 + g31*v1 + g32*v2), c*(g30*v0 + g31*v1 + g32*v2), c*(g30*v0 + g31*v1 + g32*v2)],
        ])
    }

    /// Performs the g-t-g integration with constant tensor field
    #[rustfmt::skip]
    pub fn integ_gtg(&self, tt: &Tensor2) -> Matrix {
        let c = self.volume;
        let mat = tt.to_matrix();
        let (a00, a01, a02) = (mat[0][0], mat[0][1], mat[0][2]);
        let (a10, a11, a12) = (mat[1][0], mat[1][1], mat[1][2]);
        let (a20, a21, a22) = (mat[2][0], mat[2][1], mat[2][2]);
        let (g00, g01, g02) = (self.gg[0][0], self.gg[0][1], self.gg[0][2]);
        let (g10, g11, g12) = (self.gg[1][0], self.gg[1][1], self.gg[1][2]);
        let (g20, g21, g22) = (self.gg[2][0], self.gg[2][1], self.gg[2][2]);
        let (g30, g31, g32) = (self.gg[3][0], self.gg[3][1], self.gg[3][2]);
        Matrix::from(&[
            [c*g00*(a00*g00 + a10*g01 + a20*g02) + c*g01*(a01*g00 + a11*g01 + a21*g02) + c*g02*(a02*g00 + a12*g01 + a22*g02), c*g10*(a00*g00 + a10*g01 + a20*g02) + c*g11*(a01*g00 + a11*g01 + a21*g02) + c*g12*(a02*g00 + a12*g01 + a22*g02), c*g20*(a00*g00 + a10*g01 + a20*g02) + c*g21*(a01*g00 + a11*g01 + a21*g02) + c*g22*(a02*g00 + a12*g01 + a22*g02), c*g30*(a00*g00 + a10*g01 + a20*g02) + c*g31*(a01*g00 + a11*g01 + a21*g02) + c*g32*(a02*g00 + a12*g01 + a22*g02)],
            [c*g00*(a00*g10 + a10*g11 + a20*g12) + c*g01*(a01*g10 + a11*g11 + a21*g12) + c*g02*(a02*g10 + a12*g11 + a22*g12), c*g10*(a00*g10 + a10*g11 + a20*g12) + c*g11*(a01*g10 + a11*g11 + a21*g12) + c*g12*(a02*g10 + a12*g11 + a22*g12), c*g20*(a00*g10 + a10*g11 + a20*g12) + c*g21*(a01*g10 + a11*g11 + a21*g12) + c*g22*(a02*g10 + a12*g11 + a22*g12), c*g30*(a00*g10 + a10*g11 + a20*g12) + c*g31*(a01*g10 + a11*g11 + a21*g12) + c*g32*(a02*g10 + a12*g11 + a22*g12)],
            [c*g00*(a00*g20 + a10*g21 + a20*g22) + c*g01*(a01*g20 + a11*g21 + a21*g22) + c*g02*(a02*g20 + a12*g21 + a22*g22), c*g10*(a00*g20 + a10*g21 + a20*g22) + c*g11*(a01*g20 + a11*g21 + a21*g22) + c*g12*(a02*g20 + a12*g21 + a22*g22), c*g20*(a00*g20 + a10*g21 + a20*g22) + c*g21*(a01*g20 + a11*g21 + a21*g22) + c*g22*(a02*g20 + a12*g21 + a22*g22), c*g30*(a00*g20 + a10*g21 + a20*g22) + c*g31*(a01*g20 + a11*g21 + a21*g22) + c*g32*(a02*g20 + a12*g21 + a22*g22)],
            [c*g00*(a00*g30 + a10*g31 + a20*g32) + c*g01*(a01*g30 + a11*g31 + a21*g32) + c*g02*(a02*g30 + a12*g31 + a22*g32), c*g10*(a00*g30 + a10*g31 + a20*g32) + c*g11*(a01*g30 + a11*g31 + a21*g32) + c*g12*(a02*g30 + a12*g31 + a22*g32), c*g20*(a00*g30 + a10*g31 + a20*g32) + c*g21*(a01*g30 + a11*g31 + a21*g32) + c*g22*(a02*g30 + a12*g31 + a22*g32), c*g30*(a00*g30 + a10*g31 + a20*g32) + c*g31*(a01*g30 + a11*g31 + a21*g32) + c*g32*(a02*g30 + a12*g31 + a22*g32)],
        ])
    }

    /// Performs the n-v-g integration with constant vector field
    #[rustfmt::skip]
    pub fn integ_nvg(&self, v0: f64, v1: f64, v2: f64) -> Matrix {
        let c = self.volume / 4.0;
        let (g00, g01, g02) = (self.gg[0][0], self.gg[0][1], self.gg[0][2]);
        let (g10, g11, g12) = (self.gg[1][0], self.gg[1][1], self.gg[1][2]);
        let (g20, g21, g22) = (self.gg[2][0], self.gg[2][1], self.gg[2][2]);
        let (g30, g31, g32) = (self.gg[3][0], self.gg[3][1], self.gg[3][2]);
        Matrix::from(&[
            [c*g00*v0, c*g01*v0, c*g02*v0, c*g10*v0, c*g11*v0, c*g12*v0, c*g20*v0, c*g21*v0, c*g22*v0, c*g30*v0, c*g31*v0, c*g32*v0],
            [c*g00*v1, c*g01*v1, c*g02*v1, c*g10*v1, c*g11*v1, c*g12*v1, c*g20*v1, c*g21*v1, c*g22*v1, c*g30*v1, c*g31*v1, c*g32*v1],
            [c*g00*v2, c*g01*v2, c*g02*v2, c*g10*v2, c*g11*v2, c*g12*v2, c*g20*v2, c*g21*v2, c*g22*v2, c*g30*v2, c*g31*v2, c*g32*v2],
            [c*g00*v0, c*g01*v0, c*g02*v0, c*g10*v0, c*g11*v0, c*g12*v0, c*g20*v0, c*g21*v0, c*g22*v0, c*g30*v0, c*g31*v0, c*g32*v0],
            [c*g00*v1, c*g01*v1, c*g02*v1, c*g10*v1, c*g11*v1, c*g12*v1, c*g20*v1, c*g21*v1, c*g22*v1, c*g30*v1, c*g31*v1, c*g32*v1],
            [c*g00*v2, c*g01*v2, c*g02*v2, c*g10*v2, c*g11*v2, c*g12*v2, c*g20*v2, c*g21*v2, c*g22*v2, c*g30*v2, c*g31*v2, c*g32*v2],
            [c*g00*v0, c*g01*v0, c*g02*v0, c*g10*v0, c*g11*v0, c*g12*v0, c*g20*v0, c*g21*v0, c*g22*v0, c*g30*v0, c*g31*v0, c*g32*v0],
            [c*g00*v1, c*g01*v1, c*g02*v1, c*g10*v1, c*g11*v1, c*g12*v1, c*g20*v1, c*g21*v1, c*g22*v1, c*g30*v1, c*g31*v1, c*g32*v1],
            [c*g00*v2, c*g01*v2, c*g02*v2, c*g10*v2, c*g11*v2, c*g12*v2, c*g20*v2, c*g21*v2, c*g22*v2, c*g30*v2, c*g31*v2, c*g32*v2],
            [c*g00*v0, c*g01*v0, c*g02*v0, c*g10*v0, c*g11*v0, c*g12*v0, c*g20*v0, c*g21*v0, c*g22*v0, c*g30*v0, c*g31*v0, c*g32*v0],
            [c*g00*v1, c*g01*v1, c*g02*v1, c*g10*v1, c*g11*v1, c*g12*v1, c*g20*v1, c*g21*v1, c*g22*v1, c*g30*v1, c*g31*v1, c*g32*v1],
            [c*g00*v2, c*g01*v2, c*g02*v2, c*g10*v2, c*g11*v2, c*g12*v2, c*g20*v2, c*g21*v2, c*g22*v2, c*g30*v2, c*g31*v2, c*g32*v2],
        ])
    }

    /// Performs the g-t-g integration with constant tensor field
    #[rustfmt::skip]
    pub fn integ_ntn(&self, sig: &Tensor2) -> Matrix {
        let vv = self.volume;
        let mat = sig.to_matrix();
        let (a00, a01, a02) = (mat[0][0], mat[0][1], mat[0][2]);
        let (a10, a11, a12) = (mat[1][0], mat[1][1], mat[1][2]);
        let (a20, a21, a22) = (mat[2][0], mat[2][1], mat[2][2]);
        Matrix::from(&[
            [a00*vv/10.0, a01*vv/10.0, a02*vv/10.0, a00*vv/20.0, a01*vv/20.0, a02*vv/20.0, a00*vv/20.0, a01*vv/20.0, a02*vv/20.0, a00*vv/20.0, a01*vv/20.0, a02*vv/20.0],
            [a10*vv/10.0, a11*vv/10.0, a12*vv/10.0, a10*vv/20.0, a11*vv/20.0, a12*vv/20.0, a10*vv/20.0, a11*vv/20.0, a12*vv/20.0, a10*vv/20.0, a11*vv/20.0, a12*vv/20.0],
            [a20*vv/10.0, a21*vv/10.0, a22*vv/10.0, a20*vv/20.0, a21*vv/20.0, a22*vv/20.0, a20*vv/20.0, a21*vv/20.0, a22*vv/20.0, a20*vv/20.0, a21*vv/20.0, a22*vv/20.0],
            [a00*vv/20.0, a01*vv/20.0, a02*vv/20.0, a00*vv/10.0, a01*vv/10.0, a02*vv/10.0, a00*vv/20.0, a01*vv/20.0, a02*vv/20.0, a00*vv/20.0, a01*vv/20.0, a02*vv/20.0],
            [a10*vv/20.0, a11*vv/20.0, a12*vv/20.0, a10*vv/10.0, a11*vv/10.0, a12*vv/10.0, a10*vv/20.0, a11*vv/20.0, a12*vv/20.0, a10*vv/20.0, a11*vv/20.0, a12*vv/20.0],
            [a20*vv/20.0, a21*vv/20.0, a22*vv/20.0, a20*vv/10.0, a21*vv/10.0, a22*vv/10.0, a20*vv/20.0, a21*vv/20.0, a22*vv/20.0, a20*vv/20.0, a21*vv/20.0, a22*vv/20.0],
            [a00*vv/20.0, a01*vv/20.0, a02*vv/20.0, a00*vv/20.0, a01*vv/20.0, a02*vv/20.0, a00*vv/10.0, a01*vv/10.0, a02*vv/10.0, a00*vv/20.0, a01*vv/20.0, a02*vv/20.0],
            [a10*vv/20.0, a11*vv/20.0, a12*vv/20.0, a10*vv/20.0, a11*vv/20.0, a12*vv/20.0, a10*vv/10.0, a11*vv/10.0, a12*vv/10.0, a10*vv/20.0, a11*vv/20.0, a12*vv/20.0],
            [a20*vv/20.0, a21*vv/20.0, a22*vv/20.0, a20*vv/20.0, a21*vv/20.0, a22*vv/20.0, a20*vv/10.0, a21*vv/10.0, a22*vv/10.0, a20*vv/20.0, a21*vv/20.0, a22*vv/20.0],
            [a00*vv/20.0, a01*vv/20.0, a02*vv/20.0, a00*vv/20.0, a01*vv/20.0, a02*vv/20.0, a00*vv/20.0, a01*vv/20.0, a02*vv/20.0, a00*vv/10.0, a01*vv/10.0, a02*vv/10.0],
            [a10*vv/20.0, a11*vv/20.0, a12*vv/20.0, a10*vv/20.0, a11*vv/20.0, a12*vv/20.0, a10*vv/20.0, a11*vv/20.0, a12*vv/20.0, a10*vv/10.0, a11*vv/10.0, a12*vv/10.0],
            [a20*vv/20.0, a21*vv/20.0, a22*vv/20.0, a20*vv/20.0, a21*vv/20.0, a22*vv/20.0, a20*vv/20.0, a21*vv/20.0, a22*vv/20.0, a20*vv/10.0, a21*vv/10.0, a22*vv/10.0],
        ])
    }

    /// Performs the g-d-g integration with constant tensor field (calculates the stiffness matrix)
    ///
    /// solution:
    ///
    /// ```text
    /// K = Bᵀ ⋅ D ⋅ B ⋅ volume
    /// ```
    pub fn integ_gdg(&mut self, young: f64, poisson: f64) -> Result<Matrix, StrError> {
        let ela = LinElasticity::new(young, poisson, false, false);
        let dd = ela.get_modulus();
        let dim_dd = 6;
        let dim_kk = 12;
        let mut bb_t_dd = Matrix::new(dim_kk, dim_dd);
        let mut kk = Matrix::new(dim_kk, dim_kk);
        mat_t_mat_mul(&mut bb_t_dd, 1.0, &self.bb, &dd.mat)?;
        mat_mat_mul(&mut kk, self.volume, &bb_t_dd, &self.bb)?;
        Ok(kk)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::AnalyticalTet4;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Matrix;

    #[test]
    fn analytical_tet4_works() {
        // unit tet4
        let space_ndim = 3;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Tet4).unwrap();
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(0, 2, 0.0);
        pad.set_xx(1, 0, 1.0);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(1, 2, 0.0);
        pad.set_xx(2, 0, 0.0);
        pad.set_xx(2, 1, 1.0);
        pad.set_xx(2, 2, 0.0);
        pad.set_xx(3, 0, 0.0);
        pad.set_xx(3, 1, 0.0);
        pad.set_xx(3, 2, 1.0);
        let mut tet = AnalyticalTet4::new(&pad);
        pad.calc_gradient(&[0.1, 0.1, 0.1]).unwrap();
        assert_eq!(tet.volume, 1.0 / 6.0);
        // println!("gg=\n{}", tet.gg);
        // println!("gradient=\n{}", state.gradient);
        assert_vec_approx_eq!(tet.gg.as_data(), pad.gradient.as_data(), 1e-15);
        let ee = 480.0;
        let nu = 1.0 / 3.0;
        let eb = ee / (12.0 * (1.0 - 2.0 * nu) * (1.0 + nu));
        let nt = eb * (1.0 - 2.0 * nu);
        let fsn = eb * (4.0 - 6.0 * nu);
        let tnu = eb * 2.0 * nu;
        let tnh = eb * 2.0 * (1.0 - nu);
        #[rustfmt::skip]
        let kk_correct = Matrix::from(&[
            [ fsn,   eb,   eb, -tnh, -nt, -nt, -nt, -tnu, 0.0, -nt, 0.0, -tnu],
            [  eb,  fsn,   eb, -tnu, -nt, 0.0, -nt, -tnh, -nt, 0.0, -nt, -tnu],
            [  eb,   eb,  fsn, -tnu, 0.0, -nt, 0.0, -tnu, -nt, -nt, -nt, -tnh],
            [-tnh, -tnu, -tnu,  tnh, 0.0, 0.0, 0.0,  tnu, 0.0, 0.0, 0.0,  tnu],
            [ -nt,  -nt,  0.0,  0.0,  nt, 0.0,  nt,  0.0, 0.0, 0.0, 0.0,  0.0],
            [ -nt,  0.0,  -nt,  0.0, 0.0,  nt, 0.0,  0.0, 0.0,  nt, 0.0,  0.0],
            [ -nt,  -nt,  0.0,  0.0,  nt, 0.0,  nt,  0.0, 0.0, 0.0, 0.0,  0.0],
            [-tnu, -tnh, -tnu,  tnu, 0.0, 0.0, 0.0,  tnh, 0.0, 0.0, 0.0,  tnu],
            [ 0.0,  -nt,  -nt,  0.0, 0.0, 0.0, 0.0,  0.0,  nt, 0.0,  nt,  0.0],
            [ -nt,  0.0,  -nt,  0.0, 0.0,  nt, 0.0,  0.0, 0.0,  nt, 0.0,  0.0],
            [ 0.0,  -nt,  -nt,  0.0, 0.0, 0.0, 0.0,  0.0,  nt, 0.0,  nt,  0.0],
            [-tnu, -tnu, -tnh,  tnu, 0.0, 0.0, 0.0,  tnu, 0.0, 0.0, 0.0,  tnh],
        ]);
        let kk = tet.integ_gdg(ee, nu).unwrap();
        assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), 1e-14);

        // non-right-angles tet4
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Tet4).unwrap();
        pad.set_xx(0, 0, 2.0);
        pad.set_xx(0, 1, 3.0);
        pad.set_xx(0, 2, 4.0);
        pad.set_xx(1, 0, 6.0);
        pad.set_xx(1, 1, 3.0);
        pad.set_xx(1, 2, 2.0);
        pad.set_xx(2, 0, 2.0);
        pad.set_xx(2, 1, 5.0);
        pad.set_xx(2, 2, 1.0);
        pad.set_xx(3, 0, 4.0);
        pad.set_xx(3, 1, 3.0);
        pad.set_xx(3, 2, 6.0);
        let mut tet = AnalyticalTet4::new(&pad);
        pad.calc_gradient(&[0.1, 0.2, 0.3]).unwrap();
        assert_eq!(tet.volume, 4.0);
        // println!("gg=\n{}", tet.gg);
        // println!("gradient=\n{}", state.gradient);
        assert_vec_approx_eq!(tet.gg.as_data(), pad.gradient.as_data(), 1e-15);
        let kk = tet.integ_gdg(ee, nu).unwrap();
        #[rustfmt::skip]
        let kk_correct = Matrix::from(&[
            [ 745.0,  540.0, 120.0,  -5.0,  30.0,  60.0,-270.0, -240.0,   0.0,-470.0, -330.0,-180.0],
            [ 540.0, 1720.0, 270.0,-120.0, 520.0, 210.0,-120.0,-1080.0, -60.0,-300.0,-1160.0,-420.0],
            [ 120.0,  270.0, 565.0,   0.0, 150.0, 175.0,   0.0, -120.0,-270.0,-120.0, -300.0,-470.0],
            [  -5.0, -120.0,   0.0, 145.0, -90.0, -60.0, -90.0,  120.0,   0.0, -50.0,   90.0,  60.0],
            [  30.0,  520.0, 150.0, -90.0, 220.0,  90.0,  60.0, -360.0, -60.0,   0.0, -380.0,-180.0],
            [  60.0,  210.0, 175.0, -60.0,  90.0, 145.0,   0.0, -120.0, -90.0,   0.0, -180.0,-230.0],
            [-270.0, -120.0,   0.0, -90.0,  60.0,   0.0, 180.0,    0.0,   0.0, 180.0,   60.0,   0.0],
            [-240.0,-1080.0,-120.0, 120.0,-360.0,-120.0,   0.0,  720.0,   0.0, 120.0,  720.0, 240.0],
            [   0.0,  -60.0,-270.0,   0.0, -60.0, -90.0,   0.0,    0.0, 180.0,   0.0,  120.0, 180.0],
            [-470.0, -300.0,-120.0, -50.0,   0.0,   0.0, 180.0,  120.0,   0.0, 340.0,  180.0, 120.0],
            [-330.0,-1160.0,-300.0,  90.0,-380.0,-180.0,  60.0,  720.0, 120.0, 180.0,  820.0, 360.0],
            [-180.0, -420.0,-470.0,  60.0,-180.0,-230.0,   0.0,  240.0, 180.0, 120.0,  360.0, 520.0],
        ]);
        assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), 1e-12);
    }
}
