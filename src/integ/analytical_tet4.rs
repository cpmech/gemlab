use crate::shapes::{GeoKind, Shape, StateOfShape};
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::{mat_mat_mul, mat_t_mat_mul, Matrix};
use russell_tensor::LinElasticity;

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
    /// Organized as the G matrix (nnode=3, space_ndim=2)
    pub gg: Matrix,

    /// Holds the B-matrix (6, 12)
    pub bb: Matrix,
}

impl AnalyticalTet4 {
    pub fn new(shape: &Shape, state: &StateOfShape) -> Self {
        assert_eq!(shape.kind, GeoKind::Tet4);

        let x1 = state.coords_transp[0][0];
        let x2 = state.coords_transp[0][1];
        let x3 = state.coords_transp[0][2];
        let x4 = state.coords_transp[0][3];

        let y1 = state.coords_transp[1][0];
        let y2 = state.coords_transp[1][1];
        let y3 = state.coords_transp[1][2];
        let y4 = state.coords_transp[1][3];

        let z1 = state.coords_transp[2][0];
        let z2 = state.coords_transp[2][1];
        let z3 = state.coords_transp[2][2];
        let z4 = state.coords_transp[2][3];

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
    /// * `shape` -- The same shape used in `new` because we need the nodal coordinates here.
    ///              Do not change the coordinates, otherwise the values will be wrong.
    pub fn integ_vec_a_linear_along_z(&self, state: &StateOfShape) -> Vec<f64> {
        let (z1, z2, z3, z4) = (
            state.coords_transp[2][0],
            state.coords_transp[2][1],
            state.coords_transp[2][2],
            state.coords_transp[2][3],
        );
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
    pub fn integ_vec_b_constant(&self, v0: f64, v1: f64, v2: f64) -> Vec<f64> {
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
    pub fn integ_vec_c_constant(&self, w0: f64, w1: f64, w2: f64) -> Vec<f64> {
        vec![
            (w0 * self.gg[0][0] + w1 * self.gg[0][1] + w2 * self.gg[0][2]) * self.volume,
            (w0 * self.gg[1][0] + w1 * self.gg[1][1] + w2 * self.gg[1][2]) * self.volume,
            (w0 * self.gg[2][0] + w1 * self.gg[2][1] + w2 * self.gg[2][2]) * self.volume,
            (w0 * self.gg[3][0] + w1 * self.gg[3][1] + w2 * self.gg[3][2]) * self.volume,
        ]
    }

    /// Integrates tensor dot gradient with constant tensor function σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2, σ₁₂√2, σ₀₂√2}
    ///
    /// Solution:
    ///
    /// ```text
    ///    dᵐ₀ = (σ₀₀ Gᵐ₀ + σ₀₁ Gᵐ₁ + σ₀₂ Gᵐ₂) V
    ///    dᵐ₁ = (σ₁₀ Gᵐ₀ + σ₁₁ Gᵐ₁ + σ₁₂ Gᵐ₂) V
    ///    dᵐ₂ = (σ₂₀ Gᵐ₀ + σ₂₁ Gᵐ₁ + σ₂₂ Gᵐ₂) V
    /// ```
    pub fn integ_vec_d_constant(&self, s00: f64, s11: f64, s22: f64, s01: f64, s12: f64, s02: f64) -> Vec<f64> {
        vec![
            (s00 * self.gg[0][0] + s01 * self.gg[0][1] + s02 * self.gg[0][2]) * self.volume,
            (s01 * self.gg[0][0] + s11 * self.gg[0][1] + s12 * self.gg[0][2]) * self.volume,
            (s02 * self.gg[0][0] + s12 * self.gg[0][1] + s22 * self.gg[0][2]) * self.volume,
            (s00 * self.gg[1][0] + s01 * self.gg[1][1] + s02 * self.gg[1][2]) * self.volume,
            (s01 * self.gg[1][0] + s11 * self.gg[1][1] + s12 * self.gg[1][2]) * self.volume,
            (s02 * self.gg[1][0] + s12 * self.gg[1][1] + s22 * self.gg[1][2]) * self.volume,
            (s00 * self.gg[2][0] + s01 * self.gg[2][1] + s02 * self.gg[2][2]) * self.volume,
            (s01 * self.gg[2][0] + s11 * self.gg[2][1] + s12 * self.gg[2][2]) * self.volume,
            (s02 * self.gg[2][0] + s12 * self.gg[2][1] + s22 * self.gg[2][2]) * self.volume,
            (s00 * self.gg[3][0] + s01 * self.gg[3][1] + s02 * self.gg[3][2]) * self.volume,
            (s01 * self.gg[3][0] + s11 * self.gg[3][1] + s12 * self.gg[3][2]) * self.volume,
            (s02 * self.gg[3][0] + s12 * self.gg[3][1] + s22 * self.gg[3][2]) * self.volume,
        ]
    }

    /// Calculates the stiffness matrix
    ///
    /// solution:
    ///
    /// ```text
    /// K = Bᵀ ⋅ D ⋅ B ⋅ volume
    /// ```
    pub fn integ_stiffness(&mut self, young: f64, poisson: f64) -> Result<Matrix, StrError> {
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
    use crate::shapes::{GeoKind, Shape, StateOfShape};
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Matrix;

    #[test]
    fn analytical_tet4_works() -> Result<(), StrError> {
        // unit tet4
        let shape = Shape::new(GeoKind::Tet4);
        let mut state = StateOfShape::new(
            GeoKind::Tet4,
            &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        )
        .unwrap();
        let mut tet = AnalyticalTet4::new(&shape, &mut state);
        shape.calc_gradient(&mut state, &[0.1, 0.1, 0.1])?;
        assert_eq!(tet.volume, 1.0 / 6.0);
        // println!("gg=\n{}", tet.gg);
        // println!("gradient=\n{}", state.gradient);
        assert_vec_approx_eq!(tet.gg.as_data(), state.gradient.as_data(), 1e-15);
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
        let kk = tet.integ_stiffness(ee, nu)?;
        assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), 1e-14);

        // random tet4
        let mut state = StateOfShape::new(
            GeoKind::Tet4,
            &[[2.0, 3.0, 4.0], [6.0, 3.0, 2.0], [2.0, 5.0, 1.0], [4.0, 3.0, 6.0]],
        )
        .unwrap();
        let mut tet = AnalyticalTet4::new(&shape, &mut state);
        shape.calc_gradient(&mut state, &[0.1, 0.2, 0.3])?;
        assert_eq!(tet.volume, 4.0);
        // println!("gg=\n{}", tet.gg);
        // println!("gradient=\n{}", state.gradient);
        assert_vec_approx_eq!(tet.gg.as_data(), state.gradient.as_data(), 1e-15);
        let kk = tet.integ_stiffness(ee, nu)?;
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
        Ok(())
    }
}
