use crate::shapes::{GeoKind, Shape, StateOfShape};
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
    pub fn new(shape: &Shape, state: &StateOfShape) -> Self {
        assert_eq!(shape.kind, GeoKind::Tri3);

        // coefficients
        let (x0, y0) = (state.coords_transp[0][0], state.coords_transp[1][0]);
        let (x1, y1) = (state.coords_transp[0][1], state.coords_transp[1][1]);
        let (x2, y2) = (state.coords_transp[0][2], state.coords_transp[1][2]);
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
    /// * `shape` -- The same shape used in `new` because we need the nodal coordinates here
    ///              Do not change the coordinates, otherwise the values will be wrong.
    pub fn integ_vec_c_bilinear(&self, state: &StateOfShape) -> Vec<f64> {
        let (x0, x1, x2) = (
            state.coords_transp[0][0],
            state.coords_transp[0][1],
            state.coords_transp[0][2],
        );
        let (y0, y1, y2) = (
            state.coords_transp[1][0],
            state.coords_transp[1][1],
            state.coords_transp[1][2],
        );
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
    use crate::shapes::{Shape, StateOfShape};
    use crate::StrError;
    use russell_chk::assert_approx_eq;

    #[test]
    fn new_works() -> Result<(), StrError> {
        // Tri3 # 1 from Figure 1.18, page 29 of [@bhatti]
        //
        // [@bhatti] Bhatti, M.A. (2005) Fundamental Finite Element Analysis
        //          and Applications, Wiley, 700p.
        //
        let shape = Shape::new(2, 2, 3)?;
        let mut state = StateOfShape::new(shape.geo_ndim, &[[0.0, 0.0], [0.2, 0.0], [0.1, 0.1]])?;
        let ana = AnalyticalTri3::new(&shape, &mut state);
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
