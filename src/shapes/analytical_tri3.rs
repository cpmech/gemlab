use crate::shapes::{GeoKind, Shape, StateOfShape, IP_TRI_INTERNAL_1};
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::{mat_mat_mul, mat_t_mat_mul, Matrix};
use russell_tensor::LinElasticity;

/// Performs analytical integrations on a Tri3
pub struct AnalyticalTri3 {
    // Holds the x-coordinates of nodes
    pub x: [f64; 3],

    // Holds the y-coordinates of nodes
    pub y: [f64; 3],

    // Holds the b-coefficients
    pub b: [f64; 3],

    // Holds the c-coefficients
    pub c: [f64; 3],

    // Holds the area of the triangle
    pub area: f64,

    /// Holds the B-matrix
    pub bb: Matrix,
}

impl AnalyticalTri3 {
    pub fn new(shape: &Shape, state: &mut StateOfShape) -> Self {
        assert_eq!(shape.kind, GeoKind::Tri3);

        // coefficients
        let (x0, y0) = (state.coords_transp[0][0], state.coords_transp[1][0]);
        let (x1, y1) = (state.coords_transp[0][1], state.coords_transp[1][1]);
        let (x2, y2) = (state.coords_transp[0][2], state.coords_transp[1][2]);
        let (b0, b1, b2) = (y1 - y2, y2 - y0, y0 - y1);
        let (c0, c1, c2) = (x2 - x1, x0 - x2, x1 - x0);
        let (f0, f1, f2) = (x1 * y2 - x2 * y1, x2 * y0 - x0 * y2, x0 * y1 - x1 * y0);

        // area
        let area = (f0 + f1 + f2) / 2.0;

        // check gradients
        let gg = Matrix::from(&[
            [b0 / (2.0 * area), c0 / (2.0 * area)],
            [b1 / (2.0 * area), c1 / (2.0 * area)],
            [b2 / (2.0 * area), c2 / (2.0 * area)],
        ]);
        shape.calc_gradient(state, &IP_TRI_INTERNAL_1[0]).unwrap();
        assert_eq!(state.gradient.as_data(), gg.as_data());

        // B-matrix
        let r = 2.0 * area;
        let s = r * SQRT_2;
        #[rustfmt::skip]
        let bb = Matrix::from(&[
            [b0/r,  0.0, b1/r,  0.0, b2/r,  0.0],
            [ 0.0, c0/r,  0.0, c1/r,  0.0, c2/r],
            [ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
            [c0/s, b0/s, c1/s, b1/s, c2/s, b2/s],
        ]);

        // results
        AnalyticalTri3 {
            x: [x0, x1, x2],
            y: [y0, y1, y2],
            b: [b0, b1, b2],
            c: [c0, c1, c2],
            area,
            bb,
        }
    }

    /// Performs the integ_vec_c operation given constant vector function w(x) = {w₀, w₁}
    ///
    /// solution:
    /// ```text
    /// cᵐ = ½ (w₀ bₘ + w₁ cₘ)
    /// ```
    ///
    /// # Input
    ///
    /// * `(w0,w1)` -- components of a constant vector function: w(x) = {w₀, w₁}
    pub fn integ_vec_c_constant(&mut self, w0: f64, w1: f64) -> Vec<f64> {
        vec![
            (w0 * self.b[0] + w1 * self.c[0]) / 2.0,
            (w0 * self.b[1] + w1 * self.c[1]) / 2.0,
            (w0 * self.b[2] + w1 * self.c[2]) / 2.0,
        ]
    }

    /// Performs the integ_vec_c operation given bilinear vector function w(x) = {x, y}
    ///
    /// solution:
    /// ```text
    /// cᵐ = ⅙ bₘ (x₀+x₁+x₂) + ⅙ cₘ (y₀+y₁+y₂)
    /// ```
    pub fn integ_vec_c_bilinear(&mut self) -> Vec<f64> {
        let ana = self;
        vec![
            (ana.x[0] + ana.x[1] + ana.x[2]) * ana.b[0] / 6.0 + (ana.y[0] + ana.y[1] + ana.y[2]) * ana.c[0] / 6.0,
            (ana.x[0] + ana.x[1] + ana.x[2]) * ana.b[1] / 6.0 + (ana.y[0] + ana.y[1] + ana.y[2]) * ana.c[1] / 6.0,
            (ana.x[0] + ana.x[1] + ana.x[2]) * ana.b[2] / 6.0 + (ana.y[0] + ana.y[1] + ana.y[2]) * ana.c[2] / 6.0,
        ]
    }

    /// Performs the integ_vec_d operation given constant tensor function σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2}
    ///
    /// solution:
    /// ```text
    /// dᵐ₀ = ½ (σ₀₀ bₘ + σ₀₁ cₘ)
    /// dᵐ₁ = ½ (σ₁₀ bₘ + σ₁₁ cₘ)
    /// ```
    ///
    /// # Input
    ///
    /// * `s₀₀, s₁₁, s₀₁` -- components of a constant tensor function: σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2}
    pub fn integ_vec_d_constant(&mut self, s00: f64, s11: f64, s01: f64) -> Vec<f64> {
        let ana = self;
        vec![
            (s00 * ana.b[0] + s01 * ana.c[0]) / 2.0,
            (s01 * ana.b[0] + s11 * ana.c[0]) / 2.0,
            (s00 * ana.b[1] + s01 * ana.c[1]) / 2.0,
            (s01 * ana.b[1] + s11 * ana.c[1]) / 2.0,
            (s00 * ana.b[2] + s01 * ana.c[2]) / 2.0,
            (s01 * ana.b[2] + s11 * ana.c[2]) / 2.0,
        ]
    }

    /// Calculates the stiffness matrix
    ///
    /// solution:
    ///
    /// ```text
    /// K = Bᵀ ⋅ D ⋅ B ⋅ th ⋅ area
    /// ```
    pub fn integ_stiffness(
        &mut self,
        young: f64,
        poisson: f64,
        plane_stress: bool,
        th: f64,
    ) -> Result<Matrix, StrError> {
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
    use russell_chk::{assert_approx_eq, assert_vec_approx_eq};

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
        assert_eq!(ana.x, [0.0, 0.2, 0.1]);
        assert_eq!(ana.y, [0.0, 0.0, 0.1]);
        assert_approx_eq!(ana.area, 0.01, 1e-15);
        assert_vec_approx_eq!(ana.b, [-0.1, 0.1, 0.0], 1e-15);
        assert_vec_approx_eq!(ana.c, [-0.1, -0.1, 0.2], 1e-15);
        Ok(())
    }
}
