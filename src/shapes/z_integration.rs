use super::{Shape, ShapeState};
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_tensor::{Tensor2, Tensor4};

impl Shape {
    /// Implements the shape(N)-scalar(S) integration case
    ///
    /// Interpolation functions times scalar field:
    ///
    /// ```text
    ///      ⌠    → →     →
    /// aᵐ = │ Nᵐ(x(ξ)) s(x) tₕ dΩ
    ///      ⌡
    ///      Ωₑ
    /// ```
    ///
    /// or, for lines in multi-dimensions:
    ///
    /// ```text
    ///      ⌠
    /// aᵐ = │ Nᵐ(ℓ(ξ)) s(ℓ) tₕ dℓ
    ///      ⌡
    ///      Γₑ
    /// ```
    ///
    /// The numerical integration is:
    ///
    /// ```text
    ///      nip-1     →     →          →
    /// aᵐ ≈   Σ    Nᵐ(ιᵖ) s(ιᵖ) tₕ |J|(ιᵖ) wᵖ
    ///       p=0
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
    ///          sequentially placed as shown above. `m` is the index of the node.
    ///          The length of `a` must be equal to `nnode`.
    ///
    /// # Input
    ///
    /// * `fn_s(index: usize) -> f64` -- s(x(ξ)) or s(ℓ) scalar function, however written as
    ///                                  a function of the index of the integration point.
    /// * `th` -- the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
    pub fn integ_vec_a_ns<F>(&self, state: &mut ShapeState, a: &mut Vector, fn_s: F, th: f64) -> Result<(), StrError>
    where
        F: Fn(usize) -> f64,
    {
        // check
        if a.dim() != self.nnode {
            return Err("the length of vector 'a' must be equal to nnode");
        }

        // clear output vector
        a.fill(0.0);

        // loop over integration points
        for index in 0..state.ip_data.len() {
            // ksi coordinates and weight
            let iota = &state.ip_data[index];
            let weight = state.ip_data[index][3];

            // calculate interpolation functions and Jacobian
            self.calc_interp(state, iota);
            let det_jac = self.calc_jacobian(state, iota)?;

            // calculate s
            let s = fn_s(index);

            // loop over nodes and perform sum
            let coef = th * det_jac * weight;
            for m in 0..self.nnode {
                a[m] += state.interp[m] * s * coef;
            }
        }
        Ok(())
    }

    /// Implements the the shape(N)-vector(V) integration case
    ///
    /// Interpolation functions times vector field:
    ///
    /// ```text
    /// →    ⌠    → →   → →
    /// bᵐ = │ Nᵐ(x(ξ)) v(x) tₕ dΩ
    ///      ⌡
    ///      Ωₑ
    /// ```
    ///
    /// The numerical integration is:
    ///
    /// ```text
    /// →    nip-1     →   → →          →
    /// bᵐ ≈   Σ    Nᵐ(ιᵖ) v(ιᵖ) tₕ |J|(ιᵖ) wᵖ
    ///       p=0
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
    ///     | bᵐᵢ |  ⟸  ii := i + m * space_ndim
    ///     └     ┘       
    ///
    /// m = ii / space_ndim
    /// i = ii % space_ndim
    /// ```
    ///
    /// * `b` -- A vector containing all `bᵐᵢ` values, one after another, and sequentially placed
    ///          as shown above (in 2D). `m` is the index of the node and `i` corresponds to `space_ndim`.
    ///          The length of `b` must be equal to `nnode * space_ndim`.
    ///
    /// # Input
    ///
    /// * `fn_v(v: &mut Vector, index: usize)` -- v(x(ξ)) vector function with `v.dim() == space_ndim`, however written as
    ///                                           a function of the index of the integration point.
    /// * `v_aux` -- is an auxiliary vector with size equal to `space_ndim`.
    /// * `th` -- the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
    pub fn integ_vec_b_nv<F>(
        &self,
        state: &mut ShapeState,
        b: &mut Vector,
        fn_v: F,
        v_aux: &mut Vector,
        th: f64,
    ) -> Result<(), StrError>
    where
        F: Fn(&mut Vector, usize),
    {
        // check
        if b.dim() != self.nnode * self.space_ndim {
            return Err("the length of vector 'b' must be equal to nnode * space_ndim");
        }
        if v_aux.dim() != self.space_ndim {
            return Err("the length of vector 'aux_v' must be equal to space_ndim");
        }

        // clear output vector
        b.fill(0.0);

        // loop over integration points
        for index in 0..state.ip_data.len() {
            // ksi coordinates and weight
            let iota = &state.ip_data[index];
            let weight = state.ip_data[index][3];

            // calculate interpolation functions and Jacobian
            self.calc_interp(state, iota);
            let det_jac = self.calc_jacobian(state, iota)?;

            // calculate v
            fn_v(v_aux, index);

            // add contribution to b vector
            let coef = th * det_jac * weight;
            self.add_to_vec_b(state, b, v_aux, coef);
        }
        Ok(())
    }

    /// Implements the vector(V)-gradient(G) integration case
    ///
    /// Vector dot gradient:
    ///
    /// ```text
    ///      ⌠ → →    →  → →
    /// cᵐ = │ w(x) · Gᵐ(x(ξ)) tₕ dΩ
    ///      ⌡
    ///      Ωₑ
    /// ```
    ///
    /// The numerical integration is:
    ///
    /// ```text
    ///      nip-1  → →     →  →          →
    /// cᵐ ≈   Σ    w(ιᵖ) · Gᵐ(ιᵖ) tₕ |J|(ιᵖ) wᵖ
    ///       p=0
    /// ```
    ///
    /// # Output
    ///
    /// ```text
    ///     ┌     ┐
    ///     |  c⁰ |
    ///     |  c¹ |
    /// c = |  c² |
    ///     | ··· |
    ///     |  cᵐ |
    ///     └     ┘
    /// ```
    ///
    /// * `c` -- A vector containing all `cᵐ` values, one after another, and
    ///          sequentially placed as shown above. `m` is the index of the node.
    ///          The length of `c` must be be equal to `nnode`.
    ///
    /// # Input
    ///
    /// * `fn_w(w: &mut Vector, index: usize)` -- w(x(ξ)) vector function with `w.dim() == space_ndim`, however written as
    ///                                           a function of the index of the integration point.
    /// * `w_aux` -- is an auxiliary vector with size equal to `space_ndim`.
    /// * `th` -- the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
    pub fn integ_vec_c_vg<F>(
        &self,
        state: &mut ShapeState,
        c: &mut Vector,
        fn_w: F,
        w_aux: &mut Vector,
        th: f64,
    ) -> Result<(), StrError>
    where
        F: Fn(&mut Vector, usize),
    {
        // check
        if c.dim() != self.nnode {
            return Err("the length of vector 'c' must be equal to nnode");
        }
        if w_aux.dim() != self.space_ndim {
            return Err("the length of vector 'aux_w' must be equal to space_ndim");
        }

        // clear output vector
        c.fill(0.0);

        // loop over integration points
        for index in 0..state.ip_data.len() {
            // ksi coordinates and weight
            let iota = &state.ip_data[index];
            let weight = state.ip_data[index][3];

            // calculate Jacobian and Gradient
            let det_jac = self.calc_gradient(state, iota)?;

            // calculate w
            fn_w(w_aux, index);

            // add contribution to c vector
            let coef = th * det_jac * weight;
            self.add_to_vec_c(state, c, w_aux, coef);
        }
        Ok(())
    }

    /// Implements the tensor(T)-gradient(G) integration case
    ///
    /// Tensor dot gradient:
    ///
    /// ```text
    /// →    ⌠   →    →  → →
    /// dᵐ = │ σ(x) · Gᵐ(x(ξ)) tₕ dΩ
    ///      ⌡ ▔
    ///      Ωₑ
    /// ```
    ///
    /// The numerical integration is:
    ///
    /// ```text
    /// →    nip-1    →     →  →          →
    /// dᵐ ≈   Σ    σ(ιᵖ) · Gᵐ(ιᵖ) tₕ |J|(ιᵖ) wᵖ
    ///       p=0   ▔
    /// ```
    ///
    /// # Output
    ///
    /// ```text
    ///     ┌     ┐
    ///     | d⁰₀ |
    ///     | d⁰₁ |
    ///     | d¹₀ |
    /// d = | d¹₁ |
    ///     | d²₀ |
    ///     | d²₁ |
    ///     | ··· |
    ///     | dᵐᵢ |  ⟸  ii := i + m * space_ndim
    ///     └     ┘
    ///
    /// m = ii / space_ndim
    /// i = ii % space_ndim
    /// ```
    ///
    /// * `d` -- A vector containing all `dᵐᵢ` values, one after another, and sequentially placed
    ///          as shown above (in 2D). `m` is the index of the node and `i` corresponds to `space_ndim`.
    ///          The length of `d` must be equal to `nnode * space_ndim`.
    ///
    /// # Input
    ///
    /// * `fn_sig(sig: &mut Tensor2, index: usize)` -- σ(x(ξ)) tensor function, however written as
    ///                                                a function of the index of the integration point.
    /// * `sig_aux` -- is an auxiliary Tensor2 with 4 components in 2D and 6 components in 3D.
    /// * `th` -- the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
    ///
    /// # Note
    ///
    /// This function is only available for space_ndim = 2D or 3D.
    pub fn integ_vec_d_tg<F>(
        &self,
        state: &mut ShapeState,
        d: &mut Vector,
        fn_sig: F,
        sig_aux: &mut Tensor2,
        th: f64,
    ) -> Result<(), StrError>
    where
        F: Fn(&mut Tensor2, usize),
    {
        // check
        if self.space_ndim == 1 {
            return Err("space_ndim must be 2 or 3");
        }
        if d.dim() != self.nnode * self.space_ndim {
            return Err("the length of vector 'd' must be equal to nnode * space_ndim");
        }
        if sig_aux.vec.dim() != 2 * self.space_ndim {
            return Err("'aux_sig' must be symmetric with dim equal to 4 in 2D or 6 in 3D");
        }

        // clear output vector
        d.fill(0.0);

        // loop over integration points
        for index in 0..state.ip_data.len() {
            // ksi coordinates and weight
            let iota = &state.ip_data[index];
            let weight = state.ip_data[index][3];

            // calculate Jacobian and Gradient
            let det_jac = self.calc_gradient(state, iota)?;

            // calculate σ
            fn_sig(sig_aux, index);

            // add contribution to d vector
            let coef = th * det_jac * weight;
            self.add_to_vec_d(state, d, sig_aux, coef);
        }
        Ok(())
    }

    /// Implements the shape(N)-scalar(S)-shape(N) integration case (e.g., diffusion matrix)
    pub fn integ_mat_1_nsn(&self, _: &mut ShapeState) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the gradient(G)-vector(V)-shape(N) integration case (e.g., compressibility matrix)
    pub fn integ_mat_2_gvn(&self, _: &mut ShapeState) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the gradient(G)-tensor(T)-gradient(G) integration case (e.g., conductivity matrix)
    pub fn integ_mat_3_gtg(&self, _: &mut ShapeState) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the shape(N)-scalar(S)-gradient(G) integration case with different shapes (e.g., coupling matrix)
    pub fn integ_mat_4_nsg(&self, _: &mut ShapeState, _: &mut Shape) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the gradient(G)-tensor(T)-shape(N) integration case with different shapes (e.g., coupling matrix)
    pub fn integ_mat_5_gtn(&self, _: &mut ShapeState, _: &mut Shape) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the shape(N)-vector(V)-shape(N) integration case with different shapes (e.g., coupling matrix)
    pub fn integ_mat_6_nvn(&self, _: &mut ShapeState, _: &mut Shape) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the gradient(G)-scalar(S)-shape(N) integration case with different shapes (e.g., coupling matrix)
    pub fn integ_mat_7_gsn(&self, _: &mut ShapeState, _: &mut Shape) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the shape(N)-tensor(T)-shape(N) integration case (e.g., mass matrix)
    pub fn integ_mat_8_ntn(&self, _: &mut ShapeState) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the shape(n)-vector(v)-gradient(g) integration case (e.g., variable density matrix)
    pub fn integ_mat_9_nvg(&self, _: &mut ShapeState) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the gradient(g)-4th-tensor(d)-gradient(g) integration case (e.g., stiffness matrix)
    ///
    /// Stiffness tensors:
    ///
    /// ```text
    ///       ⌠               →    →
    /// Kᵐⁿ = │ Gᵐₖ Dᵢₖⱼₗ Gⁿₗ eᵢ ⊗ eⱼ tₕ dΩ
    /// ▔     ⌡
    ///       Ωₑ
    /// ```
    ///
    /// The numerical integration is:
    ///
    /// ```text
    ///         nip-1     →         →       →          →
    /// Kᵐⁿᵢⱼ ≈   Σ   Gᵐₖ(ιᵖ) Dᵢₖⱼₗ(ιᵖ) Gⁿₗ(ιᵖ) tₕ |J|(ιᵖ) wᵖ
    ///          p=0
    /// ```
    ///
    /// # Output
    ///
    /// ```text
    ///     ┌                                               ┐
    ///     | K⁰⁰₀₀ K⁰⁰₀₁ K⁰¹₀₀ K⁰¹₀₁ K⁰²₀₀ K⁰²₀₁ ··· K⁰ⁿ₀ⱼ |
    ///     | K⁰⁰₁₀ K⁰⁰₁₁ K⁰¹₁₀ K⁰¹₁₁ K⁰²₁₀ K⁰²₁₁ ··· K⁰ⁿ₁ⱼ |
    ///     | K¹⁰₀₀ K¹⁰₀₁ K¹¹₀₀ K¹¹₀₁ K¹²₀₀ K¹²₀₁ ··· K¹ⁿ₀ⱼ |
    /// K = | K¹⁰₁₀ K¹⁰₁₁ K¹¹₁₀ K¹¹₁₁ K¹²₁₀ K¹²₁₁ ··· K¹ⁿ₁ⱼ |
    ///     | K²⁰₀₀ K²⁰₀₁ K²¹₀₀ K²¹₀₁ K²²₀₀ K²²₀₁ ··· K²ⁿ₀ⱼ |
    ///     | K²⁰₁₀ K²⁰₁₁ K²¹₁₀ K²¹₁₁ K²²₁₀ K²²₁₁ ··· K²ⁿ₁ⱼ |
    ///     |  ···   ···   ···   ···   ···   ···  ···  ···  |
    ///     | Kᵐ⁰ᵢ₀ Kᵐ⁰ᵢ₁ Kᵐ¹ᵢ₀ Kᵐ¹ᵢ₁ Kᵐ²ᵢ₀ Kᵐ²ᵢ₁ ··· Kᵐⁿᵢⱼ |  ⟸  ii := i + m * space_ndim
    ///     └                                               ┘
    ///                                                 ⇑
    ///                                                 jj := j + n * space_ndim
    ///
    /// m = ii / space_ndim    n = jj / space_ndim
    /// i = ii % space_ndim    j = jj % space_ndim
    /// ```
    ///
    /// * `kk` -- A matrix containing all `Kᵐⁿᵢⱼ` values, one after another, and sequentially placed
    ///           as shown above (in 2D). `m` and `n` are the indices of the node and `i` and `j`
    ///           correspond to `space_ndim`. The dimension of `K` must be equal to
    ///           (`nnode * space_ndim`, `nnode * space_ndim`).
    ///
    /// # Input
    ///
    /// * `fn_dd(dd: &mut Tensor4, index: usize)` -- D(x(ξ)) constitutive modulus function, given as
    ///                                              a function of the index of the integration point.
    /// * `dd_aux` -- is an auxiliary Tensor4 (minor-symmetric in 2D or 3D with 4 or 6 components, respectively).
    /// * `th` -- the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
    ///
    /// # Note
    ///
    /// This function is only available for space_ndim = 2D or 3D.
    pub fn integ_mat_10_gdg<F>(
        &self,
        state: &mut ShapeState,
        kk: &mut Matrix,
        fn_dd: F,
        dd_aux: &mut Tensor4,
        th: f64,
    ) -> Result<(), StrError>
    where
        F: Fn(&mut Tensor4, usize),
    {
        // check
        let (nrow_kk, ncol_kk) = kk.dims();
        let (nrow_dd, ncol_dd) = dd_aux.mat.dims();
        if self.space_ndim == 1 {
            return Err("space_ndim must be 2 or 3");
        }
        if nrow_kk != ncol_kk || nrow_kk != self.nnode * self.space_ndim {
            return Err("'K' matrix must be square with dim equal to nnode * space_ndim");
        }
        if nrow_dd != ncol_dd || nrow_dd != 2 * self.space_ndim {
            return Err("'D' tensor must be symmetric with dim equal to 4 in 2D or 6 in 3D");
        }

        // clear output matrix
        kk.fill(0.0);

        // loop over integration points
        for index in 0..state.ip_data.len() {
            // ksi coordinates and weight
            let iota = &state.ip_data[index];
            let weight = state.ip_data[index][3];

            // calculate Jacobian and Gradient
            let det_jac = self.calc_gradient(state, iota)?;

            // calculate constitutive modulus
            fn_dd(dd_aux, index);

            // add contribution to K matrix
            let coef = det_jac * weight * th;
            self.add_to_mat_kk(state, kk, dd_aux, coef);
        }
        Ok(())
    }

    /// Adds contribution to the b-vector in integ_vec_2_nv
    #[inline]
    fn add_to_vec_b(&self, state: &mut ShapeState, b: &mut Vector, v: &Vector, coef: f64) {
        if self.space_ndim == 1 {
            for m in 0..self.nnode {
                b[m] += coef * state.interp[m] * v[0];
            }
        } else if self.space_ndim == 2 {
            for m in 0..self.nnode {
                b[0 + m * 2] += coef * state.interp[m] * v[0];
                b[1 + m * 2] += coef * state.interp[m] * v[1];
            }
        } else {
            for m in 0..self.nnode {
                b[0 + m * 3] += coef * state.interp[m] * v[0];
                b[1 + m * 3] += coef * state.interp[m] * v[1];
                b[2 + m * 3] += coef * state.interp[m] * v[2];
            }
        }
    }

    /// Adds contribution to the c-vector in integ_vec_3_vg
    #[inline]
    fn add_to_vec_c(&self, state: &mut ShapeState, c: &mut Vector, w: &Vector, coef: f64) {
        let g = &state.gradient;
        if self.space_ndim == 1 {
            for m in 0..self.nnode {
                c[m] += coef * w[0] * g[m][0];
            }
        } else if self.space_ndim == 2 {
            for m in 0..self.nnode {
                c[m] += coef * (w[0] * g[m][0] + w[1] * g[m][1]);
            }
        } else {
            for m in 0..self.nnode {
                c[m] += coef * (w[0] * g[m][0] + w[1] * g[m][1] + w[2] * g[m][2]);
            }
        }
    }

    /// Adds contribution to the d-vector in integ_vec_4_tg
    #[inline]
    fn add_to_vec_d(&self, state: &mut ShapeState, d: &mut Vector, sig: &Tensor2, coef: f64) {
        let s = SQRT_2;
        let g = &state.gradient;
        let t = &sig.vec;
        if self.space_ndim == 2 {
            for m in 0..self.nnode {
                d[0 + m * 2] += coef * (t[0] * g[m][0] + t[3] * g[m][1] / s);
                d[1 + m * 2] += coef * (t[3] * g[m][0] / s + t[1] * g[m][1]);
            }
        } else {
            for m in 0..self.nnode {
                d[0 + m * 3] += coef * (t[0] * g[m][0] + t[3] * g[m][1] / s + t[5] * g[m][2] / s);
                d[1 + m * 3] += coef * (t[3] * g[m][0] / s + t[1] * g[m][1] + t[4] * g[m][2] / s);
                d[2 + m * 3] += coef * (t[5] * g[m][0] / s + t[4] * g[m][1] / s + t[2] * g[m][2]);
            }
        }
    }

    /// Adds contribution to the K-matrix in integ_mat_10_gdg
    #[inline]
    #[rustfmt::skip]
    fn add_to_mat_kk(&self, state: &mut ShapeState, kk: &mut Matrix, dd: &Tensor4, c: f64) {
        let s = SQRT_2;
        let g = &state.gradient;
        let d = &dd.mat;
        if self.space_ndim == 2 {
            for m in 0..self.nnode {
                for n in 0..self.nnode {
                    kk[0+m*2][0+n*2] += c * (g[m][1]*g[n][1]*d[3][3] + s*g[m][1]*g[n][0]*d[3][0] + s*g[m][0]*g[n][1]*d[0][3] + 2.0*g[m][0]*g[n][0]*d[0][0]) / 2.0;
                    kk[0+m*2][1+n*2] += c * (g[m][1]*g[n][0]*d[3][3] + s*g[m][1]*g[n][1]*d[3][1] + s*g[m][0]*g[n][0]*d[0][3] + 2.0*g[m][0]*g[n][1]*d[0][1]) / 2.0;
                    kk[1+m*2][0+n*2] += c * (g[m][0]*g[n][1]*d[3][3] + s*g[m][0]*g[n][0]*d[3][0] + s*g[m][1]*g[n][1]*d[1][3] + 2.0*g[m][1]*g[n][0]*d[1][0]) / 2.0;
                    kk[1+m*2][1+n*2] += c * (g[m][0]*g[n][0]*d[3][3] + s*g[m][0]*g[n][1]*d[3][1] + s*g[m][1]*g[n][0]*d[1][3] + 2.0*g[m][1]*g[n][1]*d[1][1]) / 2.0;
                }
            }
        } else {
            for m in 0..self.nnode {
                for n in 0..self.nnode {
                    kk[0+m*3][0+n*3] += c * (g[m][2]*g[n][2]*d[5][5] + g[m][2]*g[n][1]*d[5][3] + s*g[m][2]*g[n][0]*d[5][0] + g[m][1]*g[n][2]*d[3][5] + g[m][1]*g[n][1]*d[3][3] + s*g[m][1]*g[n][0]*d[3][0] + s*g[m][0]*g[n][2]*d[0][5] + s*g[m][0]*g[n][1]*d[0][3] + 2.0*g[m][0]*g[n][0]*d[0][0]) / 2.0;
                    kk[0+m*3][1+n*3] += c * (g[m][2]*g[n][2]*d[5][4] + g[m][2]*g[n][0]*d[5][3] + s*g[m][2]*g[n][1]*d[5][1] + g[m][1]*g[n][2]*d[3][4] + g[m][1]*g[n][0]*d[3][3] + s*g[m][1]*g[n][1]*d[3][1] + s*g[m][0]*g[n][2]*d[0][4] + s*g[m][0]*g[n][0]*d[0][3] + 2.0*g[m][0]*g[n][1]*d[0][1]) / 2.0;
                    kk[0+m*3][2+n*3] += c * (g[m][2]*g[n][0]*d[5][5] + g[m][2]*g[n][1]*d[5][4] + s*g[m][2]*g[n][2]*d[5][2] + g[m][1]*g[n][0]*d[3][5] + g[m][1]*g[n][1]*d[3][4] + s*g[m][1]*g[n][2]*d[3][2] + s*g[m][0]*g[n][0]*d[0][5] + s*g[m][0]*g[n][1]*d[0][4] + 2.0*g[m][0]*g[n][2]*d[0][2]) / 2.0;
                    kk[1+m*3][0+n*3] += c * (g[m][2]*g[n][2]*d[4][5] + g[m][2]*g[n][1]*d[4][3] + s*g[m][2]*g[n][0]*d[4][0] + g[m][0]*g[n][2]*d[3][5] + g[m][0]*g[n][1]*d[3][3] + s*g[m][0]*g[n][0]*d[3][0] + s*g[m][1]*g[n][2]*d[1][5] + s*g[m][1]*g[n][1]*d[1][3] + 2.0*g[m][1]*g[n][0]*d[1][0]) / 2.0;
                    kk[1+m*3][1+n*3] += c * (g[m][2]*g[n][2]*d[4][4] + g[m][2]*g[n][0]*d[4][3] + s*g[m][2]*g[n][1]*d[4][1] + g[m][0]*g[n][2]*d[3][4] + g[m][0]*g[n][0]*d[3][3] + s*g[m][0]*g[n][1]*d[3][1] + s*g[m][1]*g[n][2]*d[1][4] + s*g[m][1]*g[n][0]*d[1][3] + 2.0*g[m][1]*g[n][1]*d[1][1]) / 2.0;
                    kk[1+m*3][2+n*3] += c * (g[m][2]*g[n][0]*d[4][5] + g[m][2]*g[n][1]*d[4][4] + s*g[m][2]*g[n][2]*d[4][2] + g[m][0]*g[n][0]*d[3][5] + g[m][0]*g[n][1]*d[3][4] + s*g[m][0]*g[n][2]*d[3][2] + s*g[m][1]*g[n][0]*d[1][5] + s*g[m][1]*g[n][1]*d[1][4] + 2.0*g[m][1]*g[n][2]*d[1][2]) / 2.0;
                    kk[2+m*3][0+n*3] += c * (g[m][0]*g[n][2]*d[5][5] + g[m][0]*g[n][1]*d[5][3] + s*g[m][0]*g[n][0]*d[5][0] + g[m][1]*g[n][2]*d[4][5] + g[m][1]*g[n][1]*d[4][3] + s*g[m][1]*g[n][0]*d[4][0] + s*g[m][2]*g[n][2]*d[2][5] + s*g[m][2]*g[n][1]*d[2][3] + 2.0*g[m][2]*g[n][0]*d[2][0]) / 2.0;
                    kk[2+m*3][1+n*3] += c * (g[m][0]*g[n][2]*d[5][4] + g[m][0]*g[n][0]*d[5][3] + s*g[m][0]*g[n][1]*d[5][1] + g[m][1]*g[n][2]*d[4][4] + g[m][1]*g[n][0]*d[4][3] + s*g[m][1]*g[n][1]*d[4][1] + s*g[m][2]*g[n][2]*d[2][4] + s*g[m][2]*g[n][0]*d[2][3] + 2.0*g[m][2]*g[n][1]*d[2][1]) / 2.0;
                    kk[2+m*3][2+n*3] += c * (g[m][0]*g[n][0]*d[5][5] + g[m][0]*g[n][1]*d[5][4] + s*g[m][0]*g[n][2]*d[5][2] + g[m][1]*g[n][0]*d[4][5] + g[m][1]*g[n][1]*d[4][4] + s*g[m][1]*g[n][2]*d[4][2] + s*g[m][2]*g[n][0]*d[2][5] + s*g[m][2]*g[n][1]*d[2][4] + 2.0*g[m][2]*g[n][2]*d[2][2]) / 2.0;
                }
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
    use russell_lab::{copy_matrix, mat_mat_mul, mat_t_mat_mul, Matrix};
    use russell_tensor::LinElasticity;

    // to test if variables are cleared before sum
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
        shape.set_node(0, 0, xmin).unwrap();
        shape.set_node(0, 1, ymin).unwrap();
        shape.set_node(1, 0, xmin + l).unwrap();
        shape.set_node(1, 1, ymin).unwrap();
        shape.set_node(2, 0, xmin + l / 2.0).unwrap();
        shape.set_node(2, 1, ymin + h).unwrap();
        (shape, area)
    }

    // line segment from xa to xb (geo_ndim == space_ndim)
    fn gen_lin2() -> (Shape, f64, f64) {
        let mut shape = Shape::new(1, 1, 2).unwrap();
        let (xa, xb) = (3.0, 9.0);
        shape.set_node(0, 0, xa).unwrap();
        shape.set_node(1, 0, xb).unwrap();
        (shape, xa, xb)
    }

    struct AnalyticalTri3 {
        x: [f64; 3], // node x-coordinates
        y: [f64; 3], // node y-coordinates
        b: [f64; 3], // b-coefficients
        c: [f64; 3], // c-coefficients
        area: f64,   // area
    }

    impl AnalyticalTri3 {
        pub fn new(tri3: &Shape) -> Self {
            // coefficients
            let (x0, y0) = (tri3.coords_transp[0][0], tri3.coords_transp[1][0]);
            let (x1, y1) = (tri3.coords_transp[0][1], tri3.coords_transp[1][1]);
            let (x2, y2) = (tri3.coords_transp[0][2], tri3.coords_transp[1][2]);
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
            let mut state = ShapeState::new(&tri3);
            let ksi = &state.ip_data[0];
            tri3.calc_gradient(&mut state, ksi).unwrap();
            assert_eq!(state.gradient.as_data(), gg.as_data());

            // results
            AnalyticalTri3 {
                x: [x0, x1, x2],
                y: [y0, y1, y2],
                b: [b0, b1, b2],
                c: [c0, c1, c2],
                area,
            }
        }
    }

    #[test]
    fn integ_vec_a_works() -> Result<(), StrError> {
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
        let (tri3, area) = gen_tri3();
        let mut state = ShapeState::new(&tri3);
        const CS: f64 = 3.0;
        let fn_s = |_| CS;
        let mut a = Vector::filled(tri3.nnode, NOISE);
        tri3.integ_vec_a_ns(&mut state, &mut a, fn_s, 1.0)?;
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
        let (lin2, xa, xb) = gen_lin2();
        let mut state = ShapeState::new(&lin2);
        let all_int_points = lin2.calc_int_points_coords(&mut state)?;
        let fn_s = |index: usize| all_int_points[index][0];
        let mut a = Vector::new(lin2.nnode);
        lin2.integ_vec_a_ns(&mut state, &mut a, fn_s, 1.0)?;
        let cf = (xb - xa) / 6.0;
        let a_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];
        assert_vec_approx_eq!(a.as_data(), a_correct, 1e-15);
        Ok(())
    }

    #[test]
    fn integ_vec_b_works() -> Result<(), StrError> {
        // This test is similar to the case_a with tri3, however using a vector
        // So, each component of `b` equals `Fₛ`
        let (tri3, area) = gen_tri3();
        let mut state = ShapeState::new(&tri3);
        const CS: f64 = 3.0;
        let fn_v = |v: &mut Vector, _: usize| v.fill(CS);
        let mut b = Vector::filled(tri3.nnode * tri3.space_ndim, NOISE);
        let mut v_aux = Vector::new(tri3.space_ndim);
        tri3.integ_vec_b_nv(&mut state, &mut b, fn_v, &mut v_aux, 1.0)?;
        let cf = CS * area / 3.0;
        let b_correct = &[cf, cf, cf, cf, cf, cf];
        assert_vec_approx_eq!(b.as_data(), b_correct, 1e-14);

        // Likewise, this test is similar to case_a with lin2, however using a vector
        // with a single component. So, each component of `b` equals `Fₛ`
        let (lin2, xa, xb) = gen_lin2();
        let mut state = ShapeState::new(&lin2);
        let all_int_points = lin2.calc_int_points_coords(&mut state)?;
        let fn_v = |v: &mut Vector, index: usize| {
            v.fill(all_int_points[index][0]);
        };
        let mut b = Vector::filled(lin2.nnode * lin2.space_ndim, NOISE);
        let mut v_aux = Vector::new(lin2.space_ndim);
        lin2.integ_vec_b_nv(&mut state, &mut b, fn_v, &mut v_aux, 1.0)?;
        let cf = (xb - xa) / 6.0;
        let b_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];
        assert_vec_approx_eq!(b.as_data(), b_correct, 1e-15);
        Ok(())
    }

    #[test]
    fn integ_vec_c_works() -> Result<(), StrError> {
        // shape and analytical gradient
        let (tri3, area) = gen_tri3();
        let mut state = ShapeState::new(&tri3);
        let ana = AnalyticalTri3::new(&tri3);
        assert_eq!(area, ana.area);

        // constant vector function: w(x) = {w₀, w₁}
        // solution:
        //    cᵐ = ½ (w₀ bₘ + w₁ cₘ)
        const W0: f64 = 2.0;
        const W1: f64 = 3.0;
        let fn_w = |w: &mut Vector, _: usize| {
            w[0] = W0;
            w[1] = W1;
        };
        let c_correct = &[
            (W0 * ana.b[0] + W1 * ana.c[0]) / 2.0,
            (W0 * ana.b[1] + W1 * ana.c[1]) / 2.0,
            (W0 * ana.b[2] + W1 * ana.c[2]) / 2.0,
        ];
        let mut c = Vector::filled(tri3.nnode, NOISE);
        let mut w_aux = Vector::new(tri3.space_ndim);
        tri3.integ_vec_c_vg(&mut state, &mut c, fn_w, &mut w_aux, 1.0)?;
        assert_vec_approx_eq!(c.as_data(), c_correct, 1e-15);

        // bilinear vector function: w(x) = {x, y}
        // solution:
        //    cᵐ = ⅙ bₘ (x₀+x₁+x₂) + ⅙ cₘ (y₀+y₁+y₂)
        let all_int_points = tri3.calc_int_points_coords(&mut state)?;
        let fn_w = |w: &mut Vector, index: usize| {
            w[0] = all_int_points[index][0];
            w[1] = all_int_points[index][1];
        };
        let c_correct = &[
            (ana.x[0] + ana.x[1] + ana.x[2]) * ana.b[0] / 6.0 + (ana.y[0] + ana.y[1] + ana.y[2]) * ana.c[0] / 6.0,
            (ana.x[0] + ana.x[1] + ana.x[2]) * ana.b[1] / 6.0 + (ana.y[0] + ana.y[1] + ana.y[2]) * ana.c[1] / 6.0,
            (ana.x[0] + ana.x[1] + ana.x[2]) * ana.b[2] / 6.0 + (ana.y[0] + ana.y[1] + ana.y[2]) * ana.c[2] / 6.0,
        ];
        let mut c = Vector::filled(tri3.nnode, NOISE);
        let mut w_aux = Vector::new(tri3.space_ndim);
        tri3.integ_vec_c_vg(&mut state, &mut c, fn_w, &mut w_aux, 1.0)?;
        assert_vec_approx_eq!(c.as_data(), c_correct, 1e-14);
        Ok(())
    }

    #[test]
    fn integ_vec_d_works() -> Result<(), StrError> {
        // shape and analytical gradient
        let (tri3, _) = gen_tri3();
        let mut state = ShapeState::new(&tri3);
        let ana = AnalyticalTri3::new(&tri3);

        // constant tensor function: σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2}
        // solution:
        //    dᵐ₀ = ½ (σ₀₀ bₘ + σ₀₁ cₘ)
        //    dᵐ₁ = ½ (σ₁₀ bₘ + σ₁₁ cₘ)
        const S00: f64 = 2.0;
        const S11: f64 = 3.0;
        const S22: f64 = 4.0;
        const S01: f64 = 5.0;
        let fn_sig = |sig: &mut Tensor2, _: usize| {
            sig.sym_set(0, 0, S00);
            sig.sym_set(1, 1, S11);
            sig.sym_set(2, 2, S22);
            sig.sym_set(0, 1, S01);
        };
        let d_correct = &[
            (S00 * ana.b[0] + S01 * ana.c[0]) / 2.0,
            (S01 * ana.b[0] + S11 * ana.c[0]) / 2.0,
            (S00 * ana.b[1] + S01 * ana.c[1]) / 2.0,
            (S01 * ana.b[1] + S11 * ana.c[1]) / 2.0,
            (S00 * ana.b[2] + S01 * ana.c[2]) / 2.0,
            (S01 * ana.b[2] + S11 * ana.c[2]) / 2.0,
        ];
        let mut d = Vector::filled(tri3.nnode * tri3.space_ndim, NOISE);
        let mut sig_aux = Tensor2::new(true, true);
        tri3.integ_vec_d_tg(&mut state, &mut d, fn_sig, &mut sig_aux, 1.0)?;
        assert_vec_approx_eq!(d.as_data(), d_correct, 1e-15);
        Ok(())
    }

    #[test]
    fn integ_mat_10_gdg_works() -> Result<(), StrError> {
        /* Element # 0 from example 1.6 from [@bhatti] page 32

         Solid bracket with thickness = 0.25

                     1     -10                connectivity:
        y=2.0 (-100) o'-,__                    eid : vertices
                     |     '-,__ 3   -10         0 :  0, 2, 3
        y=1.5 - - -  |        ,'o-,__            1 :  3, 1, 0
                     |  1   ,'  |    '-,__ 5     2 :  2, 4, 5
                     |    ,'    |  3   ,-'o      3 :  5, 3, 2
                     |  ,'  0   |   ,-'   |
                     |,'        |,-'   2  |   constraints:
        y=0.0 (-100) o----------o---------o    -100 : fixed on x and y
                     0          2         4
                    x=0.0     x=2.0     x=4.0

        # References

        [@bhatti] Bhatti, M.A. (2005) Fundamental Finite Element Analysis
                  and Applications, Wiley, 700p.
        */

        // shape and analytical gradient
        let mut tri3 = Shape::new(2, 2, 3).unwrap();
        let mut state = ShapeState::new(&tri3);
        tri3.set_node(0, 0, 0.0).unwrap();
        tri3.set_node(0, 1, 0.0).unwrap();
        tri3.set_node(1, 0, 2.0).unwrap();
        tri3.set_node(1, 1, 0.0).unwrap();
        tri3.set_node(2, 0, 2.0).unwrap();
        tri3.set_node(2, 1, 1.5).unwrap();
        let ana = AnalyticalTri3::new(&mut tri3);

        // elasticity modulus and function
        let ela = LinElasticity::new(10000.0, 0.2, false, true);
        let dd_ela = ela.get_modulus();
        let fn_dd = |dd: &mut Tensor4, _: usize| {
            copy_matrix(&mut dd.mat, &dd_ela.mat).unwrap();
        };

        // constants
        let th = 0.25; // thickness
        let dim_dd = 2 * tri3.space_ndim;
        let dim_kk = tri3.nnode * tri3.space_ndim;

        // compute B-matrix (dim_dd,dim_kk)
        let r = 2.0 * ana.area;
        let s = r * SQRT_2;
        #[rustfmt::skip]
        let bb = Matrix::from(&[
            [ana.b[0]/r,        0.0, ana.b[1]/r,        0.0, ana.b[2]/r,        0.0],
            [       0.0, ana.c[0]/r,        0.0, ana.c[1]/r,        0.0, ana.c[2]/r],
            [       0.0,        0.0,        0.0,        0.0,        0.0,        0.0],
            [ana.c[0]/s, ana.b[0]/s, ana.c[1]/s, ana.b[1]/s, ana.c[2]/s, ana.b[2]/s],
        ]);
        assert_eq!(bb.dims(), (dim_dd, dim_kk));

        // compute K = Bᵀ ⋅ D ⋅ B
        let mut bb_t_dd = Matrix::new(dim_kk, dim_dd);
        let mut kk_correct = Matrix::new(dim_kk, dim_kk);
        mat_t_mat_mul(&mut bb_t_dd, 1.0, &bb, &dd_ela.mat)?;
        mat_mat_mul(&mut kk_correct, th * ana.area, &bb_t_dd, &bb)?;

        // perform integration
        let mut kk = Matrix::new(dim_kk, dim_kk);
        let mut dd_aux = Tensor4::new(true, true);
        tri3.integ_mat_10_gdg(&mut state, &mut kk, fn_dd, &mut dd_aux, th)?;

        // results from Bhatti's book
        #[rustfmt::skip]
        let kk_bhatti = Matrix::from( &[
            [  9.765625000000001e+02,  0.000000000000000e+00, -9.765625000000001e+02,  2.604166666666667e+02,  0.000000000000000e+00, -2.604166666666667e+02],
            [  0.000000000000000e+00,  3.906250000000000e+02,  5.208333333333334e+02, -3.906250000000000e+02, -5.208333333333334e+02,  0.000000000000000e+00],
            [ -9.765625000000001e+02,  5.208333333333334e+02,  1.671006944444445e+03, -7.812500000000000e+02, -6.944444444444445e+02,  2.604166666666667e+02],
            [  2.604166666666667e+02, -3.906250000000000e+02, -7.812500000000000e+02,  2.126736111111111e+03,  5.208333333333334e+02, -1.736111111111111e+03],
            [  0.000000000000000e+00, -5.208333333333334e+02, -6.944444444444445e+02,  5.208333333333334e+02,  6.944444444444445e+02,  0.000000000000000e+00],
            [ -2.604166666666667e+02,  0.000000000000000e+00,  2.604166666666667e+02, -1.736111111111111e+03,  0.000000000000000e+00,  1.736111111111111e+03],
        ]);

        // check
        assert_vec_approx_eq!(kk_correct.as_data(), kk_bhatti.as_data(), 1e-12);
        assert_vec_approx_eq!(kk_correct.as_data(), kk.as_data(), 1e-12);
        Ok(())
    }
}
