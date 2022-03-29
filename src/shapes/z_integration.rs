use super::z_integ_points::*;
use super::{GeoClass, Shape};
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_tensor::{Tensor2, Tensor4};

impl Shape {
    /// Selects integrations points and weights
    ///
    /// # Options
    ///
    /// ## n_integ_point for Lin class
    ///
    /// * `1` -- Conventional Legendre integration points and weights
    /// * `2` -- Conventional Legendre integration points and weights
    /// * `3` -- Conventional Legendre integration points and weights
    /// * `4` -- Conventional Legendre integration points and weights
    /// * `5` -- Conventional Legendre integration points and weights
    ///
    /// ## n_integ_point for Tri class
    ///
    /// * `1` -- Internal integration points and weights
    /// * `3` -- Internal integration points and weights
    /// * `1_003` -- Edge integration points and weights
    /// * `4` -- Internal integration points and weights
    /// * `12` -- Internal integration points and weights
    /// * `16` -- Internal integration points and weights
    ///
    /// ## n_integ_point for Qua class
    ///
    /// * `1` -- Conventional Legendre integration points and weights
    /// * `4` -- Conventional Legendre integration points and weights
    /// * `5` -- Wilson's integration points and weights. "Corner" version
    /// * `1_005` -- 5 points. Wilson's integration points and weights. "Stable" version version with w0=0.004 and wa=0.999 to mimic 4-point rule
    /// * `8` -- Wilson's integration points and weights.
    /// * `9` -- Conventional Legendre integration points and weights
    /// * `16` -- Conventional Legendre integration points and weights
    ///
    /// ## n_integ_point for Tet class
    ///
    /// * `1` -- Internal integration points and weights
    /// * `4` -- Internal integration points and weights
    /// * `5` -- Internal integration points and weights
    /// * `6` -- Internal integration points and weights
    ///
    /// ## n_integ_point for Hex class
    ///
    /// * `6` -- Iron's integration points and weights
    /// * `8` -- Conventional Legendre integration points and weights
    /// * `9` -- Wilson's integration points and weights. "Corner" version
    /// * `1_009` -- Wilson's integration points and weights. "Stable" version
    /// * `14` -- Iron's integration points and weights
    /// * `27` -- Conventional Legendre integration points and weights
    pub fn select_integ_points(&mut self, n_integ_point: usize) -> Result<(), StrError> {
        self.integ_points = match self.class {
            // Lin
            GeoClass::Lin => match n_integ_point {
                1 => &IP_LIN_LEGENDRE_1,
                2 => &IP_LIN_LEGENDRE_2,
                3 => &IP_LIN_LEGENDRE_3,
                4 => &IP_LIN_LEGENDRE_4,
                5 => &IP_LIN_LEGENDRE_5,
                _ => return Err("number of integration points is not available for Lin class"),
            },
            // Tri
            GeoClass::Tri => match n_integ_point {
                1 => &IP_TRI_INTERNAL_1,
                3 => &IP_TRI_INTERNAL_3,
                1_003 => &IP_TRI_EDGE_3,
                4 => &IP_TRI_INTERNAL_4,
                12 => &IP_TRI_INTERNAL_12,
                16 => &IP_TRI_INTERNAL_16,
                _ => return Err("number of integration points is not available for Tri class"),
            },
            // Qua
            GeoClass::Qua => match n_integ_point {
                1 => &IP_QUA_LEGENDRE_1,
                4 => &IP_QUA_LEGENDRE_4,
                5 => &IP_QUA_WILSON_CORNER_5,
                1_005 => &IP_QUA_WILSON_STABLE_5,
                8 => &IP_QUA_WILSON_8,
                9 => &IP_QUA_LEGENDRE_9,
                16 => &IP_QUA_LEGENDRE_16,
                _ => return Err("number of integration points is not available for Qua class"),
            },
            // Tet
            GeoClass::Tet => match n_integ_point {
                1 => &IP_TET_INTERNAL_1,
                4 => &IP_TET_INTERNAL_4,
                5 => &IP_TET_INTERNAL_5,
                6 => &IP_TET_INTERNAL_6,
                _ => return Err("number of integration points is not available for Tet class"),
            },
            // Hex
            GeoClass::Hex => match n_integ_point {
                6 => &IP_HEX_IRONS_6,
                8 => &IP_HEX_LEGENDRE_8,
                9 => &IP_HEX_WILSON_CORNER_9,
                1_009 => &IP_HEX_WILSON_STABLE_9,
                14 => &IP_HEX_IRONS_14,
                27 => &IP_HEX_LEGENDRE_27,
                _ => return Err("number of integration points is not available for Hex class"),
            },
        };
        Ok(())
    }

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
    /// * `thickness` -- tₕ the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
    /// * `fn_s(index_ip: usize) -> f64` -- s(x(ξ)) or s(ℓ) scalar function,
    ///   however written as a function of the index of the integration point.
    pub fn integ_vec_a_ns<F>(&mut self, a: &mut Vector, thickness: f64, fn_s: F) -> Result<(), StrError>
    where
        F: Fn(usize) -> Result<f64, StrError>,
    {
        // check
        if a.dim() != self.nnode {
            return Err("the length of vector 'a' must be equal to nnode");
        }

        // clear output vector
        a.fill(0.0);

        // loop over integration points
        for index in 0..self.integ_points.len() {
            // ksi coordinates and weight
            let iota = &self.integ_points[index];
            let weight = self.integ_points[index][3];

            // calculate interpolation functions and Jacobian
            self.calc_interp(iota);
            let det_jac = self.calc_jacobian(iota)?;

            // calculate s
            let s = fn_s(index)?;

            // loop over nodes and perform sum
            let coef = thickness * det_jac * weight;
            for m in 0..self.nnode {
                a[m] += self.temp_interp[m] * s * coef;
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
    /// * `thickness` -- tₕ the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
    /// * `fn_v(v: &mut Vector, index_ip: usize)` -- v(x(ξ)) vector function with `v.dim() == space_ndim`,
    ///    however written as function of the index of the integration point.
    pub fn integ_vec_b_nv<F>(&mut self, b: &mut Vector, thickness: f64, fn_v: F) -> Result<(), StrError>
    where
        F: Fn(&mut Vector, usize) -> Result<(), StrError>,
    {
        // check
        if b.dim() != self.nnode * self.space_ndim {
            return Err("the length of vector 'b' must be equal to nnode * space_ndim");
        }

        // allocate auxiliary vector
        let mut v = Vector::new(self.space_ndim);

        // clear output vector
        b.fill(0.0);

        // loop over integration points
        for index in 0..self.integ_points.len() {
            // ksi coordinates and weight
            let iota = &self.integ_points[index];
            let weight = self.integ_points[index][3];

            // calculate interpolation functions and Jacobian
            self.calc_interp(iota);
            let det_jac = self.calc_jacobian(iota)?;

            // calculate v
            fn_v(&mut v, index)?;

            // add contribution to b vector
            let coef = thickness * det_jac * weight;
            self.add_to_vec_b(b, &v, coef);
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
    /// * `thickness` -- tₕ the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
    /// * `fn_w(w: &mut Vector, index_ip: usize)` -- w(x(ξ)) vector function with `w.dim() == space_ndim`,
    ///   however written as a function of the index of the integration point.
    pub fn integ_vec_c_vg<F>(&mut self, c: &mut Vector, thickness: f64, fn_w: F) -> Result<(), StrError>
    where
        F: Fn(&mut Vector, usize) -> Result<(), StrError>,
    {
        // check
        if c.dim() != self.nnode {
            return Err("the length of vector 'c' must be equal to nnode");
        }

        // allocate auxiliary vector
        let mut w = Vector::new(self.space_ndim);

        // clear output vector
        c.fill(0.0);

        // loop over integration points
        for index in 0..self.integ_points.len() {
            // ksi coordinates and weight
            let iota = &self.integ_points[index];
            let weight = self.integ_points[index][3];

            // calculate Jacobian and Gradient
            let det_jac = self.calc_gradient(iota)?;

            // calculate w
            fn_w(&mut w, index)?;

            // add contribution to c vector
            let coef = thickness * det_jac * weight;
            self.add_to_vec_c(c, &w, coef);
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
    /// * `thickness` -- tₕ the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
    /// * `calc_sig` -- calculates the σ tensor at given integration point using `calc_sig(sig,index_ip)`
    ///
    /// # Note
    ///
    /// This function is only available for space_ndim = 2D or 3D.
    pub fn integ_vec_d_tg<F>(&mut self, d: &mut Vector, thickness: f64, calc_sig: F) -> Result<(), StrError>
    where
        F: Fn(&mut Tensor2, usize) -> Result<(), StrError>,
    {
        // check
        if self.space_ndim == 1 {
            return Err("space_ndim must be 2 or 3");
        }
        if d.dim() != self.nnode * self.space_ndim {
            return Err("the length of vector 'd' must be equal to nnode * space_ndim");
        }

        // allocate auxiliary tensor
        let mut sig = Tensor2::new(true, self.space_ndim == 2);

        // clear output vector
        d.fill(0.0);

        // loop over integration points
        for index in 0..self.integ_points.len() {
            // ksi coordinates and weight
            let iota = &self.integ_points[index];
            let weight = self.integ_points[index][3];

            // calculate Jacobian and Gradient
            let det_jac = self.calc_gradient(iota)?;

            // calculate σ tensor
            calc_sig(&mut sig, index)?;

            // add contribution to d vector
            let coef = thickness * det_jac * weight;
            self.add_to_vec_d(d, &sig, coef);
        }
        Ok(())
    }

    /// Implements the shape(N)-scalar(S)-shape(N) integration case (e.g., diffusion matrix)
    pub fn integ_mat_1_nsn(&self) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the gradient(G)-vector(V)-shape(N) integration case (e.g., compressibility matrix)
    pub fn integ_mat_2_gvn(&self) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the gradient(G)-tensor(T)-gradient(G) integration case (e.g., conductivity matrix)
    pub fn integ_mat_3_gtg(&self) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the shape(N)-scalar(S)-gradient(G) integration case with different shapes (e.g., coupling matrix)
    pub fn integ_mat_4_nsg(&self) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the gradient(G)-tensor(T)-shape(N) integration case with different shapes (e.g., coupling matrix)
    pub fn integ_mat_5_gtn(&self) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the shape(N)-vector(V)-shape(N) integration case with different shapes (e.g., coupling matrix)
    pub fn integ_mat_6_nvn(&self) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the gradient(G)-scalar(S)-shape(N) integration case with different shapes (e.g., coupling matrix)
    pub fn integ_mat_7_gsn(&self) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the shape(N)-tensor(T)-shape(N) integration case (e.g., mass matrix)
    pub fn integ_mat_8_ntn(&self) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the shape(n)-vector(v)-gradient(g) integration case (e.g., variable density matrix)
    pub fn integ_mat_9_nvg(&self) -> Result<(), StrError> {
        Ok(())
    }

    /// Implements the gradient(G)-4th-tensor(D)-gradient(G) integration case (e.g., stiffness matrix)
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
    /// * `thickness` -- tₕ the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
    /// * `calc_dd` -- calculates the D tensor at given integration point using `calc_dd(dd,index_ip)`
    ///
    /// # Note
    ///
    /// This function is only available for space_ndim = 2D or 3D.
    pub fn integ_mat_10_gdg<F>(&mut self, kk: &mut Matrix, thickness: f64, calc_dd: F) -> Result<(), StrError>
    where
        F: Fn(&mut Tensor4, usize) -> Result<(), StrError>,
    {
        // check
        let (nrow_kk, ncol_kk) = kk.dims();
        if self.space_ndim == 1 {
            return Err("space_ndim must be 2 or 3");
        }
        if nrow_kk != ncol_kk || nrow_kk != self.nnode * self.space_ndim {
            return Err("K matrix must be square with dim equal to nnode * space_ndim");
        }

        // allocate auxiliary tensor
        let mut dd = Tensor4::new(true, self.space_ndim == 2);

        // clear output matrix
        kk.fill(0.0);

        // loop over integration points
        for index in 0..self.integ_points.len() {
            // ksi coordinates and weight
            let iota = &self.integ_points[index];
            let weight = self.integ_points[index][3];

            // calculate Jacobian and Gradient
            let det_jac = self.calc_gradient(iota)?;

            // calculate D tensor
            calc_dd(&mut dd, index)?;

            // add contribution to K matrix
            let coef = det_jac * weight * thickness;
            self.add_to_mat_kk(kk, &dd, coef);
        }
        Ok(())
    }

    /// Adds contribution to the b-vector in integ_vec_2_nv
    #[inline]
    fn add_to_vec_b(&self, b: &mut Vector, v: &Vector, coef: f64) {
        let nn = &self.temp_interp;
        if self.space_ndim == 1 {
            for m in 0..self.nnode {
                b[m] += coef * nn[m] * v[0];
            }
        } else if self.space_ndim == 2 {
            for m in 0..self.nnode {
                b[0 + m * 2] += coef * nn[m] * v[0];
                b[1 + m * 2] += coef * nn[m] * v[1];
            }
        } else {
            for m in 0..self.nnode {
                b[0 + m * 3] += coef * nn[m] * v[0];
                b[1 + m * 3] += coef * nn[m] * v[1];
                b[2 + m * 3] += coef * nn[m] * v[2];
            }
        }
    }

    /// Adds contribution to the c-vector in integ_vec_3_vg
    #[inline]
    fn add_to_vec_c(&self, c: &mut Vector, w: &Vector, coef: f64) {
        let g = &self.temp_gradient;
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
    fn add_to_vec_d(&self, d: &mut Vector, sig: &Tensor2, coef: f64) {
        let s = SQRT_2;
        let g = &self.temp_gradient;
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
    fn add_to_mat_kk(&self,  kk: &mut Matrix, dd: &Tensor4, c: f64) {
        let s = SQRT_2;
        let g = &self.temp_gradient;
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
    use crate::{shapes::AnalyticalTri3, util::SQRT_3};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::{copy_matrix, copy_vector, Matrix};
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
        shape.set_node(0, 0, 0, xmin).unwrap();
        shape.set_node(0, 0, 1, ymin).unwrap();
        shape.set_node(1, 1, 0, xmin + l).unwrap();
        shape.set_node(1, 1, 1, ymin).unwrap();
        shape.set_node(2, 2, 0, xmin + l / 2.0).unwrap();
        shape.set_node(2, 2, 1, ymin + h).unwrap();
        (shape, area)
    }

    // line segment from xa to xb (geo_ndim == space_ndim)
    fn gen_lin2() -> (Shape, f64, f64) {
        let mut shape = Shape::new(1, 1, 2).unwrap();
        let (xa, xb) = (3.0, 9.0);
        shape.set_node(0, 0, 0, xa).unwrap();
        shape.set_node(1, 1, 0, xb).unwrap();
        (shape, xa, xb)
    }

    #[test]
    fn select_integ_points_works() -> Result<(), StrError> {
        // Lin
        let mut shape = Shape::new(1, 1, 2)?;
        for n_integ_point in [1, 2, 3, 4, 5] {
            shape.select_integ_points(n_integ_point)?;
            assert_eq!(shape.integ_points.len(), n_integ_point);
        }
        assert_eq!(
            shape.select_integ_points(100).err(),
            Some("number of integration points is not available for Lin class")
        );

        // Tri
        let mut shape = Shape::new(2, 2, 3)?;
        for n_integ_point in [1, 3, 4, 12, 16] {
            shape.select_integ_points(n_integ_point)?;
            assert_eq!(shape.integ_points.len(), n_integ_point);
        }
        shape.select_integ_points(1_003)?;
        assert_eq!(shape.integ_points.len(), 3);
        assert_eq!(
            shape.select_integ_points(100).err(),
            Some("number of integration points is not available for Tri class")
        );

        // Qua
        let mut shape = Shape::new(2, 2, 4)?;
        for n_integ_point in [1, 4, 5, 8, 9, 16] {
            shape.select_integ_points(n_integ_point)?;
            assert_eq!(shape.integ_points.len(), n_integ_point);
        }
        shape.select_integ_points(1_005)?;
        assert_eq!(shape.integ_points.len(), 5);
        assert_eq!(
            shape.select_integ_points(100).err(),
            Some("number of integration points is not available for Qua class")
        );

        // Tet
        let mut shape = Shape::new(3, 3, 4)?;
        for n_integ_point in [1, 4, 5, 6] {
            shape.select_integ_points(n_integ_point)?;
            assert_eq!(shape.integ_points.len(), n_integ_point);
        }
        assert_eq!(
            shape.select_integ_points(100).err(),
            Some("number of integration points is not available for Tet class")
        );

        // Hex
        let mut shape = Shape::new(3, 3, 8)?;
        for n_integ_point in [6, 8, 9, 14, 27] {
            shape.select_integ_points(n_integ_point)?;
            assert_eq!(shape.integ_points.len(), n_integ_point);
        }
        shape.select_integ_points(1_009)?;
        assert_eq!(shape.integ_points.len(), 9);
        assert_eq!(
            shape.select_integ_points(100).err(),
            Some("number of integration points is not available for Hex class")
        );
        Ok(())
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
        let (mut tri3, area) = gen_tri3();
        const CS: f64 = 3.0;
        let mut a = Vector::filled(tri3.nnode, NOISE);
        tri3.integ_vec_a_ns(&mut a, 1.0, |_| Ok(CS))?;
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
        let all_integ_points = lin2.calc_integ_points_coords()?;
        let mut a = Vector::new(lin2.nnode);
        lin2.integ_vec_a_ns(&mut a, 1.0, |index_ip: usize| Ok(all_integ_points[index_ip][0]))?;
        let cf = (xb - xa) / 6.0;
        let a_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];
        assert_vec_approx_eq!(a.as_data(), a_correct, 1e-15);
        Ok(())
    }

    #[test]
    fn integ_vec_b_works() -> Result<(), StrError> {
        // This test is similar to the case_a with tri3, however using a vector
        // So, each component of `b` equals `Fₛ`
        let (mut tri3, area) = gen_tri3();
        const CS: f64 = 3.0;
        let mut b = Vector::filled(tri3.nnode * tri3.space_ndim, NOISE);
        tri3.integ_vec_b_nv(&mut b, 1.0, |v: &mut Vector, _: usize| {
            v.fill(CS);
            Ok(())
        })?;
        let cf = CS * area / 3.0;
        let b_correct = &[cf, cf, cf, cf, cf, cf];
        assert_vec_approx_eq!(b.as_data(), b_correct, 1e-14);

        // Likewise, this test is similar to case_a with lin2, however using a vector
        // with a single component. So, each component of `b` equals `Fₛ`
        let (mut lin2, xa, xb) = gen_lin2();
        let all_integ_points = lin2.calc_integ_points_coords()?;
        let mut b = Vector::filled(lin2.nnode * lin2.space_ndim, NOISE);
        lin2.integ_vec_b_nv(&mut b, 1.0, |v: &mut Vector, index: usize| {
            v.fill(all_integ_points[index][0]);
            Ok(())
        })?;
        let cf = (xb - xa) / 6.0;
        let b_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];
        assert_vec_approx_eq!(b.as_data(), b_correct, 1e-15);
        Ok(())
    }

    #[test]
    fn integ_vec_c_works() -> Result<(), StrError> {
        // shape and analytical gradient
        let (mut tri3, area) = gen_tri3();
        let mut ana = AnalyticalTri3::new(&mut tri3);
        assert_eq!(area, ana.area);

        // constant vector function: w(x) = {w₀, w₁}
        // solution:
        //    cᵐ = ½ (w₀ bₘ + w₁ cₘ)
        const W0: f64 = 2.0;
        const W1: f64 = 3.0;
        let c_correct = ana.integ_vec_c_constant(W0, W1);
        let mut c = Vector::filled(tri3.nnode, NOISE);
        tri3.integ_vec_c_vg(&mut c, 1.0, |w: &mut Vector, _: usize| {
            w[0] = W0;
            w[1] = W1;
            Ok(())
        })?;
        assert_vec_approx_eq!(c.as_data(), c_correct, 1e-15);

        // bilinear vector function: w(x) = {x, y}
        // solution:
        //    cᵐ = ⅙ bₘ (x₀+x₁+x₂) + ⅙ cₘ (y₀+y₁+y₂)
        let all_integ_points = tri3.calc_integ_points_coords()?;
        let c_correct = ana.integ_vec_c_bilinear();
        let mut c = Vector::filled(tri3.nnode, NOISE);
        tri3.integ_vec_c_vg(&mut c, 1.0, |w: &mut Vector, index: usize| {
            w[0] = all_integ_points[index][0];
            w[1] = all_integ_points[index][1];
            Ok(())
        })?;
        assert_vec_approx_eq!(c.as_data(), c_correct, 1e-14);
        Ok(())
    }

    #[derive(Clone)]
    pub struct StateStress {
        pub sigma: Tensor2,            // total or effective stress
        pub internal_values: Vec<f64>, // internal values
    }

    pub struct StateElement {
        pub stress: Vec<StateStress>, // (n_integ_point)
    }

    struct LinElastModel2d {
        lin_elast: LinElasticity,
    }

    impl LinElastModel2d {
        pub fn new(young: f64, poisson: f64, plane_stress: bool) -> Self {
            LinElastModel2d {
                lin_elast: LinElasticity::new(young, poisson, true, plane_stress),
            }
        }
        pub fn consistent_modulus(&self, dd: &mut Tensor4, _state: &StateStress) -> Result<(), StrError> {
            let dd_ela = self.lin_elast.get_modulus();
            copy_matrix(&mut dd.mat, &dd_ela.mat)
        }
    }

    struct ElemElast2d {
        shape: Shape,
        model: LinElastModel2d,
        thickness: f64,
        pub f: Vector,
        pub kk: Matrix,
    }

    impl ElemElast2d {
        pub fn new(shape: Shape, young: f64, poisson: f64, plane_stress: bool, thickness: f64) -> Self {
            let neq = shape.nnode * shape.space_ndim;
            ElemElast2d {
                shape,
                model: LinElastModel2d::new(young, poisson, plane_stress),
                thickness,
                f: Vector::new(neq),
                kk: Matrix::new(neq, neq),
            }
        }
        pub fn calculate_f(&mut self, state: &StateElement) -> Result<(), StrError> {
            self.shape.integ_vec_d_tg(&mut self.f, self.thickness, |sig, index_ip| {
                copy_vector(&mut sig.vec, &state.stress[index_ip].sigma.vec)
            })
        }
        pub fn calculate_kk(&mut self, state: &StateElement) -> Result<(), StrError> {
            self.shape
                .integ_mat_10_gdg(&mut self.kk, self.thickness, |dd, index_ip| {
                    self.model.consistent_modulus(dd, &state.stress[index_ip])
                })
        }
    }

    #[test]
    fn integ_vec_d_works() -> Result<(), StrError> {
        // shape and analytical gradient
        let (mut tri3, _) = gen_tri3();
        let mut ana = AnalyticalTri3::new(&mut tri3);

        // constant tensor function: σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2}
        // solution:
        //    dᵐ₀ = ½ (σ₀₀ bₘ + σ₀₁ cₘ)
        //    dᵐ₁ = ½ (σ₁₀ bₘ + σ₁₁ cₘ)
        const S00: f64 = 2.0;
        const S11: f64 = 3.0;
        const S22: f64 = 4.0;
        const S01: f64 = 5.0;
        let f_vec_correct = ana.integ_vec_d_constant(S00, S11, S01);

        // constants
        let young = 10_000.0;
        let poisson = 0.2;
        let thickness = 1.0;

        // element instance with callback function calc_sig for integration
        let n_integ_point = tri3.integ_points.len();
        let mut element = ElemElast2d::new(tri3, young, poisson, true, thickness);

        // state at all integration points
        #[rustfmt::skip]
        let sigma = Tensor2::from_matrix(&[
            [S00, S01, 0.0],
            [S01, S11, 0.0],
            [0.0, 0.0, S22],
        ],true,true)?;
        let state = StateElement {
            stress: vec![
                StateStress {
                    sigma,
                    internal_values: Vec::new()
                };
                n_integ_point
            ],
        };

        // test integration
        element.calculate_f(&state)?;
        assert_vec_approx_eq!(element.f.as_data(), f_vec_correct, 1e-14);
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
        tri3.set_node(0, 0, 0, 0.0).unwrap();
        tri3.set_node(0, 0, 1, 0.0).unwrap();
        tri3.set_node(2, 1, 0, 2.0).unwrap();
        tri3.set_node(2, 1, 1, 0.0).unwrap();
        tri3.set_node(3, 2, 0, 2.0).unwrap();
        tri3.set_node(3, 2, 1, 1.5).unwrap();
        let mut ana = AnalyticalTri3::new(&mut tri3);

        // constants
        let young = 10_000.0;
        let poisson = 0.2;
        let thickness = 0.25; // thickness

        // element instance with callback function calc_dd for integration
        let n_integ_point = tri3.integ_points.len();
        let mut element = ElemElast2d::new(tri3, young, poisson, true, thickness);

        // state at all integration points
        #[rustfmt::skip]
        let sigma = Tensor2::from_matrix(&[
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ],true,true)?;
        let state = StateElement {
            stress: vec![
                StateStress {
                    sigma,
                    internal_values: Vec::new()
                };
                n_integ_point
            ],
        };

        // test integration
        element.calculate_kk(&state)?;

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
        let kk_correct = ana.integ_stiffness(young, poisson, thickness)?;
        assert_vec_approx_eq!(kk_correct.as_data(), kk_bhatti.as_data(), 1e-12);
        assert_vec_approx_eq!(kk_correct.as_data(), element.kk.as_data(), 1e-12);
        Ok(())
    }
}
