use super::*;
use crate::shapes::{GeoClass, Shape};
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_tensor::{Tensor2, Tensor4};

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
    ///      nip-1     →     →       →
    /// aᵐ ≈   Σ    Nᵐ(ιᵖ) s(ιᵖ) |J|(ιᵖ) wᵖ
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
    pub fn integ_case_a<F>(&mut self, a: &mut Vector, fn_s: F) -> Result<(), StrError>
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
        for index in 0..self.ip_data.len() {
            // ksi coordinates and weight
            let iota = &self.ip_data[index];
            let weight = self.ip_data[index][3];

            // calculate interpolation functions and Jacobian
            self.calc_interp(iota);
            let det_jac = self.calc_jacobian(iota)?;

            // calculate s
            let s = fn_s(index);

            // loop over nodes and perform summation
            for m in 0..self.nnode {
                a[m] += self.interp[m] * s * det_jac * weight;
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
    /// →    nip-1     →   → →       →
    /// bᵐ ≈   Σ    Nᵐ(ιᵖ) v(ιᵖ) |J|(ιᵖ) wᵖ
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
    /// * `aux_v` -- is an auxiliary vector with size equal to `space_ndim`.
    pub fn integ_case_b<F>(&mut self, b: &mut Vector, fn_v: F, aux_v: &mut Vector) -> Result<(), StrError>
    where
        F: Fn(&mut Vector, usize),
    {
        // check
        if b.dim() != self.nnode * self.space_ndim {
            return Err("the length of vector 'b' must be equal to nnode * space_ndim");
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

            // calculate v
            fn_v(aux_v, index);

            // loop over nodes and perform summation
            for m in 0..self.nnode {
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
    ///      nip-1  → →     →  →       →
    /// cᵐ ≈   Σ    w(ιᵖ) · Gᵐ(ιᵖ) |J|(ιᵖ) wᵖ
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
    /// * `aux_w` -- is an auxiliary vector with size equal to `space_ndim`.
    pub fn integ_case_c<F>(&mut self, c: &mut Vector, fn_w: F, aux_w: &mut Vector) -> Result<(), StrError>
    where
        F: Fn(&mut Vector, usize),
    {
        // check
        if c.dim() != self.nnode {
            return Err("the length of vector 'c' must be equal to nnode");
        }
        if aux_w.dim() != self.space_ndim {
            return Err("the length of vector 'aux_w' must be equal to space_ndim");
        }

        // clear output vector
        c.fill(0.0);

        // loop over integration points
        for index in 0..self.ip_data.len() {
            // ksi coordinates and weight
            let iota = &self.ip_data[index];
            let weight = self.ip_data[index][3];

            // calculate Jacobian and Gradient
            let det_jac = self.calc_gradient(iota)?;

            // calculate w
            fn_w(aux_w, index);

            // loop over nodes and perform summation
            for m in 0..self.nnode {
                let w_dot_grad = self.vec_dot_grad(m, aux_w);
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
    /// →    nip-1    →     →  →       →
    /// dᵐ ≈   Σ    σ(ιᵖ) · Gᵐ(ιᵖ) |J|(ιᵖ) wᵖ
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
    /// * `aux_sig` -- is an auxiliary Tensor2.
    /// * `aux_vec` -- is an auxiliary Vector with size equal to `space_ndim`.
    pub fn integ_case_d<F>(
        &mut self,
        d: &mut Vector,
        fn_sig: F,
        aux_sig: &mut Tensor2,
        aux_vec: &mut Vector,
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
        if aux_sig.vec.dim() != 2 * self.space_ndim {
            return Err("'aux_sig' must be symmetric with dim equal to 4 in 2D or 6 in 3D");
        }
        if aux_vec.dim() != self.space_ndim {
            return Err("the length of vector 'aux_vec' must be equal to space_ndim");
        }

        // clear output vector
        d.fill(0.0);

        // loop over integration points
        for index in 0..self.ip_data.len() {
            // ksi coordinates and weight
            let iota = &self.ip_data[index];
            let weight = self.ip_data[index][3];

            // calculate Jacobian and Gradient
            let det_jac = self.calc_gradient(iota)?;

            // calculate σ
            fn_sig(aux_sig, index);

            // loop over nodes and perform summation
            for m in 0..self.nnode {
                // aux_vec := σ · G
                self.tensor_dot_grad(aux_vec, m, &aux_sig);
                for i in 0..self.space_ndim {
                    let ii = i + m * self.space_ndim;
                    d[ii] += aux_vec[i] * det_jac * weight;
                }
            }
        }
        Ok(())
    }

    /// Computes the Nᵐ s Nⁿ integral
    ///
    /// ```text
    /// →    ⌠   →    →  → →   →
    /// dᵐ = │ σ(x) · Gᵐ(x(ξ)) dΩ
    ///      ⌡ ▔
    ///      Ω
    /// ```
    pub fn integ_mat_nsn() {}
    pub fn integ_mat_gvn() {}
    pub fn integ_mat_gtg() {}
    pub fn integ_mat_ntn() {}
    pub fn integ_mat_nvg() {}

    /// Computes the Gᵐₖ Dᵢₖⱼₗ Gⁿₗ integral (stiffness)
    ///
    /// Stiffness tensors:
    ///
    /// ```text
    ///       ⌠               →    →
    /// Kᵐⁿ = │ Gᵐₖ Dᵢₖⱼₗ Gⁿₗ eᵢ ⊗ eⱼ dΩ
    /// ▔     ⌡
    ///       Ωₑ
    /// ```
    ///
    /// The numerical integration is:
    ///
    /// ```text
    ///         nip-1     →         →       →       →
    /// Kᵐⁿᵢⱼ ≈   Σ   Gᵐₖ(ιᵖ) Dᵢₖⱼₗ(ιᵖ) Gⁿₗ(ιᵖ) |J|(ιᵖ) wᵖ
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
    /// * `aux_dd` -- is an auxiliary Tensor4 (minor-symmetric in 2D or 3D with 4 or 6 components, respectively).
    /// * `thickness` -- thickness of the 2D domain (e.g., for plane-stress problems).
    pub fn integ_mat_gdg<F>(
        &mut self,
        kk: &mut Matrix,
        fn_dd: F,
        aux_dd: &mut Tensor4,
        thickness: f64,
    ) -> Result<(), StrError>
    where
        F: Fn(&mut Tensor4, usize),
    {
        // check
        let (nrow_kk, ncol_kk) = kk.dims();
        let (nrow_dd, ncol_dd) = aux_dd.mat.dims();
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
        for index in 0..self.ip_data.len() {
            // ksi coordinates and weight
            let iota = &self.ip_data[index];
            let weight = self.ip_data[index][3];

            // calculate Jacobian and Gradient
            let det_jac = self.calc_gradient(iota)?;

            // calculate constitutive modulus
            fn_dd(aux_dd, index);

            // add contribution to K matrix
            let c = det_jac * weight * thickness;
            self.stiffness_contribution(kk, c, aux_dd);
        }
        Ok(())
    }

    /// Computes vector dot the gradient at node m
    fn vec_dot_grad(&self, m: usize, w: &Vector) -> f64 {
        let mut res = 0.0;
        for i in 0..self.space_ndim {
            res += w[i] * self.gradient[m][i];
        }
        res
    }

    /// Computes tensor dot the gradient at node m
    fn tensor_dot_grad(&self, res: &mut Vector, m: usize, sig: &Tensor2) {
        for i in 0..self.space_ndim {
            res[i] = 0.0;
            for j in 0..self.space_ndim {
                res[i] += sig.get(i, j) * self.gradient[m][j];
            }
        }
    }

    /// Computes the contribution to the stiffness matrix
    #[rustfmt::skip]
    fn stiffness_contribution(&self, kk: &mut Matrix, c: f64, dd: &Tensor4) {
        let s = SQRT_2;
        let g = &self.gradient;
        let d = &dd.mat;
        if self.space_ndim == 3 {
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
        } else {
            for m in 0..self.nnode {
                for n in 0..self.nnode {
                    kk[0+m*2][0+n*2] += c * (g[m][1]*g[n][1]*d[3][3] + s*g[m][1]*g[n][0]*d[3][0] + s*g[m][0]*g[n][1]*d[0][3] + 2.0*g[m][0]*g[n][0]*d[0][0]) / 2.0;
                    kk[0+m*2][1+n*2] += c * (g[m][1]*g[n][0]*d[3][3] + s*g[m][1]*g[n][1]*d[3][1] + s*g[m][0]*g[n][0]*d[0][3] + 2.0*g[m][0]*g[n][1]*d[0][1]) / 2.0;
                    kk[1+m*2][0+n*2] += c * (g[m][0]*g[n][1]*d[3][3] + s*g[m][0]*g[n][0]*d[3][0] + s*g[m][1]*g[n][1]*d[1][3] + 2.0*g[m][1]*g[n][0]*d[1][0]) / 2.0;
                    kk[1+m*2][1+n*2] += c * (g[m][0]*g[n][0]*d[3][3] + s*g[m][0]*g[n][1]*d[3][1] + s*g[m][1]*g[n][0]*d[1][3] + 2.0*g[m][1]*g[n][1]*d[1][1]) / 2.0;
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
    use russell_lab::Matrix;

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
        let mut a = Vector::filled(tri3.nnode, NOISE);
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
        let mut a = Vector::new(lin2.nnode);
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
        let mut b = Vector::filled(tri3.nnode * tri3.space_ndim, NOISE);
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
        let mut b = Vector::filled(lin2.nnode * lin2.space_ndim, NOISE);
        let mut aux_v = Vector::new(lin2.space_ndim);
        lin2.integ_case_b(&mut b, fn_v, &mut aux_v)?;
        let cf = (xb - xa) / 6.0;
        let b_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];
        assert_vec_approx_eq!(b.as_data(), b_correct, 1e-15);
        Ok(())
    }

    struct AnalyticalTri3 {
        x: [f64; 3],
        y: [f64; 3],
        b: [f64; 3],
        c: [f64; 3],
    }

    fn analytical_tri3(area: f64, tri3: &mut Shape) -> AnalyticalTri3 {
        let (x0, y0) = (tri3.coords_transp[0][0], tri3.coords_transp[1][0]);
        let (x1, y1) = (tri3.coords_transp[0][1], tri3.coords_transp[1][1]);
        let (x2, y2) = (tri3.coords_transp[0][2], tri3.coords_transp[1][2]);
        let (b0, b1, b2) = (y1 - y2, y2 - y0, y0 - y1);
        let (c0, c1, c2) = (x2 - x1, x0 - x2, x1 - x0);
        let (f0, f1, f2) = (x1 * y2 - x2 * y1, x2 * y0 - x0 * y2, x0 * y1 - x1 * y0);
        let aa = (f0 + f1 + f2) / 2.0;
        assert_eq!(area, aa);
        let gg = Matrix::from(&[
            [b0 / (2.0 * aa), c0 / (2.0 * aa)],
            [b1 / (2.0 * aa), c1 / (2.0 * aa)],
            [b2 / (2.0 * aa), c2 / (2.0 * aa)],
        ]);
        tri3.calc_gradient(&tri3.ip_data[0]).unwrap();
        assert_eq!(tri3.gradient.as_data(), gg.as_data());
        AnalyticalTri3 {
            x: [x0, x1, x2],
            y: [y0, y1, y2],
            b: [b0, b1, b2],
            c: [c0, c1, c2],
        }
    }

    #[test]
    fn integ_case_c_works() -> Result<(), StrError> {
        // shape and analytical gradient
        let (mut tri3, area) = gen_tri3();
        let ana = analytical_tri3(area, &mut tri3);

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
        let mut aux_w = Vector::new(tri3.space_ndim);
        tri3.integ_case_c(&mut c, fn_w, &mut aux_w)?;
        assert_vec_approx_eq!(c.as_data(), c_correct, 1e-15);

        // bilinear vector function: w(x) = {x, y}
        // solution:
        //    cᵐ = ⅙ bₘ (x₀+x₁+x₂) + ⅙ cₘ (y₀+y₁+y₂)
        let all_int_points = tri3.calc_int_points_coords()?;
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
        let mut aux_w = Vector::new(tri3.space_ndim);
        tri3.integ_case_c(&mut c, fn_w, &mut aux_w)?;
        assert_vec_approx_eq!(c.as_data(), c_correct, 1e-14);
        Ok(())
    }

    #[test]
    fn integ_case_d_works() -> Result<(), StrError> {
        // shape and analytical gradient
        let (mut tri3, area) = gen_tri3();
        let ana = analytical_tri3(area, &mut tri3);

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
        let mut aux_sig = Tensor2::new(true, true);
        let mut aux_vec = Vector::new(tri3.space_ndim);
        tri3.integ_case_d(&mut d, fn_sig, &mut aux_sig, &mut aux_vec)?;
        assert_vec_approx_eq!(d.as_data(), d_correct, 1e-15);
        Ok(())
    }
}
