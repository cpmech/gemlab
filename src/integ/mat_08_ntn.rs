use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::math::SQRT_2;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the shape(N) times tensor(T) times shape(N) integration case 08 (e.g., mass matrix)
///
/// Callback function: `α ← f(T, p, N, G)`
///
/// Mass coefficients:
///
/// ```text
///       ⌠
/// Kᵐⁿ = │ Nᵐ T Nⁿ α dΩ
/// ▔     ⌡    ▔
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///         nip-1    →       →      →       →
/// Kᵐⁿᵢⱼ ≈   Σ   Nᵐ(ιᵖ) Tᵢⱼ(ιᵖ) Nⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
///          p=0
/// ```
///
/// # Output
///
/// ```text
///     ┌                                               ┐
///     | K⁰⁰₀₀ K⁰⁰₀₁ K⁰¹₀₀ K⁰¹₀₁ K⁰²₀₀ K⁰²₀₁ ··· K⁰ⁿ₀ⱼ |  ⟸  ii0
///     | K⁰⁰₁₀ K⁰⁰₁₁ K⁰¹₁₀ K⁰¹₁₁ K⁰²₁₀ K⁰²₁₁ ··· K⁰ⁿ₁ⱼ |
///     | K¹⁰₀₀ K¹⁰₀₁ K¹¹₀₀ K¹¹₀₁ K¹²₀₀ K¹²₀₁ ··· K¹ⁿ₀ⱼ |
/// K = | K¹⁰₁₀ K¹⁰₁₁ K¹¹₁₀ K¹¹₁₁ K¹²₁₀ K¹²₁₁ ··· K¹ⁿ₁ⱼ |
///     | K²⁰₀₀ K²⁰₀₁ K²¹₀₀ K²¹₀₁ K²²₀₀ K²²₀₁ ··· K²ⁿ₀ⱼ |
///     | K²⁰₁₀ K²⁰₁₁ K²¹₁₀ K²¹₁₁ K²²₁₀ K²²₁₁ ··· K²ⁿ₁ⱼ |
///     |  ···   ···   ···   ···   ···   ···  ···  ···  |
///     | Kᵐ⁰ᵢ₀ Kᵐ⁰ᵢ₁ Kᵐ¹ᵢ₀ Kᵐ¹ᵢ₁ Kᵐ²ᵢ₀ Kᵐ²ᵢ₁ ··· Kᵐⁿᵢⱼ |  ⟸  ii := i + m ⋅ space_ndim
///     └                                               ┘
///        ⇑                                        ⇑
///       jj0                                       jj := j + n ⋅ space_ndim
///
/// m = ii / space_ndim    n = jj / space_ndim
/// i = ii % space_ndim    j = jj % space_ndim
/// ```
///
/// * `kk` -- A matrix containing all `Kᵐⁿᵢⱼ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `i` and `j` correspond to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode ⋅ space_ndim` and `ncol(K) ≥ jj0 + nnode ⋅ space_ndim`.
/// * `pad` -- Some members of the scratchpad will be modified.
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output matrix where to add components.
/// * `jj0` -- Stride marking the first column in the output matrix where to add components.
/// * `clear_kk` -- Fills `kk` matrix with zeros, otherwise accumulate values into `kk`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_tt` -- Function `f(T,p,N,G)→α` that computes `T(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions N(ιᵖ), and gradients G(ιᵖ). `T` is set for `space_ndim`.
///   `fn_tt` returns α that can accommodate plane-strain or axisymmetric simulations.
pub fn mat_08_ntn<F>(
    kk: &mut Matrix,
    pad: &mut Scratchpad,
    ii0: usize,
    jj0: usize,
    clear_kk: bool,
    ips: IntegPointData,
    mut fn_tt: F,
) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor2, usize, &Vector, &Matrix) -> Result<f64, StrError>,
{
    // check
    let (space_ndim, nnode) = pad.xxt.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    if nrow_kk < ii0 + nnode * space_ndim {
        return Err("nrow(K) must be ≥ ii0 + nnode ⋅ space_ndim");
    }
    if ncol_kk < jj0 + nnode * space_ndim {
        return Err("ncol(K) must be ≥ jj0 + nnode ⋅ space_ndim");
    }

    // allocate auxiliary tensor
    let mut tt = Tensor2::new(true, space_ndim == 2);

    // clear output matrix
    if clear_kk {
        kk.fill(0.0);
    }

    // loop over integration points
    let s = SQRT_2;
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate interpolation functions and Jacobian
        (pad.fn_interp)(&mut pad.interp, iota); // N
        let det_jac = pad.calc_gradient(iota)?; // G

        // calculate T tensor
        let nn = &pad.interp;
        let gg = &pad.gradient;
        let alpha = fn_tt(&mut tt, p, nn, gg)?;

        // add contribution to K matrix
        let c = alpha * det_jac * weight;
        let t = &tt.vec;
        if space_ndim == 2 {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + 0 + m * 2][jj0 + 0 + n * 2] += c * nn[m] * t[0] * nn[n];
                    kk[ii0 + 0 + m * 2][jj0 + 1 + n * 2] += c * nn[m] * t[3] * nn[n] / s;

                    kk[ii0 + 1 + m * 2][jj0 + 0 + n * 2] += c * nn[m] * t[3] * nn[n] / s;
                    kk[ii0 + 1 + m * 2][jj0 + 1 + n * 2] += c * nn[m] * t[1] * nn[n];
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + 0 + m * 3][jj0 + 0 + n * 3] += c * nn[m] * t[0] * nn[n];
                    kk[ii0 + 0 + m * 3][jj0 + 1 + n * 3] += c * nn[m] * t[3] * nn[n] / s;
                    kk[ii0 + 0 + m * 3][jj0 + 2 + n * 3] += c * nn[m] * t[5] * nn[n] / s;

                    kk[ii0 + 1 + m * 3][jj0 + 0 + n * 3] += c * nn[m] * t[3] * nn[n] / s;
                    kk[ii0 + 1 + m * 3][jj0 + 1 + n * 3] += c * nn[m] * t[1] * nn[n];
                    kk[ii0 + 1 + m * 3][jj0 + 2 + n * 3] += c * nn[m] * t[4] * nn[n] / s;

                    kk[ii0 + 2 + m * 3][jj0 + 0 + n * 3] += c * nn[m] * t[5] * nn[n] / s;
                    kk[ii0 + 2 + m * 3][jj0 + 1 + n * 3] += c * nn[m] * t[4] * nn[n] / s;
                    kk[ii0 + 2 + m * 3][jj0 + 2 + n * 3] += c * nn[m] * t[2] * nn[n];
                }
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, AnalyticalTet4, AnalyticalTri3, IP_LIN_LEGENDRE_1, IP_TRI_INTERNAL_1};
    use russell_chk::vec_approx_eq;
    use russell_lab::{Matrix, Vector};
    use russell_tensor::{copy_tensor2, Tensor2};

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(4, 4);
        let mut tt = Tensor2::new(true, true);
        let nn = Vector::new(0);
        let gg = Matrix::new(0, 0);
        let f = |_tt: &mut Tensor2, _p: usize, _nn: &Vector, _gg: &Matrix| Ok(1.0);
        assert_eq!(f(&mut tt, 0, &nn, &gg).unwrap(), 1.0);
        assert_eq!(
            integ::mat_08_ntn(&mut kk, &mut pad, 1, 0, false, &[], f).err(),
            Some("nrow(K) must be ≥ ii0 + nnode ⋅ space_ndim")
        );
        assert_eq!(
            integ::mat_08_ntn(&mut kk, &mut pad, 0, 1, false, &[], f).err(),
            Some("ncol(K) must be ≥ jj0 + nnode ⋅ space_ndim")
        );
        // more errors
        assert_eq!(
            integ::mat_08_ntn(&mut kk, &mut pad, 0, 0, false, &IP_LIN_LEGENDRE_1, f).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(6, 6);
        assert_eq!(
            integ::mat_08_ntn(&mut kk, &mut pad, 0, 0, false, &IP_TRI_INTERNAL_1, |_, _, _, _| Err(
                "stop"
            ))
            .err(),
            Some("stop")
        );
    }

    #[test]
    fn tri3_works() {
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3 * 2, 3 * 2);
        let ana = AnalyticalTri3::new(&pad);
        let (rho, th) = (2.7, 1.0);
        let kk_correct = ana.mat_08_ntn(rho, th);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-15, 1e-15];
        let selection: Vec<_> = [3, 6].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_08_ntn(&mut kk, &mut pad, 0, 0, true, ips, |tt, _, _, _| {
                tt.sym_set(0, 0, rho);
                tt.sym_set(1, 1, rho);
                Ok(1.0)
            })
            .unwrap();
            // println!("{}", kk);
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn tet4_works() {
        let mut pad = aux::gen_pad_tet4();
        let mut kk = Matrix::new(4 * 3, 4 * 3);
        let ana = AnalyticalTet4::new(&pad);
        #[rustfmt::skip]
        let sig = Tensor2::from_matrix(&[
            [1.1, 1.2, 1.3],
            [1.2, 2.2, 2.3],
            [1.3, 2.3, 3.3]],
        true, false).unwrap();
        let kk_correct = ana.mat_08_ntn(&sig);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_08_ntn(&mut kk, &mut pad, 0, 0, true, ips, |tt, _, _, _| {
                copy_tensor2(tt, &sig)?;
                Ok(1.0)
            })
            .unwrap();
            // println!("{}", kk);
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }
}
