use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::Matrix;

/// Implements the gradient(G) time scalar(S) time shape(Nb) integration case 07 (e.g., coupling matrix)
///
/// **Notes:**
///
/// * `m` ranges over the number of nodes of the *driver* shape specified by `pad` (for `Gᵐ`)
/// * `n` ranges over the number of nodes of the *lower-order* shape specified by `pad_b` (for `Nbⁿ`)
/// * For example, `1 ≤ m ≤ 8` for a `pad→Qua8` and `1 ≤ n ≤ 4` for `pad_b→Qua4`
/// * The determinant of the Jacobian is calculated for `pad` (`pad` is the driver of the calculations)
/// * The number of integration points must consider the nodes of `pad` and the expected order of the whole integrand
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠ →
/// Kᵐⁿ = │ Gᵐ s Nbⁿ dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///        nip-1     →     →       →       →
/// Kᵐⁿᵢ ≈   Σ   Gᵐᵢ(ιᵖ) s(ιᵖ) Nbⁿ(ιᵖ) |J|(ιᵖ) wᵖ
///         p=0
/// ```
///
/// # Output
///
/// ```text
///     ┌                         ┐
///     | K⁰⁰₀ K⁰¹₀ K⁰²₀ ··· K⁰ⁿ₀ |  ⟸  ii0
///     | K⁰⁰₁ K⁰¹₁ K⁰²₁ ··· K⁰ⁿ₁ |
///     | K¹⁰₀ K¹¹₀ K¹²₀ ··· K¹ⁿ₀ |
/// K = | K¹⁰₁ K¹¹₁ K¹²₁ ··· K¹ⁿ₁ |
///     | K²⁰₀ K²¹₀ K²²₀ ··· K²ⁿ₀ |
///     | K²⁰₁ K²¹₁ K²²₁ ··· K²ⁿ₁ |
///     |  ···  ···  ··· ···  ··· |     (pad)
///     | Kᵐ⁰ᵢ Kᵐ¹ᵢ Kᵐ²ᵢ ··· Kᵐⁿᵢ |  ⟸  ii := i + m ⋅ space_ndim
///     └                         ┘
///        ⇑                  ⇑
///       jj0         (pad_b) jj
///
/// m = ii / space_ndim
/// i = ii % space_ndim
/// ```
///
/// * `kk` -- A matrix containing all `Kᵐⁿᵢ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `i` corresponds to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode ⋅ space_ndim` and `ncol(K) ≥ jj0 + pad_b.nnode`.
/// * `pad` -- Driver scratchpad (modified) to compute G
/// * `pad_b` -- Lower-order scratchpad (modified) to compute Nb
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output matrix where to add components.
/// * `jj0` -- Stride marking the first column in the output matrix where to add components.
/// * `clear_kk` -- Fills `kk` matrix with zeros, otherwise accumulate values into `kk`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_s` -- Function `f(p)` that computes `s(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
///
/// # Warning
///
/// The two [crate::shapes::Scratchpad]s mut be compatible, otherwise **calculation errors may occur**.
/// Therefore, `pad_b` must be either the lower-version of `pad` or have the same shape as `pad`.
pub fn mat_07_gsn<F>(
    kk: &mut Matrix,
    pad: &mut Scratchpad,
    pad_b: &mut Scratchpad,
    ii0: usize,
    jj0: usize,
    clear_kk: bool,
    ips: IntegPointData,
    fn_s: F,
) -> Result<(), StrError>
where
    F: Fn(usize) -> Result<f64, StrError>,
{
    // check
    let nnode_b = pad_b.interp.dim();
    let (space_ndim, nnode) = pad.xxt.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    if nrow_kk < ii0 + nnode * space_ndim {
        return Err("nrow(K) must be ≥ ii0 + pad.nnode ⋅ space_ndim");
    }
    if ncol_kk < jj0 + nnode_b {
        return Err("ncol(K) must be ≥ jj0 + pad_b.nnode");
    }

    // clear output matrix
    if clear_kk {
        kk.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate interpolation functions and Jacobian
        (pad_b.fn_interp)(&mut pad_b.interp, iota);
        let det_jac = pad.calc_gradient(iota)?;

        // calculate s
        let s = fn_s(p)?;

        // add contribution to K matrix
        let val = s * det_jac * weight;
        let g = &pad.gradient;
        let nnb = &pad_b.interp;
        if space_ndim == 2 {
            for m in 0..nnode {
                for n in 0..nnode_b {
                    kk[ii0 + 0 + m * 2][jj0 + n] += g[m][0] * val * nnb[n];
                    kk[ii0 + 1 + m * 2][jj0 + n] += g[m][1] * val * nnb[n];
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode_b {
                    kk[ii0 + 0 + m * 3][jj0 + n] += g[m][0] * val * nnb[n];
                    kk[ii0 + 1 + m * 3][jj0 + n] += g[m][1] * val * nnb[n];
                    kk[ii0 + 2 + m * 3][jj0 + n] += g[m][2] * val * nnb[n];
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
    use crate::integ::{self, AnalyticalQua8};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Matrix;

    #[test]
    fn capture_some_errors() {
        let (a, b) = (2.0, 3.0);
        let mut pad_b = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut kk = Matrix::new(8 * 2, 4);
        assert_eq!(
            integ::mat_07_gsn(&mut kk, &mut pad, &mut pad_b, 1, 0, false, &[], |_| Ok(0.0)).err(),
            Some("nrow(K) must be ≥ ii0 + pad.nnode ⋅ space_ndim")
        );
        assert_eq!(
            integ::mat_07_gsn(&mut kk, &mut pad, &mut pad_b, 0, 1, false, &[], |_| Ok(0.0)).err(),
            Some("ncol(K) must be ≥ jj0 + pad_b.nnode")
        );
    }

    #[test]
    fn mat_07_gsn_works() {
        let (a, b) = (2.0, 3.0);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_b = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut kk = Matrix::new(8 * 2, 4);
        let ana = AnalyticalQua8::new(a, b);
        let s = 9.0;
        let kk_correct = ana.integ_gsnb(s);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [4, 9].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_07_gsn(&mut kk, &mut pad, &mut pad_b, 0, 0, true, ips, |_| Ok(s)).unwrap();
            // println!("{:.2}", kk);
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }
}