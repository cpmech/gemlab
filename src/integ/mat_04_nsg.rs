use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::Matrix;

/// Implements the shape(Nb) time scalar(S) time gradient(G) integration case 04 (e.g., coupling matrix)
///
/// **Notes:**
///
/// * `m` ranges over the number of nodes of the *lower-order* shape specified by `pad_b` (for `Nbᵐ`)
/// * `n` ranges over the number of nodes of the *driver* shape specified by `pad` (for `Gⁿ`)
/// * For example, `1 ≤ m ≤ 4` for a `pad_b→Qua4` and `1 ≤ n ≤ 8` for `pad→Qua8`
/// * The determinant of the Jacobian is calculated for `pad` (`pad` is the driver of the calculations)
/// * The number of integration points must consider the nodes of `pad` and the expected order of the whole integrand
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠       →
/// Kᵐⁿ = │ Nbᵐ s Gⁿ dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///        nip-1     →     →       →       →
/// Kᵐⁿⱼ ≈   Σ   Nbᵐ(ιᵖ) s(ιᵖ) Gⁿⱼ(ιᵖ) |J|(ιᵖ) wᵖ
///         p=0
/// ```
///
/// # Output
///
/// ```text
///     ┌                                        ┐
///     | K⁰⁰₀ K⁰⁰₁ K⁰¹₀ K⁰¹₁ K⁰²₀ K⁰²₁ ··· K⁰ⁿⱼ |  ⟸  ii0
///     | K¹⁰₀ K¹⁰₁ K¹¹₀ K¹¹₁ K¹²₀ K¹²₁ ··· K¹ⁿⱼ |
/// K = | K²⁰₀ K²⁰₁ K²¹₀ K²¹₁ K²²₀ K²²₁ ··· K²ⁿⱼ |
///     |  ··   ··   ··   ··   ··   ··  ···  ··  |    (pad_b)
///     | Kᵐ⁰₀ Kᵐ⁰₁ Kᵐ¹₀ Kᵐ¹₁ Kᵐ²₀ Kᵐ²₁ ··· Kᵐⁿⱼ |  ⟸  ii
///     └                                        ┘
///        ⇑                                  ⇑
///       jj0                          (pad)  jj := j + n ⋅ space_ndim
///
/// n = jj / space_ndim
/// j = jj % space_ndim
/// ```
///
/// * `kk` -- A matrix containing all `Kᵐⁿⱼ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `j` corresponds to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + pad_b.nnode` and `ncol(K) ≥ jj0 + pad.nnode ⋅ space_ndim`
/// * `pad_b` -- Lower-order scratchpad (modified) to compute Nb
/// * `pad` -- Driver scratchpad (modified) to compute G
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
pub fn mat_04_nsg<F>(
    kk: &mut Matrix,
    pad_b: &mut Scratchpad,
    pad: &mut Scratchpad,
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
    if nrow_kk < ii0 + nnode_b {
        return Err("nrow(K) must be ≥ ii0 + pad_b.nnode");
    }
    if ncol_kk < jj0 + nnode * space_ndim {
        return Err("ncol(K) must be ≥ jj0 + pad.nnode ⋅ space_ndim");
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
        let nnb = &pad_b.interp;
        let g = &pad.gradient;
        if space_ndim == 2 {
            for m in 0..nnode_b {
                for n in 0..nnode {
                    kk[ii0 + m][jj0 + 0 + n * 2] += nnb[m] * val * g[n][0];
                    kk[ii0 + m][jj0 + 1 + n * 2] += nnb[m] * val * g[n][1];
                }
            }
        } else {
            for m in 0..nnode_b {
                for n in 0..nnode {
                    kk[ii0 + m][jj0 + 0 + n * 3] += nnb[m] * val * g[n][0];
                    kk[ii0 + m][jj0 + 1 + n * 3] += nnb[m] * val * g[n][1];
                    kk[ii0 + m][jj0 + 2 + n * 3] += nnb[m] * val * g[n][2];
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
        let mut kk = Matrix::new(4, 8 * 2);
        assert_eq!(
            integ::mat_04_nsg(&mut kk, &mut pad_b, &mut pad, 1, 0, false, &[], |_| Ok(0.0)).err(),
            Some("nrow(K) must be ≥ ii0 + pad_b.nnode")
        );
        assert_eq!(
            integ::mat_04_nsg(&mut kk, &mut pad_b, &mut pad, 0, 1, false, &[], |_| Ok(0.0)).err(),
            Some("ncol(K) must be ≥ jj0 + pad.nnode ⋅ space_ndim")
        );
    }

    #[test]
    fn mat_04_nsg_works() {
        let (a, b) = (2.0, 3.0);
        let mut pad_b = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut kk = Matrix::new(4, 8 * 2);
        let ana = AnalyticalQua8::new(a, b);
        let s = 9.0;
        let kk_correct = ana.integ_nbsg(s);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [4, 9].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_04_nsg(&mut kk, &mut pad_b, &mut pad, 0, 0, true, ips, |_| Ok(s)).unwrap();
            // println!("{:.2}", kk);
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }
}
