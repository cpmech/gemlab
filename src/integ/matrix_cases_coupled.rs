use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::Matrix;

/// Implements the shape(Nb) time scalar(S) time gradient(G) integration case with different shapes (e.g., coupling matrix)
///
/// **Note:** `m` ranges over the number of nodes of the lower-order shape specified by `pad_b`,
/// corresponding to `Nbᵐ`, and `n` ranges over the number of the "driver" shape specified by `pad`,
/// corresponding to `Gⁿ`. For example, `m ∈ [1,4]` of a `Qua4` and `n ∈ [1,8]` of `Qua8`.
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
/// * `pad` -- "Driver" scratchpad (modified) to compute G
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output matrix where to add components.
/// * `jj0` -- Stride marking the first column in the output matrix where to add components.
/// * `clear_kk` -- Fills `kk` matrix with zeros, otherwise accumulate values into `kk`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_s` -- Function `f(p)` that computes `s(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
///
/// # Warning and additional notes
///
/// The two [crate::shapes::Scratchpad]s mut be compatible with `pad_b` being the lower-order version of `pad`.
/// For example, if `pad` corresponds to `Qua8`, then `pad_b` must be a `Qua4`. Otherwise, **calculation errors may occur**.
///
/// Note: The determinant of the Jacobian is calculated for `pad`; i.e., `pad` is the **driver** of the calculations.
/// Also, the number of integration points must be consider the number of nodes of `pad` and the expected order of the whole integrand.
pub fn mat_coupling_nbsg<F>(
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
    let nnode = pad.interp.dim();
    let space_ndim = pad.xmax.len();
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

/// Implements the gradient(G) time tensor(T) time shape(N) integration case with different shapes (e.g., coupling matrix)
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠ →
/// Kᵐⁿ = │ GBᵐ ⋅ T Nⁿ dΩ
///       ⌡       ▔
///       Ωₑ
/// ```
pub fn mat_coupling_gtn() -> Result<(), StrError> {
    Err("mat_coupling_gtn: TODO")
}

/// Implements the shape(N) time vector(V) time shape(N) integration case with different shapes (e.g., coupling matrix)
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠    →
/// Kᵐⁿ = │ Nᵐ v NBⁿ dΩ
///       ⌡
///       Ωₑ
/// ```
pub fn mat_coupling_nvn() -> Result<(), StrError> {
    Err("mat_coupling_nvn: TODO")
}

/// Implements the gradient(G) time scalar(S) time shape(Nb) integration case with different shapes (e.g., coupling matrix)
///
/// **Note:** `m` ranges over the number of nodes of the "driver" shape specified by `pad`,
/// corresponding to `Gᵐ`, and `n` ranges over the number of the lower-order shape specified by `pad_b`,
/// corresponding to `Nbⁿ`. For example, `m ∈ [1,8]` of a `Qua8` and `n ∈ [1,4]` of `Qua4`.
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
/// * `pad` -- "Driver" scratchpad (modified) to compute G
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
/// # Warning and additional notes
///
/// The two [crate::shapes::Scratchpad]s mut be compatible with `pad_b` being the lower-order version of `pad`.
/// For example, if `pad` corresponds to `Qua8`, then `pad_b` must be a `Qua4`. Otherwise, **calculation errors may occur**.
///
/// Note: The determinant of the Jacobian is calculated for `pad`; i.e., `pad` is the **driver** of the calculations.
/// Also, the number of integration points must be consider the number of nodes of `pad` and the expected order of the whole integrand.
pub fn mat_coupling_gsnb<F>(
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
    let nnode = pad.interp.dim();
    let nnode_b = pad_b.interp.dim();
    let space_ndim = pad.xmax.len();
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
        let mut kk = Matrix::new(4, 8 * 2);
        assert_eq!(
            integ::mat_coupling_nbsg(&mut kk, &mut pad_b, &mut pad, 1, 0, false, &[], |_| Ok(0.0)).err(),
            Some("nrow(K) must be ≥ ii0 + pad_b.nnode")
        );
        assert_eq!(
            integ::mat_coupling_nbsg(&mut kk, &mut pad_b, &mut pad, 0, 1, false, &[], |_| Ok(0.0)).err(),
            Some("ncol(K) must be ≥ jj0 + pad.nnode ⋅ space_ndim")
        );
        let mut kk = Matrix::new(8 * 2, 4);
        assert_eq!(
            integ::mat_coupling_gsnb(&mut kk, &mut pad, &mut pad_b, 1, 0, false, &[], |_| Ok(0.0)).err(),
            Some("nrow(K) must be ≥ ii0 + pad.nnode ⋅ space_ndim")
        );
        assert_eq!(
            integ::mat_coupling_gsnb(&mut kk, &mut pad, &mut pad_b, 0, 1, false, &[], |_| Ok(0.0)).err(),
            Some("ncol(K) must be ≥ jj0 + pad_b.nnode")
        );
    }
}
