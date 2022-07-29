use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::Matrix;

/// Implements the shape(Nb) time scalar(S) time gradient(G) integration case with different shapes (e.g., coupling matrix)
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
///     |  ··   ··   ··   ··   ··   ··  ···  ··  |
///     | Kᵐ⁰₀ Kᵐ⁰₁ Kᵐ¹₀ Kᵐ¹₁ Kᵐ²₀ Kᵐ²₁ ··· Kᵐⁿⱼ |  ⟸  ii
///     └                                        ┘
///        ⇑                                  ⇑
///       jj0                                 jj := j + n ⋅ space_ndim
///
/// n = jj / space_ndim
/// j = jj % space_ndim
/// ```
///
/// * `kk` -- A matrix containing all `Kᵐⁿⱼ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `j` corresponds to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode` and `ncol(K) ≥ jj0 + nnode ⋅ space_ndim`.
/// * `pad_b` -- (to compute Nb). Some members of the scratchpad will be modified.
/// * `pad` -- Some members of the scratchpad will be modified.
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output matrix where to add components.
/// * `jj0` -- Stride marking the first column in the output matrix where to add components.
/// * `clear_kk` -- Fills `kk` matrix with zeros, otherwise accumulate values into `kk`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_s` -- Function `f(p)` that computes `s(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
///
/// # Warning and notes
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
    let nnode = pad.interp.dim();
    let space_ndim = pad.xmax.len();
    let (nrow_kk, ncol_kk) = kk.dims();
    if nrow_kk < ii0 + nnode {
        return Err("nrow(K) must be ≥ ii0 + nnode");
    }
    if ncol_kk < jj0 + nnode * space_ndim {
        return Err("ncol(K) must be ≥ jj0 + nnode ⋅ space_ndim");
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
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + m][jj0 + 0 + n * 2] += nnb[m] * val * g[n][0];
                    kk[ii0 + m][jj0 + 1 + n * 2] += nnb[m] * val * g[n][1];
                }
            }
        } else {
            for m in 0..nnode {
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
///     |  ···  ···  ··· ···  ··· |
///     | Kᵐ⁰ᵢ Kᵐ¹ᵢ Kᵐ²ᵢ ··· Kᵐⁿᵢ |  ⟸  ii := i + m ⋅ space_ndim
///     └                         ┘
///        ⇑                  ⇑
///       jj0                 jj
///
/// m = ii / space_ndim
/// i = ii % space_ndim
/// ```
///
/// * `kk` -- A matrix containing all `Kᵐⁿᵢ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `i` corresponds to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode ⋅ space_ndim` and `ncol(K) ≥ jj0 + nnode`.
/// * `pad` -- Some members of the scratchpad will be modified.
/// * `pad_b` -- (to compute Nb). Some members of the scratchpad will be modified.
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output matrix where to add components.
/// * `jj0` -- Stride marking the first column in the output matrix where to add components.
/// * `clear_kk` -- Fills `kk` matrix with zeros, otherwise accumulate values into `kk`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_s` -- Function `f(p)` that computes `s(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
///
/// # Warning and notes
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
    let space_ndim = pad.xmax.len();
    let (nrow_kk, ncol_kk) = kk.dims();
    if nrow_kk < ii0 + nnode * space_ndim {
        return Err("nrow(K) must be ≥ ii0 + nnode ⋅ space_ndim");
    }
    if ncol_kk < jj0 + nnode {
        return Err("ncol(K) must be ≥ jj0 + nnode");
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
                for n in 0..nnode {
                    kk[ii0 + 0 + m * 2][jj0 + n] += g[m][0] * val * nnb[n];
                    kk[ii0 + 1 + m * 2][jj0 + n] += g[m][1] * val * nnb[n];
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode {
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
    use super::{mat_coupling_gsnb, mat_coupling_gtn, mat_coupling_nbsg, mat_coupling_nvn};

    #[test]
    fn functions_return_todo() {
        assert_eq!(mat_coupling_gtn().err(), Some("mat_coupling_gtn: TODO"));
        assert_eq!(mat_coupling_nvn().err(), Some("mat_coupling_nvn: TODO"));
    }
}
