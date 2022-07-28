use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the gradient(G) dot vector(V) times shape(N) integration case (e.g., compressibility matrix)
///
/// Compressibility coefficients:
///
/// ```text
///       ⌠ →    →
/// Kᵐⁿ = │ Gᵐ ⋅ v Nⁿ dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///       nip-1 →  →     → →      →       →
/// Kᵐⁿ ≈   Σ   Gᵐ(ιᵖ) ⋅ v(ιᵖ) Nⁿ(ιᵖ) |J|(ιᵖ) wᵖ
///        p=0
/// ```
///
/// # Output
///
/// ```text
///     ┌                     ┐
///     | K⁰⁰ K⁰¹ K⁰² ··· K⁰ⁿ |  ⟸  ii0
///     | K¹⁰ K¹¹ K¹² ··· K¹ⁿ |
/// K = | K²⁰ K²¹ K²² ··· K²ⁿ |
///     |  ··  ··  ·· ···  ·· |
///     | Kᵐ⁰ Kᵐ¹ Kᵐ² ··· Kᵐⁿ |  ⟸  ii
///     └                     ┘
///        ⇑                ⇑
///       jj0               jj
/// ```
///
/// * `kk` -- A matrix containing all `Kᵐⁿ` values, one after another, and
///   sequentially placed as shown above. `m` and `n` are the indices of the nodes.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode` and `ncol(K) ≥ jj0 + nnode`
/// * `pad` -- Some members of the scratchpad will be modified.
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output matrix where to add components.
/// * `jj0` -- Stride marking the first column in the output matrix where to add components.
/// * `clear_kk` -- Fills `kk` matrix with zeros, otherwise accumulate values into `kk`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_v` -- Function `f(v,p)` that computes `v(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`.
///   The dim of `v` is equal to `space_ndim`.
pub fn mat_gvn<F>(
    kk: &mut Matrix,
    pad: &mut Scratchpad,
    ii0: usize,
    jj0: usize,
    clear_kk: bool,
    ips: IntegPointData,
    fn_v: F,
) -> Result<(), StrError>
where
    F: Fn(&mut Vector, usize) -> Result<(), StrError>,
{
    // check
    let nnode = pad.interp.dim();
    let (nrow_kk, ncol_kk) = kk.dims();
    if nrow_kk < ii0 + nnode {
        return Err("nrow(K) must be ≥ ii0 + nnode");
    }
    if ncol_kk < jj0 + nnode {
        return Err("ncol(K) must be ≥ jj0 + nnode");
    }

    // allocate auxiliary vector
    let space_ndim = pad.xmax.len();
    let mut v = Vector::new(space_ndim);

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
        (pad.fn_interp)(&mut pad.interp, iota);
        let det_jac = pad.calc_gradient(iota)?;

        // calculate v
        fn_v(&mut v, p)?;

        // add contribution to K matrix
        let c = det_jac * weight;
        let nn = &pad.interp;
        let g = &pad.gradient;
        if space_ndim == 2 {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + m][jj0 + n] += c * (g[m][0] * v[0] + g[m][1] * v[1]) * nn[n];
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + m][jj0 + n] += c * (g[m][0] * v[0] + g[m][1] * v[1] + g[m][2] * v[2]) * nn[n];
                }
            }
        }
    }
    Ok(())
}

/// Implements the shape(N) times vector(V) dot gradient(G) integration case (e.g., variable density matrix)
///
/// Variable density coefficients:
///
/// ```text
///       ⌠    →   →
/// Kᵐⁿ = │ Nᵐ v ⊗ Gⁿ dΩ
/// ▔     ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///         nip-1    →   →  →   →   →       →
/// Kᵐⁿᵢⱼ ≈   Σ   Nᵐ(ιᵖ) vᵢ(ιᵖ) Gⁿⱼ(ιᵖ) |J|(ιᵖ) wᵖ
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
/// * `fn_v` -- Function `f(v,p)` that computes `v(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`.
///   The dim of `v` is equal to `space_ndim`.
pub fn mat_nvg<F>(
    kk: &mut Matrix,
    pad: &mut Scratchpad,
    ii0: usize,
    jj0: usize,
    clear_kk: bool,
    ips: IntegPointData,
    fn_v: F,
) -> Result<(), StrError>
where
    F: Fn(&mut Vector, usize) -> Result<(), StrError>,
{
    // check
    let nnode = pad.interp.dim();
    let space_ndim = pad.xmax.len();
    let (nrow_kk, ncol_kk) = kk.dims();
    if nrow_kk < ii0 + nnode * space_ndim {
        return Err("nrow(K) must be ≥ ii0 + nnode ⋅ space_ndim");
    }
    if ncol_kk < jj0 + nnode * space_ndim {
        return Err("ncol(K) must be ≥ jj0 + nnode ⋅ space_ndim");
    }

    // allocate auxiliary vector
    let mut v = Vector::new(space_ndim);

    // clear output matrix
    if clear_kk {
        kk.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate interpolation functions, Jacobian and gradient
        (pad.fn_interp)(&mut pad.interp, iota);
        let det_jac = pad.calc_gradient(iota)?;

        // calculate v
        fn_v(&mut v, p)?;

        // add contribution to K matrix
        let c = det_jac * weight;
        let nn = &pad.interp;
        let g = &pad.gradient;
        if space_ndim == 2 {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + 0 + m * 2][jj0 + 0 + n * 2] += c * nn[m] * v[0] * g[n][0];
                    kk[ii0 + 0 + m * 2][jj0 + 1 + n * 2] += c * nn[m] * v[0] * g[n][1];

                    kk[ii0 + 1 + m * 2][jj0 + 0 + n * 2] += c * nn[m] * v[1] * g[n][0];
                    kk[ii0 + 1 + m * 2][jj0 + 1 + n * 2] += c * nn[m] * v[1] * g[n][1];
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + 0 + m * 3][jj0 + 0 + n * 3] += c * nn[m] * v[0] * g[n][0];
                    kk[ii0 + 0 + m * 3][jj0 + 1 + n * 3] += c * nn[m] * v[0] * g[n][1];
                    kk[ii0 + 0 + m * 3][jj0 + 2 + n * 3] += c * nn[m] * v[0] * g[n][2];

                    kk[ii0 + 1 + m * 3][jj0 + 0 + n * 3] += c * nn[m] * v[1] * g[n][0];
                    kk[ii0 + 1 + m * 3][jj0 + 1 + n * 3] += c * nn[m] * v[1] * g[n][1];
                    kk[ii0 + 1 + m * 3][jj0 + 2 + n * 3] += c * nn[m] * v[1] * g[n][2];

                    kk[ii0 + 2 + m * 3][jj0 + 0 + n * 3] += c * nn[m] * v[2] * g[n][0];
                    kk[ii0 + 2 + m * 3][jj0 + 1 + n * 3] += c * nn[m] * v[2] * g[n][1];
                    kk[ii0 + 2 + m * 3][jj0 + 2 + n * 3] += c * nn[m] * v[2] * g[n][2];
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
    use crate::integ::{self, AnalyticalTri3};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Matrix;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(2, 2);
        assert_eq!(
            integ::mat_gvn(&mut kk, &mut pad, 1, 0, false, &[], |_, _| Ok(())).err(),
            Some("nrow(K) must be ≥ ii0 + nnode")
        );
        assert_eq!(
            integ::mat_gvn(&mut kk, &mut pad, 0, 1, false, &[], |_, _| Ok(())).err(),
            Some("ncol(K) must be ≥ jj0 + nnode")
        );
        let mut kk = Matrix::new(4, 4);
        assert_eq!(
            integ::mat_nvg(&mut kk, &mut pad, 1, 0, false, &[], |_, _| Ok(())).err(),
            Some("nrow(K) must be ≥ ii0 + nnode ⋅ space_ndim")
        );
        assert_eq!(
            integ::mat_nvg(&mut kk, &mut pad, 0, 1, false, &[], |_, _| Ok(())).err(),
            Some("ncol(K) must be ≥ jj0 + nnode ⋅ space_ndim")
        );
    }

    #[test]
    fn mat_gvn_tri3_works() {
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let ana = AnalyticalTri3::new(&pad);
        // constant
        let (vx, vy) = (2.0, 3.0);
        let kk_correct = ana.integ_gvn_constant(vx, vy);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [3].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_gvn(&mut kk, &mut pad, 0, 0, true, ips, |v, _| {
                v[0] = vx;
                v[1] = vy;
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
        // bilinear
        let kk_correct = ana.integ_gvn_bilinear(&pad);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-15];
        let selection: Vec<_> = [3, 6].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = integ::points_coords(&mut pad, ips).unwrap();
            integ::mat_gvn(&mut kk, &mut pad, 0, 0, true, ips, |v, p| {
                v[0] = x_ips[p][0];
                v[1] = x_ips[p][1];
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }
}
