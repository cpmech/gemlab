use super::CommonArgs;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the shape(N) times vector(V) dot gradient(B) integration case 09 (e.g., variable density matrix)
///
/// Callback function: `f(v, p, N, B)`
///
/// Variable density coefficients:
///
/// ```text
///       ⌠    →   →
/// Kᵐⁿ = │ Nᵐ v ⊗ Bⁿ α dΩ
/// ▔     ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///         nip-1    →   →  →   →   →       →
/// Kᵐⁿᵢⱼ ≈   Σ   Nᵐ(ιᵖ) vᵢ(ιᵖ) Bⁿⱼ(ιᵖ) |J|(ιᵖ) wᵖ α
///          p=0
/// ```
///
/// # Results
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
/// # Arguments
///
/// * `kk` -- A matrix containing all `Kᵐⁿᵢⱼ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `i` and `j` correspond to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode ⋅ space_ndim` and `ncol(K) ≥ jj0 + nnode ⋅ space_ndim`.
/// * `args` --- Common arguments
/// * `fn_v` -- Function `f(v,p,N,B)` that computes `v(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions N(ιᵖ), and gradients B(ιᵖ). `v.dim() = space_ndim`.
pub fn mat_09_nvb<F>(kk: &mut Matrix, args: &mut CommonArgs, mut fn_v: F) -> Result<(), StrError>
where
    F: FnMut(&mut Vector, usize, &Vector, &Matrix) -> Result<(), StrError>,
{
    // check
    let (space_ndim, nnode) = args.pad.xxt.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if nrow_kk < ii0 + nnode * space_ndim {
        return Err("nrow(K) must be ≥ ii0 + nnode ⋅ space_ndim");
    }
    if ncol_kk < jj0 + nnode * space_ndim {
        return Err("ncol(K) must be ≥ jj0 + nnode ⋅ space_ndim");
    }

    // allocate auxiliary vector
    let mut v = Vector::new(space_ndim);

    // clear output matrix
    if args.clear {
        kk.fill(0.0);
    }

    // loop over integration points
    for p in 0..args.ips.len() {
        // ksi coordinates and weight
        let iota = &args.ips[p];
        let weight = args.ips[p][3];

        // calculate interpolation functions, Jacobian and gradient
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B

        // calculate v
        let nn = &args.pad.interp;
        let gg = &args.pad.gradient;
        fn_v(&mut v, p, nn, gg)?;

        // calculate coefficient
        let c = if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt[0][m];
            }
            det_jac * weight * args.alpha * r
        } else {
            det_jac * weight * args.alpha
        };

        // add contribution to K matrix
        if space_ndim == 2 {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + 0 + m * 2][jj0 + 0 + n * 2] += c * nn[m] * v[0] * gg[n][0];
                    kk[ii0 + 0 + m * 2][jj0 + 1 + n * 2] += c * nn[m] * v[0] * gg[n][1];

                    kk[ii0 + 1 + m * 2][jj0 + 0 + n * 2] += c * nn[m] * v[1] * gg[n][0];
                    kk[ii0 + 1 + m * 2][jj0 + 1 + n * 2] += c * nn[m] * v[1] * gg[n][1];
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + 0 + m * 3][jj0 + 0 + n * 3] += c * nn[m] * v[0] * gg[n][0];
                    kk[ii0 + 0 + m * 3][jj0 + 1 + n * 3] += c * nn[m] * v[0] * gg[n][1];
                    kk[ii0 + 0 + m * 3][jj0 + 2 + n * 3] += c * nn[m] * v[0] * gg[n][2];

                    kk[ii0 + 1 + m * 3][jj0 + 0 + n * 3] += c * nn[m] * v[1] * gg[n][0];
                    kk[ii0 + 1 + m * 3][jj0 + 1 + n * 3] += c * nn[m] * v[1] * gg[n][1];
                    kk[ii0 + 1 + m * 3][jj0 + 2 + n * 3] += c * nn[m] * v[1] * gg[n][2];

                    kk[ii0 + 2 + m * 3][jj0 + 0 + n * 3] += c * nn[m] * v[2] * gg[n][0];
                    kk[ii0 + 2 + m * 3][jj0 + 1 + n * 3] += c * nn[m] * v[2] * gg[n][1];
                    kk[ii0 + 2 + m * 3][jj0 + 2 + n * 3] += c * nn[m] * v[2] * gg[n][2];
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
    use crate::integ::{self, AnalyticalTet4, AnalyticalTri3, CommonArgs, IP_LIN_LEGENDRE_1, IP_TRI_INTERNAL_1};
    use russell_chk::vec_approx_eq;
    use russell_lab::{Matrix, Vector};

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(4, 4);
        let mut vv = Vector::new(0);
        let nn = Vector::new(0);
        let gg = Matrix::new(0, 0);
        let f = |_: &mut Vector, _: usize, _: &Vector, _: &Matrix| Ok(());
        f(&mut vv, 0, &nn, &gg).unwrap();
        let mut args = CommonArgs::new(&mut pad, &[]);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_09_nvb(&mut kk, &mut args, f).err(),
            Some("nrow(K) must be ≥ ii0 + nnode ⋅ space_ndim")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_09_nvb(&mut kk, &mut args, f).err(),
            Some("ncol(K) must be ≥ jj0 + nnode ⋅ space_ndim")
        );
        args.jj0 = 0;
        // more errors
        args.ips = &IP_LIN_LEGENDRE_1;
        assert_eq!(
            integ::mat_09_nvb(&mut kk, &mut args, f).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(6, 6);
        let mut args = CommonArgs::new(&mut pad, &IP_TRI_INTERNAL_1);
        assert_eq!(
            integ::mat_09_nvb(&mut kk, &mut args, |_, _, _, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn mat_09_nvb_tri3_works() {
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3 * 2, 3 * 2);
        let ana = AnalyticalTri3::new(&pad);
        // constant
        let (v0, v1) = (2.0, 3.0);
        let kk_correct = ana.mat_09_nvb(v0, v1);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [3].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_09_nvb(&mut kk, &mut args, |v, _, _, _| {
                v[0] = v0;
                v[1] = v1;
                Ok(())
            })
            .unwrap();
            // println!("{}", kk);
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn mat_09_nvb_tet4_works() {
        let mut pad = aux::gen_pad_tet4();
        let mut kk = Matrix::new(4 * 3, 4 * 3);
        let ana = AnalyticalTet4::new(&pad);
        // constant
        let (v0, v1, v2) = (2.0, 3.0, 4.0);
        let kk_correct = ana.mat_09_nvb(v0, v1, v2);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_09_nvb(&mut kk, &mut args, |v, _, _, _| {
                v[0] = v0;
                v[1] = v1;
                v[2] = v2;
                Ok(())
            })
            .unwrap();
            // println!("{}", kk);
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }
}
