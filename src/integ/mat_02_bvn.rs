use super::CommonArgs;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the gradient(B) dot vector(V) times shape(N) integration case 02 (e.g., compressibility matrix)
///
/// Callback function: `f(v, p, N, B)`
///
/// Compressibility coefficients:
///
/// ```text
///       ⌠ →    →
/// Kᵐⁿ = │ Bᵐ ⋅ v Nⁿ α dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///       nip-1 →  →     → →      →       →
/// Kᵐⁿ ≈   Σ   Bᵐ(ιᵖ) ⋅ v(ιᵖ) Nⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
///        p=0
/// ```
///
/// # Results
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
/// # Arguments
///
/// * `kk` -- A matrix containing all `Kᵐⁿ` values, one after another, and
///   sequentially placed as shown above. `m` and `n` are the indices of the nodes.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode` and `ncol(K) ≥ jj0 + nnode`
/// * `args` --- Common arguments
/// * `fn_v` -- Function `f(v,p,N,B)` that computes `v(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions N(ιᵖ), and gradients B(ιᵖ). `v.dim() = space_ndim`.
pub fn mat_02_bvn<F>(kk: &mut Matrix, args: &mut CommonArgs, mut fn_v: F) -> Result<(), StrError>
where
    F: FnMut(&mut Vector, usize, &Vector, &Matrix) -> Result<(), StrError>,
{
    // check
    let (space_ndim, nnode) = args.pad.xxt.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if nrow_kk < ii0 + nnode {
        return Err("nrow(K) must be ≥ ii0 + nnode");
    }
    if ncol_kk < jj0 + nnode {
        return Err("ncol(K) must be ≥ jj0 + nnode");
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

        // calculate interpolation functions and Jacobian
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B

        // calculate v
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        fn_v(&mut v, p, nn, bb)?;

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
                    kk[ii0 + m][jj0 + n] += c * (bb[m][0] * v[0] + bb[m][1] * v[1]) * nn[n];
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + m][jj0 + n] += c * (bb[m][0] * v[0] + bb[m][1] * v[1] + bb[m][2] * v[2]) * nn[n];
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
        let mut kk = Matrix::new(2, 2);
        let mut v = Vector::new(0);
        let nn = Vector::new(0);
        let gg = Matrix::new(0, 0);
        let f = |_v: &mut Vector, _p: usize, _nn: &Vector, _gg: &Matrix| Ok(());
        f(&mut v, 0, &nn, &gg).unwrap();
        let mut args = CommonArgs::new(&mut pad, &[]);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_02_bvn(&mut kk, &mut args, f).err(),
            Some("nrow(K) must be ≥ ii0 + nnode")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_02_bvn(&mut kk, &mut args, f).err(),
            Some("ncol(K) must be ≥ jj0 + nnode")
        );
        args.jj0 = 0;
        // more errors
        args.ips = &IP_LIN_LEGENDRE_1;
        assert_eq!(
            integ::mat_02_bvn(&mut kk, &mut args, f).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_qua4(0.0, 0.0, 1.0, 1.0);
        let mut kk = Matrix::new(4, 4);
        let mut args = CommonArgs::new(&mut pad, &IP_TRI_INTERNAL_1);
        assert_eq!(
            integ::mat_02_bvn(&mut kk, &mut args, |_, _, _, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn tri3_works() {
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let ana = AnalyticalTri3::new(&pad);
        // constant
        let (v0, v1) = (2.0, 3.0);
        let kk_correct = ana.mat_02_bvn(v0, v1);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [3].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_02_bvn(&mut kk, &mut args, |v, _, _, _| {
                v[0] = v0;
                v[1] = v1;
                Ok(())
            })
            .unwrap();
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
        // bilinear
        let kk_correct = ana.mat_02_bvn_bilinear(&pad);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-15];
        let selection: Vec<_> = [3, 6].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            let x_ips = integ::points_coords(&mut args.pad, ips).unwrap();
            integ::mat_02_bvn(&mut kk, &mut args, |v, p, _, _| {
                v[0] = x_ips[p][0];
                v[1] = x_ips[p][1];
                Ok(())
            })
            .unwrap();
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn tet4_works() {
        let mut pad = aux::gen_pad_tet4();
        let mut kk = Matrix::new(4, 4);
        let ana = AnalyticalTet4::new(&pad);
        let (v0, v1, v2) = (2.0, 3.0, 4.0);
        let kk_correct = ana.mat_02_bvn(v0, v1, v2);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_02_bvn(&mut kk, &mut args, |v, _, _, _| {
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
