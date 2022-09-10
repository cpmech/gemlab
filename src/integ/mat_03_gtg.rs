use super::CommonArgs;
use crate::StrError;
use russell_lab::math::SQRT_2;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the gradient(G) dot tensor(T) dot gradient(G) integration case 03 (e.g., conductivity matrix)
///
/// Callback function: `f(T, p, N, G)`
///
/// Conductivity coefficients:
///
/// ```text
///       ⌠ →        →
/// Kᵐⁿ = │ Gᵐ ⋅ T ⋅ Gⁿ α dΩ
///       ⌡      ▔
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///       nip-1 →  →       →     →  →       →
/// Kᵐⁿ ≈   Σ   Gᵐ(ιᵖ) ⋅ T(ιᵖ) ⋅ Gⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
///        p=0           ▔
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
/// * `fn_tt` -- Function `f(T,p,N,G)` that computes `T(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions N(ιᵖ), and gradients G(ιᵖ). `T` is set for `space_ndim`.
pub fn mat_03_gtg<F>(kk: &mut Matrix, args: &mut CommonArgs, mut fn_tt: F) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor2, usize, &Vector, &Matrix) -> Result<(), StrError>,
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

    // allocate auxiliary tensor
    let mut tt = Tensor2::new(true, space_ndim == 2);

    // clear output matrix
    if args.clear {
        kk.fill(0.0);
    }

    // loop over integration points
    let s = SQRT_2;
    for p in 0..args.ips.len() {
        // ksi coordinates and weight
        let iota = &args.ips[p];
        let weight = args.ips[p][3];

        // calculate Jacobian and gradient
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // G

        // calculate T tensor
        let nn = &args.pad.interp;
        let gg = &args.pad.gradient;
        fn_tt(&mut tt, p, nn, gg)?;

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
        let t = &tt.vec;
        if space_ndim == 2 {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + m][jj0 + n] += c
                        * (gg[n][1] * (t[1] * gg[m][1] + (t[3] * gg[m][0]) / s)
                            + gg[n][0] * (t[0] * gg[m][0] + (t[3] * gg[m][1]) / s));
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + m][jj0 + n] += c
                        * (gg[n][2] * (t[2] * gg[m][2] + (t[5] * gg[m][0]) / s + (t[4] * gg[m][1]) / s)
                            + gg[n][1] * (t[1] * gg[m][1] + (t[3] * gg[m][0]) / s + (t[4] * gg[m][2]) / s)
                            + gg[n][0] * (t[0] * gg[m][0] + (t[3] * gg[m][1]) / s + (t[5] * gg[m][2]) / s));
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
    use crate::integ::{
        self, AnalyticalQua4, AnalyticalQua8, AnalyticalTet4, AnalyticalTri3, CommonArgs, IP_LIN_LEGENDRE_1,
        IP_TRI_INTERNAL_1,
    };
    use russell_chk::vec_approx_eq;
    use russell_lab::{Matrix, Vector};
    use russell_tensor::{copy_tensor2, Tensor2};

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(2, 2);
        let mut tt = Tensor2::new(true, true);
        let nn = Vector::new(0);
        let gg = Matrix::new(0, 0);
        let f = |_tt: &mut Tensor2, _p: usize, _nn: &Vector, _gg: &Matrix| Ok(());
        f(&mut tt, 0, &nn, &gg).unwrap();
        let mut args = CommonArgs::new(&mut pad, &[]);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_03_gtg(&mut kk, &mut args, f).err(),
            Some("nrow(K) must be ≥ ii0 + nnode")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_03_gtg(&mut kk, &mut args, f).err(),
            Some("ncol(K) must be ≥ jj0 + nnode")
        );
        args.jj0 = 0;
        // more errors
        args.ips = &IP_LIN_LEGENDRE_1;
        assert_eq!(
            integ::mat_03_gtg(&mut kk, &mut args, f).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_qua4(0.0, 0.0, 1.0, 1.0);
        let mut kk = Matrix::new(4, 4);
        let mut args = CommonArgs::new(&mut pad, &IP_TRI_INTERNAL_1);
        assert_eq!(
            integ::mat_03_gtg(&mut kk, &mut args, |_, _, _, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn tri3_works() {
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let ana = AnalyticalTri3::new(&pad);
        let (kx, ky) = (2.5, 3.8);
        let kk_correct = ana.mat_03_gtg(kx, ky);
        let class = pad.kind.class();
        let tolerances = [1e-15, 1e-15, 1e-15];
        let selection: Vec<_> = [1, 3, 7].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_03_gtg(&mut kk, &mut args, |tt, _, _, _| {
                tt.sym_set(0, 0, kx);
                tt.sym_set(1, 1, ky);
                Ok(())
            })
            .unwrap();
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn qua4_works() {
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua4(2.0, 1.0, a, b);
        let mut kk = Matrix::new(4, 4);
        let ana = AnalyticalQua4::new(a, b);
        let (kx, ky) = (2.5, 3.8);
        let kk_correct = ana.mat_03_gtg(kx, ky);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_03_gtg(&mut kk, &mut args, |tt, _, _, _| {
                tt.sym_set(0, 0, kx);
                tt.sym_set(1, 1, ky);
                Ok(())
            })
            .unwrap();
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn qua8_works() {
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua8(2.0, 1.0, a, b);
        let mut kk = Matrix::new(8, 8);
        let ana = AnalyticalQua8::new(a, b);
        let (kx, ky) = (2.5, 3.8);
        let kk_correct = ana.mat_03_gtg(kx, ky);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [9, 16].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_03_gtg(&mut kk, &mut args, |tt, _, _, _| {
                tt.sym_set(0, 0, kx);
                tt.sym_set(1, 1, ky);
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
        #[rustfmt::skip]
        let sig = Tensor2::from_matrix(&[
            [1.1, 1.2, 1.3],
            [1.2, 2.2, 2.3],
            [1.3, 2.3, 3.3]],
        true, false).unwrap();
        let kk_correct = ana.mat_03_gtg(&sig);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_03_gtg(&mut kk, &mut args, |tt, _, _, _| {
                copy_tensor2(tt, &sig)?;
                Ok(())
            })
            .unwrap();
            // println!("{}", kk);
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }
}
