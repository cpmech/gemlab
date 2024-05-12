use super::CommonArgs;
use crate::StrError;
use russell_lab::math::SQRT_2;
use russell_lab::{Matrix, Vector};
use russell_tensor::{Mandel, Tensor2};

/// Implements the shape(N) times tensor(T) times shape(N) integration case 08 (e.g., mass matrix)
///
/// Callback function: `f(T, p, N, B)`
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
/// * `fn_tt` -- Function `f(T,p,N,B)` that computes `T(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions N(ιᵖ), and gradients B(ιᵖ). `T` is set for `space_ndim`.
pub fn mat_08_ntn<F>(kk: &mut Matrix, args: &mut CommonArgs, mut fn_tt: F) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor2, usize, &Vector, &Matrix) -> Result<(), StrError>,
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

    // allocate auxiliary tensor
    let mut tt = Tensor2::new(Mandel::new(2 * space_ndim));

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

        // calculate interpolation functions and Jacobian
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B

        // calculate T tensor
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        fn_tt(&mut tt, p, nn, bb)?;

        // calculate coefficient
        let c = if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt.get(0, m);
            }
            det_jac * weight * args.alpha * r
        } else {
            det_jac * weight * args.alpha
        };

        // add contribution to K matrix
        let t = tt.vector();
        if space_ndim == 2 {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk.add(ii0 + 0 + m * 2, jj0 + 0 + n * 2, c * nn[m] * t[0] * nn[n]);
                    kk.add(ii0 + 0 + m * 2, jj0 + 1 + n * 2, c * nn[m] * t[3] * nn[n] / s);

                    kk.add(ii0 + 1 + m * 2, jj0 + 0 + n * 2, c * nn[m] * t[3] * nn[n] / s);
                    kk.add(ii0 + 1 + m * 2, jj0 + 1 + n * 2, c * nn[m] * t[1] * nn[n]);
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk.add(ii0 + 0 + m * 3, jj0 + 0 + n * 3, c * nn[m] * t[0] * nn[n]);
                    kk.add(ii0 + 0 + m * 3, jj0 + 1 + n * 3, c * nn[m] * t[3] * nn[n] / s);
                    kk.add(ii0 + 0 + m * 3, jj0 + 2 + n * 3, c * nn[m] * t[5] * nn[n] / s);

                    kk.add(ii0 + 1 + m * 3, jj0 + 0 + n * 3, c * nn[m] * t[3] * nn[n] / s);
                    kk.add(ii0 + 1 + m * 3, jj0 + 1 + n * 3, c * nn[m] * t[1] * nn[n]);
                    kk.add(ii0 + 1 + m * 3, jj0 + 2 + n * 3, c * nn[m] * t[4] * nn[n] / s);

                    kk.add(ii0 + 2 + m * 3, jj0 + 0 + n * 3, c * nn[m] * t[5] * nn[n] / s);
                    kk.add(ii0 + 2 + m * 3, jj0 + 1 + n * 3, c * nn[m] * t[4] * nn[n] / s);
                    kk.add(ii0 + 2 + m * 3, jj0 + 2 + n * 3, c * nn[m] * t[2] * nn[n]);
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
    use russell_lab::{mat_approx_eq, Matrix, Vector};
    use russell_tensor::{Mandel, Tensor2};

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(4, 4);
        let mut tt = Tensor2::new(Mandel::Symmetric2D);
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let f = |_tt: &mut Tensor2, _p: usize, _nn: &Vector, _bb: &Matrix| Ok(());
        f(&mut tt, 0, &nn, &bb).unwrap();
        let mut args = CommonArgs::new(&mut pad, &[]);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_08_ntn(&mut kk, &mut args, f).err(),
            Some("nrow(K) must be ≥ ii0 + nnode ⋅ space_ndim")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_08_ntn(&mut kk, &mut args, f).err(),
            Some("ncol(K) must be ≥ jj0 + nnode ⋅ space_ndim")
        );
        args.jj0 = 0;
        // more errors
        args.ips = &IP_LIN_LEGENDRE_1;
        assert_eq!(
            integ::mat_08_ntn(&mut kk, &mut args, f).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(6, 6);
        let mut args = CommonArgs::new(&mut pad, &IP_TRI_INTERNAL_1);
        assert_eq!(
            integ::mat_08_ntn(&mut kk, &mut args, |_, _, _, _| Err("stop")).err(),
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
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_08_ntn(&mut kk, &mut args, |tt, _, _, _| {
                tt.sym_set(0, 0, rho);
                tt.sym_set(1, 1, rho);
                Ok(())
            })
            .unwrap();
            // println!("{}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
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
        Mandel::Symmetric).unwrap();
        let kk_correct = ana.mat_08_ntn(&sig);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_08_ntn(&mut kk, &mut args, |tt, _, _, _| {
                tt.mirror(&sig);
                Ok(())
            })
            .unwrap();
            // println!("{}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }
}
