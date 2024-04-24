use super::CommonArgs;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the shape(N) times scalar(S) times shape(N) integration case 01 (e.g., diffusion matrix)
///
/// Callback function: `s ← f(p, N, B)`
///
/// Diffusion coefficients:
///
/// ```text
///       ⌠
/// Kᵐⁿ = │ Nᵐ s Nⁿ α dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///       nip-1    →     →      →       →
/// Kᵐⁿ ≈   Σ   Nᵐ(ιᵖ) s(ιᵖ) Nⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
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
/// * `fn_s` -- Function `f(p,N,B)→s` that computes `s(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions N(ιᵖ), and gradients B(ιᵖ).
pub fn mat_01_nsn<F>(kk: &mut Matrix, args: &mut CommonArgs, mut fn_s: F) -> Result<(), StrError>
where
    F: FnMut(usize, &Vector, &Matrix) -> Result<f64, StrError>,
{
    // check
    let nnode = args.pad.interp.dim();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if nrow_kk < ii0 + nnode {
        return Err("nrow(K) must be ≥ ii0 + nnode");
    }
    if ncol_kk < jj0 + nnode {
        return Err("ncol(K) must be ≥ jj0 + nnode");
    }

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

        // calculate s
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        let s = fn_s(p, nn, bb)?;

        // calculate coefficient
        let c = if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt.get(0, m);
            }
            s * det_jac * weight * args.alpha * r
        } else {
            s * det_jac * weight * args.alpha
        };

        // add contribution to K matrix
        for m in 0..nnode {
            for n in 0..nnode {
                kk.add(ii0 + m, jj0 + n, nn[m] * c * nn[n]);
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
    use russell_lab::{mat_approx_eq, Matrix, Vector};

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(2, 2);
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let f = |_p: usize, _nn: &Vector, _bb: &Matrix| Ok(0.0);
        assert_eq!(f(0, &nn, &bb).unwrap(), 0.0);
        let mut args = CommonArgs::new(&mut pad, &[]);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_01_nsn(&mut kk, &mut args, f).err(),
            Some("nrow(K) must be ≥ ii0 + nnode")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_01_nsn(&mut kk, &mut args, f).err(),
            Some("ncol(K) must be ≥ jj0 + nnode")
        );
        args.jj0 = 0;
        // more errors
        args.ips = &IP_LIN_LEGENDRE_1;
        assert_eq!(
            integ::mat_01_nsn(&mut kk, &mut args, f).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let mut args = CommonArgs::new(&mut pad, &IP_TRI_INTERNAL_1);
        assert_eq!(
            integ::mat_01_nsn(&mut kk, &mut args, |_, _, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn tri3_works() {
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let s = 12.0;
        let ana = AnalyticalTri3::new(&pad);
        let kk_correct = ana.mat_01_nsn(s, 1.0);
        let class = pad.kind.class();
        let tolerances = [8.34, 1e-14, 1e-14, 1e-14, 1e-12, 1e-13]; // note how bad rule-1 integ is here
        let selection: Vec<_> = [1, 3, 6, 7, 12, 16]
            .iter()
            .map(|n| integ::points(class, *n).unwrap())
            .collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_01_nsn(&mut kk, &mut args, |_, _, _| Ok(s)).unwrap();
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn qua4_works() {
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua4(2.0, 1.0, a, b);
        let mut kk = Matrix::new(4, 4);
        let s = 12.0;
        let ana = AnalyticalQua4::new(a, b);
        let kk_correct = ana.mat_01_nsn(s, 1.0);
        let class = pad.kind.class();
        let tolerances = [7.01, 1e-14, 1e-14, 1e-14]; // note how bad rule-1 integ is here
        let selection: Vec<_> = [1, 4, 9, 16]
            .iter()
            .map(|n| integ::points(class, *n).unwrap())
            .collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_01_nsn(&mut kk, &mut args, |_, _, _| Ok(s)).unwrap();
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn qua8_works() {
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua8(2.0, 1.0, a, b);
        let mut kk = Matrix::new(8, 8);
        let s = 3.0;
        let ana = AnalyticalQua8::new(a, b);
        let kk_correct = ana.mat_01_nsn(s, 1.0);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [9, 16].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_01_nsn(&mut kk, &mut args, |_, _, _| Ok(s)).unwrap();
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn tet4_works() {
        let mut pad = aux::gen_pad_tet4();
        let mut kk = Matrix::new(4, 4);
        let s = 3.0;
        let ana = AnalyticalTet4::new(&pad);
        let kk_correct = ana.mat_01_nsn(s);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_01_nsn(&mut kk, &mut args, |_, _, _| Ok(s)).unwrap();
            // println!("{}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }
}
