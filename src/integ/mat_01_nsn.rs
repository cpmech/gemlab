use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the shape(N) times scalar(S) times shape(N) integration case 01 (e.g., diffusion matrix)
///
/// Diffusion coefficients:
///
/// ```text
///       ⌠
/// Kᵐⁿ = │ Nᵐ s Nⁿ dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///       nip-1    →     →      →       →
/// Kᵐⁿ ≈   Σ   Nᵐ(ιᵖ) s(ιᵖ) Nⁿ(ιᵖ) |J|(ιᵖ) wᵖ
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
/// * `fn_s` -- Function `f(p,N,G)` that computes `s(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions N(ιᵖ), and gradients G(ιᵖ).
pub fn mat_01_nsn<F>(
    kk: &mut Matrix,
    pad: &mut Scratchpad,
    ii0: usize,
    jj0: usize,
    clear_kk: bool,
    ips: IntegPointData,
    mut fn_s: F,
) -> Result<(), StrError>
where
    F: FnMut(usize, &Vector, &Matrix) -> Result<f64, StrError>,
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
        (pad.fn_interp)(&mut pad.interp, iota); // N
        let det_jac = pad.calc_gradient(iota)?; // G

        // calculate s
        let nn = &pad.interp;
        let gg = &pad.gradient;
        let s = fn_s(p, nn, gg)?;

        // add contribution to K matrix
        let val = s * det_jac * weight;
        for m in 0..nnode {
            for n in 0..nnode {
                kk[ii0 + m][jj0 + n] += nn[m] * val * nn[n];
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
        self, AnalyticalQua4, AnalyticalQua8, AnalyticalTet4, AnalyticalTri3, IP_LIN_LEGENDRE_1, IP_TRI_INTERNAL_1,
    };
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Matrix;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(2, 2);
        assert_eq!(
            integ::mat_01_nsn(&mut kk, &mut pad, 1, 0, false, &[], |_, _, _| Ok(0.0)).err(),
            Some("nrow(K) must be ≥ ii0 + nnode")
        );
        assert_eq!(
            integ::mat_01_nsn(&mut kk, &mut pad, 0, 1, false, &[], |_, _, _| Ok(0.0)).err(),
            Some("ncol(K) must be ≥ jj0 + nnode")
        );
        // more errors
        assert_eq!(
            integ::mat_01_nsn(&mut kk, &mut pad, 0, 0, false, &IP_LIN_LEGENDRE_1, |_, _, _| Ok(0.0)).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        assert_eq!(
            integ::mat_01_nsn(&mut kk, &mut pad, 0, 0, false, &IP_TRI_INTERNAL_1, |_, _, _| Err(
                "stop"
            ))
            .err(),
            Some("stop")
        );
    }

    #[test]
    fn mat_01_nsn_tri3_works() {
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
            integ::mat_01_nsn(&mut kk, &mut pad, 0, 0, true, ips, |_, _, _| Ok(s)).unwrap();
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn mat_01_nsn_qua4_works() {
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
            integ::mat_01_nsn(&mut kk, &mut pad, 0, 0, true, ips, |_, _, _| Ok(s)).unwrap();
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn mat_01_nsn_qua8_works() {
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
            integ::mat_01_nsn(&mut kk, &mut pad, 0, 0, true, ips, |_, _, _| Ok(s)).unwrap();
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn mat_01_nsn_tet4_works() {
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
            integ::mat_01_nsn(&mut kk, &mut pad, 0, 0, true, ips, |_, _, _| Ok(s)).unwrap();
            // println!("{}", kk);
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }
}
