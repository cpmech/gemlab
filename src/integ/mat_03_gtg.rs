use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::Matrix;
use russell_tensor::Tensor2;

/// Implements the gradient(G) dot tensor(T) dot gradient(G) integration case 03 (e.g., conductivity matrix)
///
/// Conductivity coefficients:
///
/// ```text
///       ⌠ →        →
/// Kᵐⁿ = │ Gᵐ ⋅ T ⋅ Gⁿ dΩ
///       ⌡      ▔
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///       nip-1 →  →       →     →  →       →
/// Kᵐⁿ ≈   Σ   Gᵐ(ιᵖ) ⋅ T(ιᵖ) ⋅ Gⁿ(ιᵖ) |J|(ιᵖ) wᵖ
///        p=0           ▔
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
/// * `fn_tt` -- Function `f(T,p)` that computes `T(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
pub fn mat_03_gtg<F>(
    kk: &mut Matrix,
    pad: &mut Scratchpad,
    ii0: usize,
    jj0: usize,
    clear_kk: bool,
    ips: IntegPointData,
    fn_tt: F,
) -> Result<(), StrError>
where
    F: Fn(&mut Tensor2, usize) -> Result<(), StrError>,
{
    // check
    let (space_ndim, nnode) = pad.xxt.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    if nrow_kk < ii0 + nnode {
        return Err("nrow(K) must be ≥ ii0 + nnode");
    }
    if ncol_kk < jj0 + nnode {
        return Err("ncol(K) must be ≥ jj0 + nnode");
    }

    // allocate auxiliary tensor
    let mut tt = Tensor2::new(true, space_ndim == 2);

    // clear output matrix
    if clear_kk {
        kk.fill(0.0);
    }

    // loop over integration points
    let s = SQRT_2;
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate Jacobian and gradient
        let det_jac = pad.calc_gradient(iota)?;

        // calculate T tensor
        fn_tt(&mut tt, p)?;

        // add contribution to K matrix
        let c = det_jac * weight;
        let g = &pad.gradient;
        let t = &tt.vec;
        if space_ndim == 2 {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + m][jj0 + n] += c
                        * (g[n][1] * (t[1] * g[m][1] + (t[3] * g[m][0]) / s)
                            + g[n][0] * (t[0] * g[m][0] + (t[3] * g[m][1]) / s));
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + m][jj0 + n] += c
                        * (g[n][2] * (t[2] * g[m][2] + (t[5] * g[m][0]) / s + (t[4] * g[m][1]) / s)
                            + g[n][1] * (t[1] * g[m][1] + (t[3] * g[m][0]) / s + (t[4] * g[m][2]) / s)
                            + g[n][0] * (t[0] * g[m][0] + (t[3] * g[m][1]) / s + (t[5] * g[m][2]) / s));
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
        self, AnalyticalQua4, AnalyticalQua8, AnalyticalTet4, AnalyticalTri3, IP_LIN_LEGENDRE_1, IP_TRI_INTERNAL_1,
    };
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::{copy_vector, Matrix};
    use russell_tensor::Tensor2;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(2, 2);
        assert_eq!(
            integ::mat_03_gtg(&mut kk, &mut pad, 1, 0, false, &[], |_, _| Ok(())).err(),
            Some("nrow(K) must be ≥ ii0 + nnode")
        );
        assert_eq!(
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 1, false, &[], |_, _| Ok(())).err(),
            Some("ncol(K) must be ≥ jj0 + nnode")
        );
        // more errors
        assert_eq!(
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 0, false, &IP_LIN_LEGENDRE_1, |_, _| Ok(())).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_qua4(0.0, 0.0, 1.0, 1.0);
        let mut kk = Matrix::new(4, 4);
        assert_eq!(
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 0, false, &IP_TRI_INTERNAL_1, |_, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn mat_03_gtg_tri3_works() {
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let ana = AnalyticalTri3::new(&pad);
        let (kx, ky) = (2.5, 3.8);
        let kk_correct = ana.integ_gtg(kx, ky);
        let class = pad.kind.class();
        let tolerances = [1e-15, 1e-15, 1e-15];
        let selection: Vec<_> = [1, 3, 7].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 0, true, ips, |tt, _| {
                tt.sym_set(0, 0, kx);
                tt.sym_set(1, 1, ky);
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn mat_03_gtg_qua4_works() {
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua4(2.0, 1.0, a, b);
        let mut kk = Matrix::new(4, 4);
        let ana = AnalyticalQua4::new(a, b);
        let (kx, ky) = (2.5, 3.8);
        let kk_correct = ana.integ_gtg(kx, ky);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 0, true, ips, |tt, _| {
                tt.sym_set(0, 0, kx);
                tt.sym_set(1, 1, ky);
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn mat_03_gtg_qua8_works() {
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua8(2.0, 1.0, a, b);
        let mut kk = Matrix::new(8, 8);
        let ana = AnalyticalQua8::new(a, b);
        let (kx, ky) = (2.5, 3.8);
        let kk_correct = ana.integ_gtg(kx, ky);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [9, 16].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 0, true, ips, |tt, _| {
                tt.sym_set(0, 0, kx);
                tt.sym_set(1, 1, ky);
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn mat_03_gtg_tet4_works() {
        let mut pad = aux::gen_pad_tet4();
        let mut kk = Matrix::new(4, 4);
        let ana = AnalyticalTet4::new(&pad);
        #[rustfmt::skip]
        let sig = Tensor2::from_matrix(&[
            [1.1, 1.2, 1.3],
            [1.2, 2.2, 2.3],
            [1.3, 2.3, 3.3]],
        true, false).unwrap();
        let kk_correct = ana.integ_gtg(&sig);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 0, true, ips, |tt, _| {
                copy_vector(&mut tt.vec, &sig.vec)
            })
            .unwrap();
            // println!("{}", kk);
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }
}
