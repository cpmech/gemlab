use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the gradient(G) dot vector(V) times shape(N) integration case 02 (e.g., compressibility matrix)
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
/// * `fn_v` -- Function `f(v,p,N,G)` that computes `v(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`.
///   shape functions N(ιᵖ), and gradients G(ιᵖ). `v.dim() = space_ndim`.
pub fn mat_02_gvn<F>(
    kk: &mut Matrix,
    pad: &mut Scratchpad,
    ii0: usize,
    jj0: usize,
    clear_kk: bool,
    ips: IntegPointData,
    mut fn_v: F,
) -> Result<(), StrError>
where
    F: FnMut(&mut Vector, usize, &Vector, &Matrix) -> Result<(), StrError>,
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

        // calculate interpolation functions and Jacobian
        (pad.fn_interp)(&mut pad.interp, iota); // N
        let det_jac = pad.calc_gradient(iota)?; // G

        // calculate v
        let nn = &pad.interp;
        let gg = &pad.gradient;
        fn_v(&mut v, p, nn, gg)?;

        // add contribution to K matrix
        let c = det_jac * weight;
        if space_ndim == 2 {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + m][jj0 + n] += c * (gg[m][0] * v[0] + gg[m][1] * v[1]) * nn[n];
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk[ii0 + m][jj0 + n] += c * (gg[m][0] * v[0] + gg[m][1] * v[1] + gg[m][2] * v[2]) * nn[n];
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
    use crate::integ::{self, AnalyticalTet4, AnalyticalTri3, IP_LIN_LEGENDRE_1, IP_TRI_INTERNAL_1};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Matrix;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(2, 2);
        assert_eq!(
            integ::mat_02_gvn(&mut kk, &mut pad, 1, 0, false, &[], |_, _, _, _| Ok(())).err(),
            Some("nrow(K) must be ≥ ii0 + nnode")
        );
        assert_eq!(
            integ::mat_02_gvn(&mut kk, &mut pad, 0, 1, false, &[], |_, _, _, _| Ok(())).err(),
            Some("ncol(K) must be ≥ jj0 + nnode")
        );
        // more errors
        assert_eq!(
            integ::mat_02_gvn(&mut kk, &mut pad, 0, 0, false, &IP_LIN_LEGENDRE_1, |_, _, _, _| Ok(())).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_qua4(0.0, 0.0, 1.0, 1.0);
        let mut kk = Matrix::new(4, 4);
        assert_eq!(
            integ::mat_02_gvn(&mut kk, &mut pad, 0, 0, false, &IP_TRI_INTERNAL_1, |_, _, _, _| Err(
                "stop"
            ))
            .err(),
            Some("stop")
        );
    }

    #[test]
    fn mat_02_gvn_tri3_works() {
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let ana = AnalyticalTri3::new(&pad);
        // constant
        let (v0, v1) = (2.0, 3.0);
        let kk_correct = ana.mat_02_gvn(v0, v1);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [3].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_02_gvn(&mut kk, &mut pad, 0, 0, true, ips, |v, _, _, _| {
                v[0] = v0;
                v[1] = v1;
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
        // bilinear
        let kk_correct = ana.mat_02_gvn_bilinear(&pad);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-15];
        let selection: Vec<_> = [3, 6].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = integ::points_coords(&mut pad, ips).unwrap();
            integ::mat_02_gvn(&mut kk, &mut pad, 0, 0, true, ips, |v, p, _, _| {
                v[0] = x_ips[p][0];
                v[1] = x_ips[p][1];
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn mat_02_gvn_tet4_works() {
        let mut pad = aux::gen_pad_tet4();
        let mut kk = Matrix::new(4, 4);
        let ana = AnalyticalTet4::new(&pad);
        let (v0, v1, v2) = (2.0, 3.0, 4.0);
        let kk_correct = ana.mat_02_gvn(v0, v1, v2);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_02_gvn(&mut kk, &mut pad, 0, 0, true, ips, |v, _, _, _| {
                v[0] = v0;
                v[1] = v1;
                v[2] = v2;
                Ok(())
            })
            .unwrap();
            // println!("{}", kk);
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol);
        });
    }
}
