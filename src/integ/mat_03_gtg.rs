use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::math::SQRT_2;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the gradient(G) dot tensor(T) dot gradient(G) integration case 03 (e.g., conductivity matrix)
///
/// Callback function: `α ← f(T, p, N, G)`
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
/// * `fn_tt` -- Function `f(T,p,N,G)→α` that computes `T(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions N(ιᵖ), and gradients G(ιᵖ). `T` is set for `space_ndim`.
///   `fn_tt` returns α that can accommodate plane-strain simulations.
///   **NOTE:** the value α is ignored if axisymmetric = true, because the radius is calculated and used instead.
pub fn mat_03_gtg<F>(
    kk: &mut Matrix,
    pad: &mut Scratchpad,
    ii0: usize,
    jj0: usize,
    clear_kk: bool,
    axisymmetric: bool,
    ips: IntegPointData,
    mut fn_tt: F,
) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor2, usize, &Vector, &Matrix) -> Result<f64, StrError>,
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
        (pad.fn_interp)(&mut pad.interp, iota); // N
        let det_jac = pad.calc_gradient(iota)?; // G

        // calculate T tensor
        let nn = &pad.interp;
        let gg = &pad.gradient;
        let alpha = fn_tt(&mut tt, p, nn, gg)?;

        // calculate coefficient
        let c = if axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * pad.xxt[0][m];
            }
            r * det_jac * weight
        } else {
            alpha * det_jac * weight
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
        self, AnalyticalQua4, AnalyticalQua8, AnalyticalTet4, AnalyticalTri3, IP_LIN_LEGENDRE_1, IP_TRI_INTERNAL_1,
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
        let f = |_tt: &mut Tensor2, _p: usize, _nn: &Vector, _gg: &Matrix| Ok(1.0);
        assert_eq!(f(&mut tt, 0, &nn, &gg).unwrap(), 1.0);
        let (clear, axis) = (true, false);
        assert_eq!(
            integ::mat_03_gtg(&mut kk, &mut pad, 1, 0, clear, axis, &[], f).err(),
            Some("nrow(K) must be ≥ ii0 + nnode")
        );
        assert_eq!(
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 1, clear, axis, &[], f).err(),
            Some("ncol(K) must be ≥ jj0 + nnode")
        );
        // more errors
        let ips = &IP_LIN_LEGENDRE_1;
        assert_eq!(
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 0, clear, axis, ips, f).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_qua4(0.0, 0.0, 1.0, 1.0);
        let mut kk = Matrix::new(4, 4);
        let ips = &IP_TRI_INTERNAL_1;
        assert_eq!(
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 0, clear, axis, ips, |_, _, _, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn tri3_works() {
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let (clear, axis) = (true, false);
        let ana = AnalyticalTri3::new(&pad);
        let (kx, ky) = (2.5, 3.8);
        let kk_correct = ana.mat_03_gtg(kx, ky);
        let class = pad.kind.class();
        let tolerances = [1e-15, 1e-15, 1e-15];
        let selection: Vec<_> = [1, 3, 7].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 0, clear, axis, ips, |tt, _, _, _| {
                tt.sym_set(0, 0, kx);
                tt.sym_set(1, 1, ky);
                Ok(1.0)
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
        let (clear, axis) = (true, false);
        let ana = AnalyticalQua4::new(a, b);
        let (kx, ky) = (2.5, 3.8);
        let kk_correct = ana.mat_03_gtg(kx, ky);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 0, clear, axis, ips, |tt, _, _, _| {
                tt.sym_set(0, 0, kx);
                tt.sym_set(1, 1, ky);
                Ok(1.0)
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
        let (clear, axis) = (true, false);
        let ana = AnalyticalQua8::new(a, b);
        let (kx, ky) = (2.5, 3.8);
        let kk_correct = ana.mat_03_gtg(kx, ky);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [9, 16].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 0, clear, axis, ips, |tt, _, _, _| {
                tt.sym_set(0, 0, kx);
                tt.sym_set(1, 1, ky);
                Ok(1.0)
            })
            .unwrap();
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn tet4_works() {
        let mut pad = aux::gen_pad_tet4();
        let mut kk = Matrix::new(4, 4);
        let (clear, axis) = (true, false);
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
            integ::mat_03_gtg(&mut kk, &mut pad, 0, 0, clear, axis, ips, |tt, _, _, _| {
                copy_tensor2(tt, &sig)?;
                Ok(1.0)
            })
            .unwrap();
            // println!("{}", kk);
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }
}
