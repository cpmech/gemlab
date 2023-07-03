use super::CommonArgs;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the shape(Nb) times scalar(S) times gradient(B) integration case 04 (e.g., coupling matrix)
///
/// Callback function: `s ← f(p, Nb, N, B)`
///
/// **Notes:**
///
/// * `m` ranges over the number of nodes of the *lower-order* shape specified by `pad_b` (for `Nbᵐ`)
/// * `n` ranges over the number of nodes of the *driver* shape specified by `pad` (for `Bⁿ`)
/// * For example, `1 ≤ m ≤ 4` for a `pad_b→Qua4` and `1 ≤ n ≤ 8` for `pad→Qua8`
/// * The determinant of the Jacobian is calculated for `pad` (`pad` is the driver of the calculations)
/// * The number of integration points must consider the nodes of `pad` and the expected order of the whole integrand
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠       →
/// Kᵐⁿ = │ Nbᵐ s Bⁿ α dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///        nip-1     →     →       →       →
/// Kᵐⁿⱼ ≈   Σ   Nbᵐ(ιᵖ) s(ιᵖ) Bⁿⱼ(ιᵖ) |J|(ιᵖ) wᵖ α
///         p=0
/// ```
///
/// # Results
///
/// ```text
///     ┌                                        ┐
///     | K⁰⁰₀ K⁰⁰₁ K⁰¹₀ K⁰¹₁ K⁰²₀ K⁰²₁ ··· K⁰ⁿⱼ |  ⟸  ii0
///     | K¹⁰₀ K¹⁰₁ K¹¹₀ K¹¹₁ K¹²₀ K¹²₁ ··· K¹ⁿⱼ |
/// K = | K²⁰₀ K²⁰₁ K²¹₀ K²¹₁ K²²₀ K²²₁ ··· K²ⁿⱼ |
///     |  ··   ··   ··   ··   ··   ··  ···  ··  |    (pad_b)
///     | Kᵐ⁰₀ Kᵐ⁰₁ Kᵐ¹₀ Kᵐ¹₁ Kᵐ²₀ Kᵐ²₁ ··· Kᵐⁿⱼ |  ⟸  ii
///     └                                        ┘
///        ⇑                                  ⇑
///       jj0                          (pad)  jj := j + n ⋅ space_ndim
///
/// n = jj / space_ndim
/// j = jj % space_ndim
/// ```
///
/// # Arguments
///
/// * `kk` -- A matrix containing all `Kᵐⁿⱼ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `j` corresponds to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + pad_b.nnode` and `ncol(K) ≥ jj0 + pad.nnode ⋅ space_ndim`
/// * `pad_b` -- Lower-order scratchpad (modified) to compute Nb
/// * `args` --- Common arguments (`pad` is the Driver scratchpad (modified) to compute B)
/// * `fn_s` -- Function `f(p,Nb,N,B)→s` that computes `s(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions Nb(ιᵖ), shape functions N(ιᵖ), and gradients B(ιᵖ).
///
/// # Warning
///
/// The two [crate::shapes::Scratchpad]s mut be compatible, otherwise **calculation errors may occur**.
/// Therefore, `pad_b` must be either the lower-version of `pad` or have the same shape as `pad`.
pub fn mat_04_nsb<F>(
    kk: &mut Matrix,
    pad_b: &mut Scratchpad,
    args: &mut CommonArgs,
    mut fn_s: F,
) -> Result<(), StrError>
where
    F: FnMut(usize, &Vector, &Vector, &Matrix) -> Result<f64, StrError>,
{
    // check
    let nnode_b = pad_b.interp.dim();
    let (space_ndim, nnode) = args.pad.xxt.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if nrow_kk < ii0 + nnode_b {
        return Err("nrow(K) must be ≥ ii0 + pad_b.nnode");
    }
    if ncol_kk < jj0 + nnode * space_ndim {
        return Err("ncol(K) must be ≥ jj0 + pad.nnode ⋅ space_ndim");
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
        (pad_b.fn_interp)(&mut pad_b.interp, iota); // Nb
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B

        // calculate s
        let nnb = &pad_b.interp;
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        let s = fn_s(p, nnb, nn, bb)?;

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
        if space_ndim == 2 {
            for m in 0..nnode_b {
                for n in 0..nnode {
                    kk.add(ii0 + m, jj0 + 0 + n * 2, nnb[m] * c * bb.get(n, 0));
                    kk.add(ii0 + m, jj0 + 1 + n * 2, nnb[m] * c * bb.get(n, 1));
                }
            }
        } else {
            for m in 0..nnode_b {
                for n in 0..nnode {
                    kk.add(ii0 + m, jj0 + 0 + n * 3, nnb[m] * c * bb.get(n, 0));
                    kk.add(ii0 + m, jj0 + 1 + n * 3, nnb[m] * c * bb.get(n, 1));
                    kk.add(ii0 + m, jj0 + 2 + n * 3, nnb[m] * c * bb.get(n, 2));
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
    use crate::integ::{self, AnalyticalQua8, AnalyticalTet4, CommonArgs};
    use russell_chk::vec_approx_eq;
    use russell_lab::{Matrix, Vector};

    #[test]
    fn capture_some_errors() {
        let (a, b) = (2.0, 3.0);
        let mut pad_b = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut kk = Matrix::new(4, 8 * 2);
        let nnb = Vector::new(0);
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let f = |_p: usize, _nnb: &Vector, _nn: &Vector, _bb: &Matrix| Ok(0.0);
        assert_eq!(f(0, &nnb, &nn, &bb).unwrap(), 0.0);
        let mut args = CommonArgs::new(&mut pad, &[]);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_04_nsb(&mut kk, &mut pad_b, &mut args, f).err(),
            Some("nrow(K) must be ≥ ii0 + pad_b.nnode")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_04_nsb(&mut kk, &mut pad_b, &mut args, f).err(),
            Some("ncol(K) must be ≥ jj0 + pad.nnode ⋅ space_ndim")
        );
    }

    #[test]
    fn qua4_qua8_works() {
        let (a, b) = (2.0, 3.0);
        let mut pad_b = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut kk = Matrix::new(4, 8 * 2);
        let ana = AnalyticalQua8::new(a, b);
        let s = 9.0;
        let kk_correct = ana.mat_04_nsb(s);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [4, 9].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_04_nsb(&mut kk, &mut pad_b, &mut args, |_, _, _, _| Ok(s)).unwrap();
            // println!("{:.2}", kk);
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn tet4_tet4_works() {
        let mut pad_b = aux::gen_pad_tet4();
        let mut pad = pad_b.clone();
        let mut kk = Matrix::new(4, 4 * 3);
        let ana = AnalyticalTet4::new(&pad);
        let s = 9.0;
        let kk_correct = ana.mat_04_nsb(s);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_04_nsb(&mut kk, &mut pad_b, &mut args, |_, _, _, _| Ok(s)).unwrap();
            // println!("{:.2}", kk);
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }
}
