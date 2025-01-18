use super::CommonArgs;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the gradient(B) times scalar(S) times shape(Nb) integration case 07 (e.g., coupling matrix)
///
/// Callback function: `s ← f(p, N, B, Nb)`
///
/// **Notes:**
///
/// * `m` ranges over the number of nodes of the *driver* shape specified by `pad` (for `Bᵐ`)
/// * `n` ranges over the number of nodes of the *lower-order* shape specified by `pad_b` (for `Nbⁿ`)
/// * For example, `1 ≤ m ≤ 8` for a `pad→Qua8` and `1 ≤ n ≤ 4` for `pad_b→Qua4`
/// * The determinant of the Jacobian is calculated for `pad` (`pad` is the driver of the calculations)
/// * The number of integration points must consider the nodes of `pad` and the expected order of the whole integrand
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠ →
/// Kᵐⁿ = │ Bᵐ s Nbⁿ α dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///        nip-1     →     →       →       →
/// Kᵐⁿᵢ ≈   Σ   Bᵐᵢ(ιᵖ) s(ιᵖ) Nbⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
///         p=0
/// ```
///
/// # Results
///
/// ```text
///     ┌                         ┐
///     | K⁰⁰₀ K⁰¹₀ K⁰²₀ ··· K⁰ⁿ₀ |  ⟸  ii0
///     | K⁰⁰₁ K⁰¹₁ K⁰²₁ ··· K⁰ⁿ₁ |
///     | K¹⁰₀ K¹¹₀ K¹²₀ ··· K¹ⁿ₀ |
/// K = | K¹⁰₁ K¹¹₁ K¹²₁ ··· K¹ⁿ₁ |
///     | K²⁰₀ K²¹₀ K²²₀ ··· K²ⁿ₀ |
///     | K²⁰₁ K²¹₁ K²²₁ ··· K²ⁿ₁ |
///     |  ···  ···  ··· ···  ··· |     (pad)
///     | Kᵐ⁰ᵢ Kᵐ¹ᵢ Kᵐ²ᵢ ··· Kᵐⁿᵢ |  ⟸  ii := i + m ⋅ space_ndim
///     └                         ┘
///        ⇑                  ⇑
///       jj0         (pad_b) jj
///
/// m = ii / space_ndim
/// i = ii % space_ndim
/// ```
///
/// # Arguments
///
/// * `kk` -- A matrix containing all `Kᵐⁿᵢ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `i` corresponds to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode ⋅ space_ndim` and `ncol(K) ≥ jj0 + pad_b.nnode`.
/// * `args` --- Common arguments (`pad` is the Driver scratchpad (modified) to compute B)
/// * `pad_b` -- Lower-order scratchpad (modified) to compute Nb
/// * `fn_s` -- Function `f(p,N,B,Nb)→s` that computes `s(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions N(ιᵖ), gradients B(ιᵖ), and shape functions Nb(ιᵖ).
///
/// # Warning
///
/// The two [crate::shapes::Scratchpad]s mut be compatible, otherwise **calculation errors may occur**.
/// Therefore, `pad_b` must be either the lower-version of `pad` or have the same shape as `pad`.
pub fn mat_07_bsn<F>(
    kk: &mut Matrix,
    args: &mut CommonArgs,
    pad_b: &mut Scratchpad,
    mut fn_s: F,
) -> Result<(), StrError>
where
    F: FnMut(usize, &Vector, &Matrix, &Vector) -> Result<f64, StrError>,
{
    // check
    let nnode_b = pad_b.interp.dim();
    let (space_ndim, nnode) = args.pad.xxt.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if nrow_kk < ii0 + nnode * space_ndim {
        return Err("nrow(K) must be ≥ ii0 + pad.nnode ⋅ space_ndim");
    }
    if ncol_kk < jj0 + nnode_b {
        return Err("ncol(K) must be ≥ jj0 + pad_b.nnode");
    }

    // clear output matrix
    if args.clear {
        kk.fill(0.0);
    }

    // loop over integration points
    for p in 0..args.gauss.npoint() {
        // ksi coordinates and weight
        let iota = args.gauss.coords(p);
        let weight = args.gauss.weight(p);

        // calculate interpolation functions and Jacobian
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B
        (pad_b.fn_interp)(&mut pad_b.interp, iota); // Nb

        // calculate s
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        let nnb = &pad_b.interp;
        let s = fn_s(p, nn, bb, nnb)?;

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
            for m in 0..nnode {
                for n in 0..nnode_b {
                    kk.add(ii0 + 0 + m * 2, jj0 + n, bb.get(m, 0) * c * nnb[n]);
                    kk.add(ii0 + 1 + m * 2, jj0 + n, bb.get(m, 1) * c * nnb[n]);
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode_b {
                    kk.add(ii0 + 0 + m * 3, jj0 + n, bb.get(m, 0) * c * nnb[n]);
                    kk.add(ii0 + 1 + m * 3, jj0 + n, bb.get(m, 1) * c * nnb[n]);
                    kk.add(ii0 + 2 + m * 3, jj0 + n, bb.get(m, 2) * c * nnb[n]);
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
    use crate::integ::{self, AnalyticalQua8, AnalyticalTet4, CommonArgs, Gauss};
    use russell_lab::{mat_approx_eq, Matrix, Vector};

    #[test]
    fn capture_some_errors() {
        let (a, b) = (2.0, 3.0);
        let mut pad_b = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut kk = Matrix::new(8 * 2, 4);
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let nnb = Vector::new(0);
        let f = |_p: usize, _nn: &Vector, _bb: &Matrix, _nnb: &Vector| Ok(0.0);
        assert_eq!(f(0, &nn, &bb, &nnb).unwrap(), 0.0);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_07_bsn(&mut kk, &mut args, &mut pad_b, f).err(),
            Some("nrow(K) must be ≥ ii0 + pad.nnode ⋅ space_ndim")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_07_bsn(&mut kk, &mut args, &mut pad_b, f).err(),
            Some("ncol(K) must be ≥ jj0 + pad_b.nnode")
        );
    }

    #[test]
    fn qua4_qua8_works() {
        let (a, b) = (2.0, 3.0);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_b = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut kk = Matrix::new(8 * 2, 4);
        let ana = AnalyticalQua8::new(a, b);
        let s = 9.0;
        let kk_correct = ana.mat_07_bsn(s);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [4, 9].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_07_bsn(&mut kk, &mut args, &mut pad_b, |_, _, _, _| Ok(s)).unwrap();
            // println!("{:.2}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn tet4_tet8_works() {
        let mut pad_b = aux::gen_pad_tet4();
        let mut pad = pad_b.clone();
        let mut kk = Matrix::new(4 * 3, 4);
        let ana = AnalyticalTet4::new(&pad);
        let s = 9.0;
        let kk_correct = ana.mat_07_bsn(s);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_07_bsn(&mut kk, &mut args, &mut pad_b, |_, _, _, _| Ok(s)).unwrap();
            // println!("{}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }
}
