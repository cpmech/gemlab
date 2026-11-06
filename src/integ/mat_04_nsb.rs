use super::CommonArgs;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the shape(Na) times scalar(S) times gradient(B) integration case 04 (e.g., coupling matrix)
///
/// Callback function: `s ← f(p, Na, N, B)`
///
/// **Notes:**
///
/// * `pad` (in `args`) is the **main** pad: it carries the geometry, the integration measure `|J|`,
///   and its own shape functions `N` and gradients `B`
/// * `pad_alt` is the **alternative** pad: only its shape functions are used, never its coordinates,
///   and it may be of any order — lower, equal to, or higher than that of `pad`
/// * `m` ranges over the nodes of `pad_alt` (for `Naᵐ`) and `n` over the nodes of `pad` (for `Bⁿ`),
///   as in the result layout below
/// * For example, `1 ≤ m ≤ 4` with a `pad_alt→Qua4` and `pad→Qua8`, or `1 ≤ m ≤ 8` the other way round
/// * The number of integration points must consider the nodes of both pads and the expected order of
///   the whole integrand
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠       →
/// Kᵐⁿ = │ Nᵃᵐ s Bⁿ α dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///        nip-1     →     →       →       →
/// Kᵐⁿⱼ ≈   Σ   Nᵃᵐ(ιᵖ) s(ιᵖ) Bⁿⱼ(ιᵖ) |J|(ιᵖ) wᵖ α
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
///     |  ··   ··   ··   ··   ··   ··  ···  ··  |    (pad_alt)
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
///   The dimensions must be `nrow(K) ≥ ii0 + pad_alt.nnode` and `ncol(K) ≥ jj0 + pad.nnode ⋅ space_ndim`
/// * `args` --- Common arguments (`pad` is the **main** scratchpad (modified), which carries the
///   geometry)
/// * `pad_alt` -- Scratchpad of the alternative shape (modified) to compute `Na`; its
///   coordinates are not used
/// * `fn_s` -- Function `f(p,Na,N,B)→s` that computes `s(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`,
///   shape functions Na(ιᵖ), shape functions N(ιᵖ), and gradients B(ιᵖ).
///
/// # Warning
///
/// Both [crate::shapes::Scratchpad]s must be functions on the same reference element, otherwise
/// **calculation errors will occur**: one `ιᵖ` must mean the same point in both. Their classes and
/// geometric dimensions are checked here; the orders are free, and `pad_alt` may be of lower, equal,
/// or higher order than `pad`.
pub fn mat_04_nsb<F>(
    kk: &mut Matrix,
    args: &mut CommonArgs,
    pad_alt: &mut Scratchpad,
    mut fn_s: F,
) -> Result<(), StrError>
where
    F: FnMut(usize, &Vector, &Vector, &Matrix) -> Result<f64, StrError>,
{
    // check
    let nnode_alt = pad_alt.interp.dim();
    let (space_ndim, _, nnode) = args.pad.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if pad_alt.kind.class() != args.pad.kind.class() || pad_alt.kind.ndim() != args.pad.kind.ndim() {
        return Err("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)");
    }
    if nrow_kk < ii0 + nnode_alt {
        return Err("nrow(K) must be ≥ ii0 + pad_alt.nnode");
    }
    if ncol_kk < jj0 + nnode * space_ndim {
        return Err("ncol(K) must be ≥ jj0 + pad.nnode ⋅ space_ndim");
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
        (pad_alt.fn_interp)(&mut pad_alt.interp, iota); // Nb
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B

        // calculate s
        let nn_alt = &pad_alt.interp;
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        let s = fn_s(p, nn_alt, nn, bb)?;

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
            for m in 0..nnode_alt {
                for n in 0..nnode {
                    kk.add(ii0 + m, jj0 + 0 + n * 2, nn_alt[m] * c * bb.get(n, 0));
                    kk.add(ii0 + m, jj0 + 1 + n * 2, nn_alt[m] * c * bb.get(n, 1));
                }
            }
        } else {
            for m in 0..nnode_alt {
                for n in 0..nnode {
                    kk.add(ii0 + m, jj0 + 0 + n * 3, nn_alt[m] * c * bb.get(n, 0));
                    kk.add(ii0 + m, jj0 + 1 + n * 3, nn_alt[m] * c * bb.get(n, 1));
                    kk.add(ii0 + m, jj0 + 2 + n * 3, nn_alt[m] * c * bb.get(n, 2));
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
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::{mat_approx_eq, Matrix, Vector};

    #[test]
    fn capture_some_errors() {
        let (a, b) = (2.0, 3.0);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut kk = Matrix::new(4, 8 * 2);
        let nn_alt = Vector::new(0);
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let f = |_p: usize, _nn_alt: &Vector, _nn: &Vector, _bb: &Matrix| Ok(0.0);
        assert_eq!(f(0, &nn_alt, &nn, &bb).unwrap(), 0.0);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_04_nsb(&mut kk, &mut args, &mut pad_alt, f).err(),
            Some("nrow(K) must be ≥ ii0 + pad_alt.nnode")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_04_nsb(&mut kk, &mut args, &mut pad_alt, f).err(),
            Some("ncol(K) must be ≥ jj0 + pad.nnode ⋅ space_ndim")
        );
    }

    #[test]
    fn qua4_qua8_works() {
        let (a, b) = (2.0, 3.0);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut kk = Matrix::new(4, 8 * 2);
        let ana = AnalyticalQua8::new(a, b);
        let s = 9.0;
        let kk_correct = ana.mat_04_nsb(s);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [4, 9].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_04_nsb(&mut kk, &mut args, &mut pad_alt, |_, _, _, _| Ok(s)).unwrap();
            // println!("{:.2}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn tet4_tet4_works() {
        let mut pad_alt = aux::gen_pad_tet4();
        let mut pad = pad_alt.clone();
        let mut kk = Matrix::new(4, 4 * 3);
        let ana = AnalyticalTet4::new(&pad);
        let s = 9.0;
        let kk_correct = ana.mat_04_nsb(s);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14];
        let selection: Vec<_> = [4].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_04_nsb(&mut kk, &mut args, &mut pad_alt, |_, _, _, _| Ok(s)).unwrap();
            // println!("{:.2}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn works_with_any_order_combination() {
        // the main pad carries the geometry; the alternative may be of lower, equal, or higher order
        #[rustfmt::skip]
        let xx = [
            [0.00, 0.00], [2.00, 0.10], [1.70, 1.30], [-0.20, 1.00],
            [1.00, 0.35], [2.05, 0.72], [0.75, 1.45], [-0.10, 0.55],
        ];
        let cases = [
            (GeoKind::Qua4, 4, GeoKind::Qua8, 8),
            (GeoKind::Qua8, 8, GeoKind::Qua4, 4),
            (GeoKind::Qua4, 4, GeoKind::Qua4, 4),
        ];

        let s = 2.5;
        for (main_kind, n_main, alt_kind, n_alt) in cases {
            let mut pad = Scratchpad::new(2, main_kind).unwrap();
            for m in 0..n_main {
                pad.set_xx(m, 0, xx[m][0]);
                pad.set_xx(m, 1, xx[m][1]);
            }
            let mut pad_alt = Scratchpad::new(2, alt_kind).unwrap();
            for m in 0..n_alt.min(n_main) {
                pad_alt.set_xx(m, 0, xx[m][0]);
                pad_alt.set_xx(m, 1, xx[m][1]);
            }
            let gauss = Gauss::new_or_sized(main_kind, Some(9)).unwrap();

            let mut kk_ref = Matrix::new(n_alt, n_main * 2);
            let mut pad_calc = pad.clone();
            let mut pad_alt_calc = Scratchpad::new(2, alt_kind).unwrap();
            for p in 0..gauss.npoint() {
                let iota = gauss.coords(p);
                let det_jac = pad_calc.calc_gradient(iota).unwrap();
                (pad_alt_calc.fn_interp)(&mut pad_alt_calc.interp, iota);
                let c = s * det_jac * gauss.weight(p);
                for m in 0..n_alt {
                    for n in 0..n_main {
                        for j in 0..2 {
                            kk_ref.add(m, n * 2 + j, pad_alt_calc.interp[m] * c * pad_calc.gradient.get(n, j));
                        }
                    }
                }
            }

            let mut kk = Matrix::new(n_alt, n_main * 2);
            let mut args = CommonArgs::new(&mut pad, &gauss);
            integ::mat_04_nsb(&mut kk, &mut args, &mut pad_alt, |_, _, _, _| Ok(s)).unwrap();
            mat_approx_eq(&kk, &kk_ref, 1e-13);
        }
    }
}
