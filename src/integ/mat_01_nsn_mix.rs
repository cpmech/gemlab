use super::CommonArgs;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the shape(Na) times scalar(S) times shape(Na) integration case 01 (mixed shapes)
///
/// Callback function: `s ← f(p, N, B, Na)`
///
/// **Notes:**
///
/// * `pad` (in `args`) is the **main** pad: it carries the geometry, the integration measure `|J|`,
///   and its own shape functions `N` and gradients `B`
/// * `pad_alt` is the **alternative** pad: only its shape functions `Na` are used, never its
///   coordinates, and it may be of any order — lower, equal to, or higher than that of `pad`
/// * `m` and `n` range over the nodes of `pad_alt` (for `Naᵐ` and `Naⁿ`), as in the result layout
///   below
/// * For example, `0 ≤ m,n < 4` with a `pad_alt→Qua4` and `pad→Qua8`, or `0 ≤ m,n < 8` the other
///   way around
/// * The number of integration points must consider the nodes of both pads and the expected order
///   of the whole integrand
///
/// Coefficients:
///
/// ```text
///       ⌠
/// Kᵐⁿ = │ Naᵐ s Naⁿ α dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///        nip-1     →     →      →       →
/// Kᵐⁿ ≈   Σ   Naᵐ(ιᵖ) s(ιᵖ) Naⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
///         p=0
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
/// * `kk` -- A matrix containing all `Kᵐⁿ` values, one after another, and sequentially placed as
///   shown above. `m` and `n` are the indices of the nodes of `pad_alt`. The dimensions must be
///   `nrow(K) ≥ ii0 + pad_alt.nnode` and `ncol(K) ≥ jj0 + pad_alt.nnode`
/// * `args` --- Common arguments (`pad` is the **main** scratchpad (modified), which carries the
///   geometry)
/// * `pad_alt` -- Scratchpad of the alternative shape (modified) to compute `Na`; its coordinates
///   are not used
/// * `fn_s` -- Function `f(p,N,B,Na)` that computes `s(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`, the shape
///   functions `N(ιᵖ)` of the main pad, the gradients `B(ιᵖ)` of the main pad, and the shape
///   functions `Na(ιᵖ)` of the alternative pad
///
/// # Warning
///
/// Both [crate::shapes::Scratchpad]s must be functions on the same reference element, otherwise
/// **calculation errors will occur**: one `ιᵖ` must mean the same point in both. Their classes and
/// geometric dimensions are checked here; the orders are free, and `pad_alt` may be of lower, equal,
/// or higher order than `pad`.
///
/// # Examples
///
/// ```
/// use gemlab::integ::{self, CommonArgs, Gauss};
/// use gemlab::shapes::{GeoKind, Scratchpad};
/// use gemlab::StrError;
/// use russell_lab::{mat_approx_eq, Matrix};
///
/// fn main() -> Result<(), StrError> {
///     let space_ndim = 2;
///     let mut pad = Scratchpad::new(space_ndim, GeoKind::Tri3)?;
///     pad.set_xx(0, 0, 2.0);
///     pad.set_xx(0, 1, 3.0);
///     pad.set_xx(1, 0, 6.0);
///     pad.set_xx(1, 1, 3.0);
///     pad.set_xx(2, 0, 2.0);
///     pad.set_xx(2, 1, 6.0);
///     let mut pad_alt = Scratchpad::new(space_ndim, GeoKind::Tri3)?;
///     pad_alt.set_xx(0, 0, 2.0);
///     pad_alt.set_xx(0, 1, 3.0);
///     pad_alt.set_xx(1, 0, 6.0);
///     pad_alt.set_xx(1, 1, 3.0);
///     pad_alt.set_xx(2, 0, 2.0);
///     pad_alt.set_xx(2, 1, 6.0);
///     let gauss = Gauss::new(pad.kind);
///     let mut kk = Matrix::new(3, 3);
///     let mut args = CommonArgs::new(&mut pad, &gauss);
///     integ::mat_01_nsn_mix(&mut kk, &mut args, &mut pad_alt, |_, _, _, _| Ok(12.0))?;
///     // solution (s = 12, A = 6):
///     // ∫ NᵢNⱼ dΩ = A/6 on the diagonal and A/12 off it
///     #[rustfmt::skip]
///     let kk_correct = Matrix::from(&[
///         [12.0, 6.0, 6.0],
///         [6.0, 12.0, 6.0],
///         [6.0, 6.0, 12.0],
///     ]);
///     mat_approx_eq(&kk, &kk_correct, 1e-14);
///     Ok(())
/// }
/// ```
pub fn mat_01_nsn_mix<F>(
    kk: &mut Matrix,
    args: &mut CommonArgs,
    pad_alt: &mut Scratchpad,
    mut fn_s: F,
) -> Result<(), StrError>
where
    F: FnMut(usize, &Vector, &Matrix, &Vector) -> Result<f64, StrError>,
{
    // check
    let nnode_alt = pad_alt.interp.dim();
    let (_, _, nnode) = args.pad.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if pad_alt.kind.class() != args.pad.kind.class() || pad_alt.kind.ndim() != args.pad.kind.ndim() {
        return Err("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)");
    }
    if nrow_kk < ii0 + nnode_alt {
        return Err("nrow(K) must be ≥ ii0 + pad_alt.nnode");
    }
    if ncol_kk < jj0 + nnode_alt {
        return Err("ncol(K) must be ≥ jj0 + pad_alt.nnode");
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

        // interpolation functions and Jacobian of the main pad
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B

        // interpolation functions of the alternative pad
        (pad_alt.fn_interp)(&mut pad_alt.interp, iota); // Na

        // calculate s
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        let na = &pad_alt.interp;
        let s = fn_s(p, nn, bb, na)?;

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
        for m in 0..nnode_alt {
            for n in 0..nnode_alt {
                kk.add(ii0 + m, jj0 + n, na[m] * c * na[n]);
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, CommonArgs, Gauss};
    use crate::recovery;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::{approx_eq, mat_approx_eq, Matrix};
    use russell_tensor::Tensor1;

    /// Compares mat_01_nsn_mix with mat_01_nsn using two identically populated pads
    ///
    /// The same `scalar_field` closure is called from both integration functions.
    fn check_equal_order(pad: &Scratchpad, axisymmetric: bool, ngauss: usize, scalar_field: impl Fn(&Tensor1) -> f64) {
        let mut pad_mix = pad.clone();
        let mut pad_base = pad.clone();
        let mut pad_alt = pad.clone();
        let gauss = Gauss::new_or_sized(pad.kind, Some(ngauss)).unwrap();
        let nnode = pad.nnode();
        let mut kk_mix = Matrix::filled(nnode, nnode, aux::NOISE);
        let mut kk_base = Matrix::filled(nnode, nnode, aux::NOISE);
        let mut args_mix = CommonArgs::new(&mut pad_mix, &gauss);
        let mut args_base = CommonArgs::new(&mut pad_base, &gauss);
        args_mix.axisymmetric = axisymmetric;
        args_base.axisymmetric = axisymmetric;

        // the scalar field is evaluated at the integration points of the main (geometric) pad
        let x_ips = recovery::get_points_coords(args_base.pad, &gauss).unwrap();
        integ::mat_01_nsn(&mut kk_base, &mut args_base, |p, _, _| Ok(scalar_field(&x_ips[p]))).unwrap();
        integ::mat_01_nsn_mix(&mut kk_mix, &mut args_mix, &mut pad_alt, |p, _, _, _| {
            Ok(scalar_field(&x_ips[p]))
        })
        .unwrap();
        mat_approx_eq(&kk_mix, &kk_base, 1e-14);
    }

    #[test]
    fn capture_some_errors() {
        // class/geo_ndim mismatch: main Qua8 (2D) with alt Tri3 (2D)
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, 2.0, 3.0);
        let mut pad_alt = aux::gen_pad_tri3();
        let mut kk = Matrix::new(8, 8);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::mat_01_nsn_mix(&mut kk, &mut args, &mut pad_alt, |_, _, _, _| Ok(1.0)).err(),
            Some("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)")
        );

        // class/geo_ndim mismatch: main Tri3 (2D) with alt Tet4 (3D)
        let mut pad = aux::gen_pad_tri3();
        let mut pad_alt = aux::gen_pad_tet4();
        let mut kk = Matrix::new(3, 3);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::mat_01_nsn_mix(&mut kk, &mut args, &mut pad_alt, |_, _, _, _| Ok(1.0)).err(),
            Some("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)")
        );

        // nrow(K) < ii0 + pad_alt.nnode
        let mut pad = aux::gen_pad_tri3();
        let mut pad_alt = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_01_nsn_mix(&mut kk, &mut args, &mut pad_alt, |_, _, _, _| Ok(1.0)).err(),
            Some("nrow(K) must be ≥ ii0 + pad_alt.nnode")
        );

        // ncol(K) < jj0 + pad_alt.nnode
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_01_nsn_mix(&mut kk, &mut args, &mut pad_alt, |_, _, _, _| Ok(1.0)).err(),
            Some("ncol(K) must be ≥ jj0 + pad_alt.nnode")
        );

        // the callback error is propagated
        args.jj0 = 0;
        assert_eq!(
            integ::mat_01_nsn_mix(&mut kk, &mut args, &mut pad_alt, |_, _, _, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn equal_order_matches_mat_01_nsn() {
        // spatially varying scalar field: s(x) = 2 + x₀ - 3 x₁ + 2 x₂ (x₂ = 0 in 2D)
        let scalar_field = |x: &Tensor1| 2.0 + x.get(0) - 3.0 * x.get(1) + 2.0 * x.get(2);

        // 2D: Tri3 + Tri3
        let pad = aux::gen_pad_tri3();
        for ngauss in [1, 3, 7] {
            check_equal_order(&pad, false, ngauss, scalar_field);
        }

        // axisymmetric: the radius is r > 0 on the chosen triangle
        for ngauss in [1, 3, 7] {
            check_equal_order(&pad, true, ngauss, scalar_field);
        }

        // 3D: Tet4 + Tet4
        let pad = aux::gen_pad_tet4();
        for ngauss in [1, 4] {
            check_equal_order(&pad, false, ngauss, scalar_field);
        }
    }

    #[test]
    fn curved_element_uses_the_main_jacobian() {
        // Curved 8-node element: the midside nodes are pushed off the straight edges, so the
        // element's map and the bilinear (4-node) map are different domains. The reference below
        // is assembled here from the main pad's |J| and the alternative shape functions,
        // independently of mat_01_nsn_mix.
        #[rustfmt::skip]
        let xx = [
            [0.00, 0.00], [2.00, 0.10], [1.70, 1.30], [-0.20, 1.00],
            [1.00, 0.35], [2.05, 0.72], [0.75, 1.45], [-0.10, 0.55],
        ];
        let mut pad = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for m in 0..8 {
            pad.set_xx(m, 0, xx[m][0]);
            pad.set_xx(m, 1, xx[m][1]);
        }
        let mut pad_alt = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad_alt.set_xx(m, 0, xx[m][0]);
            pad_alt.set_xx(m, 1, xx[m][1]);
        }
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();

        // coordinates of the integration points from the main shape (used by the scalar field)
        let mut pad_x = pad.clone();
        let x_ips = recovery::get_points_coords(&mut pad_x, &gauss).unwrap();

        // spatially varying scalar field, shared by the reference loop and the callback
        let scalar_field = |x: &Tensor1| 1.0 + 2.0 * x.get(0) - 3.0 * x.get(1);

        // kk_ref: assembled from the main pad's |J| and the alternative shape functions
        let mut kk_ref = Matrix::new(4, 4);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_alt_calc.fn_interp)(&mut pad_alt_calc.interp, iota);
            let c = scalar_field(&x_ips[p]) * det_jac * gauss.weight(p);
            for m in 0..4 {
                for n in 0..4 {
                    kk_ref.add(m, n, pad_alt_calc.interp[m] * c * pad_alt_calc.interp[n]);
                }
            }
        }

        // the mix function must match the reference
        let mut kk = Matrix::filled(4, 4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_01_nsn_mix(&mut kk, &mut args, &mut pad_alt, |p, _, _, _| {
            Ok(scalar_field(&x_ips[p]))
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the coordinates of pad_alt are not used: a pad without coordinates gives the same matrix
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let mut kk_nocoords = Matrix::filled(4, 4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_01_nsn_mix(&mut kk_nocoords, &mut args, &mut pad_alt_nocoords, |p, _, _, _| {
            Ok(scalar_field(&x_ips[p]))
        })
        .unwrap();
        mat_approx_eq(&kk_nocoords, &kk, 1e-15);

        // the wrong variant uses |J| from the alternative pad's own bilinear map
        let mut kk_wrong = Matrix::new(4, 4);
        let mut pad_alt_wrong = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad_alt_wrong.set_xx(m, 0, xx[m][0]);
            pad_alt_wrong.set_xx(m, 1, xx[m][1]);
        }
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_alt_wrong.calc_gradient(iota).unwrap();
            let c = scalar_field(&x_ips[p]) * det_jac * gauss.weight(p);
            for m in 0..4 {
                for n in 0..4 {
                    kk_wrong.add(m, n, pad_alt_wrong.interp[m] * c * pad_alt_wrong.interp[n]);
                }
            }
        }
        let mut max_diff = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                max_diff = f64::max(max_diff, f64::abs(kk.get(i, j) - kk_wrong.get(i, j)));
            }
        }
        println!("max |K - K_wrong| = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn reversed_order_qua4_qua8() {
        // pad = Qua4 carries the geometry and pad_alt = Qua8 is evaluated in the bilinear map
        // (a subparametric element). The reference is assembled here independently of the function.
        #[rustfmt::skip]
        let xx = [
            [0.00, 0.00], [2.00, 0.10], [1.70, 1.30], [-0.20, 1.00],
            [1.00, 0.35], [2.05, 0.72], [0.75, 1.45], [-0.10, 0.55],
        ];
        let mut pad = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad.set_xx(m, 0, xx[m][0]);
            pad.set_xx(m, 1, xx[m][1]);
        }
        let mut pad_alt = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for m in 0..8 {
            pad_alt.set_xx(m, 0, xx[m][0]);
            pad_alt.set_xx(m, 1, xx[m][1]);
        }
        let gauss = Gauss::new_or_sized(GeoKind::Qua4, Some(9)).unwrap();

        // coordinates of the integration points from the main shape (used by the scalar field)
        let mut pad_x = pad.clone();
        let x_ips = recovery::get_points_coords(&mut pad_x, &gauss).unwrap();

        // spatially varying scalar field, shared by the reference loop and the callback
        let scalar_field = |x: &Tensor1| 1.0 + 2.0 * x.get(0) - 3.0 * x.get(1);

        // kk_ref: assembled from the main pad's |J| and the alternative shape functions
        let mut kk_ref = Matrix::new(8, 8);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_alt_calc.fn_interp)(&mut pad_alt_calc.interp, iota);
            let c = scalar_field(&x_ips[p]) * det_jac * gauss.weight(p);
            for m in 0..8 {
                for n in 0..8 {
                    kk_ref.add(m, n, pad_alt_calc.interp[m] * c * pad_alt_calc.interp[n]);
                }
            }
        }

        // the mix function must match the reference
        let mut kk = Matrix::filled(8, 8, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_01_nsn_mix(&mut kk, &mut args, &mut pad_alt, |p, _, _, _| {
            Ok(scalar_field(&x_ips[p]))
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the coordinates of pad_alt are not used: a pad without coordinates gives the same matrix
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        let mut kk_nocoords = Matrix::filled(8, 8, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_01_nsn_mix(&mut kk_nocoords, &mut args, &mut pad_alt_nocoords, |p, _, _, _| {
            Ok(scalar_field(&x_ips[p]))
        })
        .unwrap();
        mat_approx_eq(&kk_nocoords, &kk, 1e-15);

        // the wrong variant uses |J| from the alternative pad's own quadratic map
        let mut kk_wrong = Matrix::new(8, 8);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_alt.calc_gradient(iota).unwrap();
            let c = scalar_field(&x_ips[p]) * det_jac * gauss.weight(p);
            for m in 0..8 {
                for n in 0..8 {
                    kk_wrong.add(m, n, pad_alt.interp[m] * c * pad_alt.interp[n]);
                }
            }
        }
        let mut max_diff = 0.0;
        for i in 0..8 {
            for j in 0..8 {
                max_diff = f64::max(max_diff, f64::abs(kk.get(i, j) - kk_wrong.get(i, j)));
            }
        }
        println!("max |K - K_wrong| = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn offset_and_accumulation_work() {
        // single run used as the reference values
        let mut pad = aux::gen_pad_tri3();
        let mut pad_alt = aux::gen_pad_tri3();
        let gauss = Gauss::new(pad.kind);
        let mut kk_single = Matrix::filled(3, 3, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_01_nsn_mix(&mut kk_single, &mut args, &mut pad_alt, |_, _, _, _| Ok(2.5)).unwrap();

        // two runs accumulating at (ii0, jj0) = (1, 2); the matrix is prefilled with 0.5 to prove
        // that clear = false leaves the border entries untouched
        let mut kk = Matrix::filled(5, 6, 0.5);
        for i in 1..4 {
            for j in 2..5 {
                kk.set(i, j, 0.0); // start the accumulated block from zero
            }
        }
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.clear = false;
        args.ii0 = 1;
        args.jj0 = 2;
        for _ in 0..2 {
            integ::mat_01_nsn_mix(&mut kk, &mut args, &mut pad_alt, |_, _, _, _| Ok(2.5)).unwrap();
        }
        for i in 0..5 {
            for j in 0..6 {
                if i == 0 || i > 3 || j < 2 || j > 4 {
                    approx_eq(kk.get(i, j), 0.5, 1e-15);
                } else {
                    approx_eq(kk.get(i, j), 2.0 * kk_single.get(i - 1, j - 2), 1e-14);
                }
            }
        }
    }
}
