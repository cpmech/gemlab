use super::CommonArgs;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::Vector;

/// Implements the shape(Na) times scalar(S) integration case 01 (mixed shapes)
///
/// Callback function: `s ← f(p, N, Na)`
///
/// **Notes:**
///
/// * `pad` (in `args`) is the **main** pad: it carries the geometry, the integration measure `|J|`,
///   and its own shape functions `N`
/// * `pad_alt` is the **alternative** pad: only its shape functions `Na` are used, never its
///   coordinates, and it may be of any order — lower, equal to, or higher than that of `pad`
/// * `m` ranges over the nodes of `pad_alt` (for `Naᵐ`), as in the result layout below
/// * For example, `0 ≤ m < 4` with a `pad_alt→Qua4` and `pad→Qua8`, or `0 ≤ m < 8` the other way round
/// * The number of integration points must consider the nodes of both pads and the expected order of
///   the whole integrand
///
/// Interpolation functions times scalar field:
///
/// ```text
///      ⌠    → →     →
/// aᵐ = │ Naᵐ(x(ξ)) s(x) α dΩ
///      ⌡
///      Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///      nip-1     →      →       →
/// aᵐ ≈   Σ    Naᵐ(ιᵖ) s(ιᵖ) |J|(ιᵖ) wᵖ α
///       p=0
/// ```
///
/// # Results
///
/// ```text
///     ┌     ┐
///     |  a⁰ |  ⟸  ii0 = 0
///     |  a¹ |
/// a = |  a² |
///     | ··· |
///     |  aᵐ |  ⟸  ii
///     └     ┘
/// ```
///
/// # Arguments
///
/// * `a` -- A vector containing all `aᵐ` values, one after another, and sequentially placed as shown
///   above. `m` is the index of the node of `pad_alt`. The length must be
///   `a.len() ≥ ii0 + pad_alt.nnode`
/// * `args` --- Common arguments (`pad` is the **main** scratchpad (modified), which carries the
///   geometry)
/// * `pad_alt` -- Scratchpad of the alternative shape (modified) to compute `Na`; its coordinates are
///   not used
/// * `fn_s` -- Function `f(p,N,Na)` that computes `s(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`, the shape
///   functions `N(ιᵖ)` of the main pad, and the shape functions `Na(ιᵖ)` of the alternative pad
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
/// use russell_lab::{vec_approx_eq, Vector};
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
///     let mut a = Vector::filled(pad.nnode(), 0.0);
///     let mut args = CommonArgs::new(&mut pad, &gauss);
///     integ::vec_01_ns_mix(&mut a, &mut args, &mut pad_alt, |_, _, _| Ok(5.0))?;
///     // solution (cₛ = 5, A = 6):
///     // aᵐ = cₛ A / 3 = 10
///     vec_approx_eq(&a, &[10.0, 10.0, 10.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_01_ns_mix<F>(
    a: &mut Vector,
    args: &mut CommonArgs,
    pad_alt: &mut Scratchpad,
    mut fn_s: F,
) -> Result<(), StrError>
where
    F: FnMut(usize, &Vector, &Vector) -> Result<f64, StrError>,
{
    // check
    let nnode_alt = pad_alt.interp.dim();
    let (_, _, nnode) = args.pad.dims();
    let ii0 = args.ii0;
    if pad_alt.kind.class() != args.pad.kind.class() || pad_alt.kind.ndim() != args.pad.kind.ndim() {
        return Err("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)");
    }
    if a.dim() < ii0 + nnode_alt {
        return Err("a.len() must be ≥ ii0 + pad_alt.nnode");
    }

    // clear output vector
    if args.clear {
        a.fill(0.0);
    }

    // loop over integration points
    for p in 0..args.gauss.npoint() {
        // ksi coordinates and weight
        let iota = args.gauss.coords(p);
        let weight = args.gauss.weight(p);

        // interpolation functions and Jacobian of the main pad
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_jacobian(iota)?;

        // interpolation functions of the alternative pad
        (pad_alt.fn_interp)(&mut pad_alt.interp, iota); // Na

        // calculate s
        let nn = &args.pad.interp;
        let na = &pad_alt.interp;
        let s = fn_s(p, nn, na)?;

        // calculate coefficient
        let coef = if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt.get(0, m);
            }
            s * det_jac * weight * args.alpha * r
        } else {
            s * det_jac * weight * args.alpha
        };

        // loop over the nodes of the alternative pad and perform sum
        for m in 0..nnode_alt {
            a[ii0 + m] += na[m] * coef;
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
    use russell_lab::{approx_eq, vec_approx_eq, Vector};
    use russell_tensor::Tensor1;

    /// Spatially varying scalar field: s(x) = 2 + 3 x₀ - x₁
    fn scalar_field(x: &Tensor1) -> f64 {
        2.0 + 3.0 * x.get(0) - x.get(1)
    }

    /// Compares vec_01_ns_mix with vec_01_ns using two identically populated pads
    fn check_equal_order(pad: &Scratchpad, axisymmetric: bool, ngauss: usize) {
        let mut pad_mix = pad.clone();
        let mut pad_base = pad.clone();
        let mut pad_alt = pad.clone();
        let gauss = Gauss::new_or_sized(pad.kind, Some(ngauss)).unwrap();
        let mut a_mix = Vector::filled(pad.nnode(), aux::NOISE);
        let mut a_base = Vector::filled(pad.nnode(), aux::NOISE);
        let mut args_mix = CommonArgs::new(&mut pad_mix, &gauss);
        let mut args_base = CommonArgs::new(&mut pad_base, &gauss);
        args_mix.axisymmetric = axisymmetric;
        args_base.axisymmetric = axisymmetric;

        // the source term is evaluated at the integration points of the main (geometric) pad
        let x_ips = recovery::get_points_coords(args_base.pad, &gauss).unwrap();
        integ::vec_01_ns(&mut a_base, &mut args_base, |p, _| Ok(scalar_field(&x_ips[p]))).unwrap();
        integ::vec_01_ns_mix(&mut a_mix, &mut args_mix, &mut pad_alt, |p, _, _| {
            Ok(scalar_field(&x_ips[p]))
        })
        .unwrap();
        vec_approx_eq(&a_mix, &a_base, 1e-14);
    }

    #[test]
    fn capture_some_errors() {
        // class/geo_ndim mismatch: main Qua8 (2D) with alt Tri3 (2D)
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, 2.0, 3.0);
        let mut pad_alt = aux::gen_pad_tri3();
        let mut a = Vector::filled(8, aux::NOISE);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::vec_01_ns_mix(&mut a, &mut args, &mut pad_alt, |_, _, _| Ok(1.0)).err(),
            Some("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)")
        );

        // class/geo_ndim mismatch: main Tri3 (2D) with alt Tet4 (3D)
        let mut pad = aux::gen_pad_tri3();
        let mut pad_alt = aux::gen_pad_tet4();
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::vec_01_ns_mix(&mut a, &mut args, &mut pad_alt, |_, _, _| Ok(1.0)).err(),
            Some("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)")
        );

        // a.dim() < ii0 + nnode_alt
        let mut pad = aux::gen_pad_tri3();
        let mut pad_alt = aux::gen_pad_tri3();
        let mut a = Vector::new(3);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::vec_01_ns_mix(&mut a, &mut args, &mut pad_alt, |_, _, _| Ok(1.0)).err(),
            Some("a.len() must be ≥ ii0 + pad_alt.nnode")
        );

        // the callback error is propagated
        args.ii0 = 0;
        assert_eq!(
            integ::vec_01_ns_mix(&mut a, &mut args, &mut pad_alt, |_, _, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn equal_order_matches_vec_01_ns() {
        // 2D: Tri3 + Tri3
        let pad = aux::gen_pad_tri3();
        for ngauss in [1, 3, 7] {
            check_equal_order(&pad, false, ngauss);
        }

        // axisymmetric: the radius is r > 0 on the chosen triangle
        for ngauss in [1, 3, 7] {
            check_equal_order(&pad, true, ngauss);
        }

        // 3D: Tet4 + Tet4
        let pad = aux::gen_pad_tet4();
        for ngauss in [1, 4] {
            check_equal_order(&pad, false, ngauss);
        }
    }

    #[test]
    fn curved_element_uses_the_main_jacobian() {
        // Curved 8-node element: the midside nodes are pushed off the straight edges, so the
        // element's map and the bilinear (4-node) map are different domains. The reference below is
        // assembled here from the main pad's |J|, independently of vec_01_ns_mix.
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

        // coordinates of the integration points from the main shape (used by the source term)
        let mut pad_x = pad.clone();
        let x_ips = recovery::get_points_coords(&mut pad_x, &gauss).unwrap();
        let source = |p: usize| 1.0 + 2.0 * x_ips[p].get(0) - x_ips[p].get(1);

        // a_ref: assembled from the main pad's calc_jacobian and the alternative shape functions
        let mut a_ref = Vector::new(4);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_jacobian(iota).unwrap();
            (pad_alt_calc.fn_interp)(&mut pad_alt_calc.interp, iota);
            let coef = source(p) * det_jac * gauss.weight(p);
            for m in 0..4 {
                a_ref[m] += pad_alt_calc.interp[m] * coef;
            }
        }

        // the mix function must match the reference
        let mut a = Vector::filled(4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_01_ns_mix(&mut a, &mut args, &mut pad_alt, |p, _, _| Ok(source(p))).unwrap();
        vec_approx_eq(&a, &a_ref, 1e-13);

        // the coordinates of pad_alt are not used: a pad without coordinates gives the same vector
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let mut a_nocoords = Vector::filled(4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_01_ns_mix(&mut a_nocoords, &mut args, &mut pad_alt_nocoords, |p, _, _| {
            Ok(source(p))
        })
        .unwrap();
        vec_approx_eq(&a_nocoords, &a, 1e-15);

        // the wrong variant uses det_jac from the alternative pad's own bilinear map
        let mut a_wrong = Vector::new(4);
        let mut pad_alt_wrong = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad_alt_wrong.set_xx(m, 0, xx[m][0]);
            pad_alt_wrong.set_xx(m, 1, xx[m][1]);
        }
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_alt_wrong.calc_jacobian(iota).unwrap();
            (pad_alt_wrong.fn_interp)(&mut pad_alt_wrong.interp, iota);
            let coef = source(p) * det_jac * gauss.weight(p);
            for m in 0..4 {
                a_wrong[m] += pad_alt_wrong.interp[m] * coef;
            }
        }
        let mut max_diff = 0.0;
        for m in 0..4 {
            max_diff = f64::max(max_diff, f64::abs(a[m] - a_wrong[m]));
        }
        println!("max |a - a_wrong| = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn reversed_order_qua4_qua8() {
        // pad = Qua4 carries the geometry and pad_alt = Qua8 is differentiated in the bilinear map
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

        // coordinates of the integration points from the main shape (used by the source term)
        let mut pad_x = pad.clone();
        let x_ips = recovery::get_points_coords(&mut pad_x, &gauss).unwrap();
        let source = |p: usize| 1.0 + 2.0 * x_ips[p].get(0) - x_ips[p].get(1);

        // a_ref: assembled from the main pad's calc_jacobian and the alternative shape functions
        let mut a_ref = Vector::new(8);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_jacobian(iota).unwrap();
            (pad_alt_calc.fn_interp)(&mut pad_alt_calc.interp, iota);
            let coef = source(p) * det_jac * gauss.weight(p);
            for m in 0..8 {
                a_ref[m] += pad_alt_calc.interp[m] * coef;
            }
        }

        // the mix function must match the reference
        let mut a = Vector::filled(8, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_01_ns_mix(&mut a, &mut args, &mut pad_alt, |p, _, _| Ok(source(p))).unwrap();
        vec_approx_eq(&a, &a_ref, 1e-13);

        // the coordinates of pad_alt are not used: a pad without coordinates gives the same vector
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        let mut a_nocoords = Vector::filled(8, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_01_ns_mix(&mut a_nocoords, &mut args, &mut pad_alt_nocoords, |p, _, _| {
            Ok(source(p))
        })
        .unwrap();
        vec_approx_eq(&a_nocoords, &a, 1e-15);

        // the wrong variant uses det_jac from the alternative pad's own quadratic map
        let mut a_wrong = Vector::new(8);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_alt.calc_jacobian(iota).unwrap();
            (pad_alt.fn_interp)(&mut pad_alt.interp, iota);
            let coef = source(p) * det_jac * gauss.weight(p);
            for m in 0..8 {
                a_wrong[m] += pad_alt.interp[m] * coef;
            }
        }
        let mut max_diff = 0.0;
        for m in 0..8 {
            max_diff = f64::max(max_diff, f64::abs(a[m] - a_wrong[m]));
        }
        println!("max |a - a_wrong| = {:.3e}", max_diff);
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
        let mut a_single = Vector::filled(3, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_01_ns_mix(&mut a_single, &mut args, &mut pad_alt, |_, _, _| Ok(2.5)).unwrap();

        // two runs accumulating at ii0 = 2; the vector is prefilled with 0.5 to prove that
        // clear = false leaves the leading entries untouched
        let mut a = Vector::filled(8, 0.5);
        for m in 0..3 {
            a[2 + m] = 0.0; // start the accumulated slots from zero
        }
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.clear = false;
        args.ii0 = 2;
        for _ in 0..2 {
            integ::vec_01_ns_mix(&mut a, &mut args, &mut pad_alt, |_, _, _| Ok(2.5)).unwrap();
        }
        for i in 0..2 {
            approx_eq(a[i], 0.5, 1e-15);
        }
        for m in 0..3 {
            approx_eq(a[2 + m], 2.0 * a_single[m], 1e-14);
        }
    }
}
