use super::CommonArgs;
use crate::StrError;
use russell_lab::Vector;

/// Implements the shape(N) times scalar(S) integration case 01
///
/// Callback function: `s ← f(p, N)`
///
/// Interpolation functions times scalar field:
///
/// ```text
///      ⌠    → →     →
/// aᵐ = │ Nᵐ(x(ξ)) s(x) α dΩ
///      ⌡
///      Ωₑ
/// ```
///
/// or, for lines in multi-dimensions:
///
/// ```text
///      ⌠
/// aᵐ = │ Nᵐ(ℓ(ξ)) s(ℓ) α dℓ
///      ⌡
///      Γₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///      nip-1     →     →       →
/// aᵐ ≈   Σ    Nᵐ(ιᵖ) s(ιᵖ) |J|(ιᵖ) wᵖ α
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
/// * `a` -- A vector containing all `aᵐ` values, one after another, and
///   sequentially placed as shown above. `m` is the index of the node.
///   The length must be `a.len() ≥ ii0 + nnode`
/// * `args` --- Common arguments
/// * `fn_s` -- Function `f(p,N)` that computes `s(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`,
///   and shape functions N(ιᵖ).
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
///     let gauss = Gauss::new(pad.kind);
///     let mut a = Vector::filled(pad.nnode(), 0.0);
///     let mut args = CommonArgs::new(&mut pad, &gauss);
///     integ::vec_01_ns(&mut a, &mut args, |_, _| Ok(5.0))?;
///     // solution (cₛ = 5, A = 6):
///     // aᵐ = cₛ A / 3 = 10
///     vec_approx_eq(&a, &[10.0, 10.0, 10.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_01_ns<F>(a: &mut Vector, args: &mut CommonArgs, mut fn_s: F) -> Result<(), StrError>
where
    F: FnMut(usize, &Vector) -> Result<f64, StrError>,
{
    // check
    let nnode = args.pad.interp.dim();
    let ii0 = args.ii0;
    if a.dim() < ii0 + nnode {
        return Err("a.len() must be ≥ ii0 + nnode");
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

        // calculate interpolation functions and Jacobian
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_jacobian(iota)?;

        // radius @ x(ιᵖ) from the geometry
        let radius = if args.axisymmetric {
            args.pad.calc_radius(iota)?
        } else {
            0.0
        };

        // calculate s
        let nn = &args.pad.interp;
        let s = fn_s(p, nn)?;

        // calculate coefficient
        let coef = if args.axisymmetric {
            s * det_jac * weight * args.alpha * radius
        } else {
            s * det_jac * weight * args.alpha
        };

        // loop over nodes and perform sum
        for m in 0..nnode {
            a[ii0 + m] += nn[m] * coef;
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, AnalyticalTet4, AnalyticalTri3, CommonArgs, Gauss};
    use crate::recovery;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::{approx_eq, vec_approx_eq, Vector};

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut a = Vector::new(2);
        let nn = Vector::new(0);
        let f = |_: usize, _: &Vector| Ok(0.0);
        assert_eq!(f(0, &nn).unwrap(), 0.0);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::vec_01_ns(&mut a, &mut args, f).err(),
            Some("a.len() must be ≥ ii0 + nnode")
        );
    }

    #[test]
    fn vec_01_ns_works_lin2_linear() {
        // lin2 with linear source term:
        //
        // s(x) = x
        //
        // solution:
        //
        //       ┌           ┐
        //     L │ 2 xa + xb │
        // a = — │           │
        //     6 │ xa + 2 xb │
        //       └           ┘
        const L: f64 = 6.0;
        let mut pad = aux::gen_pad_lin2(L);

        // solution
        let cf = L / 6.0;
        let (xa, xb) = (pad.xxt.get(0, 0), pad.xxt.get(0, 1));
        let a_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-15, 1e-14, 1e-14, 1e-14];
        let selection: Vec<_> = [2, 3, 4, 5]
            .iter()
            .map(|n| Gauss::new_sized(class, *n).unwrap())
            .collect();

        // check
        let mut a = Vector::filled(pad.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            let x_ips = recovery::get_points_coords(args.pad, ips).unwrap();
            integ::vec_01_ns(&mut a, &mut args, |p, _| Ok(x_ips[p].get(0))).unwrap();
            vec_approx_eq(&a, a_correct, tol);
        });
    }

    #[test]
    fn vec_01_ns_works_tri3_constant() {
        // tri3 with a constant source term s(x) = cₛ
        let mut pad = aux::gen_pad_tri3();

        // solution
        let ana = AnalyticalTri3::new(&pad);
        const CS: f64 = 3.0;
        let a_correct = ana.vec_01_ns(CS, false);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14, 1e-14, 1e-13, 1e-14];
        let selection: Vec<_> = [1, 3, 4, 12, 16]
            .iter()
            .map(|n| Gauss::new_sized(class, *n).unwrap())
            .collect();

        // check
        let mut a = Vector::filled(pad.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_01_ns(&mut a, &mut args, |_, _| Ok(CS)).unwrap();
            vec_approx_eq(&a, &a_correct, tol);
        });
    }

    #[test]
    fn vec_01_ns_works_tet4_constant() {
        // tet 4 with a constant source term s(x) = cs
        let mut pad = aux::gen_pad_tet4();

        // solution
        const CS: f64 = 120.0;
        let ana = AnalyticalTet4::new(&pad);
        let a_correct = ana.vec_01_ns(CS);

        // integration points
        // Note that the tolerance is high for ngauss = 1
        // because the numerical integration performs poorly with few IPs
        let class = pad.kind.class();
        let tolerances = [1e-13, 1e-13];
        let selection: Vec<_> = [1, 4].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();

        // check
        let mut a = Vector::filled(pad.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_01_ns(&mut a, &mut args, |_, _| Ok(CS)).unwrap();
            vec_approx_eq(&a, &a_correct, tol);
        });
    }

    #[test]
    fn vec_01_ns_works_tet4_linear() {
        // tet 4 with a linear source term s(x) = z = x₂
        let mut pad = aux::gen_pad_tet4();

        // solution
        let ana = AnalyticalTet4::new(&pad);
        let a_correct = ana.vec_01_ns_linear_along_z(&pad);

        // integration points
        // Note that the tolerance is high for ngauss = 1
        // because the numerical integration performs poorly with few IPs
        let class = pad.kind.class();
        let tolerances = [0.56, 1e-15, 1e-14, 1e-15, 1e-15, 1e-15, 1e-15, 1e-15];
        let selection: Vec<_> = [1, 4, 5, 8, 14, 15, 24]
            .iter()
            .map(|n| Gauss::new_sized(class, *n).unwrap())
            .collect();

        // check
        let mut a = Vector::filled(pad.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            let x_ips = recovery::get_points_coords(args.pad, ips).unwrap();
            integ::vec_01_ns(&mut a, &mut args, |p, _| Ok(x_ips[p].get(2))).unwrap();
            vec_approx_eq(&a, &a_correct, tol);
        });
    }

    #[test]
    fn curved_element_uses_the_main_jacobian() {
        // A Qua4 field over a curved Qua8 geometry. The reference uses the Qua8 measure and the Qua4
        // shape functions, assembled by hand independently of vec_01_ns.
        let mut pad = aux::gen_pad_mixed_qua4_qua8(0.0);
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();

        // integration-point coordinates from the curved geometry (used by the source term)
        let mut pad_geo = aux::gen_pad_curved_qua8(0.0);
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();
        let source = |p: usize| 1.0 + 2.0 * x_ips[p].get(0) - x_ips[p].get(1);

        // reference: the Qua8 measure and the Qua4 shape functions
        let mut a_ref = Vector::new(4);
        let mut pad_calc = pad_geo.clone();
        let mut pad_field = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_jacobian(iota).unwrap();
            (pad_field.fn_interp)(&mut pad_field.interp, iota);
            let coef = source(p) * det_jac * gauss.weight(p);
            for m in 0..4 {
                a_ref[m] += pad_field.interp[m] * coef;
            }
        }

        // the single mixed pad must reproduce the reference
        let mut a = Vector::filled(4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_01_ns(&mut a, &mut args, |p, _| Ok(source(p))).unwrap();
        vec_approx_eq(&a, &a_ref, 1e-13);

        // the wrong variant uses the bilinear (Qua4) measure
        let mut pad_wrong = aux::gen_pad_bilinear_qua4(0.0);
        let mut a_wrong = Vector::new(4);
        let mut pad_field_wrong = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_wrong.calc_jacobian(iota).unwrap();
            (pad_field_wrong.fn_interp)(&mut pad_field_wrong.interp, iota);
            let coef = source(p) * det_jac * gauss.weight(p);
            for m in 0..4 {
                a_wrong[m] += pad_field_wrong.interp[m] * coef;
            }
        }
        let mut max_diff = 0.0;
        for m in 0..4 {
            max_diff = f64::max(max_diff, f64::abs(a[m] - a_wrong[m]));
        }
        println!("max |a - a_wrong| (curved) = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn reversed_order_qua4_qua8() {
        // A Qua8 field over a straight-edged Qua4 geometry. The reference uses the Qua4 measure and
        // the Qua8 shape functions, assembled by hand independently of vec_01_ns.
        let mut pad = aux::gen_pad_mixed_qua8_qua4();
        let gauss = Gauss::new_or_sized(GeoKind::Qua4, Some(9)).unwrap();

        // integration-point coordinates from the bilinear geometry (used by the source term)
        let mut pad_geo = aux::gen_pad_bilinear_qua4(0.0);
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();
        let source = |p: usize| 1.0 + 2.0 * x_ips[p].get(0) - x_ips[p].get(1);

        // reference: the Qua4 measure and the Qua8 shape functions
        let mut a_ref = Vector::new(8);
        let mut pad_calc = pad_geo.clone();
        let mut pad_field = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_jacobian(iota).unwrap();
            (pad_field.fn_interp)(&mut pad_field.interp, iota);
            let coef = source(p) * det_jac * gauss.weight(p);
            for m in 0..8 {
                a_ref[m] += pad_field.interp[m] * coef;
            }
        }

        // the single mixed pad must reproduce the reference
        let mut a = Vector::filled(8, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_01_ns(&mut a, &mut args, |p, _| Ok(source(p))).unwrap();
        vec_approx_eq(&a, &a_ref, 1e-13);

        // the wrong variant uses the quadratic (Qua8) geometry with the curved midside nodes
        let mut pad_wrong = aux::gen_pad_curved_qua8(0.0);
        let mut a_wrong = Vector::new(8);
        let mut pad_field_wrong = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_wrong.calc_jacobian(iota).unwrap();
            (pad_field_wrong.fn_interp)(&mut pad_field_wrong.interp, iota);
            let coef = source(p) * det_jac * gauss.weight(p);
            for m in 0..8 {
                a_wrong[m] += pad_field_wrong.interp[m] * coef;
            }
        }
        let mut max_diff = 0.0;
        for m in 0..8 {
            max_diff = f64::max(max_diff, f64::abs(a[m] - a_wrong[m]));
        }
        println!("max |a - a_wrong| (reversed) = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn axisymmetric_curved_element_uses_the_geometric_radius() {
        // A Qua4 field over a curved Qua8 geometry shifted to x > 0. The reference uses the Qua8
        // radius and measure, assembled by hand independently of vec_01_ns.
        let mut pad = aux::gen_pad_mixed_qua4_qua8(3.0);
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();

        let mut pad_geo = aux::gen_pad_curved_qua8(3.0);
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();
        let source = |p: usize| 1.0 + 2.0 * x_ips[p].get(0) - x_ips[p].get(1);

        // reference: the Qua8 measure and radius, and the Qua4 shape functions
        let mut a_ref = Vector::new(4);
        let mut det_jacs = vec![0.0; gauss.npoint()];
        let mut radii = vec![0.0; gauss.npoint()];
        let mut pad_calc = pad_geo.clone();
        let mut pad_field = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            det_jacs[p] = pad_calc.calc_jacobian(iota).unwrap();
            (pad_calc.fn_interp)(&mut pad_calc.interp, iota); // geometry shapes
            for k in 0..8 {
                radii[p] += pad_calc.interp[k] * pad_calc.xxt.get(0, k);
            }
            (pad_field.fn_interp)(&mut pad_field.interp, iota);
            let coef = source(p) * det_jacs[p] * gauss.weight(p) * radii[p];
            for m in 0..4 {
                a_ref[m] += pad_field.interp[m] * coef;
            }
        }

        // the mixed pad must use the geometric radius and measure
        let mut a = Vector::filled(4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.axisymmetric = true;
        integ::vec_01_ns(&mut a, &mut args, |p, _| Ok(source(p))).unwrap();
        vec_approx_eq(&a, &a_ref, 1e-13);

        // the wrong variant uses the bilinear (Qua4) radius with the correct Qua8 measure
        let pad_corners = aux::gen_pad_bilinear_qua4(3.0);
        let mut a_wrong = Vector::new(4);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            (pad_field.fn_interp)(&mut pad_field.interp, iota);
            let mut radius_wrong = 0.0;
            for m in 0..4 {
                radius_wrong += pad_field.interp[m] * pad_corners.xxt.get(0, m);
            }
            let coef = source(p) * det_jacs[p] * gauss.weight(p) * radius_wrong;
            for m in 0..4 {
                a_wrong[m] += pad_field.interp[m] * coef;
            }
        }
        let mut max_diff = 0.0;
        for m in 0..4 {
            max_diff = f64::max(max_diff, f64::abs(a[m] - a_wrong[m]));
        }
        println!("max |a - a_wrong| (axisymmetric) = {:.3e}", max_diff);
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
        let gauss = Gauss::new(pad.kind);
        let mut a_single = Vector::filled(3, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_01_ns(&mut a_single, &mut args, |_, _| Ok(2.5)).unwrap();

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
            integ::vec_01_ns(&mut a, &mut args, |_, _| Ok(2.5)).unwrap();
        }
        for i in 0..2 {
            approx_eq(a[i], 0.5, 1e-15);
        }
        for m in 0..3 {
            approx_eq(a[2 + m], 2.0 * a_single[m], 1e-14);
        }
    }

    #[test]
    fn equal_order_matches_plain_pad() {
        // Tri3 in 2D, including the axisymmetric case
        let pad = aux::gen_pad_tri3();
        for axisymmetric in [false, true] {
            for ngauss in [1, 3, 7] {
                let gauss = Gauss::new_or_sized(pad.kind, Some(ngauss)).unwrap();
                let mut pad_base = pad.clone();
                let mut pad_mixed = aux::gen_pad_equal_order_mixed(&pad);
                let mut a_base = Vector::filled(pad.nnode(), aux::NOISE);
                let mut a_mixed = Vector::filled(pad.nnode(), aux::NOISE);
                let mut args_base = CommonArgs::new(&mut pad_base, &gauss);
                let mut args_mixed = CommonArgs::new(&mut pad_mixed, &gauss);
                args_base.axisymmetric = axisymmetric;
                args_mixed.axisymmetric = axisymmetric;
                let x_ips = recovery::get_points_coords(args_base.pad, &gauss).unwrap();
                let source = |p: usize| 1.0 + 2.0 * x_ips[p].get(0) - x_ips[p].get(1);
                integ::vec_01_ns(&mut a_base, &mut args_base, |p, _| Ok(source(p))).unwrap();
                integ::vec_01_ns(&mut a_mixed, &mut args_mixed, |p, _| Ok(source(p))).unwrap();
                vec_approx_eq(&a_mixed, &a_base, 1e-14);
            }
        }

        // Tet4 in 3D
        let pad = aux::gen_pad_tet4();
        for ngauss in [1, 4] {
            let gauss = Gauss::new_or_sized(pad.kind, Some(ngauss)).unwrap();
            let mut pad_base = pad.clone();
            let mut pad_mixed = aux::gen_pad_equal_order_mixed(&pad);
            let mut a_base = Vector::filled(pad.nnode(), aux::NOISE);
            let mut a_mixed = Vector::filled(pad.nnode(), aux::NOISE);
            let mut args_base = CommonArgs::new(&mut pad_base, &gauss);
            let mut args_mixed = CommonArgs::new(&mut pad_mixed, &gauss);
            let x_ips = recovery::get_points_coords(args_base.pad, &gauss).unwrap();
            let source = |p: usize| 1.0 + 2.0 * x_ips[p].get(0) - x_ips[p].get(1);
            integ::vec_01_ns(&mut a_base, &mut args_base, |p, _| Ok(source(p))).unwrap();
            integ::vec_01_ns(&mut a_mixed, &mut args_mixed, |p, _| Ok(source(p))).unwrap();
            vec_approx_eq(&a_mixed, &a_base, 1e-14);
        }
    }
}
