use super::CommonArgs;
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor1;

/// Implements the gradient(B) dot vector(V) integration case 03
///
/// Callback function: `f(w, p, N, B)`
///
/// Vector dot gradient:
///
/// ```text
///      ⌠ →  → →     → →
/// cᵐ = │ Bᵐ(x(ξ)) · w(x) α dΩ
///      ⌡
///      Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///      nip-1  →  →     → →       →
/// cᵐ ≈   Σ    Bᵐ(ιᵖ) · w(ιᵖ) |J|(ιᵖ) wᵖ α
///       p=0
/// ```
///
/// # Results
///
/// ```text
///     ┌     ┐
///     |  c⁰ |  ⟸  ii0 = 0
///     |  c¹ |
/// c = |  c² |
///     | ··· |
///     |  cᵐ |  ⟸  ii
///     └     ┘
/// ```
///
/// # Arguments
///
/// * `c` -- A vector containing all `cᵐ` values, one after another, and sequentially placed as shown above.
///   `m` is the index of the node. The length must be `c.len() ≥ ii0 + nnode`.
/// * `args` --- Common arguments
/// * `fn_w` -- Function `f(w,p,N,B)` that computes `w(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`,
///   shape functions N(ιᵖ), and gradients B(ιᵖ). `w` has 3 components (only `space_ndim` are used).
///
/// # Examples
///
/// See also the `examples` directory.
///
/// ```
/// use gemlab::integ::{self, CommonArgs, Gauss};
/// use gemlab::shapes::{GeoKind, Scratchpad};
/// use gemlab::StrError;
/// use russell_lab::{Vector, vec_approx_eq};
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
///     let mut c = Vector::filled(pad.nnode(), 0.0);
///     let mut args = CommonArgs::new(&mut pad, &gauss);
///     integ::vec_03_bv(&mut c, &mut args, |w, _, _, _| {
///         w.set(0, 1.0);
///         w.set(1, 2.0);
///         Ok(())
///     })?;
///     // solution (A = 6):
///     // cᵐ = (w₀ Bᵐ₀ + w₁ Bᵐ₁) A
///     //     ┌       ┐
///     //     │ -¼ -⅓ │
///     // B = │  ¼  0 │
///     //     │  0  ⅓ │
///     //     └       ┘
///     vec_approx_eq(&c, &[-5.5, 1.5, 4.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_03_bv<F>(c: &mut Vector, args: &mut CommonArgs, mut fn_w: F) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor1, usize, &Vector, &Matrix) -> Result<(), StrError>,
{
    // check
    let (space_ndim, _, nnode) = args.pad.dims();
    let ii0 = args.ii0;
    if c.dim() < ii0 + nnode {
        return Err("c.len() must be ≥ ii0 + nnode");
    }

    // allocate auxiliary vector
    let mut w = Tensor1::new();

    // clear output vector
    if args.clear {
        c.fill(0.0);
    }

    // loop over integration points
    for p in 0..args.gauss.npoint() {
        // ksi coordinates and weight
        let iota = args.gauss.coords(p);
        let weight = args.gauss.weight(p);

        // calculate interpolation functions, Jacobian, and gradient
        let det_jac = args.pad.calc_gradient(iota)?; // N and B

        // radius @ x(ιᵖ) from the geometry
        let radius = if args.axisymmetric {
            args.pad.calc_radius(iota)?
        } else {
            0.0
        };

        // calculate w
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        fn_w(&mut w, p, nn, bb)?;

        // calculate coefficient
        let coef = if args.axisymmetric {
            det_jac * weight * args.alpha * radius
        } else {
            det_jac * weight * args.alpha
        };

        // add contribution to c vector
        if space_ndim == 2 {
            for m in 0..nnode {
                c[ii0 + m] += coef * (w.get(0) * bb.get(m, 0) + w.get(1) * bb.get(m, 1));
            }
        } else {
            for m in 0..nnode {
                c[ii0 + m] += coef * (w.get(0) * bb.get(m, 0) + w.get(1) * bb.get(m, 1) + w.get(2) * bb.get(m, 2));
            }
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
    use russell_lab::{approx_eq, vec_approx_eq, Matrix, Vector};
    use russell_tensor::Tensor1;

    /// Spatially varying vector field: w(x) = {2 + x₀, -1 + 3 x₁, 0.5 - x₂}
    fn vector_field(w: &mut Tensor1, x: &Tensor1) {
        w.set(0, 2.0 + x.get(0));
        w.set(1, -1.0 + 3.0 * x.get(1));
        w.set(2, 0.5 - x.get(2));
    }

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut c = Vector::new(2);
        let mut w = Tensor1::new();
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let f = |_: &mut Tensor1, _: usize, _: &Vector, _: &Matrix| Ok(());
        f(&mut w, 0, &nn, &bb).unwrap();
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::vec_03_bv(&mut c, &mut args, f).err(),
            Some("c.len() must be ≥ ii0 + nnode")
        );
    }

    #[test]
    fn tri3_constant_works() {
        // constant vector function: w(x) = {w₀, w₁}
        const W0: f64 = 2.0;
        const W1: f64 = 3.0;
        let mut pad = aux::gen_pad_tri3();

        // solution
        let ana = AnalyticalTri3::new(&pad);
        let c_correct = ana.vec_03_bv(W0, W1);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [1, 3].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();

        // check
        let mut c = Vector::filled(pad.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_03_bv(&mut c, &mut args, |w, _, _, _| {
                w.set(0, W0);
                w.set(1, W1);
                Ok(())
            })
            .unwrap();
            vec_approx_eq(&c, &c_correct, tol);
        });
    }

    #[test]
    fn tri3_bilinear_works() {
        // bilinear vector function: w(x) = {x, y}
        let mut pad = aux::gen_pad_tri3();

        // solution
        let ana = AnalyticalTri3::new(&pad);
        let c_correct = ana.vec_03_bv_bilinear(&pad);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [1, 3].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();

        // check
        let mut c = Vector::filled(pad.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            let x_ips = recovery::get_points_coords(args.pad, ips).unwrap();
            integ::vec_03_bv(&mut c, &mut args, |w, p, _, _| {
                w.set(0, x_ips[p].get(0));
                w.set(1, x_ips[p].get(1));
                Ok(())
            })
            .unwrap();
            vec_approx_eq(&c, &c_correct, tol);
        });
    }

    #[test]
    fn tet4_constant_works() {
        // tet 4 with constant vector  w(x) = {w0, w1, w2}
        let mut pad = aux::gen_pad_tet4();

        // solution
        const W0: f64 = 2.0;
        const W1: f64 = 3.0;
        const W2: f64 = 4.0;
        let ana = AnalyticalTet4::new(&pad);
        let c_correct = ana.vec_03_bv(W0, W1, W2);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14];
        let selection: Vec<_> = [1, 4, 5, 8, 14, 15, 24]
            .iter()
            .map(|n| Gauss::new_sized(class, *n).unwrap())
            .collect();

        // check
        let mut c = Vector::filled(pad.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_03_bv(&mut c, &mut args, |w, _, _, _| {
                w.set(0, W0);
                w.set(1, W1);
                w.set(2, W2);
                Ok(())
            })
            .unwrap();
            vec_approx_eq(&c, &c_correct, tol);
        });
    }

    #[test]
    fn curved_element_uses_the_main_jacobian() {
        // A Qua4 field over a curved Qua8 geometry. The reference uses dNa/dξ · J⁻¹ and the Qua8
        // measure, assembled by hand independently of vec_03_bv.
        let mut pad = aux::gen_pad_mixed_qua4_qua8(0.0);
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();

        // integration-point coordinates from the curved geometry (used by the vector field)
        let mut pad_geo = aux::gen_pad_curved_qua8(0.0);
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();

        // reference: dNa/dξ · J⁻¹ and the Qua8 measure
        let mut c_ref = Vector::new(4);
        let mut pad_calc = pad_geo.clone();
        let mut pad_field = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_field.fn_deriv)(&mut pad_field.deriv, iota);
            let (w0, w1) = (2.0 + x_ips[p].get(0), -1.0 + 3.0 * x_ips[p].get(1));
            let coef = det_jac * gauss.weight(p);
            for m in 0..4 {
                let mut ba0 = 0.0;
                let mut ba1 = 0.0;
                for k in 0..2 {
                    ba0 += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, 0);
                    ba1 += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, 1);
                }
                c_ref[m] += coef * (w0 * ba0 + w1 * ba1);
            }
        }

        // the single mixed pad must reproduce the reference
        let mut c = Vector::filled(4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_03_bv(&mut c, &mut args, |w, p, _, _| {
            vector_field(w, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        vec_approx_eq(&c, &c_ref, 1e-13);

        // the wrong variant uses Ba and detJ from the bilinear (Qua4) map
        let mut pad_wrong = aux::gen_pad_bilinear_qua4(0.0);
        let mut c_wrong = Vector::new(4);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_wrong.calc_gradient(iota).unwrap();
            let (w0, w1) = (2.0 + x_ips[p].get(0), -1.0 + 3.0 * x_ips[p].get(1));
            let coef = det_jac * gauss.weight(p);
            for m in 0..4 {
                c_wrong[m] += coef * (w0 * pad_wrong.gradient.get(m, 0) + w1 * pad_wrong.gradient.get(m, 1));
            }
        }
        let mut max_diff = 0.0;
        for m in 0..4 {
            max_diff = f64::max(max_diff, f64::abs(c[m] - c_wrong[m]));
        }
        println!("max |c - c_wrong| (curved) = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn reversed_order_qua4_qua8() {
        // A Qua8 field over a straight-edged Qua4 geometry. The reference uses dNa/dξ · J⁻¹ and
        // the Qua4 measure, assembled by hand independently of vec_03_bv.
        let mut pad = aux::gen_pad_mixed_qua8_qua4();
        let gauss = Gauss::new_or_sized(GeoKind::Qua4, Some(9)).unwrap();

        // integration-point coordinates from the bilinear geometry (used by the vector field)
        let mut pad_geo = aux::gen_pad_bilinear_qua4(0.0);
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();

        // reference: dNa/dξ · J⁻¹ and the Qua4 measure
        let mut c_ref = Vector::new(8);
        let mut pad_calc = pad_geo.clone();
        let mut pad_field = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_field.fn_deriv)(&mut pad_field.deriv, iota);
            let (w0, w1) = (2.0 + x_ips[p].get(0), -1.0 + 3.0 * x_ips[p].get(1));
            let coef = det_jac * gauss.weight(p);
            for m in 0..8 {
                let mut ba0 = 0.0;
                let mut ba1 = 0.0;
                for k in 0..2 {
                    ba0 += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, 0);
                    ba1 += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, 1);
                }
                c_ref[m] += coef * (w0 * ba0 + w1 * ba1);
            }
        }

        // the single mixed pad must reproduce the reference
        let mut c = Vector::filled(8, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_03_bv(&mut c, &mut args, |w, p, _, _| {
            vector_field(w, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        vec_approx_eq(&c, &c_ref, 1e-13);

        // the wrong variant uses Ba and detJ from the quadratic (Qua8) map
        let mut pad_wrong = aux::gen_pad_curved_qua8(0.0);
        let mut c_wrong = Vector::new(8);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_wrong.calc_gradient(iota).unwrap();
            let (w0, w1) = (2.0 + x_ips[p].get(0), -1.0 + 3.0 * x_ips[p].get(1));
            let coef = det_jac * gauss.weight(p);
            for m in 0..8 {
                c_wrong[m] += coef * (w0 * pad_wrong.gradient.get(m, 0) + w1 * pad_wrong.gradient.get(m, 1));
            }
        }
        let mut max_diff = 0.0;
        for m in 0..8 {
            max_diff = f64::max(max_diff, f64::abs(c[m] - c_wrong[m]));
        }
        println!("max |c - c_wrong| (reversed) = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn axisymmetric_curved_element_uses_the_geometric_radius() {
        // A Qua4 field over a curved Qua8 geometry shifted to x > 0. The reference uses the Qua8
        // radius and measure, assembled by hand independently of vec_03_bv.
        let mut pad = aux::gen_pad_mixed_qua4_qua8(3.0);
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();

        let mut pad_geo = aux::gen_pad_curved_qua8(3.0);
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();

        // reference: dNa/dξ · J⁻¹, the Qua8 measure, and the Qua8 radius
        let mut c_ref = Vector::new(4);
        let mut det_jacs = vec![0.0; gauss.npoint()];
        let mut radii = vec![0.0; gauss.npoint()];
        let mut pad_calc = pad_geo.clone();
        let mut pad_field = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            det_jacs[p] = pad_calc.calc_gradient(iota).unwrap();
            (pad_calc.fn_interp)(&mut pad_calc.interp, iota); // geometry shapes
            for k in 0..8 {
                radii[p] += pad_calc.interp[k] * pad_calc.xxt.get(0, k);
            }
            (pad_field.fn_deriv)(&mut pad_field.deriv, iota);
            let (w0, w1) = (2.0 + x_ips[p].get(0), -1.0 + 3.0 * x_ips[p].get(1));
            let coef = det_jacs[p] * gauss.weight(p) * radii[p];
            for m in 0..4 {
                let mut ba0 = 0.0;
                let mut ba1 = 0.0;
                for k in 0..2 {
                    ba0 += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, 0);
                    ba1 += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, 1);
                }
                c_ref[m] += coef * (w0 * ba0 + w1 * ba1);
            }
        }

        // the mixed pad must use the geometric radius and measure
        let mut c = Vector::filled(4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.axisymmetric = true;
        integ::vec_03_bv(&mut c, &mut args, |w, p, _, _| {
            vector_field(w, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        vec_approx_eq(&c, &c_ref, 1e-13);

        // the wrong variant uses the bilinear (Qua4) radius with the correct Qua8 measure and Ba
        let pad_corners = aux::gen_pad_bilinear_qua4(3.0);
        let mut c_wrong = Vector::new(4);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            pad_calc.calc_gradient(iota).unwrap();
            (pad_field.fn_interp)(&mut pad_field.interp, iota);
            let mut radius_wrong = 0.0;
            for m in 0..4 {
                radius_wrong += pad_field.interp[m] * pad_corners.xxt.get(0, m);
            }
            (pad_field.fn_deriv)(&mut pad_field.deriv, iota);
            let (w0, w1) = (2.0 + x_ips[p].get(0), -1.0 + 3.0 * x_ips[p].get(1));
            let coef = det_jacs[p] * gauss.weight(p) * radius_wrong;
            for m in 0..4 {
                let mut ba0 = 0.0;
                let mut ba1 = 0.0;
                for k in 0..2 {
                    ba0 += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, 0);
                    ba1 += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, 1);
                }
                c_wrong[m] += coef * (w0 * ba0 + w1 * ba1);
            }
        }
        let mut max_diff = 0.0;
        for m in 0..4 {
            max_diff = f64::max(max_diff, f64::abs(c[m] - c_wrong[m]));
        }
        println!("max |c - c_wrong| (axisymmetric) = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn offset_and_accumulation_work() {
        // single run used as the reference values
        let mut pad = aux::gen_pad_qua4(0.0, 0.0, 2.0, 3.0);
        let gauss = Gauss::new(pad.kind);
        let mut c_single = Vector::filled(4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_03_bv(&mut c_single, &mut args, |w, _, _, _| {
            w.set(0, 1.0);
            w.set(1, 2.0);
            Ok(())
        })
        .unwrap();

        // two runs accumulating at ii0 = 3; the vector is prefilled with 0.5 to prove that
        // clear = false leaves the leading entries untouched
        let mut c = Vector::filled(8, 0.5);
        for m in 0..4 {
            c[3 + m] = 0.0; // start the accumulated slots from zero
        }
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.clear = false;
        args.ii0 = 3;
        for _ in 0..2 {
            integ::vec_03_bv(&mut c, &mut args, |w, _, _, _| {
                w.set(0, 1.0);
                w.set(1, 2.0);
                Ok(())
            })
            .unwrap();
        }
        for i in 0..3 {
            approx_eq(c[i], 0.5, 1e-15);
        }
        for m in 0..4 {
            approx_eq(c[3 + m], 2.0 * c_single[m], 1e-14);
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
                let mut c_base = Vector::filled(pad.nnode(), aux::NOISE);
                let mut c_mixed = Vector::filled(pad.nnode(), aux::NOISE);
                let mut args_base = CommonArgs::new(&mut pad_base, &gauss);
                let mut args_mixed = CommonArgs::new(&mut pad_mixed, &gauss);
                args_base.axisymmetric = axisymmetric;
                args_mixed.axisymmetric = axisymmetric;
                let x_ips = recovery::get_points_coords(args_base.pad, &gauss).unwrap();
                integ::vec_03_bv(&mut c_base, &mut args_base, |w, p, _, _| {
                    vector_field(w, &x_ips[p]);
                    Ok(())
                })
                .unwrap();
                integ::vec_03_bv(&mut c_mixed, &mut args_mixed, |w, p, _, _| {
                    vector_field(w, &x_ips[p]);
                    Ok(())
                })
                .unwrap();
                vec_approx_eq(&c_mixed, &c_base, 1e-14);
            }
        }

        // Tet4 in 3D
        let pad = aux::gen_pad_tet4();
        for ngauss in [1, 4] {
            let gauss = Gauss::new_or_sized(pad.kind, Some(ngauss)).unwrap();
            let mut pad_base = pad.clone();
            let mut pad_mixed = aux::gen_pad_equal_order_mixed(&pad);
            let mut c_base = Vector::filled(pad.nnode(), aux::NOISE);
            let mut c_mixed = Vector::filled(pad.nnode(), aux::NOISE);
            let mut args_base = CommonArgs::new(&mut pad_base, &gauss);
            let mut args_mixed = CommonArgs::new(&mut pad_mixed, &gauss);
            let x_ips = recovery::get_points_coords(args_base.pad, &gauss).unwrap();
            integ::vec_03_bv(&mut c_base, &mut args_base, |w, p, _, _| {
                vector_field(w, &x_ips[p]);
                Ok(())
            })
            .unwrap();
            integ::vec_03_bv(&mut c_mixed, &mut args_mixed, |w, p, _, _| {
                vector_field(w, &x_ips[p]);
                Ok(())
            })
            .unwrap();
            vec_approx_eq(&c_mixed, &c_base, 1e-14);
        }
    }
}
