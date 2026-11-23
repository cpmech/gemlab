use super::CommonArgs;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the shape(N) times scalar(S) times shape(N) integration case 01 (e.g., diffusion matrix)
///
/// Callback function: `s ← f(p, N, B)`
///
/// Diffusion coefficients:
///
/// ```text
///       ⌠
/// Kᵐⁿ = │ Nᵐ s Nⁿ α dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///       nip-1    →     →      →       →
/// Kᵐⁿ ≈   Σ   Nᵐ(ιᵖ) s(ιᵖ) Nⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
///        p=0
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
/// * `kk` -- A matrix containing all `Kᵐⁿ` values, one after another, and
///   sequentially placed as shown above. `m` and `n` are the indices of the nodes.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode` and `ncol(K) ≥ jj0 + nnode`
/// * `args` --- Common arguments
/// * `fn_s` -- Function `f(p,N,B)→s` that computes `s(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`,
///   shape functions N(ιᵖ), and gradients B(ιᵖ).
pub fn mat_01_nsn<F>(kk: &mut Matrix, args: &mut CommonArgs, mut fn_s: F) -> Result<(), StrError>
where
    F: FnMut(usize, &Vector, &Matrix) -> Result<f64, StrError>,
{
    // check
    let nnode = args.pad.interp.dim();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if nrow_kk < ii0 + nnode {
        return Err("nrow(K) must be ≥ ii0 + nnode");
    }
    if ncol_kk < jj0 + nnode {
        return Err("ncol(K) must be ≥ jj0 + nnode");
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

        // calculate interpolation functions, Jacobian, and gradient
        let det_jac = args.pad.calc_gradient(iota)?; // N and B

        // radius @ x(ιᵖ) from the geometry
        let radius = if args.axisymmetric {
            args.pad.calc_radius(iota)?
        } else {
            0.0
        };

        // calculate s
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        let s = fn_s(p, nn, bb)?;

        // calculate coefficient
        let c = if args.axisymmetric {
            s * det_jac * weight * args.alpha * radius
        } else {
            s * det_jac * weight * args.alpha
        };

        // add contribution to K matrix
        for m in 0..nnode {
            for n in 0..nnode {
                kk.add(ii0 + m, jj0 + n, nn[m] * c * nn[n]);
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, AnalyticalQua4, AnalyticalQua8, AnalyticalTet4, AnalyticalTri3, CommonArgs, Gauss};
    use crate::recovery;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::{approx_eq, mat_approx_eq, Matrix, Vector};
    use russell_tensor::Tensor1;

    /// Spatially varying scalar field: s(x) = 1 + 2 x₀ - 3 x₁
    fn scalar_field(x: &Tensor1) -> f64 {
        1.0 + 2.0 * x.get(0) - 3.0 * x.get(1)
    }

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(2, 2);
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let f = |_p: usize, _nn: &Vector, _bb: &Matrix| Ok(0.0);
        assert_eq!(f(0, &nn, &bb).unwrap(), 0.0);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_01_nsn(&mut kk, &mut args, f).err(),
            Some("nrow(K) must be ≥ ii0 + nnode")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_01_nsn(&mut kk, &mut args, f).err(),
            Some("ncol(K) must be ≥ jj0 + nnode")
        );
        args.jj0 = 0;
        // more errors
        assert_eq!(
            integ::mat_01_nsn(&mut kk, &mut args, f).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::mat_01_nsn(&mut kk, &mut args, |_, _, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn tri3_works() {
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let s = 12.0;
        let ana = AnalyticalTri3::new(&pad);
        let kk_correct = ana.mat_01_nsn(s, 1.0);
        let class = pad.kind.class();
        let tolerances = [8.34, 1e-14, 1e-14, 1e-14, 1e-12, 1e-13]; // note how bad rule-1 integ is here
        let selection: Vec<_> = [1, 3, 6, 7, 12, 16]
            .iter()
            .map(|n| Gauss::new_sized(class, *n).unwrap())
            .collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_01_nsn(&mut kk, &mut args, |_, _, _| Ok(s)).unwrap();
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn qua4_works() {
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua4(2.0, 1.0, a, b);
        let mut kk = Matrix::new(4, 4);
        let s = 12.0;
        let ana = AnalyticalQua4::new(a, b);
        let kk_correct = ana.mat_01_nsn(s, 1.0);
        let class = pad.kind.class();
        let tolerances = [7.01, 1e-14, 1e-14, 1e-14]; // note how bad rule-1 integ is here
        let selection: Vec<_> = [1, 4, 9, 16]
            .iter()
            .map(|n| Gauss::new_sized(class, *n).unwrap())
            .collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_01_nsn(&mut kk, &mut args, |_, _, _| Ok(s)).unwrap();
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn qua8_works() {
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua8(2.0, 1.0, a, b);
        let mut kk = Matrix::new(8, 8);
        let s = 3.0;
        let ana = AnalyticalQua8::new(a, b);
        let kk_correct = ana.mat_01_nsn(s, 1.0);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [9, 16].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_01_nsn(&mut kk, &mut args, |_, _, _| Ok(s)).unwrap();
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn tet4_works() {
        let mut pad = aux::gen_pad_tet4();
        let mut kk = Matrix::new(4, 4);
        let s = 3.0;
        let ana = AnalyticalTet4::new(&pad);
        let kk_correct = ana.mat_01_nsn(s);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_01_nsn(&mut kk, &mut args, |_, _, _| Ok(s)).unwrap();
            // println!("{}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn curved_element_uses_the_main_jacobian() {
        // A Qua4 field over a curved Qua8 geometry. The reference uses the Qua8 measure and the
        // Qua4 shape functions, assembled by hand independently of mat_01_nsn.
        let mut pad = aux::gen_pad_mixed_qua4_qua8(0.0);
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();

        // integration-point coordinates from the curved geometry (used by the scalar field)
        let mut pad_geo = aux::gen_pad_curved_qua8(0.0);
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();

        // reference: the Qua8 measure and the Qua4 shape functions
        let mut kk_ref = Matrix::new(4, 4);
        let mut pad_calc = pad_geo.clone();
        let mut pad_field = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_jacobian(iota).unwrap();
            (pad_field.fn_interp)(&mut pad_field.interp, iota);
            let c = scalar_field(&x_ips[p]) * det_jac * gauss.weight(p);
            for m in 0..4 {
                for n in 0..4 {
                    kk_ref.add(m, n, pad_field.interp[m] * c * pad_field.interp[n]);
                }
            }
        }

        // the single mixed pad must reproduce the reference
        let mut kk = Matrix::filled(4, 4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_01_nsn(&mut kk, &mut args, |p, _, _| Ok(scalar_field(&x_ips[p]))).unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the wrong variant uses the bilinear (Qua4) measure
        let mut pad_wrong = aux::gen_pad_bilinear_qua4(0.0);
        let mut kk_wrong = Matrix::new(4, 4);
        let mut pad_field_wrong = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_wrong.calc_jacobian(iota).unwrap();
            (pad_field_wrong.fn_interp)(&mut pad_field_wrong.interp, iota);
            let c = scalar_field(&x_ips[p]) * det_jac * gauss.weight(p);
            for m in 0..4 {
                for n in 0..4 {
                    kk_wrong.add(m, n, pad_field_wrong.interp[m] * c * pad_field_wrong.interp[n]);
                }
            }
        }
        let mut max_diff = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                max_diff = f64::max(max_diff, f64::abs(kk.get(i, j) - kk_wrong.get(i, j)));
            }
        }
        println!("max |K - K_wrong| (curved) = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn reversed_order_qua4_qua8() {
        // A Qua8 field over a straight-edged Qua4 geometry. The reference uses the Qua4 measure and
        // the Qua8 shape functions, assembled by hand independently of mat_01_nsn.
        let mut pad = aux::gen_pad_mixed_qua8_qua4();
        let gauss = Gauss::new_or_sized(GeoKind::Qua4, Some(9)).unwrap();

        // integration-point coordinates from the bilinear geometry (used by the scalar field)
        let mut pad_geo = aux::gen_pad_bilinear_qua4(0.0);
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();

        // reference: the Qua4 measure and the Qua8 shape functions
        let mut kk_ref = Matrix::new(8, 8);
        let mut pad_calc = pad_geo.clone();
        let mut pad_field = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_jacobian(iota).unwrap();
            (pad_field.fn_interp)(&mut pad_field.interp, iota);
            let c = scalar_field(&x_ips[p]) * det_jac * gauss.weight(p);
            for m in 0..8 {
                for n in 0..8 {
                    kk_ref.add(m, n, pad_field.interp[m] * c * pad_field.interp[n]);
                }
            }
        }

        // the single mixed pad must reproduce the reference
        let mut kk = Matrix::filled(8, 8, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_01_nsn(&mut kk, &mut args, |p, _, _| Ok(scalar_field(&x_ips[p]))).unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the wrong variant uses the quadratic (Qua8) geometry with the curved midside nodes
        let mut pad_wrong = aux::gen_pad_curved_qua8(0.0);
        let mut kk_wrong = Matrix::new(8, 8);
        let mut pad_field_wrong = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_wrong.calc_jacobian(iota).unwrap();
            (pad_field_wrong.fn_interp)(&mut pad_field_wrong.interp, iota);
            let c = scalar_field(&x_ips[p]) * det_jac * gauss.weight(p);
            for m in 0..8 {
                for n in 0..8 {
                    kk_wrong.add(m, n, pad_field_wrong.interp[m] * c * pad_field_wrong.interp[n]);
                }
            }
        }
        let mut max_diff = 0.0;
        for i in 0..8 {
            for j in 0..8 {
                max_diff = f64::max(max_diff, f64::abs(kk.get(i, j) - kk_wrong.get(i, j)));
            }
        }
        println!("max |K - K_wrong| (reversed) = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn axisymmetric_curved_element_uses_the_geometric_radius() {
        // A Qua4 field over a curved Qua8 geometry shifted to x > 0. The reference uses the Qua8
        // radius and measure, assembled by hand independently of mat_01_nsn.
        let mut pad = aux::gen_pad_mixed_qua4_qua8(3.0);
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();

        let mut pad_geo = aux::gen_pad_curved_qua8(3.0);
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();

        // reference: the Qua8 measure and radius, and the Qua4 shape functions
        let mut kk_ref = Matrix::new(4, 4);
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
            let c = scalar_field(&x_ips[p]) * det_jacs[p] * gauss.weight(p) * radii[p];
            for m in 0..4 {
                for n in 0..4 {
                    kk_ref.add(m, n, pad_field.interp[m] * c * pad_field.interp[n]);
                }
            }
        }

        // the mixed pad must use the geometric radius and measure
        let mut kk = Matrix::filled(4, 4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.axisymmetric = true;
        integ::mat_01_nsn(&mut kk, &mut args, |p, _, _| Ok(scalar_field(&x_ips[p]))).unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the wrong variant uses the bilinear (Qua4) radius with the correct Qua8 measure
        let pad_corners = aux::gen_pad_bilinear_qua4(3.0);
        let mut kk_wrong = Matrix::new(4, 4);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            (pad_field.fn_interp)(&mut pad_field.interp, iota);
            let mut radius_wrong = 0.0;
            for m in 0..4 {
                radius_wrong += pad_field.interp[m] * pad_corners.xxt.get(0, m);
            }
            let c = scalar_field(&x_ips[p]) * det_jacs[p] * gauss.weight(p) * radius_wrong;
            for m in 0..4 {
                for n in 0..4 {
                    kk_wrong.add(m, n, pad_field.interp[m] * c * pad_field.interp[n]);
                }
            }
        }
        let mut max_diff = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                max_diff = f64::max(max_diff, f64::abs(kk.get(i, j) - kk_wrong.get(i, j)));
            }
        }
        println!("max |K - K_wrong| (axisymmetric) = {:.3e}", max_diff);
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
        let mut kk_single = Matrix::filled(3, 3, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_01_nsn(&mut kk_single, &mut args, |_, _, _| Ok(2.5)).unwrap();

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
            integ::mat_01_nsn(&mut kk, &mut args, |_, _, _| Ok(2.5)).unwrap();
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

    #[test]
    fn equal_order_matches_plain_pad() {
        // Tri3 in 2D, including the axisymmetric case
        let pad = aux::gen_pad_tri3();
        for axisymmetric in [false, true] {
            for ngauss in [1, 3, 7] {
                let gauss = Gauss::new_or_sized(pad.kind, Some(ngauss)).unwrap();
                let mut pad_base = pad.clone();
                let mut pad_mixed = aux::gen_pad_equal_order_mixed(&pad);
                let mut kk_base = Matrix::filled(pad.nnode(), pad.nnode(), aux::NOISE);
                let mut kk_mixed = Matrix::filled(pad.nnode(), pad.nnode(), aux::NOISE);
                let mut args_base = CommonArgs::new(&mut pad_base, &gauss);
                let mut args_mixed = CommonArgs::new(&mut pad_mixed, &gauss);
                args_base.axisymmetric = axisymmetric;
                args_mixed.axisymmetric = axisymmetric;
                let x_ips = recovery::get_points_coords(args_base.pad, &gauss).unwrap();
                integ::mat_01_nsn(&mut kk_base, &mut args_base, |p, _, _| Ok(scalar_field(&x_ips[p]))).unwrap();
                integ::mat_01_nsn(&mut kk_mixed, &mut args_mixed, |p, _, _| Ok(scalar_field(&x_ips[p]))).unwrap();
                mat_approx_eq(&kk_mixed, &kk_base, 1e-14);
            }
        }

        // Tet4 in 3D
        let pad = aux::gen_pad_tet4();
        for ngauss in [1, 4] {
            let gauss = Gauss::new_or_sized(pad.kind, Some(ngauss)).unwrap();
            let mut pad_base = pad.clone();
            let mut pad_mixed = aux::gen_pad_equal_order_mixed(&pad);
            let mut kk_base = Matrix::filled(pad.nnode(), pad.nnode(), aux::NOISE);
            let mut kk_mixed = Matrix::filled(pad.nnode(), pad.nnode(), aux::NOISE);
            let mut args_base = CommonArgs::new(&mut pad_base, &gauss);
            let mut args_mixed = CommonArgs::new(&mut pad_mixed, &gauss);
            let x_ips = recovery::get_points_coords(args_base.pad, &gauss).unwrap();
            integ::mat_01_nsn(&mut kk_base, &mut args_base, |p, _, _| Ok(scalar_field(&x_ips[p]))).unwrap();
            integ::mat_01_nsn(&mut kk_mixed, &mut args_mixed, |p, _, _| Ok(scalar_field(&x_ips[p]))).unwrap();
            mat_approx_eq(&kk_mixed, &kk_base, 1e-14);
        }
    }
}
