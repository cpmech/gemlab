use super::CommonArgs;
use crate::StrError;
use russell_lab::math::SQRT_2;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the gradient(B) dot tensor(T) dot gradient(B) integration case 03 (e.g., conductivity matrix)
///
/// Callback function: `f(T, p, N, B)`
///
/// Conductivity coefficients:
///
/// ```text
///       ⌠ →        →
/// Kᵐⁿ = │ Bᵐ ⋅ T ⋅ Bⁿ α dΩ
///       ⌡      ▔
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///       nip-1 →  →       →     →  →       →
/// Kᵐⁿ ≈   Σ   Bᵐ(ιᵖ) ⋅ T(ιᵖ) ⋅ Bⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
///        p=0           ▔
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
/// * `fn_tt` -- Function `f(T,p,N,B)` that computes `T(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`,
///   shape functions N(ιᵖ), and gradients B(ιᵖ). `T` is set for `space_ndim`.
pub fn mat_03_btb<const N: usize, F>(kk: &mut Matrix, args: &mut CommonArgs, mut fn_tt: F) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor2<N>, usize, &Vector, &Matrix) -> Result<(), StrError>,
{
    // check
    let (space_ndim, _, nnode) = args.pad.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if nrow_kk < ii0 + nnode {
        return Err("nrow(K) must be ≥ ii0 + nnode");
    }
    if ncol_kk < jj0 + nnode {
        return Err("ncol(K) must be ≥ jj0 + nnode");
    }
    if N != 2 * space_ndim {
        return Err("tensor dimension N must equal 2 * space_ndim");
    }

    // allocate auxiliary tensor
    let mut tt = Tensor2::<N>::new();

    // clear output matrix
    if args.clear {
        kk.fill(0.0);
    }

    // loop over integration points
    let s = SQRT_2;
    for p in 0..args.gauss.npoint() {
        // ksi coordinates and weight
        let iota = args.gauss.coords(p);
        let weight = args.gauss.weight(p);

        // calculate Jacobian and gradient
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B

        // radius @ x(ιᵖ) from the geometry
        let radius = if args.axisymmetric {
            args.pad.calc_radius(iota)?
        } else {
            0.0
        };

        // calculate T tensor
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        fn_tt(&mut tt, p, nn, bb)?;

        // calculate coefficient
        let c = if args.axisymmetric {
            det_jac * weight * args.alpha * radius
        } else {
            det_jac * weight * args.alpha
        };

        // add contribution to K matrix
        if space_ndim == 2 {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk.add(
                        ii0 + m,
                        jj0 + n,
                        c * (bb.get(n, 1) * (tt.get(1) * bb.get(m, 1) + (tt.get(3) * bb.get(m, 0)) / s)
                            + bb.get(n, 0) * (tt.get(0) * bb.get(m, 0) + (tt.get(3) * bb.get(m, 1)) / s)),
                    );
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode {
                    kk.add(
                        ii0 + m,
                        jj0 + n,
                        c * (bb.get(n, 2)
                            * (tt.get(2) * bb.get(m, 2)
                                + (tt.get(5) * bb.get(m, 0)) / s
                                + (tt.get(4) * bb.get(m, 1)) / s)
                            + bb.get(n, 1)
                                * (tt.get(1) * bb.get(m, 1)
                                    + (tt.get(3) * bb.get(m, 0)) / s
                                    + (tt.get(4) * bb.get(m, 2)) / s)
                            + bb.get(n, 0)
                                * (tt.get(0) * bb.get(m, 0)
                                    + (tt.get(3) * bb.get(m, 1)) / s
                                    + (tt.get(5) * bb.get(m, 2)) / s)),
                    );
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
    use crate::integ::{self, AnalyticalQua4, AnalyticalQua8, AnalyticalTet4, AnalyticalTri3, CommonArgs, Gauss};
    use crate::recovery;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::{approx_eq, mat_approx_eq, Matrix, Vector};
    use russell_tensor::{Tensor1, Tensor2};

    /// Spatially varying full anisotropic tensor in 2D (symmetric Mandel components)
    fn tensor_field_2d(tt: &mut Tensor2<4>, x: &Tensor1) {
        tt.sym_set_std(0, 0, 2.0 + x.get(0));
        tt.sym_set_std(1, 1, 3.0 + 0.5 * x.get(1));
        tt.sym_set_std(2, 2, 1.0 + 0.25 * (x.get(0) + x.get(1)));
        tt.sym_set_std(0, 1, 0.5 + 0.1 * x.get(0));
    }

    /// Spatially varying full anisotropic tensor in 3D (symmetric Mandel components)
    fn tensor_field_3d(tt: &mut Tensor2<6>, x: &Tensor1) {
        tt.sym_set_std(0, 0, 2.0 + x.get(0));
        tt.sym_set_std(1, 1, 3.0 + 0.5 * x.get(1));
        tt.sym_set_std(2, 2, 1.0 + 0.25 * x.get(2));
        tt.sym_set_std(0, 1, 0.5 + 0.1 * x.get(0));
        tt.sym_set_std(1, 2, -0.3 + 0.2 * x.get(1));
        tt.sym_set_std(2, 0, 0.4 - 0.15 * x.get(2));
    }

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(2, 2);
        let mut tt = Tensor2::<4>::new();
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let f = |_tt: &mut Tensor2<4>, _p: usize, _nn: &Vector, _bb: &Matrix| Ok(());
        f(&mut tt, 0, &nn, &bb).unwrap();
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_03_btb::<4, _>(&mut kk, &mut args, f).err(),
            Some("nrow(K) must be ≥ ii0 + nnode")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_03_btb::<4, _>(&mut kk, &mut args, f).err(),
            Some("ncol(K) must be ≥ jj0 + nnode")
        );
        args.jj0 = 0;
        // more errors
        assert_eq!(
            integ::mat_03_btb::<4, _>(&mut kk, &mut args, f).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_qua4(0.0, 0.0, 1.0, 1.0);
        let mut kk = Matrix::new(4, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::mat_03_btb::<4, _>(&mut kk, &mut args, |_, _, _, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn tri3_works() {
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let ana = AnalyticalTri3::new(&pad);
        let (kx, ky) = (2.5, 3.8);
        let kk_correct = ana.mat_03_btb(kx, ky, false);
        let class = pad.kind.class();
        let tolerances = [1e-15, 1e-15, 1e-15];
        let selection: Vec<_> = [1, 3, 7].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_03_btb::<4, _>(&mut kk, &mut args, |tt, _, _, _| {
                tt.sym_set_std(0, 0, kx);
                tt.sym_set_std(1, 1, ky);
                Ok(())
            })
            .unwrap();
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn qua4_works() {
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua4(2.0, 1.0, a, b);
        let mut kk = Matrix::new(4, 4);
        let ana = AnalyticalQua4::new(a, b);
        let (kx, ky) = (2.5, 3.8);
        let kk_correct = ana.mat_03_btb(kx, ky);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_03_btb::<4, _>(&mut kk, &mut args, |tt, _, _, _| {
                tt.sym_set_std(0, 0, kx);
                tt.sym_set_std(1, 1, ky);
                Ok(())
            })
            .unwrap();
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn qua8_works() {
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua8(2.0, 1.0, a, b);
        let mut kk = Matrix::new(8, 8);
        let ana = AnalyticalQua8::new(a, b);
        let (kx, ky) = (2.5, 3.8);
        let kk_correct = ana.mat_03_btb(kx, ky);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [9, 16].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_03_btb::<4, _>(&mut kk, &mut args, |tt, _, _, _| {
                tt.sym_set_std(0, 0, kx);
                tt.sym_set_std(1, 1, ky);
                Ok(())
            })
            .unwrap();
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn tet4_works() {
        let mut pad = aux::gen_pad_tet4();
        let mut kk = Matrix::new(4, 4);
        let ana = AnalyticalTet4::new(&pad);
        #[rustfmt::skip]
        let sig = Tensor2::<6>::from_std_matrix(&[
            [1.1, 1.2, 1.3],
            [1.2, 2.2, 2.3],
            [1.3, 2.3, 3.3]]).unwrap();
        let kk_correct = ana.mat_03_btb(&sig);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14];
        let selection: Vec<_> = [4].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_03_btb::<6, _>(&mut kk, &mut args, |tt, _, _, _| {
                tt.set_tensor(1.0, &sig);
                Ok(())
            })
            .unwrap();
            // println!("{}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn curved_element_uses_the_main_jacobian() {
        // A Qua4 field over a curved Qua8 geometry. The reference uses Ba = dNa/dξ · J⁻¹ and the
        // Qua8 measure, assembled by hand independently of mat_03_btb.
        let mut pad = aux::gen_pad_mixed_qua4_qua8(0.0);
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();

        // integration-point coordinates from the curved geometry (used by the tensor field)
        let mut pad_geo = aux::gen_pad_curved_qua8(0.0);
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();

        // reference: Ba = dNa/dξ · J⁻¹ and the Qua8 measure
        let mut kk_ref = Matrix::new(4, 4);
        let mut pad_calc = pad_geo.clone();
        let mut pad_field = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_field.fn_deriv)(&mut pad_field.deriv, iota);
            let mut tt = Tensor2::<4>::new();
            tensor_field_2d(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            let mut ba = Matrix::new(4, 2);
            for m in 0..4 {
                for j in 0..2 {
                    let mut g = 0.0;
                    for k in 0..2 {
                        g += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, j);
                    }
                    ba.set(m, j, g);
                }
            }
            for m in 0..4 {
                for n in 0..4 {
                    let mut sum = 0.0;
                    for i in 0..2 {
                        for j in 0..2 {
                            sum += ba.get(m, i) * tt.get_std(i, j) * ba.get(n, j);
                        }
                    }
                    kk_ref.add(m, n, c * sum);
                }
            }
        }

        // the single mixed pad must reproduce the reference
        let mut kk = Matrix::filled(4, 4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_03_btb::<4, _>(&mut kk, &mut args, |tt, p, _, _| {
            tensor_field_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the wrong variant uses Ba and detJ from the bilinear (Qua4) map
        let mut pad_wrong = aux::gen_pad_bilinear_qua4(0.0);
        let mut kk_wrong = Matrix::new(4, 4);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_wrong.calc_gradient(iota).unwrap();
            let mut tt = Tensor2::<4>::new();
            tensor_field_2d(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            for m in 0..4 {
                for n in 0..4 {
                    let mut sum = 0.0;
                    for i in 0..2 {
                        for j in 0..2 {
                            sum += pad_wrong.gradient.get(m, i) * tt.get_std(i, j) * pad_wrong.gradient.get(n, j);
                        }
                    }
                    kk_wrong.add(m, n, c * sum);
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
        // A Qua8 field over a straight-edged Qua4 geometry. The reference uses Ba = dNa/dξ · J⁻¹
        // and the Qua4 measure, assembled by hand independently of mat_03_btb.
        let mut pad = aux::gen_pad_mixed_qua8_qua4();
        let gauss = Gauss::new_or_sized(GeoKind::Qua4, Some(9)).unwrap();

        // integration-point coordinates from the bilinear geometry (used by the tensor field)
        let mut pad_geo = aux::gen_pad_bilinear_qua4(0.0);
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();

        // reference: Ba = dNa/dξ · J⁻¹ and the Qua4 measure
        let mut kk_ref = Matrix::new(8, 8);
        let mut pad_calc = pad_geo.clone();
        let mut pad_field = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_field.fn_deriv)(&mut pad_field.deriv, iota);
            let mut tt = Tensor2::<4>::new();
            tensor_field_2d(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            let mut ba = Matrix::new(8, 2);
            for m in 0..8 {
                for j in 0..2 {
                    let mut g = 0.0;
                    for k in 0..2 {
                        g += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, j);
                    }
                    ba.set(m, j, g);
                }
            }
            for m in 0..8 {
                for n in 0..8 {
                    let mut sum = 0.0;
                    for i in 0..2 {
                        for j in 0..2 {
                            sum += ba.get(m, i) * tt.get_std(i, j) * ba.get(n, j);
                        }
                    }
                    kk_ref.add(m, n, c * sum);
                }
            }
        }

        // the single mixed pad must reproduce the reference
        let mut kk = Matrix::filled(8, 8, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_03_btb::<4, _>(&mut kk, &mut args, |tt, p, _, _| {
            tensor_field_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the wrong variant uses Ba and detJ from the quadratic (Qua8) map
        let mut pad_wrong = aux::gen_pad_curved_qua8(0.0);
        let mut kk_wrong = Matrix::new(8, 8);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_wrong.calc_gradient(iota).unwrap();
            let mut tt = Tensor2::<4>::new();
            tensor_field_2d(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            for m in 0..8 {
                for n in 0..8 {
                    let mut sum = 0.0;
                    for i in 0..2 {
                        for j in 0..2 {
                            sum += pad_wrong.gradient.get(m, i) * tt.get_std(i, j) * pad_wrong.gradient.get(n, j);
                        }
                    }
                    kk_wrong.add(m, n, c * sum);
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
        // radius and measure, assembled by hand independently of mat_03_btb.
        let mut pad = aux::gen_pad_mixed_qua4_qua8(3.0);
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();

        let mut pad_geo = aux::gen_pad_curved_qua8(3.0);
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();

        // reference: Ba = dNa/dξ · J⁻¹, the Qua8 measure, and the Qua8 radius
        let mut kk_ref = Matrix::new(4, 4);
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
            let mut tt = Tensor2::<4>::new();
            tensor_field_2d(&mut tt, &x_ips[p]);
            let c = det_jacs[p] * gauss.weight(p) * radii[p];
            let mut ba = Matrix::new(4, 2);
            for m in 0..4 {
                for j in 0..2 {
                    let mut g = 0.0;
                    for k in 0..2 {
                        g += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, j);
                    }
                    ba.set(m, j, g);
                }
            }
            for m in 0..4 {
                for n in 0..4 {
                    let mut sum = 0.0;
                    for i in 0..2 {
                        for j in 0..2 {
                            sum += ba.get(m, i) * tt.get_std(i, j) * ba.get(n, j);
                        }
                    }
                    kk_ref.add(m, n, c * sum);
                }
            }
        }

        // the mixed pad must use the geometric radius and measure
        let mut kk = Matrix::filled(4, 4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.axisymmetric = true;
        integ::mat_03_btb::<4, _>(&mut kk, &mut args, |tt, p, _, _| {
            tensor_field_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the wrong variant uses the bilinear (Qua4) radius with the correct Qua8 measure and Ba
        let pad_corners = aux::gen_pad_bilinear_qua4(3.0);
        let mut kk_wrong = Matrix::new(4, 4);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            pad_calc.calc_gradient(iota).unwrap();
            (pad_field.fn_interp)(&mut pad_field.interp, iota);
            let mut radius_wrong = 0.0;
            for m in 0..4 {
                radius_wrong += pad_field.interp[m] * pad_corners.xxt.get(0, m);
            }
            (pad_field.fn_deriv)(&mut pad_field.deriv, iota);
            let mut tt = Tensor2::<4>::new();
            tensor_field_2d(&mut tt, &x_ips[p]);
            let c = det_jacs[p] * gauss.weight(p) * radius_wrong;
            let mut ba = Matrix::new(4, 2);
            for m in 0..4 {
                for j in 0..2 {
                    let mut g = 0.0;
                    for k in 0..2 {
                        g += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, j);
                    }
                    ba.set(m, j, g);
                }
            }
            for m in 0..4 {
                for n in 0..4 {
                    let mut sum = 0.0;
                    for i in 0..2 {
                        for j in 0..2 {
                            sum += ba.get(m, i) * tt.get_std(i, j) * ba.get(n, j);
                        }
                    }
                    kk_wrong.add(m, n, c * sum);
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
    fn anisotropic_3d_tet4_over_tet10() {
        // Tet10 geometry built from the Tet4 coordinates with the midside nodes at the edge
        // midpoints (the quadratic map reproduces the linear one); the field is Tet4.
        #[rustfmt::skip]
        let xx = [
            [2.0, 3.0, 4.0], [6.0, 3.0, 2.0], [2.0, 5.0, 1.0], [4.0, 3.0, 6.0],
        ];
        let mut pad = Scratchpad::new_mixed(3, GeoKind::Tet4, GeoKind::Tet10).unwrap();
        let mut pad_geo = Scratchpad::new(3, GeoKind::Tet10).unwrap();
        for m in 0..4 {
            for j in 0..3 {
                pad.set_xx(m, j, xx[m][j]);
                pad_geo.set_xx(m, j, xx[m][j]);
            }
        }
        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
        for (i, (a, b)) in edges.iter().enumerate() {
            for j in 0..3 {
                let x = 0.5 * (xx[*a][j] + xx[*b][j]);
                pad.set_xx(4 + i, j, x);
                pad_geo.set_xx(4 + i, j, x);
            }
        }
        let gauss = Gauss::new(GeoKind::Tet10);

        // integration-point coordinates from the geometry (used by the tensor field)
        let x_ips = recovery::get_points_coords(&mut pad_geo, &gauss).unwrap();

        // reference: Ba = dNa/dξ · J⁻¹ and the Tet10 measure
        let mut kk_ref = Matrix::new(4, 4);
        let mut pad_calc = pad_geo.clone();
        let mut pad_field = Scratchpad::new(3, GeoKind::Tet4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_field.fn_deriv)(&mut pad_field.deriv, iota);
            let mut tt = Tensor2::<6>::new();
            tensor_field_3d(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            let mut ba = Matrix::new(4, 3);
            for m in 0..4 {
                for j in 0..3 {
                    let mut g = 0.0;
                    for k in 0..3 {
                        g += pad_field.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, j);
                    }
                    ba.set(m, j, g);
                }
            }
            for m in 0..4 {
                for n in 0..4 {
                    let mut sum = 0.0;
                    for i in 0..3 {
                        for j in 0..3 {
                            sum += ba.get(m, i) * tt.get_std(i, j) * ba.get(n, j);
                        }
                    }
                    kk_ref.add(m, n, c * sum);
                }
            }
        }

        // the single mixed pad must reproduce the reference
        let mut kk = Matrix::filled(4, 4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_03_btb::<6, _>(&mut kk, &mut args, |tt, p, _, _| {
            tensor_field_3d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);
    }

    #[test]
    fn offset_and_accumulation_work() {
        // constant anisotropic tensor used by both runs
        let tensor_field = |tt: &mut Tensor2<4>, _x: &Tensor1| {
            tt.sym_set_std(0, 0, 6.0);
            tt.sym_set_std(1, 1, 5.5);
            tt.sym_set_std(2, 2, 3.25);
            tt.sym_set_std(0, 1, 0.9);
        };
        let x = Tensor1::new();

        // single run used as the reference values
        let mut pad = aux::gen_pad_tri3();
        let gauss = Gauss::new(pad.kind);
        let mut kk_single = Matrix::filled(3, 3, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_03_btb::<4, _>(&mut kk_single, &mut args, |tt, _, _, _| {
            tensor_field(tt, &x);
            Ok(())
        })
        .unwrap();

        // two runs accumulating at (ii0, jj0) = (1, 1); the matrix is prefilled with 0.5 to prove
        // that clear = false leaves the border entries untouched
        let mut kk = Matrix::filled(4, 4, 0.5);
        for i in 1..4 {
            for j in 1..4 {
                kk.set(i, j, 0.0); // start the accumulated block from zero
            }
        }
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.clear = false;
        args.ii0 = 1;
        args.jj0 = 1;
        for _ in 0..2 {
            integ::mat_03_btb::<4, _>(&mut kk, &mut args, |tt, _, _, _| {
                tensor_field(tt, &x);
                Ok(())
            })
            .unwrap();
        }
        for i in 0..4 {
            for j in 0..4 {
                if i == 0 || j == 0 {
                    approx_eq(kk.get(i, j), 0.5, 1e-15);
                } else {
                    approx_eq(kk.get(i, j), 2.0 * kk_single.get(i - 1, j - 1), 1e-14);
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
                integ::mat_03_btb::<4, _>(&mut kk_base, &mut args_base, |tt, p, _, _| {
                    tensor_field_2d(tt, &x_ips[p]);
                    Ok(())
                })
                .unwrap();
                integ::mat_03_btb::<4, _>(&mut kk_mixed, &mut args_mixed, |tt, p, _, _| {
                    tensor_field_2d(tt, &x_ips[p]);
                    Ok(())
                })
                .unwrap();
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
            integ::mat_03_btb::<6, _>(&mut kk_base, &mut args_base, |tt, p, _, _| {
                tensor_field_3d(tt, &x_ips[p]);
                Ok(())
            })
            .unwrap();
            integ::mat_03_btb::<6, _>(&mut kk_mixed, &mut args_mixed, |tt, p, _, _| {
                tensor_field_3d(tt, &x_ips[p]);
                Ok(())
            })
            .unwrap();
            mat_approx_eq(&kk_mixed, &kk_base, 1e-14);
        }
    }
}
