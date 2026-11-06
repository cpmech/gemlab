use super::CommonArgs;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::math::SQRT_2;
use russell_lab::{mat_mat_mul, Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the gradient(Bb) times tensor(T) times shape(N) integration case 05 (e.g., coupling matrix)
///
/// Callback function: `f(T, p, Bb, N, B)`
///
/// **Notes:**
///
/// * `m` ranges over the nodes of the shape specified by `pad_b` (for `Bbᵐ`), which is the shape
///   being differentiated
/// * `n` ranges over the nodes of the *driver* shape specified by `pad` (for `Nⁿ`); the driver carries
///   the geometry and the integration measure
/// * For example, `1 ≤ m ≤ 4` for a `pad_b→Qua4` and `1 ≤ n ≤ 8` for `pad→Qua8`. The roles may be
///   reversed (`pad_b→Qua8`, `pad→Qua4`): the geometry is then the 4-node one and the quadratic shape
///   is differentiated in that map, as in a subparametric element
/// * The determinant of the Jacobian is calculated for `pad` (`pad` is the driver of the calculations)
/// * The gradient `Bb` is taken in the driver's mapping as well: `Bb = dNb/dξ · J⁻¹` with the `J` of
///   `pad`. Thus only the shape functions of `pad_b` are used (its coordinates are not read), of any
///   order, and the two pads describe the same element when both are built from the same cell
/// * The number of integration points must consider the nodes of `pad` and the expected order of the whole integrand
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠ →
/// Kᵐⁿ = │ Bbᵐ ⋅ T Nⁿ α dΩ
///       ⌡       ▔
///       Ωₑ
/// ```
///
/// The numerical integration is (assuming an implicit sum over repeated lower indices):
///
/// ```text
///        nip-1      →       →      →       →
/// Kᵐⁿⱼ ≈   Σ   Bbᵐᵢ(ιᵖ) Tᵢⱼ(ιᵖ) Nⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
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
/// * `pad_b` -- Scratchpad of the lower-order shape (modified) to compute `Nb` and `Bb`; its
///   coordinates are not used
/// * `args` --- Common arguments (`pad` is the Driver scratchpad (modified) to compute B)
/// * `fn_tt` -- Function `f(T,p,Bb,N,B)` that computes `T(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`,
///   the gradients Bb(ιᵖ) (in the driver's mapping), shape functions N(ιᵖ), and gradients B(ιᵖ).
///   `T` is set for `space_ndim`.
///
/// # Warning
///
/// The two [crate::shapes::Scratchpad]s must be functions on the same reference element, otherwise
/// **calculation errors will occur**: the same `ιᵖ` must refer to the same point of both. Since the
/// coordinates of `pad_b` are not used, the requirement is on its shape functions, not on its
/// coordinates, and `pad_b` may be of the lower or of the higher order than `pad` (for example
/// `pad_b→Qua4` with `pad→Qua8`, or the reverse).
pub fn mat_05_btn<const N: usize, F>(
    kk: &mut Matrix,
    pad_b: &mut Scratchpad,
    args: &mut CommonArgs,
    mut fn_tt: F,
) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor2<N>, usize, &Matrix, &Vector, &Matrix) -> Result<(), StrError>,
{
    // check
    let nnode_b = pad_b.interp.dim();
    let (space_ndim, _, nnode) = args.pad.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if pad_b.dims().1 != space_ndim {
        return Err("mat_05_btn requires pad_b with geo_ndim = space_ndim");
    }
    if nrow_kk < ii0 + nnode_b {
        return Err("nrow(K) must be ≥ ii0 + pad_b.nnode");
    }
    if ncol_kk < jj0 + nnode * space_ndim {
        return Err("ncol(K) must be ≥ jj0 + pad.nnode ⋅ space_ndim");
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

        // interpolation functions of both pads and the Jacobian of the driver
        (pad_b.fn_interp)(&mut pad_b.interp, iota); // Nb
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B and J of the driver

        // gradient of the lower-order shape in the driver's mapping: Bb = dNb/dξ · J⁻¹
        (pad_b.fn_deriv)(&mut pad_b.deriv, iota);
        mat_mat_mul(&mut pad_b.gradient, 1.0, &pad_b.deriv, &args.pad.inv_jacobian, 0.0).unwrap(); // cannot fail because the dimensions are checked above

        // calculate T tensor
        let ggb = &pad_b.gradient;
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        fn_tt(&mut tt, p, ggb, nn, bb)?;

        // calculate coefficient
        let c = if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt.get(0, m);
            }
            det_jac * weight * args.alpha * r
        } else {
            det_jac * weight * args.alpha
        };

        // add contribution to K matrix
        if space_ndim == 2 {
            for m in 0..nnode_b {
                for n in 0..nnode {
                    kk.add(
                        ii0 + m,
                        jj0 + 0 + n * 2,
                        c * (ggb.get(m, 0) * tt.get(0) + ggb.get(m, 1) * tt.get(3) / s) * nn[n],
                    );
                    kk.add(
                        ii0 + m,
                        jj0 + 1 + n * 2,
                        c * (ggb.get(m, 0) * tt.get(3) / s + ggb.get(m, 1) * tt.get(1)) * nn[n],
                    );
                }
            }
        } else {
            for m in 0..nnode_b {
                for n in 0..nnode {
                    kk.add(
                        ii0 + m,
                        jj0 + 0 + n * 3,
                        c * (ggb.get(m, 0) * tt.get(0) + ggb.get(m, 1) * tt.get(3) / s + ggb.get(m, 2) * tt.get(5) / s)
                            * nn[n],
                    );
                    kk.add(
                        ii0 + m,
                        jj0 + 1 + n * 3,
                        c * (ggb.get(m, 0) * tt.get(3) / s + ggb.get(m, 1) * tt.get(1) + ggb.get(m, 2) * tt.get(4) / s)
                            * nn[n],
                    );
                    kk.add(
                        ii0 + m,
                        jj0 + 2 + n * 3,
                        c * (ggb.get(m, 0) * tt.get(5) / s + ggb.get(m, 1) * tt.get(4) / s + ggb.get(m, 2) * tt.get(2))
                            * nn[n],
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
    use crate::integ::{self, AnalyticalQua8, AnalyticalTet4, CommonArgs, Gauss};
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::{mat_approx_eq, Matrix, Vector};
    use russell_tensor::Tensor2;

    #[test]
    fn capture_some_errors() {
        let (a, b) = (2.0, 3.0);
        let mut pad_b = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut kk = Matrix::new(4, 8 * 2);
        let mut tt = Tensor2::<4>::new();
        let ggb = Matrix::new(0, 0);
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let f = |_tt: &mut Tensor2<4>, _p: usize, _ggb: &Matrix, _nn: &Vector, _bb: &Matrix| Ok(());
        f(&mut tt, 0, &ggb, &nn, &bb).unwrap();
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_05_btn::<4, _>(&mut kk, &mut pad_b, &mut args, f).err(),
            Some("nrow(K) must be ≥ ii0 + pad_b.nnode")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_05_btn::<4, _>(&mut kk, &mut pad_b, &mut args, f).err(),
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
        #[rustfmt::skip]
        let tt = Tensor2::<4>::from_std_matrix(&[
            [2.0, 5.0, 0.0],
            [5.0, 3.0, 0.0],
            [0.0, 0.0, 4.0]]).unwrap();
        let kk_correct = ana.mat_05_btn(&tt);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [4, 9].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_05_btn::<4, _>(&mut kk, &mut pad_b, &mut args, |ten, _, _, _, _| {
                ten.set_tensor(1.0, &tt);
                Ok(())
            })
            .unwrap();
            // println!("{}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn tet4_tet4_works() {
        let mut pad_b = aux::gen_pad_tet4();
        let mut pad = pad_b.clone();
        let mut kk = Matrix::new(4, 4 * 3);
        let ana = AnalyticalTet4::new(&pad);
        #[rustfmt::skip]
        let tt = Tensor2::<6>::from_std_matrix(&[
            [1.0, 4.0, 6.0],
            [4.0, 2.0, 5.0],
            [6.0, 5.0, 3.0]]).unwrap();
        let kk_correct = ana.mat_05_btn(&tt);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14];
        let selection: Vec<_> = [4].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_05_btn::<6, _>(&mut kk, &mut pad_b, &mut args, |ten, _, _, _, _| {
                ten.set_tensor(1.0, &tt);
                Ok(())
            })
            .unwrap();
            // println!("{}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }
    #[test]
    fn curved_element_uses_the_driver_jacobian() {
        // Curved 8-node element: the midside nodes are pushed off the straight edges, so the
        // element's map and the bilinear (4-node) map are different domains. The reference below
        // is assembled here from dNb/dξ and the driver's J⁻¹, independently of mat_05_btn.
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
        let mut pad_b = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad_b.set_xx(m, 0, xx[m][0]);
            pad_b.set_xx(m, 1, xx[m][1]);
        }
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();
        #[rustfmt::skip]
        let tt = Tensor2::<4>::from_std_matrix(&[
            [2.0, 0.3, 0.0],
            [0.3, 1.0, 0.0],
            [0.0, 0.0, 0.0]]).unwrap();

        // kk_ref: Bb = dNb/dξ · J⁻¹ (driver's J). kk_own: Bb from the bilinear map (old behaviour).
        let mut kk_ref = Matrix::new(4, 16);
        let mut kk_own = Matrix::new(4, 16);
        let mut pad_calc = pad.clone();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            pad_b.calc_gradient(iota).unwrap(); // bilinear map: its own gradient
            let c = det_jac * gauss.weight(p);
            for m in 0..4 {
                for n in 0..8 {
                    for j in 0..2 {
                        let mut g_ref = 0.0;
                        let mut g_own = 0.0;
                        for k in 0..2 {
                            g_ref += (pad_b.deriv.get(m, 0) * pad_calc.inv_jacobian.get(0, k)
                                + pad_b.deriv.get(m, 1) * pad_calc.inv_jacobian.get(1, k))
                                * tt.get_std(k, j);
                            g_own += pad_b.gradient.get(m, k) * tt.get_std(k, j);
                        }
                        kk_ref.add(m, n * 2 + j, g_ref * pad_calc.interp[n] * c);
                        kk_own.add(m, n * 2 + j, g_own * pad_calc.interp[n] * c);
                    }
                }
            }
        }

        let mut kk = Matrix::new(4, 16);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn::<4, _>(&mut kk, &mut pad_b, &mut args, |ten, _, _, _, _| {
            ten.set_tensor(1.0, &tt);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the coordinates of pad_b are not used: the same matrix comes out of a pad without coordinates
        let mut pad_b_nocoords = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let mut kk_nocoords = Matrix::new(4, 16);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn::<4, _>(&mut kk_nocoords, &mut pad_b_nocoords, &mut args, |ten, _, _, _, _| {
            ten.set_tensor(1.0, &tt);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk_nocoords, &kk_ref, 1e-13);

        // on this element the old convention differs by an O(0.1) amount, so the check discriminates
        let mut max_diff = 0.0;
        for i in 0..4 {
            for j in 0..16 {
                max_diff = f64::max(max_diff, f64::abs(kk.get(i, j) - kk_own.get(i, j)));
            }
        }
        println!("max |K - K_own| = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }
    #[test]
    fn reversed_order_gives_the_subparametric_element() {
        // pad_b = Qua8 (differentiated), driver = Qua4 (geometry): the quadratic shape functions are
        // differentiated in the bilinear map. The reference is assembled here from dNb/dξ · J⁻¹.
        #[rustfmt::skip]
        let xx = [
            [0.00, 0.00], [2.00, 0.10], [1.70, 1.30], [-0.20, 1.00],
            [1.00, 0.35], [2.05, 0.72], [0.75, 1.45], [-0.10, 0.55],
        ];
        let mut pad_b = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for m in 0..8 {
            pad_b.set_xx(m, 0, xx[m][0]);
            pad_b.set_xx(m, 1, xx[m][1]);
        }
        let mut pad = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad.set_xx(m, 0, xx[m][0]);
            pad.set_xx(m, 1, xx[m][1]);
        }
        let gauss = Gauss::new_or_sized(GeoKind::Qua4, Some(9)).unwrap();
        #[rustfmt::skip]
        let tt = Tensor2::<4>::from_std_matrix(&[
            [2.0, 0.3, 0.0],
            [0.3, 1.0, 0.0],
            [0.0, 0.0, 0.0]]).unwrap();

        let mut kk_ref = Matrix::new(8, 4 * 2);
        let mut pad_calc = pad.clone();
        let mut pad_deriv = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_deriv.fn_deriv)(&mut pad_deriv.deriv, iota);
            let c = det_jac * gauss.weight(p);
            for m in 0..8 {
                for n in 0..4 {
                    for j in 0..2 {
                        let mut g = 0.0;
                        for k in 0..2 {
                            g += (pad_deriv.deriv.get(m, 0) * pad_calc.inv_jacobian.get(0, k)
                                + pad_deriv.deriv.get(m, 1) * pad_calc.inv_jacobian.get(1, k))
                                * tt.get_std(k, j);
                        }
                        kk_ref.add(m, n * 2 + j, g * pad_calc.interp[n] * c);
                    }
                }
            }
        }

        let mut kk = Matrix::new(8, 4 * 2);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn::<4, _>(&mut kk, &mut pad_b, &mut args, |ten, _, _, _, _| {
            ten.set_tensor(1.0, &tt);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);
    }
}
