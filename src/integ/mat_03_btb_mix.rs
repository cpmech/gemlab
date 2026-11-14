use super::CommonArgs;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::math::SQRT_2;
use russell_lab::{mat_mat_mul, Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the gradient(Ba) dot tensor(T) dot gradient(Ba) integration case 03 (mixed shapes)
///
/// Callback function: `f(T, p, N, B, Na, Ba)`
///
/// **Notes:**
///
/// * `pad` (in `args`) is the **main** pad: it carries the geometry, the integration measure `|J|`,
///   and its own shape functions `N` and gradients `B`
/// * `pad_alt` is the **alternative** pad: only its shape functions `Na` and derivatives `dNa/dξ`
///   are used, never its coordinates, and it may be of any order — lower, equal to, or higher than
///   that of `pad`
/// * `m` and `n` range over the nodes of `pad_alt` (for `Baᵐ` and `Baⁿ`), whose gradient is taken
///   in the **main** pad's mapping, `Ba = dNa/dξ · J⁻¹`
/// * For example, `0 ≤ m,n < 4` with a `pad_alt→Qua4` and `pad→Qua8`, or the reverse, which
///   differentiates the quadratic shape in the 4-node map (a subparametric element)
/// * The number of integration points must consider the nodes of both pads and the expected order
///   of the whole integrand
///
/// Conductivity coefficients:
///
/// ```text
///       ⌠ →        →
/// Kᵐⁿ = │ Baᵐ ⋅ T ⋅ Baⁿ α dΩ
///       ⌡       ▔
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///       nip-1 →   →       →     →   →       →
/// Kᵐⁿ ≈   Σ   Baᵐ(ιᵖ) ⋅ T(ιᵖ) ⋅ Baⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
///        p=0            ▔
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
/// * `pad_alt` -- Scratchpad of the alternative shape (modified) to compute `Na` and `Ba`; its
///   coordinates are not used
/// * `fn_tt` -- Function `f(T,p,N,B,Na,Ba)` that computes `T(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`,
///   shape functions `N(ιᵖ)`, gradients `B(ιᵖ)`, shape functions `Na(ιᵖ)`, and gradients `Ba(ιᵖ)`
///   (in the main pad's mapping). `T` is set in the symmetric Mandel representation with
///   `N = 2 * space_ndim` components.
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
///     integ::mat_03_btb_mix::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, _, _, _, _, _| {
///         tt.sym_set_std(0, 0, 2.0);
///         tt.sym_set_std(1, 1, 3.0);
///         Ok(())
///     })?;
///     // solution (T = diag(2, 3), A = 6):
///     // Kᵐⁿ = A Bᵐᵀ T Bⁿ
///     #[rustfmt::skip]
///     let kk_correct = Matrix::from(&[
///         [2.75, -0.75, -2.0],
///         [-0.75, 0.75, 0.0],
///         [-2.0, 0.0, 2.0],
///     ]);
///     mat_approx_eq(&kk, &kk_correct, 1e-14);
///     Ok(())
/// }
/// ```
pub fn mat_03_btb_mix<const N: usize, F>(
    kk: &mut Matrix,
    args: &mut CommonArgs,
    pad_alt: &mut Scratchpad,
    mut fn_tt: F,
) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor2<N>, usize, &Vector, &Matrix, &Vector, &Matrix) -> Result<(), StrError>,
{
    // check
    let nnode_alt = pad_alt.interp.dim();
    let (space_ndim, _, nnode) = args.pad.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if pad_alt.kind.class() != args.pad.kind.class() || pad_alt.kind.ndim() != args.pad.kind.ndim() {
        return Err("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)");
    }
    if pad_alt.dims().1 != space_ndim {
        return Err("mat_03_btb_mix requires pad_alt with geo_ndim = space_ndim");
    }
    if nrow_kk < ii0 + nnode_alt {
        return Err("nrow(K) must be ≥ ii0 + pad_alt.nnode");
    }
    if ncol_kk < jj0 + nnode_alt {
        return Err("ncol(K) must be ≥ jj0 + pad_alt.nnode");
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

        // interpolation functions of both pads and the Jacobian of the main pad
        (pad_alt.fn_interp)(&mut pad_alt.interp, iota); // Na
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B and J of the main pad

        // gradient of the alternative shape in the main pad's mapping: Ba = dNa/dξ · J⁻¹
        (pad_alt.fn_deriv)(&mut pad_alt.deriv, iota);
        mat_mat_mul(&mut pad_alt.gradient, 1.0, &pad_alt.deriv, &args.pad.inv_jacobian, 0.0).unwrap(); // cannot fail because the dimensions are checked above

        // calculate T tensor
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        let na = &pad_alt.interp;
        let ba = &pad_alt.gradient;
        fn_tt(&mut tt, p, nn, bb, na, ba)?;

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
            for m in 0..nnode_alt {
                for n in 0..nnode_alt {
                    kk.add(
                        ii0 + m,
                        jj0 + n,
                        c * (ba.get(n, 1) * (tt.get(1) * ba.get(m, 1) + (tt.get(3) * ba.get(m, 0)) / s)
                            + ba.get(n, 0) * (tt.get(0) * ba.get(m, 0) + (tt.get(3) * ba.get(m, 1)) / s)),
                    );
                }
            }
        } else {
            for m in 0..nnode_alt {
                for n in 0..nnode_alt {
                    kk.add(
                        ii0 + m,
                        jj0 + n,
                        c * (ba.get(n, 2)
                            * (tt.get(2) * ba.get(m, 2)
                                + (tt.get(5) * ba.get(m, 0)) / s
                                + (tt.get(4) * ba.get(m, 1)) / s)
                            + ba.get(n, 1)
                                * (tt.get(1) * ba.get(m, 1)
                                    + (tt.get(3) * ba.get(m, 0)) / s
                                    + (tt.get(4) * ba.get(m, 2)) / s)
                            + ba.get(n, 0)
                                * (tt.get(0) * ba.get(m, 0)
                                    + (tt.get(3) * ba.get(m, 1)) / s
                                    + (tt.get(5) * ba.get(m, 2)) / s)),
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
    use crate::integ::{self, CommonArgs, Gauss};
    use crate::recovery;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::{approx_eq, mat_approx_eq, Matrix};
    use russell_tensor::{Tensor1, Tensor2};

    /// Compares mat_03_btb_mix with mat_03_btb using two identically populated pads
    ///
    /// The same `tensor_field` closure is called from both integration functions.
    fn check_equal_order<const N: usize>(
        pad: &Scratchpad,
        axisymmetric: bool,
        ngauss: usize,
        tensor_field: impl Fn(&mut Tensor2<N>, &Tensor1),
    ) {
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

        // the tensor field is evaluated at the integration points of the main (geometric) pad
        let x_ips = recovery::get_points_coords(args_base.pad, &gauss).unwrap();
        integ::mat_03_btb::<N, _>(&mut kk_base, &mut args_base, |tt, p, _, _| {
            tensor_field(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        integ::mat_03_btb_mix::<N, _>(&mut kk_mix, &mut args_mix, &mut pad_alt, |tt, p, _, _, _, _| {
            tensor_field(tt, &x_ips[p]);
            Ok(())
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
            integ::mat_03_btb_mix::<4, _>(&mut kk, &mut args, &mut pad_alt, |_, _, _, _, _, _| Ok(())).err(),
            Some("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)")
        );

        // pad_alt.dims().1 != space_ndim: a line (geo_ndim = 1) in a 2D space
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut pad_alt = Scratchpad::new(2, GeoKind::Lin3).unwrap();
        let mut kk = Matrix::new(2, 2);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::mat_03_btb_mix::<4, _>(&mut kk, &mut args, &mut pad_alt, |_, _, _, _, _, _| Ok(())).err(),
            Some("mat_03_btb_mix requires pad_alt with geo_ndim = space_ndim")
        );

        // nrow(K) < ii0 + pad_alt.nnode
        let mut pad = aux::gen_pad_tri3();
        let mut pad_alt = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_03_btb_mix::<4, _>(&mut kk, &mut args, &mut pad_alt, |_, _, _, _, _, _| Ok(())).err(),
            Some("nrow(K) must be ≥ ii0 + pad_alt.nnode")
        );

        // ncol(K) < jj0 + pad_alt.nnode
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_03_btb_mix::<4, _>(&mut kk, &mut args, &mut pad_alt, |_, _, _, _, _, _| Ok(())).err(),
            Some("ncol(K) must be ≥ jj0 + pad_alt.nnode")
        );

        // N != 2 * space_ndim
        args.jj0 = 0;
        assert_eq!(
            integ::mat_03_btb_mix::<6, _>(&mut kk, &mut args, &mut pad_alt, |_, _, _, _, _, _| Ok(())).err(),
            Some("tensor dimension N must equal 2 * space_ndim")
        );

        // the callback error is propagated
        assert_eq!(
            integ::mat_03_btb_mix::<4, _>(&mut kk, &mut args, &mut pad_alt, |_, _, _, _, _, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn equal_order_matches_mat_03_btb() {
        // spatially varying full anisotropic tensor in 2D (symmetric Mandel components)
        let tensor_field_2d = |tt: &mut Tensor2<4>, x: &Tensor1| {
            tt.sym_set_std(0, 0, 2.0 + x.get(0));
            tt.sym_set_std(1, 1, 3.0 + 0.5 * x.get(1));
            tt.sym_set_std(2, 2, 1.0 + 0.25 * (x.get(0) + x.get(1)));
            tt.sym_set_std(0, 1, 0.5 + 0.1 * x.get(0));
        };

        // 2D: Tri3 + Tri3
        let pad = aux::gen_pad_tri3();
        for ngauss in [1, 3, 7] {
            check_equal_order(&pad, false, ngauss, tensor_field_2d);
        }

        // axisymmetric: the radius is r > 0 on the chosen triangle and T₂₂ != 0
        for ngauss in [1, 3, 7] {
            check_equal_order(&pad, true, ngauss, tensor_field_2d);
        }

        // spatially varying full anisotropic tensor in 3D (symmetric Mandel components)
        let tensor_field_3d = |tt: &mut Tensor2<6>, x: &Tensor1| {
            tt.sym_set_std(0, 0, 2.0 + x.get(0));
            tt.sym_set_std(1, 1, 3.0 + 0.5 * x.get(1));
            tt.sym_set_std(2, 2, 1.0 + 0.25 * x.get(2));
            tt.sym_set_std(0, 1, 0.5 + 0.1 * x.get(0));
            tt.sym_set_std(1, 2, -0.3 + 0.2 * x.get(1));
            tt.sym_set_std(2, 0, 0.4 - 0.15 * x.get(2));
        };

        // 3D: Tet4 + Tet4
        let pad = aux::gen_pad_tet4();
        for ngauss in [1, 4] {
            check_equal_order(&pad, false, ngauss, tensor_field_3d);
        }
    }

    #[test]
    fn curved_element_uses_the_main_jacobian() {
        // Curved 8-node element: the midside nodes are pushed off the straight edges, so the
        // element's map and the bilinear (4-node) map are different domains. The reference below
        // is assembled here from dNa/dξ · J⁻¹ and the main pad's det(J), independently of
        // mat_03_btb_mix.
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

        // coordinates of the integration points from the main shape (used by the tensor field)
        let mut pad_x = pad.clone();
        let x_ips = recovery::get_points_coords(&mut pad_x, &gauss).unwrap();

        // spatially varying full anisotropic tensor, shared by the reference loop and the callback
        let tensor_field = |tt: &mut Tensor2<4>, x: &Tensor1| {
            tt.sym_set_std(0, 0, 2.0 + x.get(0));
            tt.sym_set_std(1, 1, 3.0 + 0.5 * x.get(1));
            tt.sym_set_std(2, 2, 1.0 + 0.25 * (x.get(0) + x.get(1)));
            tt.sym_set_std(0, 1, 0.5 + 0.1 * x.get(0));
        };

        // kk_ref: assembled from dNa/dξ · J⁻¹ and the main pad's det(J)
        let mut kk_ref = Matrix::new(4, 4);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_alt_calc.fn_deriv)(&mut pad_alt_calc.deriv, iota);
            let mut tt = Tensor2::<4>::new();
            tensor_field(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            let mut ba = Matrix::new(4, 2);
            for m in 0..4 {
                for j in 0..2 {
                    let mut g = 0.0;
                    for k in 0..2 {
                        g += pad_alt_calc.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, j);
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

        // the mix function must match the reference
        let mut kk = Matrix::filled(4, 4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_03_btb_mix::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _, _| {
            tensor_field(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the coordinates of pad_alt are not used: a pad without coordinates gives the same matrix
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let mut kk_nocoords = Matrix::filled(4, 4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_03_btb_mix::<4, _>(
            &mut kk_nocoords,
            &mut args,
            &mut pad_alt_nocoords,
            |tt, p, _, _, _, _| {
                tensor_field(tt, &x_ips[p]);
                Ok(())
            },
        )
        .unwrap();
        mat_approx_eq(&kk_nocoords, &kk, 1e-15);

        // the wrong variant uses Ba and detJ from the alternative pad's own bilinear map
        let mut kk_wrong = Matrix::new(4, 4);
        let mut pad_alt_wrong = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad_alt_wrong.set_xx(m, 0, xx[m][0]);
            pad_alt_wrong.set_xx(m, 1, xx[m][1]);
        }
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_alt_wrong.calc_gradient(iota).unwrap();
            let mut tt = Tensor2::<4>::new();
            tensor_field(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            for m in 0..4 {
                for n in 0..4 {
                    let mut sum = 0.0;
                    for i in 0..2 {
                        for j in 0..2 {
                            sum +=
                                pad_alt_wrong.gradient.get(m, i) * tt.get_std(i, j) * pad_alt_wrong.gradient.get(n, j);
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
        println!("max |K - K_wrong| = {:.3e}", max_diff);
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

        // coordinates of the integration points from the main shape (used by the tensor field)
        let mut pad_x = pad.clone();
        let x_ips = recovery::get_points_coords(&mut pad_x, &gauss).unwrap();

        // spatially varying full anisotropic tensor, shared by the reference loop and the callback
        let tensor_field = |tt: &mut Tensor2<4>, x: &Tensor1| {
            tt.sym_set_std(0, 0, 2.0 + x.get(0));
            tt.sym_set_std(1, 1, 3.0 + 0.5 * x.get(1));
            tt.sym_set_std(2, 2, 1.0 + 0.25 * (x.get(0) + x.get(1)));
            tt.sym_set_std(0, 1, 0.5 + 0.1 * x.get(0));
        };

        // kk_ref: assembled from dNa/dξ · J⁻¹ and the main pad's det(J)
        let mut kk_ref = Matrix::new(8, 8);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_alt_calc.fn_deriv)(&mut pad_alt_calc.deriv, iota);
            let mut tt = Tensor2::<4>::new();
            tensor_field(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            let mut ba = Matrix::new(8, 2);
            for m in 0..8 {
                for j in 0..2 {
                    let mut g = 0.0;
                    for k in 0..2 {
                        g += pad_alt_calc.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, j);
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

        // the mix function must match the reference
        let mut kk = Matrix::filled(8, 8, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_03_btb_mix::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _, _| {
            tensor_field(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the coordinates of pad_alt are not used: a pad without coordinates gives the same matrix
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        let mut kk_nocoords = Matrix::filled(8, 8, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_03_btb_mix::<4, _>(
            &mut kk_nocoords,
            &mut args,
            &mut pad_alt_nocoords,
            |tt, p, _, _, _, _| {
                tensor_field(tt, &x_ips[p]);
                Ok(())
            },
        )
        .unwrap();
        mat_approx_eq(&kk_nocoords, &kk, 1e-15);

        // the wrong variant uses Ba and detJ from the alternative pad's own quadratic map
        let mut kk_wrong = Matrix::new(8, 8);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_alt.calc_gradient(iota).unwrap();
            let mut tt = Tensor2::<4>::new();
            tensor_field(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            for m in 0..8 {
                for n in 0..8 {
                    let mut sum = 0.0;
                    for i in 0..2 {
                        for j in 0..2 {
                            sum += pad_alt.gradient.get(m, i) * tt.get_std(i, j) * pad_alt.gradient.get(n, j);
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
        println!("max |K - K_wrong| = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn anisotropic_3d_tet10_tet4() {
        // Tet10 main pad built from the Tet4 coordinates with the midside nodes at the edge
        // midpoints (the quadratic map reproduces the linear one). The alternative has only 4 nodes.
        #[rustfmt::skip]
        let xx = [
            [2.0, 3.0, 4.0], [6.0, 3.0, 2.0], [2.0, 5.0, 1.0], [4.0, 3.0, 6.0],
        ];
        let mut pad = Scratchpad::new(3, GeoKind::Tet10).unwrap();
        for m in 0..4 {
            for j in 0..3 {
                pad.set_xx(m, j, xx[m][j]);
            }
        }
        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
        for (i, (a, b)) in edges.iter().enumerate() {
            for j in 0..3 {
                pad.set_xx(4 + i, j, 0.5 * (xx[*a][j] + xx[*b][j]));
            }
        }
        let mut pad_alt = aux::gen_pad_tet4();
        let gauss = Gauss::new(pad.kind);

        // coordinates of the integration points from the main shape (used by the tensor field)
        let mut pad_x = pad.clone();
        let x_ips = recovery::get_points_coords(&mut pad_x, &gauss).unwrap();

        // spatially varying full anisotropic tensor with T₁₂, T₂₀ != 0, shared with the reference loop
        let tensor_field = |tt: &mut Tensor2<6>, x: &Tensor1| {
            tt.sym_set_std(0, 0, 2.0 + x.get(0));
            tt.sym_set_std(1, 1, 3.0 + 0.5 * x.get(1));
            tt.sym_set_std(2, 2, 1.0 + 0.25 * x.get(2));
            tt.sym_set_std(0, 1, 0.5 + 0.1 * x.get(0));
            tt.sym_set_std(1, 2, -0.3 + 0.2 * x.get(1));
            tt.sym_set_std(2, 0, 0.4 - 0.15 * x.get(2));
        };

        // kk_ref: assembled from dNa/dξ · J⁻¹ and the main pad's det(J)
        let mut kk_ref = Matrix::new(10, 4);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(3, GeoKind::Tet4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_alt_calc.fn_deriv)(&mut pad_alt_calc.deriv, iota);
            let mut tt = Tensor2::<6>::new();
            tensor_field(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            let mut ba = Matrix::new(4, 3);
            for m in 0..4 {
                for j in 0..3 {
                    let mut g = 0.0;
                    for k in 0..3 {
                        g += pad_alt_calc.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, j);
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

        // the mix function must match the reference
        let mut kk = Matrix::filled(10, 4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_03_btb_mix::<6, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _, _| {
            tensor_field(tt, &x_ips[p]);
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
        let mut pad_alt = aux::gen_pad_tri3();
        let gauss = Gauss::new(pad.kind);
        let mut kk_single = Matrix::filled(3, 3, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_03_btb_mix::<4, _>(&mut kk_single, &mut args, &mut pad_alt, |tt, _, _, _, _, _| {
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
            integ::mat_03_btb_mix::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, _, _, _, _, _| {
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
}
