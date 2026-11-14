use super::CommonArgs;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{mat_mat_mul, Matrix, Vector};
use russell_tensor::Tensor1;

/// Implements the gradient(Ba) dot vector(V) integration case 03 (mixed shapes)
///
/// Callback function: `f(w, p, N, B, Na, Ba)`
///
/// **Notes:**
///
/// * `pad` (in `args`) is the **main** pad: it carries the geometry, the integration measure `|J|`,
///   and its own shape functions `N` and gradients `B`
/// * `pad_alt` is the **alternative** pad: only its shape functions `Na` and derivatives `dNa/dξ` are
///   used, never its coordinates, and it may be of any order — lower, equal to, or higher than that
///   of `pad`
/// * `m` ranges over the nodes of `pad_alt` (for `Baᵐ`), whose gradient is taken in the **main** pad's
///   mapping, `Ba = dNa/dξ · J⁻¹`
/// * For example, `0 ≤ m < 4` with a `pad_alt→Qua4` and `pad→Qua8`, or the reverse, which
///   differentiates the quadratic shape in the 4-node map (a subparametric element)
/// * The number of integration points must consider the nodes of both pads and the expected order of
///   the whole integrand
///
/// Vector dot gradient:
///
/// ```text
///      ⌠ →  → →     → →
/// cᵐ = │ Baᵐ(x(ξ)) · w(x) α dΩ
///      ⌡
///      Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///      nip-1  →  →     → →       →
/// cᵐ ≈   Σ    Baᵐ(ιᵖ) · w(ιᵖ) |J|(ιᵖ) wᵖ α
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
/// * `c` -- A vector containing all `cᵐ` values, one after another, and sequentially placed as shown
///   above. `m` is the index of the node of `pad_alt`. The length must be
///   `c.len() ≥ ii0 + pad_alt.nnode`
/// * `args` --- Common arguments (`pad` is the **main** scratchpad (modified), which carries the
///   geometry)
/// * `pad_alt` -- Scratchpad of the alternative shape (modified) to compute `Na` and `Ba`; its
///   coordinates are not used
/// * `fn_w` -- Function `f(w,p,N,B,Na,Ba)` that computes `w(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`,
///   shape functions `N(ιᵖ)`, gradients `B(ιᵖ)`, shape functions `Na(ιᵖ)`, and gradients `Ba(ιᵖ)` (in
///   the main pad's mapping). `w` has 3 components (only `space_ndim` are used).
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
///     let mut c = Vector::filled(pad.nnode(), 0.0);
///     let mut args = CommonArgs::new(&mut pad, &gauss);
///     integ::vec_03_bv_mix(&mut c, &mut args, &mut pad_alt, |w, _, _, _, _, _| {
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
pub fn vec_03_bv_mix<F>(
    c: &mut Vector,
    args: &mut CommonArgs,
    pad_alt: &mut Scratchpad,
    mut fn_w: F,
) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor1, usize, &Vector, &Matrix, &Vector, &Matrix) -> Result<(), StrError>,
{
    // check
    let nnode_alt = pad_alt.interp.dim();
    let (space_ndim, _, nnode) = args.pad.dims();
    let ii0 = args.ii0;
    if pad_alt.kind.class() != args.pad.kind.class() || pad_alt.kind.ndim() != args.pad.kind.ndim() {
        return Err("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)");
    }
    if pad_alt.dims().1 != space_ndim {
        return Err("vec_03_bv_mix requires pad_alt with geo_ndim = space_ndim");
    }
    if c.dim() < ii0 + nnode_alt {
        return Err("c.len() must be ≥ ii0 + pad_alt.nnode");
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

        // interpolation functions of both pads and the Jacobian of the main pad
        (pad_alt.fn_interp)(&mut pad_alt.interp, iota); // Na
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B and J of the main pad

        // gradient of the alternative shape in the main pad's mapping: Ba = dNa/dξ · J⁻¹
        (pad_alt.fn_deriv)(&mut pad_alt.deriv, iota);
        mat_mat_mul(&mut pad_alt.gradient, 1.0, &pad_alt.deriv, &args.pad.inv_jacobian, 0.0).unwrap(); // cannot fail because the dimensions are checked above

        // calculate w
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        let na = &pad_alt.interp;
        let ba = &pad_alt.gradient;
        fn_w(&mut w, p, nn, bb, na, ba)?;

        // calculate coefficient
        let coef = if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt.get(0, m);
            }
            det_jac * weight * args.alpha * r
        } else {
            det_jac * weight * args.alpha
        };

        // add contribution to c vector
        if space_ndim == 2 {
            for m in 0..nnode_alt {
                c[ii0 + m] += coef * (w.get(0) * ba.get(m, 0) + w.get(1) * ba.get(m, 1));
            }
        } else {
            for m in 0..nnode_alt {
                c[ii0 + m] += coef * (w.get(0) * ba.get(m, 0) + w.get(1) * ba.get(m, 1) + w.get(2) * ba.get(m, 2));
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
    use russell_lab::{approx_eq, vec_approx_eq, Vector};
    use russell_tensor::Tensor1;

    /// Spatially varying vector field: w(x) = {2 + x₀, -1 + 3 x₁, 0.5 - x₂}
    fn vector_field(w: &mut Tensor1, x: &Tensor1) {
        w.set(0, 2.0 + x.get(0));
        w.set(1, -1.0 + 3.0 * x.get(1));
        w.set(2, 0.5 - x.get(2));
    }

    /// Compares vec_03_bv_mix with vec_03_bv using two identically populated pads
    fn check_equal_order(pad: &Scratchpad, axisymmetric: bool, ngauss: usize) {
        let mut pad_mix = pad.clone();
        let mut pad_base = pad.clone();
        let mut pad_alt = pad.clone();
        let gauss = Gauss::new_or_sized(pad.kind, Some(ngauss)).unwrap();
        let mut c_mix = Vector::filled(pad.nnode(), aux::NOISE);
        let mut c_base = Vector::filled(pad.nnode(), aux::NOISE);
        let mut args_mix = CommonArgs::new(&mut pad_mix, &gauss);
        let mut args_base = CommonArgs::new(&mut pad_base, &gauss);
        args_mix.axisymmetric = axisymmetric;
        args_base.axisymmetric = axisymmetric;

        // the vector field is evaluated at the integration points of the main (geometric) pad
        let x_ips = recovery::get_points_coords(args_base.pad, &gauss).unwrap();
        integ::vec_03_bv(&mut c_base, &mut args_base, |w, p, _, _| {
            vector_field(w, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        integ::vec_03_bv_mix(&mut c_mix, &mut args_mix, &mut pad_alt, |w, p, _, _, _, _| {
            vector_field(w, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        vec_approx_eq(&c_mix, &c_base, 1e-14);
    }

    #[test]
    fn capture_some_errors() {
        // class/geo_ndim mismatch: main Qua8 (2D) with alt Tri3 (2D)
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, 2.0, 3.0);
        let mut pad_alt = aux::gen_pad_tri3();
        let mut c = Vector::filled(8, aux::NOISE);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::vec_03_bv_mix(&mut c, &mut args, &mut pad_alt, |_, _, _, _, _, _| Ok(())).err(),
            Some("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)")
        );

        // pad_alt.dims().1 != space_ndim: a line (geo_ndim = 1) in a 2D space
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut pad_alt = Scratchpad::new(2, GeoKind::Lin3).unwrap();
        let mut c = Vector::filled(2, aux::NOISE);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::vec_03_bv_mix(&mut c, &mut args, &mut pad_alt, |_, _, _, _, _, _| Ok(())).err(),
            Some("vec_03_bv_mix requires pad_alt with geo_ndim = space_ndim")
        );

        // c.dim() < ii0 + nnode_alt
        let mut pad = aux::gen_pad_tri3();
        let mut pad_alt = aux::gen_pad_tri3();
        let mut c = Vector::new(3);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::vec_03_bv_mix(&mut c, &mut args, &mut pad_alt, |_, _, _, _, _, _| Ok(())).err(),
            Some("c.len() must be ≥ ii0 + pad_alt.nnode")
        );

        // the callback error is propagated
        args.ii0 = 0;
        assert_eq!(
            integ::vec_03_bv_mix(&mut c, &mut args, &mut pad_alt, |_, _, _, _, _, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn equal_order_matches_vec_03_bv() {
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
        // assembled here from dNa/dξ · J⁻¹ and the main pad's det(J), independently of vec_03_bv_mix.
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

        // coordinates of the integration points from the main shape (used by the vector field)
        let mut pad_x = pad.clone();
        let x_ips = recovery::get_points_coords(&mut pad_x, &gauss).unwrap();

        // c_ref: assembled from dNa/dξ · J⁻¹ and the main pad's det(J)
        let mut c_ref = Vector::new(4);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_alt_calc.fn_deriv)(&mut pad_alt_calc.deriv, iota);
            let (w0, w1) = (2.0 + x_ips[p].get(0), -1.0 + 3.0 * x_ips[p].get(1));
            let coef = det_jac * gauss.weight(p);
            for m in 0..4 {
                let mut ba0 = 0.0;
                let mut ba1 = 0.0;
                for k in 0..2 {
                    ba0 += pad_alt_calc.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, 0);
                    ba1 += pad_alt_calc.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, 1);
                }
                c_ref[m] += coef * (w0 * ba0 + w1 * ba1);
            }
        }

        // the mix function must match the reference
        let mut c = Vector::filled(4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_03_bv_mix(&mut c, &mut args, &mut pad_alt, |w, p, _, _, _, _| {
            vector_field(w, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        vec_approx_eq(&c, &c_ref, 1e-13);

        // the coordinates of pad_alt are not used: a pad without coordinates gives the same vector
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let mut c_nocoords = Vector::filled(4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_03_bv_mix(&mut c_nocoords, &mut args, &mut pad_alt_nocoords, |w, p, _, _, _, _| {
            vector_field(w, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        vec_approx_eq(&c_nocoords, &c, 1e-15);

        // the wrong variant uses Ba and detJ from the alternative pad's own bilinear map
        let mut c_wrong = Vector::new(4);
        let mut pad_alt_wrong = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad_alt_wrong.set_xx(m, 0, xx[m][0]);
            pad_alt_wrong.set_xx(m, 1, xx[m][1]);
        }
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_alt_wrong.calc_gradient(iota).unwrap();
            let (w0, w1) = (2.0 + x_ips[p].get(0), -1.0 + 3.0 * x_ips[p].get(1));
            let coef = det_jac * gauss.weight(p);
            for m in 0..4 {
                c_wrong[m] += coef * (w0 * pad_alt_wrong.gradient.get(m, 0) + w1 * pad_alt_wrong.gradient.get(m, 1));
            }
        }
        let mut max_diff = 0.0;
        for m in 0..4 {
            max_diff = f64::max(max_diff, f64::abs(c[m] - c_wrong[m]));
        }
        println!("max |c - c_wrong| = {:.3e}", max_diff);
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

        // coordinates of the integration points from the main shape (used by the vector field)
        let mut pad_x = pad.clone();
        let x_ips = recovery::get_points_coords(&mut pad_x, &gauss).unwrap();

        // c_ref: assembled from dNa/dξ · J⁻¹ and the main pad's det(J)
        let mut c_ref = Vector::new(8);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_alt_calc.fn_deriv)(&mut pad_alt_calc.deriv, iota);
            let (w0, w1) = (2.0 + x_ips[p].get(0), -1.0 + 3.0 * x_ips[p].get(1));
            let coef = det_jac * gauss.weight(p);
            for m in 0..8 {
                let mut ba0 = 0.0;
                let mut ba1 = 0.0;
                for k in 0..2 {
                    ba0 += pad_alt_calc.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, 0);
                    ba1 += pad_alt_calc.deriv.get(m, k) * pad_calc.inv_jacobian.get(k, 1);
                }
                c_ref[m] += coef * (w0 * ba0 + w1 * ba1);
            }
        }

        // the mix function must match the reference
        let mut c = Vector::filled(8, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_03_bv_mix(&mut c, &mut args, &mut pad_alt, |w, p, _, _, _, _| {
            vector_field(w, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        vec_approx_eq(&c, &c_ref, 1e-13);

        // the coordinates of pad_alt are not used: a pad without coordinates gives the same vector
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        let mut c_nocoords = Vector::filled(8, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_03_bv_mix(&mut c_nocoords, &mut args, &mut pad_alt_nocoords, |w, p, _, _, _, _| {
            vector_field(w, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        vec_approx_eq(&c_nocoords, &c, 1e-15);

        // the wrong variant uses Ba and detJ from the alternative pad's own quadratic map
        let mut c_wrong = Vector::new(8);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_alt.calc_gradient(iota).unwrap();
            let (w0, w1) = (2.0 + x_ips[p].get(0), -1.0 + 3.0 * x_ips[p].get(1));
            let coef = det_jac * gauss.weight(p);
            for m in 0..8 {
                c_wrong[m] += coef * (w0 * pad_alt.gradient.get(m, 0) + w1 * pad_alt.gradient.get(m, 1));
            }
        }
        let mut max_diff = 0.0;
        for m in 0..8 {
            max_diff = f64::max(max_diff, f64::abs(c[m] - c_wrong[m]));
        }
        println!("max |c - c_wrong| = {:.3e}", max_diff);
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
        let mut pad_alt = pad.clone();
        let gauss = Gauss::new(pad.kind);
        let mut c_single = Vector::filled(4, aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::vec_03_bv_mix(&mut c_single, &mut args, &mut pad_alt, |w, _, _, _, _, _| {
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
            integ::vec_03_bv_mix(&mut c, &mut args, &mut pad_alt, |w, _, _, _, _, _| {
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
}
