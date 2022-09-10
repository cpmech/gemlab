use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::Vector;

/// Implements the the shape(N) times vector(V) integration case 02 (boundary integral version)
///
/// Callback function: `α ← f(v, p, un, N)`
///
/// Interpolation functions times vector field:
///
/// ```text
/// →    ⌠    → →   → →
/// bᵐ = │ Nᵐ(x(ξ)) v(x) α dΓ
///      ⌡
///      Γₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
/// →    nip-1     →   → →     →   →
/// bᵐ ≈   Σ    Nᵐ(ιᵖ) v(ιᵖ) ||n||(ιᵖ) wᵖ α
///       p=0
/// ```
///
/// # Output
///
/// ```text
///     ┌     ┐
///     | b⁰₀ |  ⟸  ii0 = 0
///     | b⁰₁ |
///     | b¹₀ |
/// b = | b¹₁ |
///     | b²₀ |
///     | b²₁ |
///     | ··· |
///     | bᵐᵢ |  ⟸  ii := i + m ⋅ space_ndim
///     └     ┘       
///
/// m = ii / space_ndim
/// i = ii % space_ndim
/// ```
///
/// * `b` -- A vector containing all `bᵐᵢ` values, one after another, and sequentially placed
///   as shown above (in 2D). `m` is the index of the node and `i` corresponds to `space_ndim`.
///   The length must be `b.len() ≥ ii0 + nnode ⋅ space_ndim`
/// * `pad` -- Some members of the scratchpad will be modified
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output vector where to add components
/// * `clear_b` -- fills `b` vector with zeros, otherwise accumulate values into `b`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_v` -- Function `f(v,p,un,N)→α` that calculates `v(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   the **unit** normal vector `un(x(ιᵖ))`, and shape functions N(ιᵖ).
///   `v.dim() = space_ndim` and `un.dim() = space_ndim`.
///   `fn_v` returns α that can accommodate plane-strain simulations.
///   **NOTE:** the value α is ignored if axisymmetric = true, because the radius is calculated and used instead.
///
/// # Examples
///
pub fn vec_02_nv_bry<F>(
    b: &mut Vector,
    pad: &mut Scratchpad,
    ii0: usize,
    clear_b: bool,
    axisymmetric: bool,
    ips: IntegPointData,
    mut fn_v: F,
) -> Result<(), StrError>
where
    F: FnMut(&mut Vector, usize, &Vector, &Vector) -> Result<f64, StrError>,
{
    // check
    let (space_ndim, nnode) = pad.xxt.dims();
    let geo_ndim = pad.deriv.dims().1;
    if space_ndim == 2 && geo_ndim != 1 {
        return Err("in 2D, geometry ndim must be equal to 1 (a line)");
    }
    if space_ndim == 3 && geo_ndim != 2 {
        return Err("in 3D, geometry ndim must be equal to 2 (a surface)");
    }
    if b.dim() < ii0 + nnode * space_ndim {
        return Err("b.len() must be ≥ ii0 + nnode ⋅ space_ndim");
    }

    // allocate auxiliary vectors
    let mut v = Vector::new(space_ndim);
    let mut un = Vector::new(space_ndim); // unit normal vector

    // clear output vector
    if clear_b {
        b.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate interpolation functions and unit normal vector
        (pad.fn_interp)(&mut pad.interp, iota); // N
        let mag_n = pad.calc_normal_vector(&mut un, iota)?;

        // calculate t
        let nn = &pad.interp;
        let alpha = fn_v(&mut v, p, &un, nn)?;

        // calculate coefficient
        let coef = if axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * pad.xxt[0][m];
            }
            r * mag_n * weight
        } else {
            alpha * mag_n * weight
        };

        // add contribution to b vector
        if space_ndim == 2 {
            for m in 0..nnode {
                b[ii0 + 0 + m * 2] += coef * nn[m] * v[0];
                b[ii0 + 1 + m * 2] += coef * nn[m] * v[1];
            }
        } else {
            for m in 0..nnode {
                b[ii0 + 0 + m * 3] += coef * nn[m] * v[0];
                b[ii0 + 1 + m * 3] += coef * nn[m] * v[1];
                b[ii0 + 2 + m * 3] += coef * nn[m] * v[2];
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ;
    use crate::integ::testing::aux;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_chk::vec_approx_eq;
    use russell_lab::math::SQRT_2;
    use russell_lab::Vector;

    // to test if variables are cleared before sum
    const NOISE: f64 = 1234.56;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_tri3();
        let mut b = Vector::new(6);
        let mut v = Vector::new(0);
        let un = Vector::new(0);
        let nn = Vector::new(0);
        let f = |_: &mut Vector, _: usize, _: &Vector, _: &Vector| Ok(0.0);
        assert_eq!(f(&mut v, 0, &un, &nn).unwrap(), 0.0);
        let (clear, axis) = (true, false);
        assert_eq!(
            integ::vec_02_nv_bry(&mut b, &mut pad, 0, clear, axis, &[], f).err(),
            Some("in 2D, geometry ndim must be equal to 1 (a line)")
        );
        let mut pad = aux::gen_pad_tet4();
        let mut b = Vector::new(8);
        assert_eq!(
            integ::vec_02_nv_bry(&mut b, &mut pad, 0, clear, axis, &[], f).err(),
            Some("in 3D, geometry ndim must be equal to 2 (a surface)")
        );
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut b = Vector::new(4);
        assert_eq!(
            integ::vec_02_nv_bry(&mut b, &mut pad, 1, clear, axis, &[], f).err(),
            Some("b.len() must be ≥ ii0 + nnode ⋅ space_ndim")
        );
    }

    #[test]
    fn vec_02_nv_bry_works_2d() {
        // Reference:
        // * `sgm:14` -- Smith, Griffiths, Margetts (2014) Programming the Finite Element Method, 5th ed.

        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Lin2).unwrap();
        let ll = 4.0;
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, ll);
        pad.set_xx(1, 1, 0.0);
        let mut b = Vector::filled(pad.kind.nnode() * space_ndim, NOISE);
        let (clear, axis) = (true, false);
        let ips = integ::default_points(pad.kind);
        // uniform
        integ::vec_02_nv_bry(&mut b, &mut pad, 0, clear, axis, ips, |t, _, _, _| {
            t[0] = 0.0;
            t[1] = -1.0;
            Ok(1.0)
        })
        .unwrap();
        vec_approx_eq(b.as_data(), &[0.0, -2.0, 0.0, -2.0], 1e-15);
        // triangular (see @sgm:14\page{605})
        let x_ips = integ::points_coords(&mut pad, ips).unwrap();
        integ::vec_02_nv_bry(&mut b, &mut pad, 0, clear, axis, ips, |t, p, _, _| {
            let c = x_ips[p][0] / ll;
            t[0] = 0.0;
            t[1] = -c;
            Ok(1.0)
        })
        .unwrap();
        vec_approx_eq(b.as_data(), &[0.0, -ll / 6.0, 0.0, -ll / 3.0], 1e-15);

        // example from @sgm:14\page{183}
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Lin3).unwrap();
        let ll = 3.0;
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, ll);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(2, 0, ll / 2.0);
        pad.set_xx(2, 1, 0.0);
        let mut b = Vector::filled(pad.kind.nnode() * space_ndim, NOISE);
        let ips = integ::default_points(pad.kind);
        // uniform
        integ::vec_02_nv_bry(&mut b, &mut pad, 0, clear, axis, ips, |t, _, _, _| {
            t[0] = 0.0;
            t[1] = -1.0;
            Ok(1.0)
        })
        .unwrap();
        vec_approx_eq(b.as_data(), &[0.0, -0.5, 0.0, -0.5, 0.0, -2.0], 1e-15);
        // triangular (see @sgm:14\page{605})
        let x_ips = integ::points_coords(&mut pad, ips).unwrap();
        integ::vec_02_nv_bry(&mut b, &mut pad, 0, clear, axis, ips, |t, p, _, _| {
            let c = x_ips[p][0] / ll;
            t[0] = 0.0;
            t[1] = -c;
            Ok(1.0)
        })
        .unwrap();
        vec_approx_eq(b.as_data(), &[0.0, 0.0, 0.0, -ll / 6.0, 0.0, -ll / 3.0], 1e-15);

        // example from @sgm:14\page{183}
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Lin5).unwrap();
        let ll = 4.0;
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, ll);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(2, 0, ll / 2.0);
        pad.set_xx(2, 1, 0.0);
        pad.set_xx(3, 0, ll / 4.0);
        pad.set_xx(3, 1, 0.0);
        pad.set_xx(4, 0, 3.0 * ll / 4.0);
        pad.set_xx(4, 1, 0.0);
        let mut b = Vector::filled(pad.kind.nnode() * space_ndim, NOISE);
        let ips = integ::default_points(pad.kind);
        // uniform
        integ::vec_02_nv_bry(&mut b, &mut pad, 0, clear, axis, ips, |t, _, _, _| {
            t[0] = 0.0;
            t[1] = -1.0;
            Ok(1.0)
        })
        .unwrap();
        let correct = &[
            0.0,
            -ll * 7.0 / 90.0, // 0
            0.0,
            -ll * 7.0 / 90.0, // 1
            0.0,
            -ll * 2.0 / 15.0, // 2
            0.0,
            -ll * 16.0 / 45.0, // 3
            0.0,
            -ll * 16.0 / 45.0, // 4
        ];
        vec_approx_eq(b.as_data(), correct, 1e-15);
        // triangular (see @sgm:14\page{605})
        let x_ips = integ::points_coords(&mut pad, ips).unwrap();
        integ::vec_02_nv_bry(&mut b, &mut pad, 0, clear, axis, ips, |t, p, _, _| {
            let c = x_ips[p][0] / ll;
            t[0] = 0.0;
            t[1] = -c;
            Ok(1.0)
        })
        .unwrap();
        let correct = &[
            0.0,
            0.0, // 0
            0.0,
            -ll * 7.0 / 90.0, // 1
            0.0,
            -ll * 1.0 / 15.0, // 2
            0.0,
            -ll * 4.0 / 45.0, // 3
            0.0,
            -ll * 4.0 / 15.0, // 4
        ];
        vec_approx_eq(b.as_data(), correct, 1e-15);
    }

    #[test]
    fn vec_02_nv_bry_works_3d() {
        let space_ndim = 3;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Qua4).unwrap();
        let (dx, dy) = (0.5, 1.0);
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(0, 2, 0.0);
        pad.set_xx(1, 0, dx);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(1, 2, 0.0);
        pad.set_xx(2, 0, dx);
        pad.set_xx(2, 1, dy);
        pad.set_xx(2, 2, 0.0);
        pad.set_xx(3, 0, 0.0);
        pad.set_xx(3, 1, dy);
        pad.set_xx(3, 2, 0.0);
        let mut b = Vector::filled(pad.kind.nnode() * space_ndim, NOISE);
        let (clear, axis) = (true, false);
        let ips = integ::default_points(pad.kind);
        integ::vec_02_nv_bry(&mut b, &mut pad, 0, clear, axis, ips, |t, _, _, _| {
            t[0] = 0.0;
            t[1] = 0.0;
            t[2] = -1.0;
            Ok(1.0)
        })
        .unwrap();
        let aa = dx * dy;
        let correct = &[
            0.0,
            0.0,
            -aa / 4.0, // 0
            0.0,
            0.0,
            -aa / 4.0, // 1
            0.0,
            0.0,
            -aa / 4.0, // 2
            0.0,
            0.0,
            -aa / 4.0, // 3
        ];
        vec_approx_eq(b.as_data(), correct, 1e-15);

        // example from @sgm:14\page{195}
        // @sgm:14 Smith, Griffiths, Margetts (2014) Programming the Finite Element Method, 5th ed.
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Qua8).unwrap();
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(0, 2, 0.0);
        pad.set_xx(1, 0, dx);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(1, 2, 0.0);
        pad.set_xx(2, 0, dx);
        pad.set_xx(2, 1, dy);
        pad.set_xx(2, 2, 0.0);
        pad.set_xx(3, 0, 0.0);
        pad.set_xx(3, 1, dy);
        pad.set_xx(3, 2, 0.0);
        pad.set_xx(4, 0, dx / 2.0);
        pad.set_xx(4, 1, 0.0);
        pad.set_xx(4, 2, 0.0);
        pad.set_xx(5, 0, dx);
        pad.set_xx(5, 1, dy / 2.0);
        pad.set_xx(5, 2, 0.0);
        pad.set_xx(6, 0, dx / 2.0);
        pad.set_xx(6, 1, dy);
        pad.set_xx(6, 2, 0.0);
        pad.set_xx(7, 0, 0.0);
        pad.set_xx(7, 1, dy / 2.0);
        pad.set_xx(7, 2, 0.0);
        let mut b = Vector::filled(pad.kind.nnode() * space_ndim, NOISE);
        let ips = integ::default_points(pad.kind);
        integ::vec_02_nv_bry(&mut b, &mut pad, 0, clear, axis, ips, |t, _, _, _| {
            t[0] = 0.0;
            t[1] = 0.0;
            t[2] = -1.0;
            Ok(1.0)
        })
        .unwrap();
        let correct = &[
            0.0,
            0.0,
            aa / 12.0, // 0
            0.0,
            0.0,
            aa / 12.0, // 1
            0.0,
            0.0,
            aa / 12.0, // 2
            0.0,
            0.0,
            aa / 12.0, // 3
            0.0,
            0.0,
            -aa / 3.0, // 4
            0.0,
            0.0,
            -aa / 3.0, // 5
            0.0,
            0.0,
            -aa / 3.0, // 6
            0.0,
            0.0,
            -aa / 3.0, // 7
        ];
        vec_approx_eq(b.as_data(), correct, 1e-15);
    }

    #[test]
    fn vec_02_nv_bry_works_arc() {
        // @bhatti:05 Example 7.9, page 518
        // Reference: Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Lin3).unwrap();
        let r = 5.0;
        pad.set_xx(0, 0, r);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, 0.0);
        pad.set_xx(1, 1, r);
        pad.set_xx(2, 0, r * SQRT_2 / 2.0);
        pad.set_xx(2, 1, r * SQRT_2 / 2.0);
        let mut b = Vector::filled(pad.kind.nnode() * space_ndim, NOISE);
        let (clear, axis) = (true, false);
        let ips = integ::default_points(pad.kind);
        let p = -20.0;
        integ::vec_02_nv_bry(&mut b, &mut pad, 0, clear, axis, ips, |t, _, un, _| {
            t[0] = p * un[0];
            t[1] = p * un[1];
            Ok(1.0)
        })
        .unwrap();
        let correct = &[
            30.4738, 2.85955, // 0
            2.85955, 30.4738, // 1
            66.6667, 66.6667, // 2
        ];
        vec_approx_eq(b.as_data(), correct, 1e-4);
    }
}
