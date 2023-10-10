use super::CommonArgs;
use crate::StrError;
use russell_lab::Vector;

/// Implements the the shape(N) times vector(V) integration case 02 (boundary integral version)
///
/// Callback function: `f(v, p, un, N)`
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
/// # Results
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
/// # Arguments
///
/// * `b` -- A vector containing all `bᵐᵢ` values, one after another, and sequentially placed
///   as shown above (in 2D). `m` is the index of the node and `i` corresponds to `space_ndim`.
///   The length must be `b.len() ≥ ii0 + nnode ⋅ space_ndim`
/// * `args` --- Common arguments
/// * `fn_v` -- Function `f(v,p,un,N)` that calculates `v(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   the **unit** normal vector `un(x(ιᵖ))`, and shape functions N(ιᵖ).
///   `v.dim() = space_ndim` and `un.dim() = space_ndim`.
///
/// # Examples
///
pub fn vec_02_nv_bry<F>(b: &mut Vector, args: &mut CommonArgs, mut fn_v: F) -> Result<(), StrError>
where
    F: FnMut(&mut Vector, usize, &Vector, &Vector) -> Result<(), StrError>,
{
    // check
    let (space_ndim, nnode) = args.pad.xxt.dims();
    let geo_ndim = args.pad.deriv.dims().1;
    if space_ndim == 2 && geo_ndim != 1 {
        return Err("in 2D, geometry ndim must be equal to 1 (a line)");
    }
    if space_ndim == 3 && geo_ndim != 2 {
        return Err("in 3D, geometry ndim must be equal to 2 (a surface)");
    }
    let ii0 = args.ii0;
    if b.dim() < ii0 + nnode * space_ndim {
        return Err("b.len() must be ≥ ii0 + nnode ⋅ space_ndim");
    }

    // allocate auxiliary vectors
    let mut v = Vector::new(space_ndim);
    let mut un = Vector::new(space_ndim); // unit normal vector

    // clear output vector
    if args.clear {
        b.fill(0.0);
    }

    // loop over integration points
    for p in 0..args.ips.len() {
        // ksi coordinates and weight
        let iota = &args.ips[p];
        let weight = args.ips[p][3];

        // calculate interpolation functions and unit normal vector
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let mag_n = args.pad.calc_normal_vector(&mut un, iota)?;

        // calculate t
        let nn = &args.pad.interp;
        fn_v(&mut v, p, &un, nn)?;

        // calculate coefficient
        let coef = if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt.get(0, m);
            }
            mag_n * weight * args.alpha * r
        } else {
            mag_n * weight * args.alpha
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
    use crate::integ::testing::aux;
    use crate::integ::{self, CommonArgs};
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::math::SQRT_2;
    use russell_lab::{vec_approx_eq, Vector};

    // to test if variables are cleared before sum
    const NOISE: f64 = 1234.56;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_tri3();
        let mut b = Vector::new(6);
        let mut v = Vector::new(0);
        let un = Vector::new(0);
        let nn = Vector::new(0);
        let f = |_: &mut Vector, _: usize, _: &Vector, _: &Vector| Ok(());
        f(&mut v, 0, &un, &nn).unwrap();
        let mut args = CommonArgs::new(&mut pad, &[]);
        assert_eq!(
            integ::vec_02_nv_bry(&mut b, &mut args, f).err(),
            Some("in 2D, geometry ndim must be equal to 1 (a line)")
        );
        let mut pad = aux::gen_pad_tet4();
        let mut args = CommonArgs::new(&mut pad, &[]);
        let mut b = Vector::new(8);
        assert_eq!(
            integ::vec_02_nv_bry(&mut b, &mut args, f).err(),
            Some("in 3D, geometry ndim must be equal to 2 (a surface)")
        );
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut args = CommonArgs::new(&mut pad, &[]);
        let mut b = Vector::new(4);
        args.ii0 = 1;
        assert_eq!(
            integ::vec_02_nv_bry(&mut b, &mut args, f).err(),
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
        let ips = integ::default_points(pad.kind);
        // uniform
        let mut args = CommonArgs::new(&mut pad, ips);
        integ::vec_02_nv_bry(&mut b, &mut args, |t, _, _, _| {
            t[0] = 0.0;
            t[1] = -1.0;
            Ok(())
        })
        .unwrap();
        vec_approx_eq(b.as_data(), &[0.0, -2.0, 0.0, -2.0], 1e-15);
        // triangular (see @sgm:14\page{605})
        let mut args = CommonArgs::new(&mut pad, ips);
        let x_ips = integ::points_coords(&mut args.pad, ips).unwrap();
        integ::vec_02_nv_bry(&mut b, &mut args, |t, p, _, _| {
            let c = x_ips[p][0] / ll;
            t[0] = 0.0;
            t[1] = -c;
            Ok(())
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
        let mut args = CommonArgs::new(&mut pad, ips);
        // uniform
        integ::vec_02_nv_bry(&mut b, &mut args, |t, _, _, _| {
            t[0] = 0.0;
            t[1] = -1.0;
            Ok(())
        })
        .unwrap();
        vec_approx_eq(b.as_data(), &[0.0, -0.5, 0.0, -0.5, 0.0, -2.0], 1e-15);
        // triangular (see @sgm:14\page{605})
        let x_ips = integ::points_coords(&mut pad, ips).unwrap();
        let mut args = CommonArgs::new(&mut pad, ips);
        integ::vec_02_nv_bry(&mut b, &mut args, |t, p, _, _| {
            let c = x_ips[p][0] / ll;
            t[0] = 0.0;
            t[1] = -c;
            Ok(())
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
        let mut args = CommonArgs::new(&mut pad, ips);
        // uniform
        integ::vec_02_nv_bry(&mut b, &mut args, |t, _, _, _| {
            t[0] = 0.0;
            t[1] = -1.0;
            Ok(())
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
        let mut args = CommonArgs::new(&mut pad, ips);
        integ::vec_02_nv_bry(&mut b, &mut args, |t, p, _, _| {
            let c = x_ips[p][0] / ll;
            t[0] = 0.0;
            t[1] = -c;
            Ok(())
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
        let ips = integ::default_points(pad.kind);
        let mut args = CommonArgs::new(&mut pad, ips);
        integ::vec_02_nv_bry(&mut b, &mut args, |t, _, _, _| {
            t[0] = 0.0;
            t[1] = 0.0;
            t[2] = -1.0;
            Ok(())
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
        let mut args = CommonArgs::new(&mut pad, ips);
        integ::vec_02_nv_bry(&mut b, &mut args, |t, _, _, _| {
            t[0] = 0.0;
            t[1] = 0.0;
            t[2] = -1.0;
            Ok(())
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
        let ips = integ::default_points(pad.kind);
        let mut args = CommonArgs::new(&mut pad, ips);
        let p = -20.0;
        integ::vec_02_nv_bry(&mut b, &mut args, |t, _, un, _| {
            t[0] = p * un[0];
            t[1] = p * un[1];
            Ok(())
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
