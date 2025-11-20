use super::CommonArgs;
use crate::StrError;
use russell_lab::Vector;

/// Implements the shape(N) times scalar(S) integration case 01 (boundary integral version)
///
/// Callback function: `s ← f(p, un, N)`
///
/// Interpolation functions times scalar field:
///
/// ```text
///      ⌠    → →     →
/// aᵐ = │ Nᵐ(x(ξ)) s(x) α dΩ
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
/// * `fn_s` -- Function `f(p,un,N)` that computes `s(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`,
///   the **unit** normal vector `un(x(ιᵖ))`, and shape functions N(ιᵖ).
///
/// # Requirements
///
/// * In 2D, geometry ndim must be equal to 1 (a line)
/// * In 3D, geometry ndim must be equal to 2 (a surface)
///
/// # Examples
///
pub fn vec_01_ns_bry<F>(a: &mut Vector, args: &mut CommonArgs, mut fn_s: F) -> Result<(), StrError>
where
    F: FnMut(usize, &Vector, &Vector) -> Result<f64, StrError>,
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
    if a.dim() < ii0 + nnode {
        return Err("a.len() must be ≥ ii0 + nnode");
    }

    // allocate auxiliary vectors
    let mut un = Vector::new(space_ndim); // unit normal vector

    // clear output vector
    if args.clear {
        a.fill(0.0);
    }

    // loop over integration points
    for p in 0..args.gauss.npoint() {
        // ksi coordinates and weight
        let iota = args.gauss.coords(p);
        let weight = args.gauss.weight(p);

        // calculate interpolation functions and unit normal vector
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let mag_n = args.pad.calc_normal_vector(&mut un, iota)?;

        // calculate s
        let nn = &args.pad.interp;
        let s = fn_s(p, &un, nn)?;

        // calculate coefficient
        let coef = if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt.get(0, m);
            }
            s * mag_n * weight * args.alpha * r
        } else {
            s * mag_n * weight * args.alpha
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
    use crate::integ::{self, CommonArgs, Gauss};
    use crate::mesh::{GeoKind, Mesh};
    use crate::recovery;
    use crate::shapes::Scratchpad;
    use russell_lab::{vec_approx_eq, vec_inner, Vector};

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut a = Vector::new(2);
        let nn = Vector::new(0);
        let un = Vector::new(0);
        let f = |_: usize, _: &Vector, _: &Vector| Ok(0.0);
        assert_eq!(f(0, &un, &nn).unwrap(), 0.0);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::vec_01_ns_bry(&mut a, &mut args, f).err(),
            Some("a.len() must be ≥ ii0 + nnode")
        );
    }

    #[test]
    fn vec_01_ns_bry_works_lin2_ignoring_un() {
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
        let ips = Gauss::new_sized(class, 2).unwrap();

        // check
        let mut a = Vector::filled(pad.kind.nnode(), aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &ips);
        let x_ips = recovery::get_points_coords(&mut args.pad, &ips).unwrap();
        integ::vec_01_ns_bry(&mut a, &mut args, |p, _, _| Ok(x_ips[p][0])).unwrap();
        vec_approx_eq(&a, a_correct, 1e-15);
    }

    #[test]
    fn vec_01_ns_bry_works_lin2_constant_w() {
        // lin2 with constant flow vector
        //
        // s(x) = [w0, w1]ᵀ · un
        //
        let mut pad = Scratchpad::new(2, GeoKind::Lin2).unwrap();
        let (xa, ya) = (1.0, 2.0);
        let (xb, yb) = (4.0, 6.0);
        pad.set_xx(0, 0, xa);
        pad.set_xx(0, 1, ya);
        pad.set_xx(1, 0, xb);
        pad.set_xx(1, 1, yb);

        // solution
        let (w0, w1) = (3.0, 5.0);
        let dx = xb - xa;
        let dy = yb - ya;
        let a_correct = &[(w1 * dx - w0 * dy) / 2.0, (w1 * dx - w0 * dy) / 2.0];

        // integration points
        let class = pad.kind.class();
        let ips = Gauss::new_sized(class, 2).unwrap();

        // check
        let mut a = Vector::filled(pad.kind.nnode(), aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &ips);
        integ::vec_01_ns_bry(&mut a, &mut args, |_, un, _| {
            let s = w0 * un[0] + w1 * un[1];
            Ok(s)
        })
        .unwrap();
        vec_approx_eq(&a, a_correct, 1e-15);
    }

    #[test]
    fn vec_01_ns_bry_works_lin2_linear_w() {
        // lin2 with linear flow vector
        //
        // s(x) = [x+y, y-x]ᵀ · un
        //
        let mut pad = Scratchpad::new(2, GeoKind::Lin2).unwrap();
        let (xa, ya) = (1.0, 2.0);
        let (xb, yb) = (4.0, 6.0);
        pad.set_xx(0, 0, xa);
        pad.set_xx(0, 1, ya);
        pad.set_xx(1, 0, xb);
        pad.set_xx(1, 1, yb);

        // solution
        let dx = xb - xa;
        let dy = yb - ya;
        let aa = dx - dy;
        let bb = dx * dx + dy * dy;
        let cc = dx + dy;
        let mm = 3.0 * aa * ya - 3.0 * xa * cc;
        let a_correct = &[(mm - bb) / 6.0, (mm - 2.0 * bb) / 6.0];

        // integration points
        let class = pad.kind.class();
        let ips = Gauss::new_sized(class, 2).unwrap();

        // check
        let mut a = Vector::filled(pad.kind.nnode(), aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &ips);
        let x_ips = recovery::get_points_coords(&mut args.pad, &ips).unwrap();
        integ::vec_01_ns_bry(&mut a, &mut args, |p, un, _| {
            let x = x_ips[p][0];
            let y = x_ips[p][1];
            let s = (x + y) * un[0] + (y - x) * un[1];
            Ok(s)
        })
        .unwrap();
        vec_approx_eq(&a, a_correct, 1e-14);
    }

    #[test]
    fn vec_01_ns_bry_works_3d_plane() {
        // Blender:
        // Apply a Quaternion rotation of (5/6, 1/6, 1/2, 1/6) to the plane with normal (0,0,1)
        let mesh = Mesh::read("data/blender/single-plane.blend.msh").unwrap();

        // check unit normal
        // Mathematica:
        // w = ResourceFunction["Quaternion"][5/6, 1/6, 1/2, 1/6];
        // m = ResourceFunction["QuaternionToRotationMatrix"][w];
        // un = {0, 0, 1};
        // m . un
        // Output: {8/9, -(1/9), 4/9}
        let un_correct = Vector::from(&[8.0 / 9.0, -(1.0 / 9.0), 4.0 / 9.0]);

        // integration points
        let mut pad = mesh.get_pad(0);
        let class = pad.kind.class();
        let ips = Gauss::new_sized(class, 4).unwrap();

        // perform integration with constant flow vector
        let w = Vector::from(&[1.0, 2.0, 3.0]);
        let s = vec_inner(&w, &un_correct);
        let mut a = Vector::filled(pad.kind.nnode(), aux::NOISE);
        let mut args = CommonArgs::new(&mut pad, &ips);
        integ::vec_01_ns_bry(&mut a, &mut args, |_, un, _| {
            vec_approx_eq(&un, &un_correct, 1e-7); // note that blender is single precision
            Ok(vec_inner(&w, un))
        })
        .unwrap();
        vec_approx_eq(&a, &[s, s, s, s], 1e-7);
    }
}
