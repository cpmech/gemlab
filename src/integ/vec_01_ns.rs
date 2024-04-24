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
/// * `fn_s` -- Function `f(p,N)` that computes `s(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   and shape functions N(ιᵖ).
///
/// # Examples
///
/// ```
/// use gemlab::integ;
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
///     let ips = integ::default_points(pad.kind);
///     let mut a = Vector::filled(pad.kind.nnode(), 0.0);
///     let mut args = integ::CommonArgs::new(&mut pad, ips);
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
    for p in 0..args.ips.len() {
        // ksi coordinates and weight
        let iota = &args.ips[p];
        let weight = args.ips[p][3];

        // calculate interpolation functions and Jacobian
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_jacobian(iota)?;

        // calculate s
        let nn = &args.pad.interp;
        let s = fn_s(p, nn)?;

        // calculate coefficient
        let coef = if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt.get(0, m);
            }
            s * det_jac * weight * args.alpha * r
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
    use crate::integ::{self, AnalyticalTet4, AnalyticalTri3, CommonArgs};
    use russell_lab::{vec_approx_eq, Vector};

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut a = Vector::new(2);
        let nn = Vector::new(0);
        let f = |_: usize, _: &Vector| Ok(0.0);
        assert_eq!(f(0, &nn).unwrap(), 0.0);
        let mut args = CommonArgs::new(&mut pad, &[]);
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
        let selection: Vec<_> = [2, 3, 4, 5].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let mut a = Vector::filled(pad.kind.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            let x_ips = integ::points_coords(&mut args.pad, ips).unwrap();
            integ::vec_01_ns(&mut a, &mut args, |p, _| Ok(x_ips[p][0])).unwrap();
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
            .map(|n| integ::points(class, *n).unwrap())
            .collect();

        // check
        let mut a = Vector::filled(pad.kind.nnode(), aux::NOISE);
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
        // Note that the tolerance is high for n_integ_point = 1
        // because the numerical integration performs poorly with few IPs
        let class = pad.kind.class();
        let tolerances = [1e-13, 1e-13];
        let selection: Vec<_> = [1, 4].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let mut a = Vector::filled(pad.kind.nnode(), aux::NOISE);
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
        // Note that the tolerance is high for n_integ_point = 1
        // because the numerical integration performs poorly with few IPs
        let class = pad.kind.class();
        let tolerances = [0.56, 1e-15, 1e-14, 1e-15, 1e-15, 1e-15, 1e-15, 1e-15];
        let selection: Vec<_> = [1, 4, 5, 8, 14, 15, 24]
            .iter()
            .map(|n| integ::points(class, *n).unwrap())
            .collect();

        // check
        let mut a = Vector::filled(pad.kind.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            let x_ips = integ::points_coords(&mut args.pad, ips).unwrap();
            integ::vec_01_ns(&mut a, &mut args, |p, _| Ok(x_ips[p][2])).unwrap();
            vec_approx_eq(&a, &a_correct, tol);
        });
    }
}
