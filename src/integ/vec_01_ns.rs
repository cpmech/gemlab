use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::Vector;

/// Implements the shape(N) times scalar(S) integration case 01
///
/// Interpolation functions times scalar field:
///
/// ```text
///      ⌠    → →     →
/// aᵐ = │ Nᵐ(x(ξ)) s(x) dΩ
///      ⌡
///      Ωₑ
/// ```
///
/// or, for lines in multi-dimensions:
///
/// ```text
///      ⌠
/// aᵐ = │ Nᵐ(ℓ(ξ)) s(ℓ) dℓ
///      ⌡
///      Γₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///      nip-1     →     →       →
/// aᵐ ≈   Σ    Nᵐ(ιᵖ) s(ιᵖ) |J|(ιᵖ) wᵖ
///       p=0
/// ```
///
/// # Output
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
/// * `a` -- A vector containing all `aᵐ` values, one after another, and
///   sequentially placed as shown above. `m` is the index of the node.
///   The length must be `a.len() ≥ ii0 + nnode`
/// * `pad` -- Some members of the scratchpad will be modified
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output vector where to add components
/// * `clear_a` -- Fills `a` vector with zeros, otherwise accumulate values into `a`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_s` -- Function `f(p)` that computes `s(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
///
/// # Examples
///
/// ```
/// use gemlab::integ;
/// use gemlab::shapes::{GeoKind, Scratchpad};
/// use gemlab::StrError;
/// use russell_chk::assert_vec_approx_eq;
/// use russell_lab::Vector;
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
///     integ::vec_01_ns(&mut a, &mut pad, 0, true, ips, |_| Ok(5.0))?;
///     // solution (cₛ = 5, A = 6):
///     // aᵐ = cₛ A / 3 = 10
///     assert_vec_approx_eq!(a.as_data(), &[10.0, 10.0, 10.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_01_ns<F>(
    a: &mut Vector,
    pad: &mut Scratchpad,
    ii0: usize,
    clear_a: bool,
    ips: IntegPointData,
    fn_s: F,
) -> Result<(), StrError>
where
    F: Fn(usize) -> Result<f64, StrError>,
{
    // check
    let nnode = pad.interp.dim();
    if a.dim() < ii0 + nnode {
        return Err("a.len() must be ≥ ii0 + nnode");
    }

    // clear output vector
    if clear_a {
        a.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate interpolation functions and Jacobian
        (pad.fn_interp)(&mut pad.interp, iota);
        let det_jac = pad.calc_jacobian(iota)?;

        // calculate s
        let s = fn_s(p)?;

        // loop over nodes and perform sum
        let val = s * det_jac * weight;
        for m in 0..nnode {
            a[ii0 + m] += pad.interp[m] * val;
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, AnalyticalTet4, AnalyticalTri3};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut a = Vector::new(2);
        assert_eq!(
            integ::vec_01_ns(&mut a, &mut pad, 1, false, &[], |_| Ok(0.0)).err(),
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
        let (xa, xb) = (pad.xxt[0][0], pad.xxt[0][1]);
        let a_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-15, 1e-14, 1e-14, 1e-14];
        let selection: Vec<_> = [2, 3, 4, 5].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let mut a = Vector::filled(pad.kind.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = integ::points_coords(&mut pad, ips).unwrap();
            integ::vec_01_ns(&mut a, &mut pad, 0, true, ips, |p| Ok(x_ips[p][0])).unwrap();
            assert_vec_approx_eq!(a.as_data(), a_correct, tol);
        });
    }

    #[test]
    fn vec_01_ns_works_tri3_constant() {
        // tri3 with a constant source term s(x) = cₛ
        let mut pad = aux::gen_pad_tri3();

        // solution
        let ana = AnalyticalTri3::new(&pad);
        const CS: f64 = 3.0;
        let a_correct = ana.integ_vec_a(CS);

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
            integ::vec_01_ns(&mut a, &mut pad, 0, true, ips, |_| Ok(CS)).unwrap();
            assert_vec_approx_eq!(a.as_data(), a_correct, tol);
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
            let x_ips = integ::points_coords(&mut pad, ips).unwrap();
            integ::vec_01_ns(&mut a, &mut pad, 0, true, ips, |p| Ok(x_ips[p][2])).unwrap();
            assert_vec_approx_eq!(a.as_data(), a_correct, tol);
        });
    }
}
