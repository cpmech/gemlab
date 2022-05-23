use crate::shapes::{IntegPointData, Shape, StateOfShape};
use crate::util::AsArray1D;
use crate::StrError;
use russell_lab::Vector;

/// Implements the shape(N)-scalar(S) integration case
///
/// Interpolation functions times scalar field:
///
/// ```text
///      ⌠    → →     →
/// aᵐ = │ Nᵐ(x(ξ)) s(x) tₕ dΩ
///      ⌡
///      Ωₑ
/// ```
///
/// or, for lines in multi-dimensions:
///
/// ```text
///      ⌠
/// aᵐ = │ Nᵐ(ℓ(ξ)) s(ℓ) tₕ dℓ
///      ⌡
///      Γₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///      nip-1     →     →          →
/// aᵐ ≈   Σ    Nᵐ(ιᵖ) s(ιᵖ) tₕ |J|(ιᵖ) wᵖ
///       p=0
/// ```
///
/// # Output
///
/// ```text
///     ┌     ┐
///     |  a⁰ |
///     |  a¹ |
/// a = |  a² |
///     | ··· |
///     |  aᵐ |
///     └     ┘
/// ```
///
/// * `a` -- A vector containing all `aᵐ` values, one after another, and
///          sequentially placed as shown above. `m` is the index of the node.
///          The length of `a` must be equal to `nnode`.
///
/// # Updated
///
/// * `state` -- Will be updated by the Shape functions
///
/// # Input
///
/// * `shape` -- Shape functions
/// * `integ_points` -- Integration points (n_integ_point)
/// * `s` -- All values produced by `s(x(ιᵖ))` (n_integ_point)
/// * `erase_a` -- Fills `a` vector with zeros, otherwise accumulate values into `a`
/// * `th` -- The out-of-plane thickness (`tₕ`) in 2D. Use 1.0 for 3D or for plane-stress models.
pub fn shape_scalar<'a, T>(
    a: &mut Vector,
    state: &mut StateOfShape,
    shape: &Shape,
    ips: IntegPointData,
    s: &'a T,
    th: f64,
    erase_a: bool,
) -> Result<(), StrError>
where
    T: AsArray1D<'a, f64>,
{
    // check
    if a.dim() != shape.nnode {
        return Err("the length of vector 'a' must be equal to nnode");
    }

    // clear output vector
    if erase_a {
        a.fill(0.0);
    }

    // loop over integration points
    for index in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[index];
        let weight = ips[index][3];

        // calculate interpolation functions and Jacobian
        shape.calc_interp(state, iota)?;
        let det_jac = shape.calc_jacobian(state, iota)?;

        // loop over nodes and perform sum
        for m in 0..shape.nnode {
            a[m] += state.interp[m] * s.at(index) * th * det_jac * weight;
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::shape_scalar;
    use crate::shapes::{Verification, IP_LIN_LEGENDRE_2, IP_TRI_INTERNAL_1};
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;

    // to test if variables are cleared before sum
    const NOISE: f64 = 1234.56;

    #[test]
    fn shape_scalar_works_lin2() -> Result<(), StrError> {
        // lin2 with linear source term:
        //
        // s(x) = x
        //
        // solution:
        //
        //        ┌           ┐
        //      l │ 2 xa + xb │
        // Fₛ = — │           │
        //      6 │ xa + 2 xb │
        //        └           ┘
        let l = 6.0;
        let (shape, mut state) = Verification::line_segment_lin2(l);
        let ips = &IP_LIN_LEGENDRE_2;
        let mut a = Vector::filled(shape.nnode, NOISE);
        let x_ips = shape.calc_integ_points_coords(&mut state, ips)?;
        let s: Vec<f64> = x_ips.iter().map(|x| x[0]).collect();
        shape_scalar(&mut a, &mut state, &shape, ips, &s, 1.0, true)?;
        let cf = l / 6.0;
        let (xa, xb) = (state.coords_transp[0][0], state.coords_transp[0][1]);
        let a_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];
        assert_vec_approx_eq!(a.as_data(), a_correct, 1e-15);
        Ok(())
    }

    #[test]
    fn shape_scalar_works_tri3() -> Result<(), StrError> {
        // tri3 with a constant source term:
        //
        // s(x) = cₛ
        //
        // solution:
        //           ┌   ┐
        //      cₛ A │ 1 │
        // Fₛ = ———— │ 1 │
        //        3  │ 1 │
        //           └   ┘
        let l = 5.0;
        let (shape, mut state, area) = Verification::equilateral_triangle_tri3(l);
        let ips = &IP_TRI_INTERNAL_1;
        let mut a = Vector::filled(shape.nnode, NOISE);
        const CS: f64 = 3.0;
        let s: Vec<f64> = (0..ips.len()).map(|_| CS).collect();
        shape_scalar(&mut a, &mut state, &shape, ips, &s, 1.0, true)?;
        let cf = CS * area / 3.0;
        let a_correct = &[cf, cf, cf];
        assert_vec_approx_eq!(a.as_data(), a_correct, 1e-14);
        Ok(())
    }
}
