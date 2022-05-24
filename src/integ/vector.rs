use crate::shapes::{IntegPointData, Shape, StateOfShape};
use crate::StrError;
use russell_lab::Vector;

/// Implements the shape(N) times scalar(S) integration case
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
/// * `ips` -- Integration points (n_integ_point)
/// * `erase_a` -- Fills `a` vector with zeros, otherwise accumulate values into `a`
/// * `th` -- The out-of-plane thickness (`tₕ`) in 2D. Use 1.0 for 3D or for plane-stress models.
/// * `fn_s` -- Function `f(p)` corresponding to `s(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
pub fn shape_times_scalar<F>(
    a: &mut Vector,
    state: &mut StateOfShape,
    shape: &Shape,
    ips: IntegPointData,
    th: f64,
    erase_a: bool,
    fn_s: F,
) -> Result<(), StrError>
where
    F: Fn(usize) -> Result<f64, StrError>,
{
    // check
    if a.dim() != shape.nnode {
        return Err("a.len() must be equal to nnode");
    }

    // clear output vector
    if erase_a {
        a.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate interpolation functions and Jacobian
        shape.calc_interp(state, iota)?;
        let det_jac = shape.calc_jacobian(state, iota)?;

        // calculate s
        let s = fn_s(p)?;

        // loop over nodes and perform sum
        let val = s * th * det_jac * weight;
        for m in 0..shape.nnode {
            a[m] += state.interp[m] * val;
        }
    }
    Ok(())
}

/// Implements the the shape(N) times vector(V) integration case
///
/// Interpolation functions times vector field:
///
/// ```text
/// →    ⌠    → →   → →
/// bᵐ = │ Nᵐ(x(ξ)) v(x) tₕ dΩ
///      ⌡
///      Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
/// →    nip-1     →   → →          →
/// bᵐ ≈   Σ    Nᵐ(ιᵖ) v(ιᵖ) tₕ |J|(ιᵖ) wᵖ
///       p=0
/// ```
///
/// # Output
///
/// ```text
///     ┌     ┐
///     | b⁰₀ |
///     | b⁰₁ |
///     | b¹₀ |
/// b = | b¹₁ |
///     | b²₀ |
///     | b²₁ |
///     | ··· |
///     | bᵐᵢ |  ⟸  ii := i + m * space_ndim
///     └     ┘       
///
/// m = ii / space_ndim
/// i = ii % space_ndim
/// ```
///
/// * `b` -- A vector containing all `bᵐᵢ` values, one after another, and sequentially placed
///          as shown above (in 2D). `m` is the index of the node and `i` corresponds to `space_ndim`.
///          The length of `b` must be equal to `nnode * space_ndim`.
///
/// # Updated
///
/// * `state` -- Will be updated by the Shape functions
///
/// # Input
///
/// * `shape` -- Shape functions
/// * `ips` -- Integration points (n_integ_point)
/// * `v` -- All values produced by `v(x(ιᵖ))` (n_integ_point); each v has len = space_ndim
/// * `erase_b` -- fills `b` vector with zeros, otherwise accumulate values into `b`
/// * `th` -- tₕ the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
pub fn shape_times_vector(
    b: &mut Vector,
    state: &mut StateOfShape,
    shape: &Shape,
    ips: IntegPointData,
    v: &Vec<Vector>,
    th: f64,
    erase_b: bool,
) -> Result<(), StrError> {
    // check
    if b.dim() != shape.nnode * shape.space_ndim {
        return Err("b.len() must be equal to nnode * space_ndim");
    }

    // clear output vector
    if erase_b {
        b.fill(0.0);
    }

    // loop over integration points
    for index in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[index];
        let weight = ips[index][3];

        // calculate interpolation functions and Jacobian
        shape.calc_interp(state, iota)?;
        let det_jac = shape.calc_jacobian(state, iota)?;

        // add contribution to b vector
        let coef = th * det_jac * weight;
        let nn = &state.interp;
        if shape.space_ndim == 1 {
            for m in 0..shape.nnode {
                b[m] += coef * nn[m] * v[index][0];
            }
        } else if shape.space_ndim == 2 {
            for m in 0..shape.nnode {
                b[0 + m * 2] += coef * nn[m] * v[index][0];
                b[1 + m * 2] += coef * nn[m] * v[index][1];
            }
        } else {
            for m in 0..shape.nnode {
                b[0 + m * 3] += coef * nn[m] * v[index][0];
                b[1 + m * 3] += coef * nn[m] * v[index][1];
                b[2 + m * 3] += coef * nn[m] * v[index][2];
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{shape_times_scalar, shape_times_vector};
    use crate::shapes::{Shape, StateOfShape, Verification, IP_LIN_LEGENDRE_2, IP_TRI_INTERNAL_1};
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;

    // to test if variables are cleared before sum
    const NOISE: f64 = 1234.56;

    #[test]
    fn capture_some_errors() {
        let shape = Shape::new(2, 1, 2).unwrap();
        let mut state = StateOfShape::new(shape.geo_ndim, &[[0.0, 0.0], [1.0, 0.0]]).unwrap();
        let mut a = Vector::new(3);
        assert_eq!(
            shape_times_scalar(&mut a, &mut state, &shape, &[], 1.0, false, |_| Ok(0.0)).err(),
            Some("a.len() must be equal to nnode")
        );
    }

    #[test]
    fn shape_times_scalar_works_lin2() -> Result<(), StrError> {
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
        let x_ips = shape.calc_integ_points_coords(&mut state, ips)?;
        let mut a = Vector::filled(shape.nnode, NOISE);
        shape_times_scalar(&mut a, &mut state, &shape, ips, 1.0, true, |p| Ok(x_ips[p][0]))?;
        let cf = l / 6.0;
        let (xa, xb) = (state.coords_transp[0][0], state.coords_transp[0][1]);
        let a_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];
        assert_vec_approx_eq!(a.as_data(), a_correct, 1e-15);
        Ok(())
    }

    #[test]
    fn shape_times_scalar_works_tri3() -> Result<(), StrError> {
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
        const CS: f64 = 3.0;
        let mut a = Vector::filled(shape.nnode, NOISE);
        shape_times_scalar(&mut a, &mut state, &shape, ips, 1.0, true, |_| Ok(CS))?;
        let cf = CS * area / 3.0;
        let a_correct = &[cf, cf, cf];
        assert_vec_approx_eq!(a.as_data(), a_correct, 1e-14);
        Ok(())
    }

    #[test]
    fn shape_times_vector_works_lin2() -> Result<(), StrError> {
        // This test is similar to the shape_times_scalar with lin2, however using a vector
        // So, each component of `b` equals `Fₛ`
        let l = 6.0;
        let (shape, mut state) = Verification::line_segment_lin2(l);
        let ips = &IP_LIN_LEGENDRE_2;
        let x_ips = shape.calc_integ_points_coords(&mut state, ips)?;
        let v: Vec<_> = x_ips.iter().map(|x| Vector::filled(shape.space_ndim, x[0])).collect();
        let mut b = Vector::filled(shape.nnode * shape.space_ndim, NOISE);
        shape_times_vector(&mut b, &mut state, &shape, ips, &v, 1.0, true)?;
        let cf = l / 6.0;
        let (xa, xb) = (state.coords_transp[0][0], state.coords_transp[0][1]);
        let b_correct = &[
            cf * (2.0 * xa + xb),
            cf * (2.0 * xa + xb),
            cf * (xa + 2.0 * xb),
            cf * (xa + 2.0 * xb),
        ];
        assert_vec_approx_eq!(b.as_data(), b_correct, 1e-15);
        Ok(())
    }

    #[test]
    fn shape_times_vector_works_tri3() -> Result<(), StrError> {
        // This test is similar to the shape_times_scalar with tri3, however using a vector
        // So, each component of `b` equals `Fₛ`
        let l = 5.0;
        let (shape, mut state, area) = Verification::equilateral_triangle_tri3(l);
        let ips = &IP_TRI_INTERNAL_1;
        const CS: f64 = 3.0;
        let v: Vec<_> = (0..ips.len()).map(|_| Vector::filled(shape.space_ndim, CS)).collect();
        let mut b = Vector::filled(shape.nnode * shape.space_ndim, NOISE);
        shape_times_vector(&mut b, &mut state, &shape, ips, &v, 1.0, true)?;
        let cf = CS * area / 3.0;
        let b_correct = &[cf, cf, cf, cf, cf, cf];
        assert_vec_approx_eq!(b.as_data(), b_correct, 1e-14);
        Ok(())
    }
}
