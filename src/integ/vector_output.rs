use crate::shapes::{IntegPointData, Shape, StateOfShape};
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::Tensor2;

/// Implements the shape(N) times scalar(S) integration case A
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
/// * `th` -- The out-of-plane thickness (`tₕ`) in 2D. Use 1.0 for 3D or for plane-stress models.
/// * `erase_a` -- Fills `a` vector with zeros, otherwise accumulate values into `a`
/// * `fn_s` -- Function `f(p)` corresponding to `s(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
pub fn a_shape_times_scalar<F>(
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

/// Implements the the shape(N) times vector(V) integration case B
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
/// * `th` -- tₕ the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
/// * `erase_b` -- fills `b` vector with zeros, otherwise accumulate values into `b`
/// * `fn_v` -- Function `f(v,p)` corresponding to `v(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
///             The dim of `v` is equal to `space_ndim`.
pub fn b_shape_times_vector<F>(
    b: &mut Vector,
    state: &mut StateOfShape,
    shape: &Shape,
    ips: IntegPointData,
    th: f64,
    erase_b: bool,
    fn_v: F,
) -> Result<(), StrError>
where
    F: Fn(&mut Vector, usize) -> Result<(), StrError>,
{
    // check
    if b.dim() != shape.nnode * shape.space_ndim {
        return Err("b.len() must be equal to nnode * space_ndim");
    }

    // allocate auxiliary vector
    let mut v = Vector::new(shape.space_ndim);

    // clear output vector
    if erase_b {
        b.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate interpolation functions and Jacobian
        shape.calc_interp(state, iota)?;
        let det_jac = shape.calc_jacobian(state, iota)?;

        // calculate v
        fn_v(&mut v, p)?;

        // add contribution to b vector
        let coef = th * det_jac * weight;
        let nn = &state.interp;
        if shape.space_ndim == 2 {
            for m in 0..shape.nnode {
                b[0 + m * 2] += coef * nn[m] * v[0];
                b[1 + m * 2] += coef * nn[m] * v[1];
            }
        } else {
            for m in 0..shape.nnode {
                b[0 + m * 3] += coef * nn[m] * v[0];
                b[1 + m * 3] += coef * nn[m] * v[1];
                b[2 + m * 3] += coef * nn[m] * v[2];
            }
        }
    }
    Ok(())
}

/// Implements the vector(V) dot gradient(G) integration case C
///
/// Vector dot gradient:
///
/// ```text
///      ⌠ → →    →  → →
/// cᵐ = │ w(x) · Gᵐ(x(ξ)) tₕ dΩ
///      ⌡
///      Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///      nip-1  → →     →  →          →
/// cᵐ ≈   Σ    w(ιᵖ) · Gᵐ(ιᵖ) tₕ |J|(ιᵖ) wᵖ
///       p=0
/// ```
///
/// # Output
///
/// ```text
///     ┌     ┐
///     |  c⁰ |
///     |  c¹ |
/// c = |  c² |
///     | ··· |
///     |  cᵐ |
///     └     ┘
/// ```
///
/// * `c` -- A vector containing all `cᵐ` values, one after another, and
///          sequentially placed as shown above. `m` is the index of the node.
///          The length of `c` must be be equal to `nnode`.
///
/// # Updated
///
/// * `state` -- Will be updated by the Shape functions
///
/// # Input
///
/// * `shape` -- Shape functions
/// * `ips` -- Integration points (n_integ_point)
/// * `th` -- The out-of-plane thickness (`tₕ`) in 2D. Use 1.0 for 3D or for plane-stress models.
/// * `erase_c` -- Fills `c` vector with zeros, otherwise accumulate values into `c`
/// * `fn_w` -- Function `f(w,p)` corresponding to `w(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
///             The dim of `w` is equal to `space_ndim`.
pub fn c_vector_dot_gradient<F>(
    c: &mut Vector,
    state: &mut StateOfShape,
    shape: &Shape,
    ips: IntegPointData,
    th: f64,
    erase_c: bool,
    fn_w: F,
) -> Result<(), StrError>
where
    F: Fn(&mut Vector, usize) -> Result<(), StrError>,
{
    // check
    if c.dim() != shape.nnode {
        return Err("c.len() must be equal to nnode");
    }

    // allocate auxiliary vector
    let mut w = Vector::new(shape.space_ndim);

    // clear output vector
    if erase_c {
        c.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate Jacobian and Gradient
        let det_jac = shape.calc_gradient(state, iota)?;

        // calculate w
        fn_w(&mut w, p)?;

        // add contribution to c vector
        let coef = th * det_jac * weight;
        let g = &state.gradient;
        if shape.space_ndim == 2 {
            for m in 0..shape.nnode {
                c[m] += coef * (w[0] * g[m][0] + w[1] * g[m][1]);
            }
        } else {
            for m in 0..shape.nnode {
                c[m] += coef * (w[0] * g[m][0] + w[1] * g[m][1] + w[2] * g[m][2]);
            }
        }
    }
    Ok(())
}

/// Implements the tensor(T) dot gradient(G) integration case D
///
/// Tensor dot gradient:
///
/// ```text
/// →    ⌠   →    →  → →
/// dᵐ = │ σ(x) · Gᵐ(x(ξ)) tₕ dΩ
///      ⌡ ▔
///      Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
/// →    nip-1    →     →  →          →
/// dᵐ ≈   Σ    σ(ιᵖ) · Gᵐ(ιᵖ) tₕ |J|(ιᵖ) wᵖ
///       p=0   ▔
/// ```
///
/// # Output
///
/// ```text
///     ┌     ┐
///     | d⁰₀ |
///     | d⁰₁ |
///     | d¹₀ |
/// d = | d¹₁ |
///     | d²₀ |
///     | d²₁ |
///     | ··· |
///     | dᵐᵢ |  ⟸  ii := i + m * space_ndim
///     └     ┘
///
/// m = ii / space_ndim
/// i = ii % space_ndim
/// ```
///
/// * `d` -- A vector containing all `dᵐᵢ` values, one after another, and sequentially placed
///          as shown above (in 2D). `m` is the index of the node and `i` corresponds to `space_ndim`.
///          The length of `d` must be equal to `nnode * space_ndim`.
///
/// # Updated
///
/// * `state` -- Will be updated by the Shape functions
///
/// # Input
///
/// * `shape` -- Shape functions
/// * `ips` -- Integration points (n_integ_point)
/// * `th` -- The out-of-plane thickness (`tₕ`) in 2D. Use 1.0 for 3D or for plane-stress models.
/// * `erase_d` -- Fills `d` vector with zeros, otherwise accumulate values into `d`
/// * `fn_sig` -- Function `f(sig,p)` corresponding to `σ(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
pub fn d_tensor_dot_gradient<F>(
    d: &mut Vector,
    state: &mut StateOfShape,
    shape: &Shape,
    ips: IntegPointData,
    th: f64,
    erase_d: bool,
    fn_sig: F,
) -> Result<(), StrError>
where
    F: Fn(&mut Tensor2, usize) -> Result<(), StrError>,
{
    // check
    if d.dim() != shape.nnode * shape.space_ndim {
        return Err("d.len() must be equal to nnode * space_ndim");
    }

    // allocate auxiliary tensor
    let mut sig = Tensor2::new(true, shape.space_ndim == 2);

    // clear output vector
    if erase_d {
        d.fill(0.0);
    }

    // loop over integration points
    let s = SQRT_2;
    for index in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[index];
        let weight = ips[index][3];

        // calculate Jacobian and Gradient
        let det_jac = shape.calc_gradient(state, iota)?;

        // calculate σ tensor
        fn_sig(&mut sig, index)?;

        // add contribution to d vector
        let coef = th * det_jac * weight;
        let g = &state.gradient;
        let t = &sig.vec;
        if shape.space_ndim == 2 {
            for m in 0..shape.nnode {
                d[0 + m * 2] += coef * (t[0] * g[m][0] + t[3] * g[m][1] / s);
                d[1 + m * 2] += coef * (t[3] * g[m][0] / s + t[1] * g[m][1]);
            }
        } else {
            for m in 0..shape.nnode {
                d[0 + m * 3] += coef * (t[0] * g[m][0] + t[3] * g[m][1] / s + t[5] * g[m][2] / s);
                d[1 + m * 3] += coef * (t[3] * g[m][0] / s + t[1] * g[m][1] + t[4] * g[m][2] / s);
                d[2 + m * 3] += coef * (t[5] * g[m][0] / s + t[4] * g[m][1] / s + t[2] * g[m][2]);
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{a_shape_times_scalar, b_shape_times_vector, c_vector_dot_gradient, d_tensor_dot_gradient};
    use crate::shapes::{AnalyticalTri3, Shape, StateOfShape, Verification, IP_LIN_LEGENDRE_2, IP_TRI_INTERNAL_1};
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
            a_shape_times_scalar(&mut a, &mut state, &shape, &[], 1.0, false, |_| Ok(0.0)).err(),
            Some("a.len() must be equal to nnode")
        );
        let mut b = Vector::new(5);
        assert_eq!(
            b_shape_times_vector(&mut b, &mut state, &shape, &[], 1.0, false, |_, _| Ok(())).err(),
            Some("b.len() must be equal to nnode * space_ndim")
        );
        let mut c = Vector::new(3);
        assert_eq!(
            c_vector_dot_gradient(&mut c, &mut state, &shape, &[], 1.0, false, |_, _| Ok(())).err(),
            Some("c.len() must be equal to nnode")
        );
        let mut d = Vector::new(5);
        assert_eq!(
            d_tensor_dot_gradient(&mut d, &mut state, &shape, &[], 1.0, false, |_, _| Ok(())).err(),
            Some("d.len() must be equal to nnode * space_ndim")
        );
    }

    #[test]
    fn a_shape_times_scalar_works_lin2() -> Result<(), StrError> {
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
        a_shape_times_scalar(&mut a, &mut state, &shape, ips, 1.0, true, |p| Ok(x_ips[p][0]))?;
        let cf = l / 6.0;
        let (xa, xb) = (state.coords_transp[0][0], state.coords_transp[0][1]);
        let a_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];
        assert_vec_approx_eq!(a.as_data(), a_correct, 1e-15);
        Ok(())
    }

    #[test]
    fn a_shape_times_scalar_works_tri3() -> Result<(), StrError> {
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
        a_shape_times_scalar(&mut a, &mut state, &shape, ips, 1.0, true, |_| Ok(CS))?;
        let cf = CS * area / 3.0;
        let a_correct = &[cf, cf, cf];
        assert_vec_approx_eq!(a.as_data(), a_correct, 1e-14);
        Ok(())
    }

    #[test]
    fn b_shape_times_vector_works_lin2() -> Result<(), StrError> {
        // This test is similar to the shape_times_scalar with lin2, however using a vector
        // So, each component of `b` equals `Fₛ`
        let l = 6.0;
        let (shape, mut state) = Verification::line_segment_lin2(l);
        let ips = &IP_LIN_LEGENDRE_2;
        let x_ips = shape.calc_integ_points_coords(&mut state, ips)?;
        let mut b = Vector::filled(shape.nnode * shape.space_ndim, NOISE);
        b_shape_times_vector(&mut b, &mut state, &shape, ips, 1.0, true, |v, p| {
            v[0] = x_ips[p][0];
            v[1] = x_ips[p][0]; // << note use of x component here too
            Ok(())
        })?;
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
    fn b_shape_times_vector_works_tri3() -> Result<(), StrError> {
        // This test is similar to the shape_times_scalar with tri3, however using a vector
        // So, each component of `b` equals `Fₛ`
        let l = 5.0;
        let (shape, mut state, area) = Verification::equilateral_triangle_tri3(l);
        let ips = &IP_TRI_INTERNAL_1;
        const CS: f64 = 3.0;
        let mut b = Vector::filled(shape.nnode * shape.space_ndim, NOISE);
        b_shape_times_vector(&mut b, &mut state, &shape, ips, 1.0, true, |v, _| {
            v[0] = CS;
            v[1] = CS;
            Ok(())
        })?;
        let cf = CS * area / 3.0;
        let b_correct = &[cf, cf, cf, cf, cf, cf];
        assert_vec_approx_eq!(b.as_data(), b_correct, 1e-14);
        Ok(())
    }

    #[test]
    fn c_vector_dot_gradient_works_constant() -> Result<(), StrError> {
        // constant vector function: w(x) = {w₀, w₁}
        // solution:
        //    cᵐ = ½ (w₀ bₘ + w₁ cₘ)
        const W0: f64 = 2.0;
        const W1: f64 = 3.0;
        let (shape, mut state, _) = Verification::equilateral_triangle_tri3(5.0);
        let ips = &IP_TRI_INTERNAL_1;
        let mut c = Vector::filled(shape.nnode, NOISE);
        c_vector_dot_gradient(&mut c, &mut state, &shape, ips, 1.0, true, |w, _| {
            w[0] = W0;
            w[1] = W1;
            Ok(())
        })?;
        let mut ana = AnalyticalTri3::new(&shape, &mut state);
        let c_correct = ana.integ_vec_c_constant(W0, W1);
        assert_vec_approx_eq!(c.as_data(), c_correct, 1e-15);
        Ok(())
    }

    #[test]
    fn c_vector_dot_gradient_works_bilinear() -> Result<(), StrError> {
        // bilinear vector function: w(x) = {x, y}
        // solution:
        //    cᵐ = ⅙ bₘ (x₀+x₁+x₂) + ⅙ cₘ (y₀+y₁+y₂)
        let (shape, mut state, _) = Verification::equilateral_triangle_tri3(5.0);
        let ips = &IP_TRI_INTERNAL_1;
        let x_ips = shape.calc_integ_points_coords(&mut state, ips)?;
        let mut c = Vector::filled(shape.nnode, NOISE);
        c_vector_dot_gradient(&mut c, &mut state, &shape, ips, 1.0, true, |w, p| {
            w[0] = x_ips[p][0];
            w[1] = x_ips[p][1];
            Ok(())
        })?;
        let mut ana = AnalyticalTri3::new(&shape, &mut state);
        let c_correct = ana.integ_vec_c_bilinear();
        assert_vec_approx_eq!(c.as_data(), c_correct, 1e-14);
        Ok(())
    }

    #[test]
    fn d_tensor_dot_gradient_works() -> Result<(), StrError> {
        // constant tensor function: σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2}
        // solution:
        //    dᵐ₀ = ½ (σ₀₀ bₘ + σ₀₁ cₘ)
        //    dᵐ₁ = ½ (σ₁₀ bₘ + σ₁₁ cₘ)
        const S00: f64 = 2.0;
        const S11: f64 = 3.0;
        const S22: f64 = 4.0;
        const S01: f64 = 5.0;
        let (shape, mut state, _) = Verification::equilateral_triangle_tri3(5.0);
        let ips = &IP_TRI_INTERNAL_1;
        let mut d = Vector::filled(shape.nnode * shape.space_ndim, NOISE);
        d_tensor_dot_gradient(&mut d, &mut state, &shape, ips, 1.0, true, |sig, _| {
            sig.sym_set(0, 0, S00);
            sig.sym_set(1, 1, S11);
            sig.sym_set(2, 2, S22);
            sig.sym_set(0, 1, S01);
            Ok(())
        })?;
        let mut ana = AnalyticalTri3::new(&shape, &mut state);
        let d_correct = ana.integ_vec_d_constant(S00, S11, S01);
        assert_vec_approx_eq!(d.as_data(), d_correct, 1e-15);
        Ok(())
    }
}
