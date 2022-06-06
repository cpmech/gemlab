use super::IntegPointData;
use crate::shapes::{Shape, StateOfShape};
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::Tensor2;

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
/// * `th` -- The out-of-plane thickness (`tₕ`) in 2D. Use 1.0 for 3D or for plane-stress models.
/// * `erase_a` -- Fills `a` vector with zeros, otherwise accumulate values into `a`
/// * `fn_s` -- Function `f(p)` corresponding to `s(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
///
/// # Examples
///
/// See also the `examples` directory.
///
/// ```
/// use gemlab::integ::{default_integ_points, vec_a_shape_times_scalar};
/// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
/// use gemlab::StrError;
/// use russell_chk::assert_vec_approx_eq;
/// use russell_lab::Vector;
///
/// fn main() -> Result<(), StrError> {
///     let shape = Shape::new(GeoKind::Tri3);
///     let mut state = StateOfShape::new(shape.kind, &[[2.0, 3.0], [6.0, 3.0], [2.0, 6.0]])?;
///     let ips = default_integ_points(shape.kind);
///     let mut a = Vector::filled(shape.nnode, 0.0);
///     vec_a_shape_times_scalar(&mut a, &mut state, &shape, ips, 1.0, true, |_| Ok(5.0))?;
///     // solution (A = 6):
///     // aᵐ = cₛ A / 3
///     assert_vec_approx_eq!(a.as_data(), &[10.0, 10.0, 10.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_a_shape_times_scalar<F>(
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
/// * `th` -- tₕ the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
/// * `erase_b` -- fills `b` vector with zeros, otherwise accumulate values into `b`
/// * `fn_v` -- Function `f(v,p)` corresponding to `v(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
///             The dim of `v` is equal to `space_ndim`.
///
/// # Examples
///
/// See also the `examples` directory.
///
/// ```
/// use gemlab::integ::{default_integ_points, vec_b_shape_times_vector};
/// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
/// use gemlab::StrError;
/// use russell_chk::assert_vec_approx_eq;
/// use russell_lab::Vector;
///
/// fn main() -> Result<(), StrError> {
///     let space_ndim = 2;
///     let shape = Shape::new(GeoKind::Tri3);
///     let mut state = StateOfShape::new(shape.kind, &[[2.0, 3.0], [6.0, 3.0], [2.0, 6.0]])?;
///     let ips = default_integ_points(shape.kind);
///     let mut b = Vector::filled(shape.nnode * space_ndim, 0.0);
///     vec_b_shape_times_vector(&mut b, &mut state, &shape, ips, 1.0, true, |v, _| {
///         v[0] = 1.0;
///         v[1] = 2.0;
///         Ok(())
///     })?;
///     // solution (A = 6):
///     // bᵐ₀ = v₀ A / 3
///     // bᵐ₁ = v₁ A / 3
///     assert_vec_approx_eq!(b.as_data(), &[2.0, 4.0, 2.0, 4.0, 2.0, 4.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_b_shape_times_vector<F>(
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
    let space_ndim = state.coords_min.len();
    if b.dim() != shape.nnode * space_ndim {
        return Err("b.len() must be equal to nnode * space_ndim");
    }

    // allocate auxiliary vector
    let mut v = Vector::new(space_ndim);

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
        if space_ndim == 2 {
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

/// Implements the vector(V) dot gradient(G) integration case
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
///
/// # Examples
///
/// See also the `examples` directory.
///
/// ```
/// use gemlab::integ::{default_integ_points, vec_c_vector_dot_gradient};
/// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
/// use gemlab::StrError;
/// use russell_chk::assert_vec_approx_eq;
/// use russell_lab::Vector;
///
/// fn main() -> Result<(), StrError> {
///     let shape = Shape::new(GeoKind::Tri3);
///     let mut state = StateOfShape::new(shape.kind, &[[2.0, 3.0], [6.0, 3.0], [2.0, 6.0]])?;
///     let ips = default_integ_points(shape.kind);
///     let mut c = Vector::filled(shape.nnode, 0.0);
///     vec_c_vector_dot_gradient(&mut c, &mut state, &shape, ips, 1.0, true, |w, _| {
///         w[0] = 1.0;
///         w[1] = 2.0;
///         Ok(())
///     })?;
///     // solution (A = 6):
///     // cᵐ = (w₀ Gᵐ₀ + w₁ Gᵐ₁) A
///     //     ┌       ┐
///     //     │ -¼ -⅓ │
///     // G = │  ¼  0 │
///     //     │  0  ⅓ │
///     //     └       ┘
///     assert_vec_approx_eq!(c.as_data(), &[-5.5, 1.5, 4.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_c_vector_dot_gradient<F>(
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
    let space_ndim = state.coords_min.len();
    let mut w = Vector::new(space_ndim);

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
        if space_ndim == 2 {
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

/// Implements the tensor(T) dot gradient(G) integration case
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
///
/// # Examples
///
/// See also the `examples` directory.
///
/// ```
/// use gemlab::integ::{default_integ_points, vec_d_tensor_dot_gradient};
/// use gemlab::shapes::{GeoKind, Shape, StateOfShape};
/// use gemlab::StrError;
/// use russell_chk::assert_vec_approx_eq;
/// use russell_lab::Vector;
///
/// fn main() -> Result<(), StrError> {
///     let space_ndim = 2;
///     let shape = Shape::new(GeoKind::Tri3);
///     let mut state = StateOfShape::new(shape.kind, &[[2.0, 3.0], [6.0, 3.0], [2.0, 6.0]])?;
///     let ips = default_integ_points(shape.kind);
///     let mut d = Vector::filled(shape.nnode * space_ndim, 0.0);
///     vec_d_tensor_dot_gradient(&mut d, &mut state, &shape, ips, 1.0, true, |sig, _| {
///         sig.sym_set(0, 0, 1.0);
///         sig.sym_set(1, 1, 2.0);
///         sig.sym_set(0, 1, 3.0);
///         Ok(())
///     })?;
///     // solution (A = 6):
///     // dᵐ₀ = (σ₀₀ Gᵐ₀ + σ₀₁ Gᵐ₁) A
///     // dᵐ₁ = (σ₁₀ Gᵐ₀ + σ₁₁ Gᵐ₁) A
///     //     ┌       ┐
///     //     │ -¼ -⅓ │
///     // G = │  ¼  0 │
///     //     │  0  ⅓ │
///     //     └       ┘
///     assert_vec_approx_eq!(d.as_data(), &[-7.5, -8.5, 1.5, 4.5, 6.0, 4.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_d_tensor_dot_gradient<F>(
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
    let space_ndim = state.coords_min.len();
    if d.dim() != shape.nnode * space_ndim {
        return Err("d.len() must be equal to nnode * space_ndim");
    }

    // allocate auxiliary tensor
    let mut sig = Tensor2::new(true, space_ndim == 2);

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
        if space_ndim == 2 {
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
    use super::{
        vec_a_shape_times_scalar, vec_b_shape_times_vector, vec_c_vector_dot_gradient, vec_d_tensor_dot_gradient,
    };
    use crate::integ::{select_integ_points, AnalyticalTet4, AnalyticalTri3};
    use crate::shapes::{GeoKind, Shape, StateOfShape};
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;

    // to test if variables are cleared before sum
    const NOISE: f64 = 1234.56;

    #[test]
    fn capture_some_errors() {
        let shape = Shape::new(GeoKind::Lin2);
        let mut state = StateOfShape::new(shape.kind, &[[0.0, 0.0], [1.0, 0.0]]).unwrap();
        let mut a = Vector::new(3);
        assert_eq!(
            vec_a_shape_times_scalar(&mut a, &mut state, &shape, &[], 1.0, false, |_| Ok(0.0)).err(),
            Some("a.len() must be equal to nnode")
        );
        let mut b = Vector::new(5);
        assert_eq!(
            vec_b_shape_times_vector(&mut b, &mut state, &shape, &[], 1.0, false, |_, _| Ok(())).err(),
            Some("b.len() must be equal to nnode * space_ndim")
        );
        let mut c = Vector::new(3);
        assert_eq!(
            vec_c_vector_dot_gradient(&mut c, &mut state, &shape, &[], 1.0, false, |_, _| Ok(())).err(),
            Some("c.len() must be equal to nnode")
        );
        let mut d = Vector::new(5);
        assert_eq!(
            vec_d_tensor_dot_gradient(&mut d, &mut state, &shape, &[], 1.0, false, |_, _| Ok(())).err(),
            Some("d.len() must be equal to nnode * space_ndim")
        );
    }

    #[test]
    fn vec_a_shape_times_scalar_works_lin2_linear() -> Result<(), StrError> {
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
        let shape = Shape::new(GeoKind::Lin2);
        let mut state = StateOfShape::new(shape.kind, &[[3.0, 4.0], [3.0 + L, 4.0]]).unwrap();
        let cf = L / 6.0;
        let (xa, xb) = (state.coords_transp[0][0], state.coords_transp[0][1]);
        let a_correct = &[cf * (2.0 * xa + xb), cf * (xa + 2.0 * xb)];
        // integration points
        let tolerances = [1e-15, 1e-14, 1e-15, 1e-15];
        let selection: Vec<_> = [2, 3, 4, 5]
            .iter()
            .map(|n| select_integ_points(shape.kind.class(), *n).unwrap())
            .collect();
        // check
        let mut a = Vector::filled(shape.nnode, NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = shape.calc_integ_points_coords(&mut state, ips).unwrap();
            vec_a_shape_times_scalar(&mut a, &mut state, &shape, ips, 1.0, true, |p| Ok(x_ips[p][0])).unwrap();
            assert_vec_approx_eq!(a.as_data(), a_correct, tol);
        });
        Ok(())
    }

    #[test]
    fn vec_a_shape_times_scalar_works_tri3_constant() -> Result<(), StrError> {
        // tri3 with a constant source term s(x) = cₛ
        const CS: f64 = 3.0;
        let shape = Shape::new(GeoKind::Tri3);
        let mut state = StateOfShape::new(shape.kind, &[[3.0, 4.0], [8.0, 4.0], [5.0, 9.0]]).unwrap();
        let ana = AnalyticalTri3::new(&shape, &state);
        let a_correct = ana.integ_vec_a_constant(CS);
        // integration points
        let tolerances = [1e-14, 1e-14, 1e-15, 1e-14, 1e-13, 1e-14];
        let selection: Vec<_> = [1, 3, 1_003, 4, 12, 16]
            .iter()
            .map(|n| select_integ_points(shape.kind.class(), *n).unwrap())
            .collect();
        // check
        let mut a = Vector::filled(shape.nnode, NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            vec_a_shape_times_scalar(&mut a, &mut state, &shape, ips, 1.0, true, |_| Ok(CS)).unwrap();
            assert_vec_approx_eq!(a.as_data(), a_correct, tol);
        });
        Ok(())
    }

    #[test]
    fn vec_a_shape_times_scalar_works_tet4_linear() -> Result<(), StrError> {
        // tet 4 with a linear source term s(x) = z = x₂
        let shape = Shape::new(GeoKind::Tet4);
        let mut state = StateOfShape::new(
            shape.kind,
            &[[2.0, 3.0, 4.0], [6.0, 3.0, 2.0], [2.0, 5.0, 1.0], [4.0, 3.0, 6.0]],
        )?;
        let ana = AnalyticalTet4::new(&shape, &state);
        let a_correct = ana.integ_vec_a_linear_along_z(&state);
        // integration points
        // Note that the tolerance is high for n_integ_point = 1
        // because the numerical integration performs poorly with few IPs
        let tolerances = [0.56, 1e-15, 1e-14, 1e-15, 1e-15, 1e-15];
        let selection: Vec<_> = [1, 4, 5, 8, 14]
            .iter()
            .map(|n| select_integ_points(shape.kind.class(), *n).unwrap())
            .collect();
        // check
        let mut a = Vector::filled(shape.nnode, NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = shape.calc_integ_points_coords(&mut state, ips).unwrap();
            vec_a_shape_times_scalar(&mut a, &mut state, &shape, ips, 1.0, true, |p| Ok(x_ips[p][2])).unwrap();
            assert_vec_approx_eq!(a.as_data(), a_correct, tol);
        });
        Ok(())
    }

    #[test]
    fn vec_b_shape_times_vector_works_lin2_linear() -> Result<(), StrError> {
        // This test is similar to the shape_times_scalar with lin2
        const L: f64 = 6.0;
        let shape = Shape::new(GeoKind::Lin2);
        let mut state = StateOfShape::new(shape.kind, &[[3.0, 4.0], [3.0 + L, 4.0]]).unwrap();
        let cf = L / 6.0;
        let (xa, xb) = (state.coords_transp[0][0], state.coords_transp[0][1]);
        let b_correct = &[
            cf * (2.0 * xa + xb),
            cf * (2.0 * xa + xb),
            cf * (xa + 2.0 * xb),
            cf * (xa + 2.0 * xb),
        ];
        // integration points
        let tolerances = [1e-15, 1e-15];
        let selection: Vec<_> = [2, 3]
            .iter()
            .map(|n| select_integ_points(shape.kind.class(), *n).unwrap())
            .collect();
        // check
        let space_ndim = state.coords_min.len();
        let mut b = Vector::filled(shape.nnode * space_ndim, NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = shape.calc_integ_points_coords(&mut state, ips).unwrap();
            vec_b_shape_times_vector(&mut b, &mut state, &shape, ips, 1.0, true, |v, p| {
                v[0] = x_ips[p][0];
                v[1] = x_ips[p][0]; // << note use of x component here too
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(b.as_data(), b_correct, tol);
        });
        Ok(())
    }

    #[test]
    fn vec_b_shape_times_vector_works_tri3_constant() -> Result<(), StrError> {
        // This test is similar to the shape_times_scalar with tri3, however using a vector
        // So, each component of `b` equals `Fₛ`
        const V0: f64 = -3.0;
        const V1: f64 = 8.0;
        let shape = Shape::new(GeoKind::Tri3);
        let mut state = StateOfShape::new(shape.kind, &[[3.0, 4.0], [8.0, 4.0], [5.0, 9.0]]).unwrap();
        let ana = AnalyticalTri3::new(&shape, &state);
        let b_correct = ana.integ_vec_b_constant(V0, V1);
        // integration points
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [1, 3]
            .iter()
            .map(|n| select_integ_points(shape.kind.class(), *n).unwrap())
            .collect();
        // check
        let space_ndim = state.coords_min.len();
        let mut b = Vector::filled(shape.nnode * space_ndim, NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            vec_b_shape_times_vector(&mut b, &mut state, &shape, ips, 1.0, true, |v, _| {
                v[0] = V0;
                v[1] = V1;
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(b.as_data(), b_correct, tol);
        });
        Ok(())
    }

    #[test]
    fn vec_b_shape_times_vector_works_tet4_constant() -> Result<(), StrError> {
        // tet 4 with constant vector
        const V0: f64 = 2.0;
        const V1: f64 = 3.0;
        const V2: f64 = 4.0;
        let shape = Shape::new(GeoKind::Tet4);
        let mut state = StateOfShape::new(
            shape.kind,
            &[[2.0, 3.0, 4.0], [6.0, 3.0, 2.0], [2.0, 5.0, 1.0], [4.0, 3.0, 6.0]],
        )?;
        let ana = AnalyticalTet4::new(&shape, &state);
        let b_correct = ana.integ_vec_b_constant(V0, V1, V2);
        // integration points
        let tolerances = [1e-15, 1e-15];
        let selection: Vec<_> = [1, 4]
            .iter()
            .map(|n| select_integ_points(shape.kind.class(), *n).unwrap())
            .collect();
        // check
        let space_ndim = state.coords_min.len();
        let mut b = Vector::filled(shape.nnode * space_ndim, NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            vec_b_shape_times_vector(&mut b, &mut state, &shape, ips, 1.0, true, |v, _| {
                v[0] = V0;
                v[1] = V1;
                v[2] = V2;
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(b.as_data(), b_correct, tol);
        });
        Ok(())
    }

    #[test]
    fn vec_c_vector_dot_gradient_works_tri3_constant() -> Result<(), StrError> {
        // constant vector function: w(x) = {w₀, w₁}
        const W0: f64 = 2.0;
        const W1: f64 = 3.0;
        let shape = Shape::new(GeoKind::Tri3);
        let mut state = StateOfShape::new(shape.kind, &[[3.0, 4.0], [8.0, 4.0], [5.0, 9.0]]).unwrap();
        let ana = AnalyticalTri3::new(&shape, &state);
        let c_correct = ana.integ_vec_c_constant(W0, W1);
        // integration points
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [1, 3]
            .iter()
            .map(|n| select_integ_points(shape.kind.class(), *n).unwrap())
            .collect();
        // check
        let mut c = Vector::filled(shape.nnode, NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            vec_c_vector_dot_gradient(&mut c, &mut state, &shape, ips, 1.0, true, |w, _| {
                w[0] = W0;
                w[1] = W1;
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(c.as_data(), c_correct, tol);
        });
        Ok(())
    }

    #[test]
    fn vec_c_vector_dot_gradient_works_tri3_bilinear() -> Result<(), StrError> {
        // bilinear vector function: w(x) = {x, y}
        let shape = Shape::new(GeoKind::Tri3);
        let mut state = StateOfShape::new(shape.kind, &[[3.0, 4.0], [8.0, 4.0], [5.0, 9.0]]).unwrap();
        let ana = AnalyticalTri3::new(&shape, &state);
        let c_correct = ana.integ_vec_c_bilinear(&state);
        // integration points
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [1, 3]
            .iter()
            .map(|n| select_integ_points(shape.kind.class(), *n).unwrap())
            .collect();
        // check
        let mut c = Vector::filled(shape.nnode, NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = shape.calc_integ_points_coords(&mut state, ips).unwrap();
            vec_c_vector_dot_gradient(&mut c, &mut state, &shape, ips, 1.0, true, |w, p| {
                w[0] = x_ips[p][0];
                w[1] = x_ips[p][1];
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(c.as_data(), c_correct, tol);
        });
        Ok(())
    }

    #[test]
    fn vec_c_vector_dot_gradient_works_tet4_constant() -> Result<(), StrError> {
        // tet 4 with constant vector  w(x) = {w0, w1, w2}
        const W0: f64 = 2.0;
        const W1: f64 = 3.0;
        const W2: f64 = 4.0;
        let shape = Shape::new(GeoKind::Tet4);
        let mut state = StateOfShape::new(
            shape.kind,
            &[[2.0, 3.0, 4.0], [6.0, 3.0, 2.0], [2.0, 5.0, 1.0], [4.0, 3.0, 6.0]],
        )?;
        let ana = AnalyticalTet4::new(&shape, &state);
        let c_correct = ana.integ_vec_c_constant(W0, W1, W2);
        // integration points
        let tolerances = [1e-14, 1e-14, 1e-14, 1e-14, 1e-14];
        let selection: Vec<_> = [1, 4, 5, 8, 14]
            .iter()
            .map(|n| select_integ_points(shape.kind.class(), *n).unwrap())
            .collect();
        // check
        let mut c = Vector::filled(shape.nnode, NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            vec_c_vector_dot_gradient(&mut c, &mut state, &shape, ips, 1.0, true, |w, _| {
                w[0] = W0;
                w[1] = W1;
                w[2] = W2;
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(c.as_data(), c_correct, tol);
        });
        Ok(())
    }

    #[test]
    fn vec_d_tensor_dot_gradient_tri3_works_constant() -> Result<(), StrError> {
        // constant tensor function: σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2}
        // solution:
        //    dᵐ₀ = ½ (σ₀₀ bₘ + σ₀₁ cₘ)
        //    dᵐ₁ = ½ (σ₁₀ bₘ + σ₁₁ cₘ)
        const S00: f64 = 2.0;
        const S11: f64 = 3.0;
        const S22: f64 = 4.0;
        const S01: f64 = 5.0;
        let shape = Shape::new(GeoKind::Tri3);
        let mut state = StateOfShape::new(shape.kind, &[[3.0, 4.0], [8.0, 4.0], [5.0, 9.0]]).unwrap();
        let ana = AnalyticalTri3::new(&shape, &state);
        let d_correct = ana.integ_vec_d_constant(S00, S11, S01);
        // integration points
        let tolerances = [1e-14, 1e-14, 1e-14, 1e-14, 1e-13, 1e-14];
        let selection: Vec<_> = [1, 3, 1_003, 4, 12, 16]
            .iter()
            .map(|n| select_integ_points(shape.kind.class(), *n).unwrap())
            .collect();
        // check
        let space_ndim = state.coords_min.len();
        let mut d = Vector::filled(shape.nnode * space_ndim, NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            vec_d_tensor_dot_gradient(&mut d, &mut state, &shape, ips, 1.0, true, |sig, _| {
                sig.sym_set(0, 0, S00);
                sig.sym_set(1, 1, S11);
                sig.sym_set(2, 2, S22);
                sig.sym_set(0, 1, S01);
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(d.as_data(), d_correct, tol);
        });
        Ok(())
    }

    #[test]
    fn vec_d_tensor_dot_gradient_tet4_works_constant() -> Result<(), StrError> {
        // constant tensor function: σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2, σ₁₂√2, σ₀₂√2}
        let shape = Shape::new(GeoKind::Tet4);
        let mut state = StateOfShape::new(
            shape.kind,
            &[[2.0, 3.0, 4.0], [6.0, 3.0, 2.0], [2.0, 5.0, 1.0], [4.0, 3.0, 6.0]],
        )?;
        const S00: f64 = 2.0;
        const S11: f64 = 3.0;
        const S22: f64 = 4.0;
        const S01: f64 = 5.0;
        const S12: f64 = 6.0;
        const S02: f64 = 7.0;
        let ana = AnalyticalTet4::new(&shape, &state);
        let d_correct = ana.integ_vec_d_constant(S00, S11, S22, S01, S12, S02);
        // integration points
        let tolerances = [1e-14, 1e-14, 1e-13, 1e-14, 1e-14];
        let selection: Vec<_> = [1, 4, 5, 8, 14]
            .iter()
            .map(|n| select_integ_points(shape.kind.class(), *n).unwrap())
            .collect();
        // check
        let space_ndim = state.coords_min.len();
        let mut d = Vector::filled(shape.nnode * space_ndim, NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            vec_d_tensor_dot_gradient(&mut d, &mut state, &shape, ips, 1.0, true, |sig, _| {
                sig.sym_set(0, 0, S00);
                sig.sym_set(1, 1, S11);
                sig.sym_set(2, 2, S22);
                sig.sym_set(0, 1, S01);
                sig.sym_set(1, 2, S12);
                sig.sym_set(0, 2, S02);
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(d.as_data(), d_correct, tol);
        });
        Ok(())
    }
}
