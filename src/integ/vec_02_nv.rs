use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::Vector;

/// Implements the the shape(N) times vector(V) integration case 02
///
/// Interpolation functions times vector field:
///
/// ```text
/// →    ⌠    → →   → →
/// bᵐ = │ Nᵐ(x(ξ)) v(x) dΩ
///      ⌡
///      Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
/// →    nip-1     →   → →       →
/// bᵐ ≈   Σ    Nᵐ(ιᵖ) v(ιᵖ) |J|(ιᵖ) wᵖ
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
/// * `fn_v` -- Function `f(v,p)` that computes `v(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`.
///   The dim of `v` is equal to `space_ndim`.
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
///     let mut b = Vector::filled(pad.kind.nnode() * space_ndim, 0.0);
///     integ::vec_b(&mut b, &mut pad, 0, true, ips, |v, _| {
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
pub fn vec_b<F>(
    b: &mut Vector,
    pad: &mut Scratchpad,
    ii0: usize,
    clear_b: bool,
    ips: IntegPointData,
    fn_v: F,
) -> Result<(), StrError>
where
    F: Fn(&mut Vector, usize) -> Result<(), StrError>,
{
    // check
    let (space_ndim, nnode) = pad.xxt.dims();
    if b.dim() < ii0 + nnode * space_ndim {
        return Err("b.len() must be ≥ ii0 + nnode ⋅ space_ndim");
    }

    // allocate auxiliary vector
    let mut v = Vector::new(space_ndim);

    // clear output vector
    if clear_b {
        b.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate interpolation functions and Jacobian
        (pad.fn_interp)(&mut pad.interp, iota);
        let det_jac = pad.calc_jacobian(iota)?;

        // calculate v
        fn_v(&mut v, p)?;

        // add contribution to b vector
        let coef = det_jac * weight;
        let nn = &pad.interp;
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
    use crate::integ::{self, AnalyticalTet4, AnalyticalTri3};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut b = Vector::new(4);
        assert_eq!(
            integ::vec_b(&mut b, &mut pad, 1, false, &[], |_, _| Ok(())).err(),
            Some("b.len() must be ≥ ii0 + nnode ⋅ space_ndim")
        );
    }

    #[test]
    fn vec_02_nv_works_lin2_linear() {
        // This test is similar to the shape_times_scalar with lin2
        const L: f64 = 6.0;
        let mut pad = aux::gen_pad_lin2(L);

        // solution
        let cf = L / 6.0;
        let (xa, xb) = (pad.xxt[0][0], pad.xxt[0][1]);
        let b_correct = &[
            cf * (2.0 * xa + xb),
            cf * (2.0 * xa + xb),
            cf * (xa + 2.0 * xb),
            cf * (xa + 2.0 * xb),
        ];

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-15, 1e-15];
        let selection: Vec<_> = [2, 3].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let (space_ndim, nnode) = pad.xxt.dims();
        let mut b = Vector::filled(nnode * space_ndim, aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = integ::points_coords(&mut pad, ips).unwrap();
            integ::vec_b(&mut b, &mut pad, 0, true, ips, |v, p| {
                v[0] = x_ips[p][0];
                v[1] = x_ips[p][0]; // << note use of x component here too
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(b.as_data(), b_correct, tol);
        });
    }

    #[test]
    fn vec_02_nv_works_tri3_constant() {
        // This test is similar to the shape_times_scalar with tri3, however using a vector
        // So, each component of `b` equals `Fₛ`
        let mut pad = aux::gen_pad_tri3();

        // solution
        let ana = AnalyticalTri3::new(&pad);
        const V0: f64 = -3.0;
        const V1: f64 = 8.0;
        let b_correct = ana.integ_vec_b(V0, V1);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [1, 3].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let (space_ndim, nnode) = pad.xxt.dims();
        let mut b = Vector::filled(nnode * space_ndim, aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::vec_b(&mut b, &mut pad, 0, true, ips, |v, _| {
                v[0] = V0;
                v[1] = V1;
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(b.as_data(), b_correct, tol);
        });
    }

    #[test]
    fn vec_02_nv_works_tet4_constant() {
        // tet 4 with constant vector
        const V0: f64 = 2.0;
        const V1: f64 = 3.0;
        const V2: f64 = 4.0;
        let mut pad = aux::gen_pad_tet4();

        // solution
        let ana = AnalyticalTet4::new(&pad);
        let b_correct = ana.integ_vec_b(V0, V1, V2);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-15, 1e-15];
        let selection: Vec<_> = [1, 4].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let (space_ndim, nnode) = pad.xxt.dims();
        let mut b = Vector::filled(nnode * space_ndim, aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::vec_b(&mut b, &mut pad, 0, true, ips, |v, _| {
                v[0] = V0;
                v[1] = V1;
                v[2] = V2;
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(b.as_data(), b_correct, tol);
        });
    }
}
