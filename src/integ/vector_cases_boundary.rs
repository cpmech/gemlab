use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::Vector;

/// Implements the the shape(N) times vector(V) integration case (boundary integral version)
///
/// Interpolation functions times vector field:
///
/// ```text
/// →    ⌠    → →   → →
/// bᵐ = │ Nᵐ(x(ξ)) t(x) tₕ dΓ
///      ⌡
///      Γₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
/// →    nip-1     →   → →        →   →
/// bᵐ ≈   Σ    Nᵐ(ιᵖ) t(ιᵖ) tₕ ||n||(ιᵖ) wᵖ
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
/// # Input
///
/// * `pad` -- **modified** Scratchpad
/// * `ips` -- Integration points (n_integ_point)
/// * `th` -- tₕ the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
/// * `clear_b` -- fills `b` vector with zeros, otherwise accumulate values into `b`
/// * `fn_t` -- Function `f(t,p,n)` that calculates `t(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
///    and given the normal vector `n(x(ιᵖ))`. The dim of `v` and `n` is equal to `space_ndim`.
///
/// # Examples
///
pub fn vec_b_shape_times_vector_boundary<F>(
    b: &mut Vector,
    pad: &mut Scratchpad,
    ips: IntegPointData,
    th: f64,
    clear_b: bool,
    fn_t: F,
) -> Result<(), StrError>
where
    F: Fn(&mut Vector, usize, &Vector) -> Result<(), StrError>,
{
    // check
    let nnode = pad.interp.dim();
    let space_ndim = pad.xmax.len();
    if b.dim() != nnode * space_ndim {
        return Err("b.len() must be equal to nnode * space_ndim");
    }

    // allocate auxiliary vectors
    let mut t = Vector::new(space_ndim);
    let mut n = Vector::new(space_ndim); // normal vector

    // clear output vector
    if clear_b {
        b.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate interpolation functions and normal vector
        (pad.fn_interp)(&mut pad.interp, iota);
        let mag_n = pad.calc_normal_vector(&mut n, iota)?;

        // calculate t
        fn_t(&mut t, p, &n)?;

        // add contribution to b vector
        let coef = th * mag_n * weight;
        let nn = &pad.interp;
        if space_ndim == 2 {
            for m in 0..nnode {
                b[0 + m * 2] += coef * nn[m] * t[0];
                b[1 + m * 2] += coef * nn[m] * t[1];
            }
        } else {
            for m in 0..nnode {
                b[0 + m * 3] += coef * nn[m] * t[0];
                b[1 + m * 3] += coef * nn[m] * t[1];
                b[2 + m * 3] += coef * nn[m] * t[2];
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::vec_b_shape_times_vector_boundary;
    use crate::integ::default_integ_points;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;

    // to test if variables are cleared before sum
    const NOISE: f64 = 1234.56;

    #[test]
    fn vec_b_shape_times_vector_boundary_works_2d() {
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Lin2).unwrap();
        let ll = 4.0;
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, ll);
        pad.set_xx(1, 1, 0.0);
        let mut b = Vector::filled(pad.kind.nnode() * space_ndim, NOISE);
        let ips = default_integ_points(pad.kind);
        vec_b_shape_times_vector_boundary(&mut b, &mut pad, ips, 1.0, true, |t, _, _| {
            t[0] = 0.0;
            t[1] = -1.0;
            Ok(())
        })
        .unwrap();
        assert_vec_approx_eq!(b.as_data(), &[0.0, -2.0, 0.0, -2.0], 1e-15);

        // example from [@sgm:14]\page{183}
        // * `sgm:14` -- Smith, Griffiths, Margetts (2014) Programming the Finite Element Method, 5th ed.
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Lin3).unwrap();
        let ll = 3.0;
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, ll);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(2, 0, ll / 2.0);
        pad.set_xx(2, 1, 0.0);
        let mut b = Vector::filled(pad.kind.nnode() * space_ndim, NOISE);
        let ips = default_integ_points(pad.kind);
        vec_b_shape_times_vector_boundary(&mut b, &mut pad, ips, 1.0, true, |t, _, _| {
            t[0] = 0.0;
            t[1] = -1.0;
            Ok(())
        })
        .unwrap();
        assert_vec_approx_eq!(b.as_data(), &[0.0, -0.5, 0.0, -0.5, 0.0, -2.0], 1e-15);
    }

    #[test]
    fn vec_b_shape_times_vector_boundary_works_3d() {
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
        let ips = default_integ_points(pad.kind);
        vec_b_shape_times_vector_boundary(&mut b, &mut pad, ips, 1.0, true, |t, _, _| {
            t[0] = 0.0;
            t[1] = 0.0;
            t[2] = -1.0;
            Ok(())
        })
        .unwrap();
        let aa = dx * dy;
        assert_vec_approx_eq!(
            b.as_data(),
            &[
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
            ],
            1e-15
        );

        // example from [@sgm:14]\page{195}
        // * `sgm:14` -- Smith, Griffiths, Margetts (2014) Programming the Finite Element Method, 5th ed.
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
        let ips = default_integ_points(pad.kind);
        vec_b_shape_times_vector_boundary(&mut b, &mut pad, ips, 1.0, true, |t, _, _| {
            t[0] = 0.0;
            t[1] = 0.0;
            t[2] = -1.0;
            Ok(())
        })
        .unwrap();
        assert_vec_approx_eq!(
            b.as_data(),
            &[
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
            ],
            1e-15
        );
    }
}
