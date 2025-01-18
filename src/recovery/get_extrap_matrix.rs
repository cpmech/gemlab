use super::get_interp_matrix;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{mat_inverse, mat_pseudo_inverse, Matrix};

/// Calculates the extrapolation matrix (integration points to nodes)
///
/// ```text
/// u_nodal =      E       u_points
/// (nnode)   (nnode,nip)    (nip)
/// ```
///
/// # Input
///
/// * `pad` -- The scratchpad
/// * `integ_points` -- Integration points' constants (n_integ_point)
///
/// # Output
///
/// * `E` -- The (nnode,n_integ_point) extrapolation matrix; aka "inverse" of the
///    interpolation matrix `P` calculated by [get_interp_matrix()]. See note below
///    regarding the inversion problem.
///
/// This function returns the "inverse" of the interpolation matrix; however, the inverse
/// is only possible if `n_integ_point == nnode`. If there are more interpolation points
/// than nodes (`n_integ_point > node`), the problem is over-determined, and we use the
/// pseudo-inverse instead. If there are fewer interpolation points than nodes
/// (`n_integ_point < nnode`), the problem is under-determined and may even yield
/// spurious results. To avoid such results and minimize the errors, we adopt the method
/// proposed by Durand and Farias in Reference #1. In this case, a correction matrix is
/// applied to the pseudo-inverse, and a translation is applied to the result.
/// An exception exists: with a single integration point, the pseudo-inverse is returned
/// without any correction.
///
/// # Reference
///
/// 1. Durand R and Farias MM (2014) A local extrapolation method for finite elements,
///    Advances in Engineering Software, 67:1-9 <https://doi.org/10.1016/j.advengsoft.2013.07.002>
///
/// # Examples
///
/// ```
/// use gemlab::recovery::{get_extrap_matrix, get_interp_matrix};
/// use gemlab::shapes::{GeoKind, Scratchpad};
/// use gemlab::StrError;
/// use russell_lab::{mat_vec_mul, vec_approx_eq, Vector};
///
/// fn main() -> Result<(), StrError> {
///     //  6 2
///     //  5 | `.    * indicates the
///     //  4 | * `.    location of ips
///     //  3 |     `.
///     //  2 |       `.
///     //  1 | *     * `.
///     //  0 0-----------1
///     //    0 1 2 3 4 5 6
///
///     let space_ndim = 2;
///     let mut pad = Scratchpad::new(space_ndim, GeoKind::Tri3)?;
///     pad.set_xx(0, 0, 0.0);
///     pad.set_xx(0, 1, 0.0);
///     pad.set_xx(1, 0, 6.0);
///     pad.set_xx(1, 1, 0.0);
///     pad.set_xx(2, 0, 0.0);
///     pad.set_xx(2, 1, 6.0);
///
///     // the last column of the array below contains the weight
///     const IP_TRI_INTERNAL_3: [[f64; 4]; 3] = [
///         [1.0 / 6.0, 1.0 / 6.0, 0.0, 1.0 / 6.0],
///         [2.0 / 3.0, 1.0 / 6.0, 0.0, 1.0 / 6.0],
///         [1.0 / 6.0, 2.0 / 3.0, 0.0, 1.0 / 6.0],
///     ];
///
///     // nodal values
///     let u_nodal = Vector::from(&[1.0, 2.0, 3.0]);
///
///     // interpolated values
///     let mut u_points = Vector::new(IP_TRI_INTERNAL_3.len());
///     let pp = get_interp_matrix(&mut pad, &IP_TRI_INTERNAL_3);
///     mat_vec_mul(&mut u_points, 1.0, &pp, &u_nodal)?;
///
///     // extrapolated values (recovered)
///     let nnode = pad.xxt.dims().1;
///     let mut u_nodal_rec = Vector::new(nnode);
///     let ee = get_extrap_matrix(&mut pad, &IP_TRI_INTERNAL_3)?;
///     mat_vec_mul(&mut u_nodal_rec, 1.0, &ee, &u_points)?;
///
///     // check
///     vec_approx_eq(&u_nodal_rec, &u_nodal, 1e-14);
///     Ok(())
/// }
/// ```
pub fn get_extrap_matrix(pad: &mut Scratchpad, integ_points: &[[f64; 4]]) -> Result<Matrix, StrError> {
    // constants
    let (nnode, geo_ndim) = pad.deriv.dims();
    let n_integ_point = integ_points.len();

    // calculate interpolation matrix P
    let mut pp = get_interp_matrix(pad, integ_points);

    // allocate output
    let mut ee = Matrix::new(nnode, n_integ_point);

    // handle cases
    if n_integ_point == nnode {
        // E ← P⁻¹
        mat_inverse(&mut ee, &pp)?;
    } else if n_integ_point > nnode {
        // E ← P⁺
        mat_pseudo_inverse(&mut ee, &mut pp)?;
    } else if n_integ_point == 1 {
        // E ← P⁺
        mat_pseudo_inverse(&mut ee, &mut pp)?;
    } else {
        // find pseudo-inverse of interpolation matrix, P⁺
        let mut ppi = Matrix::new(nnode, n_integ_point);
        mat_pseudo_inverse(&mut ppi, &mut pp)?;

        // ξ matrix (nnode,geo_ndim+1) with natural coordinates of nodes, augmented by a column of ones.
        // From Ref #1: ξ is a matrix containing the local coordinates of nodes (Eq. (31) of Ref #1).
        let mut x = Matrix::new(nnode, geo_ndim + 1);
        for m in 0..nnode {
            let r = pad.kind.reference_coords(m);
            for d in 0..geo_ndim {
                x.set(m, d, r[d]);
            }
            x.set(m, geo_ndim, 1.0);
        }

        // ξ_hat matrix (n_integ_point,geo_ndim+1) with the natural coordinates of the integration points, augmented by a column of ones.
        // From Ref #1: ξ_hat is a matrix containing the local coordinates of the sampling (integration) points (Eq. (30) of Ref #1)
        let mut xh = Matrix::new(n_integ_point, geo_ndim + 1);
        for p in 0..n_integ_point {
            for d in 0..geo_ndim {
                xh.set(p, d, integ_points[p][d]);
            }
            xh.set(p, geo_ndim, 1.0);
        }

        // calculate ξ_hat_inv (geo_ndim+1,n_integ_point), the pseudo-inverse of ξ_hat
        let mut xhi = Matrix::new(geo_ndim + 1, n_integ_point);
        let mut copy = xh.clone();
        mat_pseudo_inverse(&mut xhi, &mut copy)?;

        // auxiliary computations
        let mut aux = Matrix::new(n_integ_point, n_integ_point);
        for d in 0..geo_ndim + 1 {
            for q in 0..n_integ_point {
                // aux ← ξ_hat * ξ_hat_inv
                for p in 0..n_integ_point {
                    aux.add(p, q, xh.get(p, d) * xhi.get(d, q));
                }
                // E ← ξ * ξ_hat_inv
                for m in 0..nnode {
                    ee.add(m, q, x.get(m, d) * xhi.get(d, q));
                }
            }
        }

        // extrapolation matrix
        for m in 0..nnode {
            for q in 0..n_integ_point {
                for p in 0..n_integ_point {
                    let c = if q == p { 1.0 } else { 0.0 };
                    ee.add(m, q, ppi.get(m, p) * (c - aux.get(p, q)));
                }
            }
        }
    }
    return Ok(ee);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::get_extrap_matrix;
    use crate::integ::{IP_QUA_LEGENDRE_1, IP_QUA_LEGENDRE_4, IP_QUA_LEGENDRE_9};
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::math::PI;
    use russell_lab::{mat_approx_eq, mat_vec_mul, vec_approx_eq, Matrix, Vector};

    fn gen_qua4(w: f64, h: f64, skew_angle: f64, save_fig: bool) -> Scratchpad {
        //  3-------------2         ξ₀   ξ₁
        //  |      ξ₁     |  node    r    s
        //  |      |      |     0 -1.0 -1.0
        //  |      +--ξ₀  |     1  1.0 -1.0
        //  |             |     2  1.0  1.0
        //  |             |     3 -1.0  1.0
        //  0-------------1
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Qua4).unwrap();
        let dx = h * f64::tan(skew_angle);
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, w);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(2, 0, w + dx);
        pad.set_xx(2, 1, h);
        pad.set_xx(3, 0, 0.0 + dx);
        pad.set_xx(3, 1, h);
        assert!(pad.ok_xxt);
        if save_fig {
            pad.draw_shape_simple("/tmp/gemlab/gen_qua4.svg").unwrap();
        }
        pad
    }

    fn gen_qua8(w: f64, h: f64, skew_angle: f64, imprecision: f64, save_fig: bool) -> Scratchpad {
        //  3------6------2         ξ₀   ξ₁
        //  |      ξ₁     |  node    r    s
        //  |      |      |     0 -1.0 -1.0   4  0.0 -1.0
        //  7      *--ξ₀  5     1  1.0 -1.0   5  1.0  0.0
        //  |             |     2  1.0  1.0   6  0.0  1.0
        //  |             |     3 -1.0  1.0   7 -1.0  0.0
        //  0------4------1
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Qua8).unwrap();
        let dx = h * f64::tan(skew_angle);
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, w);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(2, 0, w + dx);
        pad.set_xx(2, 1, h);
        pad.set_xx(3, 0, 0.0 + dx);
        pad.set_xx(3, 1, h);
        pad.set_xx(4, 0, w / 2.0);
        pad.set_xx(4, 1, 0.0 + imprecision);
        pad.set_xx(5, 0, w + dx / 2.0 - imprecision);
        pad.set_xx(5, 1, h / 2.0);
        pad.set_xx(6, 0, w / 2.0 + dx - imprecision);
        pad.set_xx(6, 1, h - imprecision);
        pad.set_xx(7, 0, 0.0 + dx / 2.0 + imprecision);
        pad.set_xx(7, 1, h / 2.0);
        assert!(pad.ok_xxt);
        if save_fig {
            pad.draw_shape_simple("/tmp/gemlab/gen_qua8.svg").unwrap();
        }
        pad
    }

    fn do_interpolate(pad: &mut Scratchpad, u_nodal: &Vector, integ_points: &[[f64; 4]]) -> Vector {
        let nnode = u_nodal.dim();
        let n_integ_point = integ_points.len();
        let mut u_point = Vector::new(n_integ_point);
        for i in 0..n_integ_point {
            (pad.fn_interp)(&mut pad.interp, &integ_points[i]);
            for m in 0..nnode {
                u_point[i] += pad.interp[m] * u_nodal[m];
            }
        }
        u_point
    }

    fn do_extrapolate(ee: &Matrix, u_point: &Vector) -> Vector {
        let nnode = ee.dims().0;
        let mut u_nodal = Vector::new(nnode);
        mat_vec_mul(&mut u_nodal, 1.0, &ee, &u_point).unwrap();
        u_nodal
    }

    #[test]
    pub fn get_extrap_matrix_works_qua4_ip4() {
        let mut pad = gen_qua4(20.0, 10.0, PI / 6.0, false);

        // extrapolation matrix
        let ips = &IP_QUA_LEGENDRE_4;
        let ee = get_extrap_matrix(&mut pad, ips).unwrap();

        // check
        let u_nodal_original = Vector::from(&[1.0, 2.0, 3.0, 4.0]); // original U values at nodes
        let u_point = do_interpolate(&mut pad, &u_nodal_original, ips); // interpolated U values @ integration points
        let u_nodal = do_extrapolate(&ee, &u_point);
        vec_approx_eq(&u_nodal, &u_nodal_original, 1e-15);
    }

    #[test]
    pub fn get_extrap_matrix_works_qua4_ip1() {
        let mut pad = gen_qua4(20.0, 10.0, PI / 6.0, false);

        // extrapolation matrix
        let ips = &IP_QUA_LEGENDRE_1;
        let ee = get_extrap_matrix(&mut pad, ips).unwrap();
        let ee_correct = Matrix::from(&[[1.0], [1.0], [1.0], [1.0]]); // the pseudo-inverse of P = [0.25, 0.25, 0.25, 0.25]
        mat_approx_eq(&ee, &ee_correct, 1e-15);

        // check
        let u_nodal_original = Vector::from(&[1.0, 2.0, 3.0, 4.0]); // original U values at nodes
        let u_point = do_interpolate(&mut pad, &u_nodal_original, ips); // interpolated U values @ integration points
        let u_nodal = do_extrapolate(&ee, &u_point);
        // println!("u_point =\n{}", u_point);
        // println!("u_nodal =\n{}", u_nodal);
        vec_approx_eq(&u_nodal, &u_nodal_original, 1.5); // we cannot get a better result with such low n_integ_point
    }

    #[test]
    pub fn get_extrap_matrix_works_qua8_ip9() {
        let mut pad = gen_qua8(20.0, 10.0, PI / 6.0, 1.0, false);

        // extrapolation matrix
        let ips = &IP_QUA_LEGENDRE_9;
        let ee = get_extrap_matrix(&mut pad, ips).unwrap();

        // check
        let u_nodal_original = Vector::from(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]); // original U values at nodes
        let u_point = do_interpolate(&mut pad, &u_nodal_original, ips); // interpolated U values @ integration points
        let u_nodal = do_extrapolate(&ee, &u_point);
        // println!("u_point =\n{}", u_point);
        // println!("u_nodal =\n{}", u_nodal);
        vec_approx_eq(&u_nodal, &u_nodal_original, 1e-13); // pretty good, with overdetermined system
    }

    #[test]
    pub fn get_extrap_matrix_works_qua8_ip4() {
        let mut pad = gen_qua8(20.0, 10.0, 0.0, 0.0, false);

        // extrapolation matrix
        let ips = &IP_QUA_LEGENDRE_4;
        let ee = get_extrap_matrix(&mut pad, ips).unwrap();

        // check
        let u_nodal_original = Vector::from(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]); // original U values at nodes
        let u_point = do_interpolate(&mut pad, &u_nodal_original, ips); // interpolated U values @ integration points
        let u_nodal = do_extrapolate(&ee, &u_point);
        // println!("u_point =\n{}", u_point);
        // println!("u_nodal =\n{}", u_nodal);
        vec_approx_eq(&u_nodal, &u_nodal_original, 12.3); // we cannot get better results with only 4 integ points
    }

    #[test]
    pub fn get_extrap_matrix_works_durand_farias_example1() {
        let mut pad = gen_qua8(1.0, 1.0, 0.0, 0.0, false);

        // extrapolation matrix
        let ips = &IP_QUA_LEGENDRE_4;
        let ee = get_extrap_matrix(&mut pad, ips).unwrap();

        // check
        let u_point = Vector::from(&[0.5, 0.5, 0.5, 0.5]);
        let u_nodal = do_extrapolate(&ee, &u_point);
        // println!("u_point =\n{}", u_point);
        // println!("u_nodal =\n{}", u_nodal);

        // let u_nodal_without_correction = Vector::from(&[-0.088, -0.088, -0.088, -0.088, 0.353, 0.353, 0.353, 0.353]);
        // vec_approx_eq(&u_nodal, &u_nodal_without_correction, 1e-3);

        let u_nodal_expected = Vector::from(&[0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]);
        vec_approx_eq(&u_nodal, &u_nodal_expected, 1e-14);
    }
}
