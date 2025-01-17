use super::get_interp_matrix;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{mat_inverse, mat_pseudo_inverse, Matrix};

/// Calculates the extrapolation matrix (integration points to nodes)
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
/// is only possible if `n_integ_points == nnode`. If there are more interpolation points
/// than nodes (`n_integ_points > node`), the problem is over-determined, and we use the
/// pseudo-inverse instead. If there are fewer interpolation points than nodes
/// (`n_integ_points < nnode`), the problem is under-determined and may even yield
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
        mat_pseudo_inverse(&mut xhi, &mut xh)?;

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
    use crate::integ::{IP_QUA_LEGENDRE_1, IP_QUA_LEGENDRE_4};
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::{mat_approx_eq, mat_vec_mul, vec_approx_eq, Matrix, Vector};

    #[test]
    pub fn get_extrap_matrix_works() {
        //  3-------------2         ξ₀   ξ₁
        //  | *    ξ₁   * |  node    r    s
        //  |      |      |     0 -1.0 -1.0
        //  |      +--ξ₀  |     1  1.0 -1.0
        //  |             |     2  1.0  1.0
        //  | *         * |     3 -1.0  1.0
        //  0-------------1

        let (w, h) = (20.0, 10.0);
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Qua4).unwrap();
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, w);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(2, 0, w);
        pad.set_xx(2, 1, h);
        pad.set_xx(3, 0, 0.0);
        pad.set_xx(3, 1, h);

        let ips = &IP_QUA_LEGENDRE_4;
        let ee = get_extrap_matrix(&mut pad, ips).unwrap();
        assert_eq!(ee.dims(), (4, 4));

        // original U values at nodes
        let u_nodes_original = Vector::from(&[1.0, 2.0, 3.0, 4.0]);

        // interpolated U values @ integration points
        let nip = ips.len();
        let nnode = 4;
        let mut u_ips = Vector::new(nip);
        for i in 0..nip {
            (pad.fn_interp)(&mut pad.interp, &ips[i]);
            for m in 0..nnode {
                u_ips[i] += pad.interp[m] * u_nodes_original[m];
            }
        }

        // extrapolated U values @ nodes points
        let mut u_nodes = Vector::new(nnode);
        mat_vec_mul(&mut u_nodes, 1.0, &ee, &u_ips).unwrap();

        // println!("E =\n{}", ee);
        // println!("U_ips =\n{}", u_ips);
        // println!("U_nodes =\n{}", u_nodes);

        // check
        vec_approx_eq(&u_nodes, &u_nodes_original, 1e-15);
    }

    #[test]
    pub fn get_extrap_matrix_works_1ip() {
        //  3-------------2         ξ₀   ξ₁
        //  |      ξ₁     |  node    r    s
        //  |      |      |     0 -1.0 -1.0
        //  |      *--ξ₀  |     1  1.0 -1.0
        //  |             |     2  1.0  1.0
        //  |             |     3 -1.0  1.0
        //  0-------------1

        let (w, h) = (20.0, 10.0);
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Qua4).unwrap();
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, w);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(2, 0, w);
        pad.set_xx(2, 1, h);
        pad.set_xx(3, 0, 0.0);
        pad.set_xx(3, 1, h);

        let ips = &IP_QUA_LEGENDRE_1;
        let ee = get_extrap_matrix(&mut pad, ips).unwrap();
        // println!("E =\n{}", ee);
        let ee_correct = Matrix::from(&[[1.0], [1.0], [1.0], [1.0]]); // the pseudo-inverse of P = [0.25, 0.25, 0.25, 0.25]
        mat_approx_eq(&ee, &ee_correct, 1e-15);
    }

    #[test]
    pub fn get_extrap_matrix_works_q8_4ip() {
        //  3------6------2         ξ₀   ξ₁
        //  |      ξ₁     |  node    r    s
        //  |      |      |     0 -1.0 -1.0   4  0.0 -1.0
        //  7      *--ξ₀  5     1  1.0 -1.0   5  1.0  0.0
        //  |             |     2  1.0  1.0   6  0.0  1.0
        //  |             |     3 -1.0  1.0   7 -1.0  0.0
        //  0------4------1

        let (w, h) = (20.0, 10.0);
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Qua8).unwrap();
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, w);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(2, 0, w);
        pad.set_xx(2, 1, h);
        pad.set_xx(3, 0, 0.0);
        pad.set_xx(3, 1, h);
        pad.set_xx(4, 0, w / 2.0);
        pad.set_xx(4, 1, 0.0);
        pad.set_xx(5, 0, w);
        pad.set_xx(5, 1, h / 2.0);
        pad.set_xx(6, 0, w / 2.0);
        pad.set_xx(6, 1, h);
        pad.set_xx(7, 0, 0.0);
        pad.set_xx(7, 1, h / 2.0);
        assert!(pad.ok_xxt);
        // pad.draw_shape_simple("/tmp/gemlab/test_get_extrap_matrix_works_q8_4ip.svg").unwrap();

        let ips = &IP_QUA_LEGENDRE_4;
        let ee = get_extrap_matrix(&mut pad, ips).unwrap();
        assert_eq!(ee.dims(), (8, 4));

        // original U values at nodes
        let u_nodes_original = Vector::from(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);

        // interpolated U values @ integration points
        let nip = ips.len();
        let nnode = 8;
        let mut u_ips = Vector::new(nip);
        for i in 0..nip {
            (pad.fn_interp)(&mut pad.interp, &ips[i]);
            for m in 0..nnode {
                u_ips[i] += pad.interp[m] * u_nodes_original[m];
            }
        }

        // extrapolated U values @ nodes points
        let mut u_nodes = Vector::new(nnode);
        mat_vec_mul(&mut u_nodes, 1.0, &ee, &u_ips).unwrap();

        println!("E =\n{}", ee);
        println!("U_ips =\n{}", u_ips);
        println!("U_nodes =\n{}", u_nodes);

        // check
        // vec_approx_eq(&u_nodes, &u_nodes_original, 1e-15);
    }
}
