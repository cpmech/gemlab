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
    let (nnode, geo_ndim) = pad.deriv.dims();
    let n_integ_point = integ_points.len();
    let mut ee = Matrix::new(nnode, n_integ_point);
    let mut pp = get_interp_matrix(pad, integ_points);
    // println!("P =\n{}", pp);
    if n_integ_point == nnode {
        mat_inverse(&mut ee, &pp)?;
    } else if n_integ_point > nnode {
        mat_pseudo_inverse(&mut ee, &mut pp)?;
    } else if n_integ_point == 1 {
        mat_pseudo_inverse(&mut ee, &mut pp)?;
    } else {
        // From Reference # 1:

        // ξ matrix (nnode,geo_ndim+1) with natural coordinates of nodes, augmented by a column with ones.
        // From Ref #1: ξ is a matrix containing the local (reference) coordinates of nodes (Eq. (31) of Ref #1).
        let mut ksi = Matrix::new(nnode, geo_ndim + 1);
        for m in 0..nnode {
            let r = pad.kind.reference_coords(m);
            for d in 0..geo_ndim {
                ksi.set(m, d, r[d]);
            }
            ksi.set(m, geo_ndim, 1.0);
        }

        // ξ_hat matrix (n_integ_point,geo_ndim+1) with the natural coordinates of the integration points, augmented by a column with ones.
        // From Ref #1: ξ_hat is a matrix containing the local (reference) coordinates of the sampling (integration) points (Eq. (30) of Ref #1)
        let mut ksi_hat = Matrix::new(n_integ_point, geo_ndim + 1);
        for i in 0..n_integ_point {
            for d in 0..geo_ndim {
                ksi_hat.set(i, d, integ_points[i][d]);
            }
            ksi_hat.set(i, geo_ndim, 1.0);
        }

        println!("ksi =\n{}", ksi);
        println!("ksi_hat =\n{}", ksi_hat);
        return Err("TODO");
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
}
