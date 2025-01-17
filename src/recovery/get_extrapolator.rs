#![allow(unused)]

use super::get_interp_matrix;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{mat_inverse, mat_pseudo_inverse, Matrix};

/// Calculates the extrapolator matrix (integration points to nodes)
///
/// # Input
///
/// * `pad` -- The scratchpad
/// * `integ_points` -- Integration points' constants (n_integ_point)
///
/// # Output
///
/// * `E` -- The (nnode,n_integ_point) extrapolator matrix
pub fn get_extrapolator(pad: &mut Scratchpad, integ_points: &[[f64; 4]]) -> Result<Matrix, StrError> {
    let nnode = pad.interp.dim();
    let n_integ_point = integ_points.len();
    let mut ee = Matrix::new(nnode, n_integ_point);
    let mut mm = get_interp_matrix(pad, integ_points);
    if n_integ_point < nnode {
        return Err("TODO");
    } else if n_integ_point == nnode {
        mat_inverse(&mut ee, &mm)?;
    } else {
        mat_pseudo_inverse(&mut ee, &mut mm)?;
    }
    return Ok(ee);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::get_extrapolator;
    use crate::integ::IP_QUA_LEGENDRE_4;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::{approx_eq, mat_approx_eq, mat_vec_mul, vec_approx_eq, Matrix, Vector};

    #[test]
    pub fn get_shape_matrix_works() {
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
        let ee = get_extrapolator(&mut pad, ips).unwrap();
        assert_eq!(ee.dims(), (4, 4));

        // U values @ integration points
        let nip = ips.len();
        let (_, nnode) = pad.xxt.dims();
        let u_ips = Vector::from(&[4.0, 4.0, 4.0, 4.0]);

        // extrapolated U values @ nodes points
        let mut u_nodes = Vector::new(nnode);
        mat_vec_mul(&mut u_nodes, 1.0, &ee, &u_ips).unwrap();

        // println!("E =\n{}", ee);
        // println!("U_ips =\n{}", u_ips);
        // println!("U_nodes =\n{}", u_nodes);

        vec_approx_eq(&u_ips, &u_nodes, 1e-14);
    }
}
