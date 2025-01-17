use crate::shapes::Scratchpad;
use russell_lab::Matrix;

/// Returns a matrix formed by the interpolation functions evaluated at all integration points
///
/// # Input
///
/// * `pad` -- The scratchpad
/// * `integ_points` -- Integration points' constants (n_integ_point)
///
/// # Output
///
/// * `M` -- The `(n_integ_point,nnode)` interpolation matrix
///
/// Possible use:
///
/// Since, for one integration point (isoparametric property):
///
/// ```text
/// → →          → →   →
/// x(ιᵖ) = Σ Nᵐ(ξ=ιᵖ) xᵐ
///         m
/// ```
///
/// Then, for all integration points:
///
/// ```text
///   X_ips    =      M           X
/// (nip,ndim)   (nip,nnode) (nnode,ndim)
/// ```
///
/// where `X_ips` are the (actual) coordinates of all integration points,
/// and `X` is a matrix with the coordinates of all vertices (`= transpose(xxt)`)
pub fn get_interp_matrix(pad: &mut Scratchpad, integ_points: &[[f64; 4]]) -> Matrix {
    let nnode = pad.interp.dim();
    let n_integ_point = integ_points.len();
    let mut mm = Matrix::new(n_integ_point, nnode);
    for i in 0..n_integ_point {
        (pad.fn_interp)(&mut pad.interp, &integ_points[i]);
        for j in 0..nnode {
            mm.set(i, j, pad.interp[j]);
        }
    }
    mm
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::get_interp_matrix;
    use crate::integ::IP_QUA_LEGENDRE_4;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::{approx_eq, Matrix};

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
        let mm = get_interp_matrix(&mut pad, ips);
        assert_eq!(mm.dims(), (4, 4));

        // For one integration point:
        //
        // → →          → →   →
        // x(ιᵖ) = Σ Nᵐ(ξ=ιᵖ) xᵐ
        //         m
        //
        // Thus, for all integration points:
        //
        //   X_ips    =      M           X
        // (nip,ndim)   (nip,nnode) (nnode,ndim)

        let nip = ips.len();
        let (ndim, nnode) = pad.xxt.dims();
        let mut xx_ips = Matrix::new(nip, ndim);
        for i in 0..nip {
            for j in 0..ndim {
                for k in 0..nnode {
                    xx_ips.add(i, j, mm.get(i, k) * pad.xxt.get(j, k));
                }
            }
        }

        // println!("{}", mm);
        // println!("{}", xx_ips);

        approx_eq(xx_ips.get(0, 0), w * (1.0 - f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(xx_ips.get(0, 1), h * (1.0 - f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(xx_ips.get(1, 0), w * (1.0 + f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(xx_ips.get(1, 1), h * (1.0 - f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(xx_ips.get(2, 0), w * (1.0 - f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(xx_ips.get(2, 1), h * (1.0 + f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(xx_ips.get(3, 0), w * (1.0 + f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(xx_ips.get(3, 1), h * (1.0 + f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
    }
}
