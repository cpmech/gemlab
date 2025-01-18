use crate::integ::Gauss;
use crate::shapes::Scratchpad;
use russell_lab::Matrix;

/// Calculates the interpolation matrix (nodes to integration points)
///
/// ```text
/// u_points =      P       u_nodal
///   (nip)    (nip,nnode)  (nnode)
/// ```
///
/// # Input
///
/// * `pad` -- The scratchpad
/// * `integ_points` -- Integration points' constants (n_integ_point)
///
/// # Output
///
/// * `P` -- The `(n_integ_point,nnode)` interpolation matrix. Returns a matrix
///   formed by the interpolation functions evaluated at all integration points.
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
///
/// # Examples
///
/// ```
/// use gemlab::integ::Gauss;
/// use gemlab::recovery::get_interp_matrix;
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
///     // nodal values
///     let u_nodal = Vector::from(&[1.0, 2.0, 3.0]);
///
///     // interpolated values
///     let gauss = Gauss::new_sized(pad.kind.class(), 3)?;
///     let mut u_points = Vector::new(gauss.data.len());
///     let pp = get_interp_matrix(&mut pad, &gauss);
///     mat_vec_mul(&mut u_points, 1.0, &pp, &u_nodal)?;
///
///     // check
///     let u_points_correct = Vector::from(&[1.5, 2.0, 2.5]);
///     vec_approx_eq(&u_points, &u_points_correct, 1e-14);
///     Ok(())
/// }
/// ```
pub fn get_interp_matrix(pad: &mut Scratchpad, gauss: &Gauss) -> Matrix {
    let nnode = pad.interp.dim();
    let n_integ_point = gauss.data.len();
    let mut pp = Matrix::new(n_integ_point, nnode);
    for i in 0..n_integ_point {
        (pad.fn_interp)(&mut pad.interp, &gauss.data[i]);
        for j in 0..nnode {
            pp.set(i, j, pad.interp[j]);
        }
    }
    pp
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::get_interp_matrix;
    use crate::integ::Gauss;
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

        let gauss = Gauss::new_sized(pad.kind.class(), 4).unwrap();
        let pp = get_interp_matrix(&mut pad, &gauss);
        assert_eq!(pp.dims(), (4, 4));

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

        let nip = gauss.data.len();
        let (ndim, nnode) = pad.xxt.dims();
        let mut xx_ips = Matrix::new(nip, ndim);
        for i in 0..nip {
            for j in 0..ndim {
                for k in 0..nnode {
                    xx_ips.add(i, j, pp.get(i, k) * pad.xxt.get(j, k));
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
