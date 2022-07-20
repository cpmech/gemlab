use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::Vector;

/// Calculates the x-y-z coordinates of all integration points
///
/// This function calls [op::calc_coords] to apply the isoparametric formula
/// to all p-th integration points `ιᵖ` according to:
///
/// ```text
/// → →          → →   →
/// x(ιᵖ) = Σ Nᵐ(ξ=ιᵖ) xᵐ
///         m         
/// ```
///
/// # Input
///
/// * `integ_points` -- Integration points' constants (n_integ_point)
///
/// # Output
///
/// * Returns an array with `n_integ_point` (number of integration points) vectors, where
///   each vector has a dimension equal to `space_ndim`.
///
/// # Examples
///
/// ```
/// use gemlab::integ::calc_ips_coords;
/// use gemlab::shapes::{GeoKind, Scratchpad};
/// use gemlab::StrError;
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
///     let x_ips = calc_ips_coords(&mut pad, &IP_TRI_INTERNAL_3)?;
///     assert_eq!(x_ips[0].as_data(), &[1.0, 1.0]);
///     assert_eq!(x_ips[1].as_data(), &[4.0, 1.0]);
///     assert_eq!(x_ips[2].as_data(), &[1.0, 4.0]);
///     Ok(())
/// }
///
/// ```
pub fn calc_ips_coords(pad: &mut Scratchpad, integ_points: &[[f64; 4]]) -> Result<Vec<Vector>, StrError> {
    let space_ndim = pad.xmax.len();
    let mut all_coords = Vec::new();
    for iota in integ_points {
        let mut x = Vector::new(space_ndim);
        pad.calc_coords(&mut x, iota)?;
        all_coords.push(x);
    }
    Ok(all_coords)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use russell_chk::assert_approx_eq;

    use super::calc_ips_coords;
    use crate::integ::IP_QUA_LEGENDRE_4;
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::StrError;

    #[test]
    pub fn calc_ips_coords_works() -> Result<(), StrError> {
        //  3-------------2         ξ₀   ξ₁
        //  | *    ξ₁   * |  node    r    s
        //  |      |      |     0 -1.0 -1.0
        //  |      +--ξ₀  |     1  1.0 -1.0
        //  |             |     2  1.0  1.0
        //  | *         * |     3 -1.0  1.0
        //  0-------------1

        let (w, h) = (20.0, 10.0);
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Qua4)?;
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, w);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(2, 0, w);
        pad.set_xx(2, 1, h);
        pad.set_xx(3, 0, 0.0);
        pad.set_xx(3, 1, h);

        let x_ips = calc_ips_coords(&mut pad, &IP_QUA_LEGENDRE_4)?;

        assert_approx_eq!(x_ips[0][0], w * (1.0 - f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        assert_approx_eq!(x_ips[0][1], h * (1.0 - f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        assert_approx_eq!(x_ips[1][0], w * (1.0 + f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        assert_approx_eq!(x_ips[1][1], h * (1.0 - f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        assert_approx_eq!(x_ips[2][0], w * (1.0 - f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        assert_approx_eq!(x_ips[2][1], h * (1.0 + f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        assert_approx_eq!(x_ips[3][0], w * (1.0 + f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        assert_approx_eq!(x_ips[3][1], h * (1.0 + f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        Ok(())
    }
}
