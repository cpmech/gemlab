use crate::integ::Gauss;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_tensor::Tensor1;

/// Calculates the x-y-z coordinates of all integration points
///
/// This function applies the isoparametric formula
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
/// * `integ_points` -- Integration points' constants (ngauss)
///
/// # Output
///
/// * Returns an array with `ngauss` (number of integration points) vectors, where
///   each vector has a dimension equal to `space_ndim`.
///
/// # Examples
///
/// ```
/// use gemlab::integ::Gauss;
/// use gemlab::recovery::get_points_coords;
/// use gemlab::shapes::{GeoKind, Scratchpad};
/// use gemlab::StrError;
/// use russell_tensor::{Tensor1, t1_approx_eq};
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
///     let gauss = Gauss::new_sized(pad.kind.class(), 3)?;
///     let x_ips = get_points_coords(&mut pad, &gauss)?;
///
///     // expected
///     let x_ref0 = Tensor1::from(&[1.0, 1.0, 0.0]);
///     let x_ref1 = Tensor1::from(&[4.0, 1.0, 0.0]);
///     let x_ref2 = Tensor1::from(&[1.0, 4.0, 0.0]);
///
///     t1_approx_eq(&x_ips[0], &x_ref0, 1e-15);
///     t1_approx_eq(&x_ips[1], &x_ref1, 1e-15);
///     t1_approx_eq(&x_ips[2], &x_ref2, 1e-15);
///     Ok(())
/// }
/// ```
pub fn get_points_coords(pad: &mut Scratchpad, gauss: &Gauss) -> Result<Vec<Tensor1>, StrError> {
    let mut all_coords = Vec::new();
    let ngauss = gauss.npoint();
    for p in 0..ngauss {
        let mut x = Tensor1::new();
        pad.calc_coords(&mut x, gauss.coords(p))?;
        all_coords.push(x);
    }
    Ok(all_coords)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::get_points_coords;
    use crate::integ::Gauss;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::approx_eq;

    #[test]
    pub fn points_coords_works() {
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
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();

        approx_eq(x_ips[0].get(0), w * (1.0 - f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(x_ips[0].get(1), h * (1.0 - f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(x_ips[1].get(0), w * (1.0 + f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(x_ips[1].get(1), h * (1.0 - f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(x_ips[2].get(0), w * (1.0 - f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(x_ips[2].get(1), h * (1.0 + f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(x_ips[3].get(0), w * (1.0 + f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
        approx_eq(x_ips[3].get(1), h * (1.0 + f64::sqrt(3.0) / 3.0) / 2.0, 1e-15);
    }
}
