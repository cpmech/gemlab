use crate::StrError;

/// Interpolates a value at a given point over a triangle
/// 
/// See file data/derivations/interpolation-over-triangle.pdf
/// 
/// # Input
/// 
/// `x,y` -- are the coordinates of the point of interest
/// `xx,yy` -- are the coordinates of the triangle vertices (len = 3)
/// `temp` -- are the known values at the triangle vertices; e.g., the temperatures (len = 3)
/// 
/// # Output
/// 
/// Returns the "temperature" `T @ x`
#[rustfmt::skip]
pub fn tri_interpolation(x: f64, y: f64, xx: &[f64], yy: &[f64], temp: &[f64]) -> Result<f64, StrError> {
    if xx.len() != 3 || yy.len() != 3 || temp.len() != 3 {
        return Err("all arrays must have len = 3");
    }
    let aa2    =  xx[0] * (yy[1] - yy[2]) + xx[1] * (yy[2] - yy[0]) + xx[2] * (yy[0] - yy[1]);
    let zeta_0 = (x     * (yy[1] - yy[2]) + xx[1] * (yy[2] - y    ) + xx[2] * (y     - yy[1])) / aa2;
    let zeta_1 = (xx[0] * (y -     yy[2]) + x     * (yy[2] - yy[0]) + xx[2] * (yy[0] - y    )) / aa2;
    let zeta_2 = (xx[0] * (yy[1] -     y) + xx[1] * (y     - yy[0]) + x     * (yy[0] - yy[1])) / aa2;
    Ok(zeta_0 * temp[0] + zeta_1 * temp[1] + zeta_2 * temp[2])
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::tri_interpolation;
    use crate::StrError;

    #[test]
    fn tri_interpolation_handles_errors() {
        assert_eq!(
            tri_interpolation(0.0, 0.0, &[0.0], &[0.0, 0.0, 0.0], &[0.0, 0.0, 0.0]).err(),
            Some("all arrays must have len = 3")
        );
        assert_eq!(
            tri_interpolation(0.0, 0.0, &[0.0, 0.0, 0.0], &[0.0], &[0.0, 0.0, 0.0]).err(),
            Some("all arrays must have len = 3")
        );
        assert_eq!(
            tri_interpolation(0.0, 0.0, &[0.0, 0.0, 0.0], &[0.0, 0.0, 0.0], &[0.0]).err(),
            Some("all arrays must have len = 3")
        );
    }

    #[test]
    fn tri_interpolation_works() -> Result<(), StrError> {
        let xx = &[0.0, 1.5, 0.0];
        let yy = &[0.0, 0.0, 1.0];
        let temp = &[3.0, 2.0, 1.0];
        assert_eq!(tri_interpolation(xx[0], yy[0], xx, yy, temp)?, temp[0]);
        assert_eq!(tri_interpolation(xx[1], yy[1], xx, yy, temp)?, temp[1]);
        assert_eq!(tri_interpolation(xx[2], yy[2], xx, yy, temp)?, temp[2]);
        assert_eq!(tri_interpolation(1.0 / 2.0, 1.0 / 4.0, xx, yy, temp)?, 13.0 / 6.0);
        Ok(())
    }
}
