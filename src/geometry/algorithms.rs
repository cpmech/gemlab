use crate::StrError;

/// Computes the point-to-point distance
pub fn point_point_distance(a: &[f64], b: &[f64]) -> Result<f64, StrError> {
    let ndim = a.len();
    if ndim < 2 || ndim > 3 {
        return Err("a.len() == ndim must be 2 or 3");
    }
    if b.len() != ndim {
        return Err("b.len() must equal a.len() == ndim");
    }
    let distance: f64;
    if ndim == 2 {
        distance = f64::sqrt((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]));
    } else {
        distance =
            f64::sqrt((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]) + (b[2] - a[2]) * (b[2] - a[2]));
    }
    Ok(distance)
}

/// Computes the distance between a point and a segment
///
/// ```text
///                 c
///                /^
/// ⇀             / | ⇀   ⇀   ⇀
/// x = vec(c-a) /  | q = x - p
///             /   |       
///            /    |
///           a----->-------------b
///             ⇀              ⇀
///             p = projection(x)
///
///           a-------------------b
///               ⇀
///               n = vec(b-c)
/// ```
///
/// ```text
/// x := c - a
/// n := b - a
/// p := ((x ⋅ n) / (n ⋅ n)) * n
/// q := x - p
/// distance := norm(q)
/// ```
pub fn point_segment_distance(a: &[f64], b: &[f64], c: &[f64]) -> Result<f64, StrError> {
    let ndim = a.len();
    if ndim < 2 || ndim > 3 {
        return Err("a.len() == ndim must be 2 or 3");
    }
    if b.len() != ndim {
        return Err("b.len() must equal a.len() == ndim");
    }
    if c.len() != ndim {
        return Err("c.len() must equal a.len() == ndim");
    }
    let distance: f64;
    if ndim == 2 {
        // n ⋅ n  with n = b - a
        let n_dot_n = (b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]);
        if n_dot_n <= f64::EPSILON {
            return Err("segment is too short");
        }
        // x ⋅ n  with x = c - a
        let x_dot_n = (c[0] - a[0]) * (b[0] - a[0]) + (c[1] - a[1]) * (b[1] - a[1]);
        // q = x - p = (c-a) - (b-a) * (x⋅n)/(n⋅n)
        let q0 = (c[0] - a[0]) - (b[0] - a[0]) * x_dot_n / n_dot_n;
        let q1 = (c[1] - a[1]) - (b[1] - a[1]) * x_dot_n / n_dot_n;
        // norm(q)
        distance = f64::sqrt(q0 * q0 + q1 * q1);
    } else {
        // n ⋅ n  with n = b - a
        let n_dot_n = (b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]) + (b[2] - a[2]) * (b[2] - a[2]);
        if n_dot_n <= f64::EPSILON {
            return Err("segment is too short");
        }
        // x ⋅ n  with x = c - a
        let x_dot_n = (c[0] - a[0]) * (b[0] - a[0]) + (c[1] - a[1]) * (b[1] - a[1]) + (c[2] - a[2]) * (b[2] - a[2]);
        // q = x - p = (c-a) - (b-a) * (x⋅n)/(n⋅n)
        let q0 = (c[0] - a[0]) - (b[0] - a[0]) * x_dot_n / n_dot_n;
        let q1 = (c[1] - a[1]) - (b[1] - a[1]) * x_dot_n / n_dot_n;
        let q2 = (c[2] - a[2]) - (b[2] - a[2]) * x_dot_n / n_dot_n;
        // norm(q)
        distance = f64::sqrt(q0 * q0 + q1 * q1 + q2 * q2);
    }
    Ok(distance)
}

/// Computes the signed distance from point to circumference (negative is inside)
///
/// # Note
///
/// This works in 2D only
pub fn point_circumference_distance(center: &[f64], radius: f64, x: &[f64]) -> Result<f64, StrError> {
    let ndim = center.len();
    if ndim != 2 {
        return Err("center.len() == ndim must be 2");
    }
    if x.len() != ndim {
        return Err("x.len() must equal center.len() == ndim");
    }
    let center_distance = point_point_distance(center, x)?;
    let distance = center_distance - radius;
    Ok(distance)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::{SQRT_2, SQRT_2_BY_3, SQRT_3};
    use russell_chk::assert_approx_eq;

    #[test]
    fn point_point_distance_fails_on_wrong_input() {
        assert_eq!(
            point_point_distance(&[0.0], &[1.0, 1.0]).err(),
            Some("a.len() == ndim must be 2 or 3")
        );
        assert_eq!(
            point_point_distance(&[0.0, 0.0], &[1.0]).err(),
            Some("b.len() must equal a.len() == ndim")
        );
    }

    #[test]
    fn point_point_distance_2d_works() -> Result<(), StrError> {
        let a = &[0.0, 0.0];
        let b = &[1.0, 0.0];
        assert_approx_eq!(point_point_distance(a, b)?, 1.0, 1e-15);

        let b = &[0.0, 2.0];
        assert_approx_eq!(point_point_distance(a, b)?, 2.0, 1e-15);

        let a = &[1.0, 1.0];
        let b = &[2.0, 2.0];
        assert_approx_eq!(point_point_distance(a, b)?, SQRT_2, 1e-15);
        Ok(())
    }

    #[test]
    fn point_point_distance_3d_works() -> Result<(), StrError> {
        let a = &[0.0, 0.0, 0.0];
        let b = &[1.0, 0.0, 0.0];
        assert_approx_eq!(point_point_distance(a, b)?, 1.0, 1e-15);

        let b = &[0.0, 2.0, 0.0];
        assert_approx_eq!(point_point_distance(a, b)?, 2.0, 1e-15);

        let b = &[0.0, 0.0, 3.0];
        assert_approx_eq!(point_point_distance(a, b)?, 3.0, 1e-15);

        let a = &[1.0, 1.0, 1.0];
        let b = &[2.0, 2.0, 2.0];
        assert_approx_eq!(point_point_distance(a, b)?, SQRT_3, 1e-15);
        Ok(())
    }

    #[test]
    fn point_segment_distance_fails_on_wrong_input() {
        assert_eq!(
            point_segment_distance(&[0.0], &[1.0, 1.0], &[2.0, 2.0]).err(),
            Some("a.len() == ndim must be 2 or 3")
        );
        assert_eq!(
            point_segment_distance(&[0.0, 0.0], &[1.0], &[2.0, 2.0]).err(),
            Some("b.len() must equal a.len() == ndim")
        );
        assert_eq!(
            point_segment_distance(&[0.0, 0.0], &[1.0, 1.0], &[2.0]).err(),
            Some("c.len() must equal a.len() == ndim")
        );
    }

    #[test]
    fn point_segment_distance_2d_fails_on_short_segment() -> Result<(), StrError> {
        let a = &[0.0, 0.0];
        let b = &[0.0, f64::EPSILON];
        let c = &[1.0, 1.0];
        assert_eq!(point_segment_distance(a, b, c).err(), Some("segment is too short"));
        Ok(())
    }

    #[test]
    fn point_segment_distance_3d_fails_on_short_segment() -> Result<(), StrError> {
        let a = &[0.0, 0.0, 0.0];
        let b = &[0.0, 0.0, f64::EPSILON];
        let c = &[1.0, 1.0, 1.0];
        assert_eq!(point_segment_distance(a, b, c).err(), Some("segment is too short"));
        Ok(())
    }

    #[test]
    fn point_segment_distance_2d_works() -> Result<(), StrError> {
        let a = &[0.0, 0.0];
        let b = &[2.0, 2.0];
        let c = &a.clone();
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 0.0);

        let c = &b.clone();
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 0.0);

        let c = &[1.0, 0.0];
        let distance = point_segment_distance(a, b, c)?;
        assert_approx_eq!(distance, SQRT_2 / 2.0, 1e-15);

        let c = &[0.0, 1.0];
        let distance = point_segment_distance(a, b, c)?;
        assert_approx_eq!(distance, SQRT_2 / 2.0, 1e-15);

        let c = &[0.5, 0.5];
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 0.0);

        let a = &[0.0, 0.0];
        let b = &[8.0, 0.0];
        let c = &[1.0, 0.5];
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 0.5);

        let a = &[0.0, 0.0];
        let b = &[0.0, 3.0];
        let c = &[0.5, 1.0];
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 0.5);

        let a = &[0.6, 0.0];
        let b = &[0.6, 1.8];
        let c = &[0.3, 1.2];
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 0.3);
        Ok(())
    }

    #[test]
    fn point_segment_distance_3d_works() -> Result<(), StrError> {
        let a = &[0.0, 0.0, 0.0];
        let b = &[2.0, 2.0, 2.0];
        let c = &a.clone();
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 0.0);

        let c = &b.clone();
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 0.0);

        let c = &[1.0, 1.0, 0.0];
        let distance = point_segment_distance(a, b, c)?;
        assert_approx_eq!(distance, SQRT_2_BY_3, 1e-15);

        let c = &[0.0, 1.0, 1.0];
        let distance = point_segment_distance(a, b, c)?;
        assert_approx_eq!(distance, SQRT_2_BY_3, 1e-15);

        let c = &[1.0, 0.0, 1.0];
        let distance = point_segment_distance(a, b, c)?;
        assert_approx_eq!(distance, SQRT_2_BY_3, 1e-15);

        let c = &[0.5, 0.5, 0.5];
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 0.0);

        let a = &[0.0, 0.0, 0.0];
        let b = &[2.0, 0.0, 0.0];
        let c = &[1.0, 0.5, 0.0];
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 0.5);

        let a = &[0.0, 0.0, 0.0];
        let b = &[0.0, 3.0, 0.0];
        let c = &[0.5, 1.0, 0.0];
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 0.5);

        let a = &[0.0, 0.0, 0.0];
        let b = &[0.0, 0.0, 4.0];
        let c = &[0.0, 0.5, 1.0];
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 0.5);

        let a = &[1.0, 2.0, 1.0];
        let b = &[1.0, 2.0, 5.0];
        let c = &[1.0, 1.0, 2.0];
        let distance = point_segment_distance(a, b, c)?;
        assert_eq!(distance, 1.0);
        Ok(())
    }

    #[test]
    fn point_circumference_distance_fails_on_wrong_input() {
        assert_eq!(
            point_circumference_distance(&[0.0], 1.0, &[2.0, 2.0]).err(),
            Some("center.len() == ndim must be 2")
        );
        assert_eq!(
            point_circumference_distance(&[0.0, 0.0], 1.0, &[2.0]).err(),
            Some("x.len() must equal center.len() == ndim")
        );
    }

    #[test]
    fn point_circumference_distance_works() -> Result<(), StrError> {
        let center = &[3.0, 4.0];
        let radius = 5.0;
        let distance = point_circumference_distance(center, radius, &[0.0, 0.0])?;
        assert_eq!(distance, 0.0);

        let distance = point_circumference_distance(center, radius, &[3.0, 4.0])?;
        assert_eq!(distance, -5.0);

        let distance = point_circumference_distance(center, radius, &[6.0, 8.0])?;
        assert_eq!(distance, 0.0);

        let distance = point_circumference_distance(center, radius, &[9.0, 12.0])?;
        assert_eq!(distance, 5.0);
        Ok(())
    }
}
