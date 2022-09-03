use crate::StrError;

/// Computes the point-to-point distance
///
/// # Input
///
/// * `a` -- first point
/// * `b` -- second point
///
/// # Output
///
/// Returns the unsigned distance between the two points.
///
/// # Examples
///
/// ```
/// use gemlab::geometry::point_point_distance;
/// use gemlab::StrError;
/// use gemlab::util::SQRT_2;
/// use russell_chk::approx_eq;
///
/// fn main() -> Result<(), StrError> {
///     let d = point_point_distance(&[0.0, 0.0], &[1.0, 1.0])?;
///     approx_eq(d, SQRT_2, 1e-15);
///     Ok(())
/// }
/// ```
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

/// Computes the distance between a point and a line
///
/// ```text
///                 c                x := c - a
///                /^                n := b - a
/// ⇀             / | ⇀   ⇀   ⇀      p := ((x ⋅ n) / (n ⋅ n)) * n
/// x = vec(c-a) /  | q = x - p      q := x - p
///             /   |                distance := norm(q)
///            /    |
///           a----->-------------b
///             ⇀              ⇀
///             p = projection(x)
///
///               ⇀
///               n = vec(b-c)
/// ```
///
/// # Input
///
/// * `a` -- first point on the line
/// * `b` -- second point on the line (different than `a`)
/// * `c` -- point whose distance is sought
///
/// # Output
///
/// Returns the unsigned distance.
///
/// # Examples
///
/// ```
/// use gemlab::geometry::point_line_distance;
/// use gemlab::StrError;
/// use gemlab::util::{SQRT_2, SQRT_3};
/// use russell_chk::approx_eq;
///
/// fn main() -> Result<(), StrError> {
///     let a = &[-10.0, -10.0, -10.0];
///     let b = &[100.0, 100.0, 100.0];
///     let d = point_line_distance(a, b, &[1.0, 1.0, 1.0])?;
///     approx_eq(d, 0.0, 1e-15);
///
///     let d = point_line_distance(a, b, &[1.0, 1.0, 0.0])?;
///     approx_eq(d, SQRT_2 / SQRT_3, 1e-15);
///
///     let d = point_line_distance(a, b, &[0.0, 1.0, 1.0])?;
///     approx_eq(d, SQRT_2 / SQRT_3, 1e-15);
///     Ok(())
/// }
/// ```
pub fn point_line_distance(a: &[f64], b: &[f64], c: &[f64]) -> Result<f64, StrError> {
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
            return Err("a-to-b segment is too short");
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
            return Err("a-to-b segment is too short");
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

/// Computes the signed distance from point to circle perimeter
///
/// # Input
///
/// `center` -- 2D circle center
/// `radius` -- circle radius
/// `p` -- 2D point
///
/// # Output
///
/// Returns the signed distance (negative is inside).
///
/// # Note
///
/// This works in 2D only.
///
/// # Examples
///
/// ```
/// use gemlab::geometry::point_circle_distance;
/// use gemlab::StrError;
/// use gemlab::util::SQRT_2;
/// use russell_chk::approx_eq;
///
/// fn main() -> Result<(), StrError> {
///     let center = &[3.0, 4.0];
///     let radius = 5.0;
///     let d = point_circle_distance(center, radius, &[3.0, 4.0])?;
///     assert_eq!(d, -radius);
///
///     let d = point_circle_distance(center, radius, &[8.0, 4.0])?;
///     assert_eq!(d, 0.0);
///
///     let d = point_circle_distance(center, radius, &[3.0, 9.0])?;
///     assert_eq!(d, 0.0);
///
///     let d = point_circle_distance(center, radius, &[8.0, 9.0])?;
///     approx_eq(d, radius * SQRT_2 - radius, 1e-15);
///     Ok(())
/// }
/// ```
pub fn point_circle_distance(center: &[f64], radius: f64, p: &[f64]) -> Result<f64, StrError> {
    let ndim = center.len();
    if ndim != 2 {
        return Err("center.len() == ndim must be 2");
    }
    if p.len() != ndim {
        return Err("x.len() must equal center.len() == ndim");
    }
    let center_distance = point_point_distance(center, p).unwrap(); // we can use unwrap here because center.len() == p.len() == 2
    let distance = center_distance - radius;
    Ok(distance)
}

/// Computes the signed distance from point to the surface of a cylinder parallel to x
///
/// # Input
///
/// `a` -- First 3D point on the cylinder axis
/// `b` -- Second 3D point on the cylinder axis (must be different than `a`)
/// `radius` -- cylinder radius
/// `p` -- 3D point
///
/// # Output
///
/// Returns the signed distance (negative is inside)
///
/// # Note
///
/// This works in 3D only.
///
/// # Examples
///
/// ```
/// use gemlab::geometry::point_cylinder_distance;
/// use gemlab::StrError;
/// use gemlab::util::{SQRT_2, SQRT_3};
/// use russell_chk::approx_eq;
///
/// fn main() -> Result<(), StrError> {
///     let a = &[-10.0, -10.0, -10.0];
///     let b = &[100.0, 100.0, 100.0];
///     let radius = 5.0;
///     let d = point_cylinder_distance(a, b, radius, &[0.0, 0.0, 0.0])?;
///     assert_eq!(d, -radius);
///
///     let d = point_cylinder_distance(a, b, radius, &[1.0, 1.0, 0.0])?;
///     approx_eq(d, SQRT_2 / SQRT_3 - radius, 1e-15);
///
///     let z = radius * SQRT_3 / SQRT_2;
///     let d = point_cylinder_distance(a, b, radius, &[0.0, 0.0, z])?;
///     approx_eq(d, 0.0, 1e-14);
///     Ok(())
/// }
/// ```
pub fn point_cylinder_distance(a: &[f64], b: &[f64], radius: f64, p: &[f64]) -> Result<f64, StrError> {
    let ndim = a.len();
    if ndim != 3 {
        return Err("axis_a.len() == ndim must be 3");
    }
    if b.len() != ndim {
        return Err("axis_b.len() must equal center.len() == ndim");
    }
    if p.len() != ndim {
        return Err("p.len() must equal center.len() == ndim");
    }
    let center_distance = point_line_distance(a, b, p)?;
    let distance = center_distance - radius;
    Ok(distance)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::{SQRT_2, SQRT_2_BY_3, SQRT_3};
    use russell_chk::approx_eq;

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
    fn point_point_distance_2d_works() {
        let a = &[0.0, 0.0];
        let b = &[1.0, 0.0];
        approx_eq(point_point_distance(a, b).unwrap(), 1.0, 1e-15);

        let b = &[0.0, 2.0];
        approx_eq(point_point_distance(a, b).unwrap(), 2.0, 1e-15);

        let a = &[1.0, 1.0];
        let b = &[2.0, 2.0];
        approx_eq(point_point_distance(a, b).unwrap(), SQRT_2, 1e-15);
    }

    #[test]
    fn point_point_distance_3d_works() {
        let a = &[0.0, 0.0, 0.0];
        let b = &[1.0, 0.0, 0.0];
        approx_eq(point_point_distance(a, b).unwrap(), 1.0, 1e-15);

        let b = &[0.0, 2.0, 0.0];
        approx_eq(point_point_distance(a, b).unwrap(), 2.0, 1e-15);

        let b = &[0.0, 0.0, 3.0];
        approx_eq(point_point_distance(a, b).unwrap(), 3.0, 1e-15);

        let a = &[1.0, 1.0, 1.0];
        let b = &[2.0, 2.0, 2.0];
        approx_eq(point_point_distance(a, b).unwrap(), SQRT_3, 1e-15);
    }

    #[test]
    fn point_line_distance_fails_on_wrong_input() {
        assert_eq!(
            point_line_distance(&[0.0], &[1.0, 1.0], &[2.0, 2.0]).err(),
            Some("a.len() == ndim must be 2 or 3")
        );
        assert_eq!(
            point_line_distance(&[0.0, 0.0], &[1.0], &[2.0, 2.0]).err(),
            Some("b.len() must equal a.len() == ndim")
        );
        assert_eq!(
            point_line_distance(&[0.0, 0.0], &[1.0, 1.0], &[2.0]).err(),
            Some("c.len() must equal a.len() == ndim")
        );
    }

    #[test]
    fn point_line_distance_2d_fails_on_short_segment() {
        let a = &[0.0, 0.0];
        let b = &[0.0, f64::EPSILON];
        let c = &[1.0, 1.0];
        assert_eq!(point_line_distance(a, b, c).err(), Some("a-to-b segment is too short"));
    }

    #[test]
    fn point_line_distance_3d_fails_on_short_segment() {
        let a = &[0.0, 0.0, 0.0];
        let b = &[0.0, 0.0, f64::EPSILON];
        let c = &[1.0, 1.0, 1.0];
        assert_eq!(point_line_distance(a, b, c).err(), Some("a-to-b segment is too short"));
    }

    #[test]
    fn point_line_distance_2d_works() {
        let a = &[0.0, 0.0];
        let b = &[2.0, 2.0];
        let c = &a.clone();
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 0.0);

        let c = &b.clone();
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 0.0);

        let c = &[1.0, 0.0];
        let distance = point_line_distance(a, b, c).unwrap();
        approx_eq(distance, SQRT_2 / 2.0, 1e-15);

        let c = &[0.0, 1.0];
        let distance = point_line_distance(a, b, c).unwrap();
        approx_eq(distance, SQRT_2 / 2.0, 1e-15);

        let c = &[0.5, 0.5];
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 0.0);

        let a = &[0.0, 0.0];
        let b = &[8.0, 0.0];
        let c = &[1.0, 0.5];
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 0.5);

        let a = &[0.0, 0.0];
        let b = &[0.0, 3.0];
        let c = &[0.5, 1.0];
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 0.5);

        let a = &[0.6, 0.0];
        let b = &[0.6, 1.8];
        let c = &[0.3, 1.2];
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 0.3);
    }

    #[test]
    fn point_line_distance_3d_works() {
        let a = &[0.0, 0.0, 0.0];
        let b = &[2.0, 2.0, 2.0];
        let c = &a.clone();
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 0.0);

        let c = &b.clone();
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 0.0);

        let c = &[1.0, 1.0, 0.0];
        let distance = point_line_distance(a, b, c).unwrap();
        approx_eq(distance, SQRT_2_BY_3, 1e-15);

        let c = &[0.0, 1.0, 1.0];
        let distance = point_line_distance(a, b, c).unwrap();
        approx_eq(distance, SQRT_2_BY_3, 1e-15);

        let c = &[1.0, 0.0, 1.0];
        let distance = point_line_distance(a, b, c).unwrap();
        approx_eq(distance, SQRT_2_BY_3, 1e-15);

        let c = &[0.5, 0.5, 0.5];
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 0.0);

        let a = &[0.0, 0.0, 0.0];
        let b = &[2.0, 0.0, 0.0];
        let c = &[1.0, 0.5, 0.0];
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 0.5);

        let a = &[0.0, 0.0, 0.0];
        let b = &[0.0, 3.0, 0.0];
        let c = &[0.5, 1.0, 0.0];
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 0.5);

        let a = &[0.0, 0.0, 0.0];
        let b = &[0.0, 0.0, 4.0];
        let c = &[0.0, 0.5, 1.0];
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 0.5);

        let a = &[1.0, 2.0, 1.0];
        let b = &[1.0, 2.0, 5.0];
        let c = &[1.0, 1.0, 2.0];
        let distance = point_line_distance(a, b, c).unwrap();
        assert_eq!(distance, 1.0);
    }

    #[test]
    fn point_circle_distance_fails_on_wrong_input() {
        assert_eq!(
            point_circle_distance(&[0.0], 1.0, &[2.0, 2.0]).err(),
            Some("center.len() == ndim must be 2")
        );
        assert_eq!(
            point_circle_distance(&[0.0, 0.0], 1.0, &[2.0]).err(),
            Some("x.len() must equal center.len() == ndim")
        );
    }

    #[test]
    fn point_circle_distance_works() {
        let center = &[3.0, 4.0];
        let radius = 5.0;
        let distance = point_circle_distance(center, radius, &[0.0, 0.0]).unwrap();
        assert_eq!(distance, 0.0);

        let distance = point_circle_distance(center, radius, &[3.0, 4.0]).unwrap();
        assert_eq!(distance, -5.0);

        let distance = point_circle_distance(center, radius, &[6.0, 8.0]).unwrap();
        assert_eq!(distance, 0.0);

        let distance = point_circle_distance(center, radius, &[9.0, 12.0]).unwrap();
        assert_eq!(distance, 5.0);
    }

    #[test]
    fn point_cylinder_distance_fails_on_wrong_input() {
        assert_eq!(
            point_cylinder_distance(&[0.0, 0.0], &[1.0, 0.0, 0.0], 1.0, &[2.0, 2.0, 2.0]).err(),
            Some("axis_a.len() == ndim must be 3")
        );
        assert_eq!(
            point_cylinder_distance(&[0.0, 0.0, 0.0], &[1.0, 0.0], 1.0, &[2.0, 2.0, 2.0]).err(),
            Some("axis_b.len() must equal center.len() == ndim")
        );
        assert_eq!(
            point_cylinder_distance(&[0.0, 0.0, 0.0], &[1.0, 0.0, 0.0], 1.0, &[2.0, 2.0]).err(),
            Some("p.len() must equal center.len() == ndim")
        );
    }

    #[test]
    fn point_cylinder_distance_works() {
        // parallel to x
        let axis_a = &[10.0, 3.0, 4.0];
        let axis_b = &[20.0, 3.0, 4.0];
        let radius = 5.0;
        let distance = point_cylinder_distance(axis_a, axis_b, radius, &[0.0, 0.0, 0.0]).unwrap();
        assert_eq!(distance, 0.0);

        let distance = point_cylinder_distance(axis_a, axis_b, radius, &[80.0, 3.0, 4.0]).unwrap();
        assert_eq!(distance, -5.0);

        let distance = point_cylinder_distance(axis_a, axis_b, radius, &[-11.0, 6.0, 8.0]).unwrap();
        assert_eq!(distance, 0.0);

        let distance = point_cylinder_distance(axis_a, axis_b, radius, &[-5.0, 9.0, 12.0]).unwrap();
        assert_eq!(distance, 5.0);

        // parallel to y
        let axis_a = &[3.0, 10.0, 4.0];
        let axis_b = &[3.0, -10.0, 4.0];
        let radius = 5.0;
        let distance = point_cylinder_distance(axis_a, axis_b, radius, &[0.0, 0.0, 0.0]).unwrap();
        assert_eq!(distance, 0.0);

        let distance = point_cylinder_distance(axis_a, axis_b, radius, &[3.0, 80.0, 4.0]).unwrap();
        assert_eq!(distance, -5.0);

        let distance = point_cylinder_distance(axis_a, axis_b, radius, &[6.0, -11.0, 8.0]).unwrap();
        assert_eq!(distance, 0.0);

        let distance = point_cylinder_distance(axis_a, axis_b, radius, &[9.0, -5.0, 12.0]).unwrap();
        assert_eq!(distance, 5.0);

        // parallel to z
        let axis_a = &[3.0, 4.0, 10.0];
        let axis_b = &[3.0, 4.0, 30.0];
        let radius = 5.0;
        let distance = point_cylinder_distance(axis_a, axis_b, radius, &[0.0, 0.0, 0.0]).unwrap();
        assert_eq!(distance, 0.0);

        let distance = point_cylinder_distance(axis_a, axis_b, radius, &[3.0, 4.0, 80.0]).unwrap();
        assert_eq!(distance, -5.0);

        let distance = point_cylinder_distance(axis_a, axis_b, radius, &[6.0, 8.0, -11.0]).unwrap();
        assert_eq!(distance, 0.0);

        let distance = point_cylinder_distance(axis_a, axis_b, radius, &[9.0, 12.0, -5.0]).unwrap();
        assert_eq!(distance, 5.0);
    }
}
