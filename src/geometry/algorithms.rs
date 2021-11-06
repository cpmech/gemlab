use super::{Point2D, Point3D, Vector2D, Vector3D};
use crate::StrError;

/// Computes the distance between a point and a segment in 2D
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
pub fn point_segment_distance_2d(a: &Point2D, b: &Point2D, c: &Point2D) -> Result<f64, StrError> {
    let n = Vector2D::new_from_points(a, b); // b - a
    let denominator = n.dot(&n); // n ⋅ n
    if denominator <= f64::EPSILON {
        return Err("segment is too short");
    }
    let x = Vector2D::new_from_points(a, c); // c - a
    let scale = x.dot(&n) / denominator; // (x ⋅ n) / (n ⋅ n)
    let p = Vector2D::new_scaled(n.u, n.v, scale);
    let q = x.get_subtracted(&p); // x - p
    Ok(q.norm())
}

/// Computes the distance between a point and a segment in 3D
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
pub fn point_segment_distance_3d(a: &Point3D, b: &Point3D, c: &Point3D) -> Result<f64, StrError> {
    let n = Vector3D::new_from_points(a, b); // b - a
    let denominator = n.dot(&n); // n ⋅ n
    if denominator <= f64::EPSILON {
        return Err("segment is too short");
    }
    let x = Vector3D::new_from_points(a, c); // c - a
    let scale = x.dot(&n) / denominator; // (x ⋅ n) / (n ⋅ n)
    let p = Vector3D::new_scaled(n.u, n.v, n.w, scale);
    let q = x.get_subtracted(&p); // x - p
    Ok(q.norm())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::{SQRT_2, SQRT_2_BY_3};
    use russell_chk::assert_approx_eq;

    #[test]
    fn point_segment_distance_2d_fails_on_short_segment() -> Result<(), StrError> {
        let a = Point2D { x: 0.0, y: 0.0 };
        let b = Point2D {
            x: f64::EPSILON,
            y: f64::EPSILON,
        };
        let c = Point2D { x: 1.0, y: 1.0 };
        assert_eq!(
            point_segment_distance_2d(&a, &b, &c).err(),
            Some("segment is too short")
        );
        Ok(())
    }

    #[test]
    fn point_segment_distance_2d_works() -> Result<(), StrError> {
        let a = Point2D { x: 0.0, y: 0.0 };
        let b = Point2D { x: 1.0, y: 1.0 };
        let c = a.clone();
        let distance = point_segment_distance_2d(&a, &b, &c)?;
        assert_eq!(distance, 0.0);

        let c = b.clone();
        let distance = point_segment_distance_2d(&a, &b, &c)?;
        assert_eq!(distance, 0.0);

        let c = Point2D { x: 1.0, y: 0.0 };
        let distance = point_segment_distance_2d(&a, &b, &c)?;
        assert_approx_eq!(distance, SQRT_2 / 2.0, 1e-15);

        let c = Point2D { x: 0.0, y: 1.0 };
        let distance = point_segment_distance_2d(&a, &b, &c)?;
        assert_approx_eq!(distance, SQRT_2 / 2.0, 1e-15);

        let c = Point2D { x: 0.5, y: 0.5 };
        let distance = point_segment_distance_2d(&a, &b, &c)?;
        assert_eq!(distance, 0.0);

        let a = Point2D { x: 0.0, y: 0.0 };
        let b = Point2D { x: 1.0, y: 0.0 };
        let c = Point2D { x: 1.0, y: 0.5 };
        let distance = point_segment_distance_2d(&a, &b, &c)?;
        assert_eq!(distance, 0.5);

        let a = Point2D { x: 0.0, y: 0.0 };
        let b = Point2D { x: 0.0, y: 1.0 };
        let c = Point2D { x: 0.5, y: 1.0 };
        let distance = point_segment_distance_2d(&a, &b, &c)?;
        assert_eq!(distance, 0.5);
        Ok(())
    }

    #[test]
    fn point_segment_distance_3d_fails_on_short_segment() -> Result<(), StrError> {
        let a = Point3D { x: 0.0, y: 0.0, z: 0.0 };
        let b = Point3D {
            x: f64::EPSILON,
            y: f64::EPSILON,
            z: f64::EPSILON,
        };
        let c = Point3D { x: 1.0, y: 1.0, z: 1.0 };
        assert_eq!(
            point_segment_distance_3d(&a, &b, &c).err(),
            Some("segment is too short")
        );
        Ok(())
    }

    #[test]
    fn point_segment_distance_3d_works() -> Result<(), StrError> {
        let a = Point3D { x: 0.0, y: 0.0, z: 0.0 };
        let b = Point3D { x: 1.0, y: 1.0, z: 1.0 };
        let c = a.clone();
        let distance = point_segment_distance_3d(&a, &b, &c)?;
        assert_eq!(distance, 0.0);

        let c = b.clone();
        let distance = point_segment_distance_3d(&a, &b, &c)?;
        assert_eq!(distance, 0.0);

        let c = Point3D { x: 1.0, y: 1.0, z: 0.0 };
        let distance = point_segment_distance_3d(&a, &b, &c)?;
        assert_approx_eq!(distance, SQRT_2_BY_3, 1e-15);

        let c = Point3D { x: 0.0, y: 1.0, z: 1.0 };
        let distance = point_segment_distance_3d(&a, &b, &c)?;
        assert_approx_eq!(distance, SQRT_2_BY_3, 1e-15);

        let c = Point3D { x: 1.0, y: 0.0, z: 1.0 };
        let distance = point_segment_distance_3d(&a, &b, &c)?;
        assert_approx_eq!(distance, SQRT_2_BY_3, 1e-15);

        let c = Point3D { x: 0.5, y: 0.5, z: 0.5 };
        let distance = point_segment_distance_3d(&a, &b, &c)?;
        assert_eq!(distance, 0.0);

        let a = Point3D { x: 0.0, y: 0.0, z: 0.0 };
        let b = Point3D { x: 1.0, y: 0.0, z: 0.0 };
        let c = Point3D { x: 1.0, y: 0.5, z: 0.0 };
        let distance = point_segment_distance_3d(&a, &b, &c)?;
        assert_eq!(distance, 0.5);

        let a = Point3D { x: 0.0, y: 0.0, z: 0.0 };
        let b = Point3D { x: 0.0, y: 1.0, z: 0.0 };
        let c = Point3D { x: 0.5, y: 1.0, z: 0.0 };
        let distance = point_segment_distance_3d(&a, &b, &c)?;
        assert_eq!(distance, 0.5);

        let a = Point3D { x: 0.0, y: 0.0, z: 0.0 };
        let b = Point3D { x: 0.0, y: 0.0, z: 1.0 };
        let c = Point3D { x: 0.0, y: 0.5, z: 1.0 };
        let distance = point_segment_distance_3d(&a, &b, &c)?;
        assert_eq!(distance, 0.5);
        Ok(())
    }
}
