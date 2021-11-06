/// Holds the Cartesian coordinates of a point in 3D
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Represents a vector in 3D
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector3D {
    pub u: f64,
    pub v: f64,
    pub w: f64,
}

/// Represents a directed segment from A to B in 3D
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Segment3D {
    pub a: Point3D,
    pub b: Point3D,
}

impl Point3D {
    /// Returns a new translated point
    #[inline]
    pub fn get_translated(&self, dx: f64, dy: f64, dz: f64) -> Point3D {
        Point3D {
            x: self.x + dx,
            y: self.y + dy,
            z: self.z + dz,
        }
    }

    /// Computes the unsigned distance between this point and another
    #[inline]
    pub fn dist_from_point(&self, another: &Point3D) -> f64 {
        f64::sqrt(
            (self.x - another.x) * (self.x - another.x)
                + (self.y - another.y) * (self.y - another.y)
                + (self.z - another.z) * (self.z - another.z),
        )
    }
}

impl Vector3D {
    /// Creates a new vector from two points (a to b)
    ///
    /// ```text
    /// vector := b - a
    /// ```
    #[inline]
    pub fn new_from_points(a: &Point3D, b: &Point3D) -> Self {
        Vector3D {
            u: b.x - a.x,
            v: b.y - a.y,
            w: b.z - a.z,
        }
    }

    /// Creates a new scaled vector
    #[inline]
    pub fn new_scaled(u: f64, v: f64, w: f64, scale: f64) -> Self {
        Vector3D {
            u: scale * u,
            v: scale * v,
            w: scale * w,
        }
    }

    /// Creates a new vector by subtracting this vector by another
    ///
    /// ```text
    /// result := this - another
    /// ```
    #[inline]
    pub fn get_subtracted(&self, another: &Vector3D) -> Vector3D {
        Vector3D {
            u: self.u - another.u,
            v: self.v - another.v,
            w: self.w - another.w,
        }
    }

    /// Computes the inner product of this vector with another
    #[inline]
    pub fn dot(&self, another: &Vector3D) -> f64 {
        self.u * another.u + self.v * another.v + self.w * another.w
    }

    /// Computes the Euclidean norm of this vector
    #[inline]
    pub fn norm(&self) -> f64 {
        f64::sqrt(self.u * self.u + self.v * self.v + self.w * self.w)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SQRT_3;
    use russell_chk::assert_approx_eq;

    #[test]
    fn point3d_get_translated_works() {
        let a = Point3D { x: 1.0, y: 2.0, z: 3.0 };
        let b = a.get_translated(1.0, 2.0, 3.0);
        assert_eq!(b.x, 2.0);
        assert_eq!(b.y, 4.0);
        assert_eq!(b.z, 6.0);
    }

    #[test]
    fn point3d_dist_from_point_works() {
        let a = Point3D { x: 0.0, y: 0.0, z: 0.0 };
        let b = Point3D { x: 1.0, y: 0.0, z: 0.0 };
        assert_approx_eq!(a.dist_from_point(&b), 1.0, 1e-15);

        let b = Point3D { x: 0.0, y: 2.0, z: 0.0 };
        assert_approx_eq!(a.dist_from_point(&b), 2.0, 1e-15);

        let b = Point3D { x: 0.0, y: 0.0, z: 3.0 };
        assert_approx_eq!(a.dist_from_point(&b), 3.0, 1e-15);

        let a = Point3D { x: 1.0, y: 1.0, z: 1.0 };
        let b = Point3D { x: 2.0, y: 2.0, z: 2.0 };
        assert_approx_eq!(a.dist_from_point(&b), SQRT_3, 1e-15);
    }

    #[test]
    fn vector3d_new_from_points_works() {
        let a = Point3D { x: 1.0, y: 2.0, z: 3.0 };
        let b = Point3D { x: 4.0, y: 5.0, z: 6.0 };
        let vector = Vector3D::new_from_points(&a, &b);
        assert_eq!(vector.u, 3.0);
        assert_eq!(vector.v, 3.0);
        assert_eq!(vector.w, 3.0);
    }

    #[test]
    fn vector3d_new_scaled_works() {
        let vector = Vector3D::new_scaled(1.0, 2.0, 3.0, 0.5);
        assert_eq!(vector.u, 0.5);
        assert_eq!(vector.v, 1.0);
        assert_eq!(vector.w, 1.5);
    }

    #[test]
    fn vector3d_get_subtracted_works() {
        let x = Vector3D { u: 1.0, v: 2.0, w: 3.0 };
        let y = Vector3D { u: 3.0, v: 2.0, w: 1.0 };
        let z = x.get_subtracted(&y);
        assert_eq!(z.u, -2.0);
        assert_eq!(z.v, 0.0);
        assert_eq!(z.w, 2.0);
    }

    #[test]
    fn vector3d_dot_works() {
        let x = Vector3D { u: 1.0, v: 2.0, w: 3.0 };
        let y = Vector3D { u: 3.0, v: 2.0, w: 1.0 };
        assert_eq!(x.dot(&y), 10.0);
    }

    #[test]
    fn vector3d_norm() {
        let x = Vector3D { u: 3.0, v: 4.0, w: 5.0 };
        assert_approx_eq!(x.norm(), f64::sqrt(50.0), 1e-15);
    }
}
