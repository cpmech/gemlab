/// Holds the Cartesian coordinates of a point in 3D
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Represents a directed segment from A to B in 3D
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Segment3D {
    pub a: Point3D,
    pub b: Point3D,
}

/// Represents a vector in 3D
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector3D {
    pub u: f64,
    pub v: f64,
    pub w: f64,
}

impl Point3D {
    /// Returns a new translated point
    pub fn get_translated(&self, dx: f64, dy: f64, dz: f64) -> Point3D {
        Point3D {
            x: self.x + dx,
            y: self.y + dy,
            z: self.z + dz,
        }
    }

    /// Computes the unsigned distance between this point and another
    pub fn dist_from_point(&self, another: &Point3D) -> f64 {
        f64::sqrt(
            (self.x - another.x) * (self.x - another.x)
                + (self.y - another.y) * (self.y - another.y)
                + (self.z - another.z) * (self.z - another.z),
        )
    }
}

impl Segment3D {
    /// Creates a new segment from A to B
    pub fn new(a: &Point3D, b: &Point3D) -> Self {
        Segment3D {
            a: a.clone(),
            b: b.clone(),
        }
    }

    /// Creates a new segment scaled by m and starting from A
    pub fn get_scaled(&self, m: f64) -> Segment3D {
        Segment3D {
            a: self.a.clone(),
            b: Point3D {
                x: self.a.x + m * (self.b.x - self.a.x),
                y: self.a.y + m * (self.b.y - self.a.y),
                z: self.a.z + m * (self.b.z - self.a.z),
            },
        }
    }

    /// Computes the length of this segment (Euclidean norm)
    pub fn length(&self) -> f64 {
        self.a.dist_from_point(&self.b)
    }

    /// Returns the vector representing this segment, scaled by m
    pub fn get_vector(&self, m: f64) -> Vector3D {
        Vector3D {
            u: m * (self.b.x - self.a.x),
            v: m * (self.b.y - self.a.y),
            w: m * (self.b.z - self.a.z),
        }
    }
}

impl Vector3D {
    /// Creates a new vector scaled by m
    pub fn new_scaled(u: f64, v: f64, w: f64, m: f64) -> Self {
        Vector3D {
            u: m * u,
            v: m * v,
            w: m * w,
        }
    }

    /// Creates a new vector by adding this vector with another
    ///
    /// ```text
    /// result := α * this + β * another
    /// ```
    pub fn get_added(&self, alpha: f64, beta: f64, another: &Vector3D) -> Vector3D {
        Vector3D {
            u: alpha * self.u + beta * another.u,
            v: alpha * self.v + beta * another.v,
            w: alpha * self.w + beta * another.w,
        }
    }

    /// Computes the inner product of this vector with another
    pub fn dot(&self, another: &Vector3D) -> f64 {
        self.u * another.u + self.v * another.v + self.w * another.w
    }

    /// Computes the Euclidean norm of this vector
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
        let a = Point3D {
            x: -1.0,
            y: 0.0,
            z: -2.0,
        };
        let b = a.get_translated(1.0, 2.0, 3.0);
        assert_eq!(b.x, 0.0);
        assert_eq!(b.y, 2.0);
        assert_eq!(b.z, 1.0);
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
    fn segment3d_new_works() {
        let a = Point3D { x: 0.0, y: 2.0, z: 3.0 };
        let b = Point3D { x: 3.0, y: 2.0, z: 1.0 };
        let segment = Segment3D::new(&a, &b);
        assert_eq!(segment.a.x, 0.0);
        assert_eq!(segment.a.y, 2.0);
        assert_eq!(segment.a.z, 3.0);
        assert_eq!(segment.b.x, 3.0);
        assert_eq!(segment.b.y, 2.0);
        assert_eq!(segment.b.z, 1.0);
    }

    #[test]
    fn segment3d_get_scaled_works() {
        let a = Point3D { x: 0.0, y: 2.0, z: 3.0 };
        let b = Point3D { x: 3.0, y: 2.0, z: 1.0 };
        let segment = Segment3D::new(&a, &b);
        let scaled = segment.get_scaled(0.5);
        assert_eq!(scaled.a.x, 0.0);
        assert_eq!(scaled.a.y, 2.0);
        assert_eq!(scaled.a.z, 3.0);
        assert_eq!(scaled.b.x, 1.5);
        assert_eq!(scaled.b.y, 2.0);
        assert_eq!(scaled.b.z, 2.0);
    }

    #[test]
    fn segment3d_length_works() {
        let a = Point3D { x: 0.0, y: 2.0, z: 3.0 };
        let b = Point3D { x: 3.0, y: 2.0, z: 1.0 };
        let segment = Segment3D::new(&a, &b);
        assert_approx_eq!(segment.length(), f64::sqrt(13.0), 1e-15)
    }

    #[test]
    fn segment3d_get_vector_works() {
        let a = Point3D { x: 0.0, y: 2.0, z: 3.0 };
        let b = Point3D { x: 3.0, y: 2.0, z: 1.0 };
        let segment = Segment3D::new(&a, &b);
        let vector = segment.get_vector(0.5);
        assert_eq!(vector.u, 1.5);
        assert_eq!(vector.v, 0.0);
        assert_eq!(vector.w, -1.0);
    }

    #[test]
    fn vector3d_new_scaled_works() {
        let vector = Vector3D::new_scaled(1.0, 2.0, 3.0, 0.5);
        assert_eq!(vector.u, 0.5);
        assert_eq!(vector.v, 1.0);
        assert_eq!(vector.w, 1.5);
    }

    #[test]
    fn vector3d_get_added_works() {
        let x = Vector3D { u: 1.0, v: 2.0, w: 3.0 };
        let y = Vector3D { u: 3.0, v: 2.0, w: 1.0 };
        let z = x.get_added(0.5, 3.0, &y);
        assert_eq!(z.u, 9.5);
        assert_eq!(z.v, 7.0);
        assert_eq!(z.w, 4.5);
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
