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
    pub fn new_added(&self, alpha: f64, beta: f64, another: &Vector3D) -> Vector3D {
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

    #[test]
    fn get_translated_works() {
        let a = Point3D {
            x: -1.0,
            y: 0.0,
            z: -2.0,
        };
        let b = a.get_translated(1.0, 2.0, 3.0);
        assert_eq!(b.x, 0.0);
    }
}
