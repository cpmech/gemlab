/// Defines a 2D point
#[derive(Debug, Clone, Copy)]
pub struct Point2d {
    pub x: f64,
    pub y: f64,
}

impl Point2d {
    /// Allocates a new instance
    pub fn new(x: f64, y: f64) -> Self {
        Point2d { x, y }
    }

    /// Allocates a new instance from a slice
    pub fn from_slice(slice: &[f64]) -> Self {
        assert!(slice.len() >= 2, "slice must have length at least 2");
        Point2d {
            x: slice[0],
            y: slice[1],
        }
    }
}

/// Defines a 2D vector
#[derive(Debug, Clone, Copy)]
pub struct Vector2d {
    pub ux: f64,
    pub uy: f64,
}

impl Vector2d {
    /// Allocates a new instance
    pub fn new(ux: f64, uy: f64) -> Self {
        Vector2d { ux, uy }
    }

    /// Allocates a new instance from a slice
    pub fn from_slice(slice: &[f64]) -> Self {
        assert!(slice.len() >= 2, "slice must have length at least 2");
        Vector2d {
            ux: slice[0],
            uy: slice[1],
        }
    }

    /// Allocates a new instance from two points
    pub fn from_points(a: &Point2d, b: &Point2d) -> Self {
        Vector2d {
            ux: b.x - a.x,
            uy: b.y - a.y,
        }
    }

    /// Calculates the cross product for 2D vectors
    ///
    /// Returns a scalar value corresponding to the magnitude of the out-of-plane vector
    pub fn cross(&self, other: &Vector2d) -> f64 {
        self.ux * other.uy - self.uy * other.ux
    }
}

/// Defines a parallelogram in 2D space
pub struct Parallelogram2d {
    pub u: Vector2d,
    pub v: Vector2d,
}

impl Parallelogram2d {
    /// Allocates a new instance from vectors
    pub fn from_vectors(u: &Vector2d, v: &Vector2d) -> Self {
        Parallelogram2d { u: *u, v: *v }
    }

    /// Allocates a new instance from vectors represented as slices
    pub fn from_slices(u: &[f64], v: &[f64]) -> Self {
        assert!(u.len() >= 2, "u slice must have length at least 2");
        assert!(v.len() >= 2, "v slice must have length at least 2");
        Parallelogram2d {
            u: Vector2d::new(u[0], u[1]),
            v: Vector2d::new(v[0], v[1]),
        }
    }

    /// Computes the signed area of the parallelogram defined by two 2D vectors
    pub fn signed_area(&self) -> f64 {
        self.u.cross(&self.v)
    }
}

/// Defines a triangle in 2D space
#[derive(Debug, Clone, Copy)]
pub struct Triangle2d {
    pub a: Point2d,
    pub b: Point2d,
    pub c: Point2d,
}

impl Triangle2d {
    /// Allocates a new instance from points
    pub fn from_points(a: &Point2d, b: &Point2d, c: &Point2d) -> Self {
        Triangle2d { a: *a, b: *b, c: *c }
    }

    /// Allocates a new instance from slices
    pub fn from_slices(a: &[f64], b: &[f64], c: &[f64]) -> Self {
        Triangle2d {
            a: Point2d::from_slice(a),
            b: Point2d::from_slice(b),
            c: Point2d::from_slice(c),
        }
    }

    /// Computes the signed area
    pub fn signed_area(&self) -> f64 {
        let vec_ab = Vector2d {
            ux: self.b.x - self.a.x,
            uy: self.b.y - self.a.y,
        };
        let vec_ac = Vector2d {
            ux: self.c.x - self.a.x,
            uy: self.c.y - self.a.y,
        };
        vec_ab.cross(&vec_ac) / 2.0
    }
}

/// Holds data defining a circle in 2D
///
/// # Examples
///
/// ```
/// use gemlab::geometry::Circle;
///
/// let circle = Circle {
///     center: [0.0, 0.0],
///     radius: 1.0,
/// };
///
/// assert_eq!(circle.center, [0.0, 0.0]);
/// assert_eq!(circle.radius, 1.0);
/// ```
#[derive(Clone, Debug)]
pub struct Circle {
    /// The center of the circle
    pub center: [f64; 2],

    /// The radius of the circle
    pub radius: f64,
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derive_methods_work() {
        let point = Point2d::new(1.0, 2.0);
        let clone = point.clone();
        let correct = "Point2d { x: 1.0, y: 2.0 }";
        assert_eq!(format!("{:?}", point), correct);
        assert_eq!(format!("{:?}", clone), correct);

        let vector = Vector2d::new(3.0, 4.0);
        let clone = vector.clone();
        let correct = "Vector2d { ux: 3.0, uy: 4.0 }";
        assert_eq!(format!("{:?}", vector), correct);
        assert_eq!(format!("{:?}", clone), correct);

        let triangle = Triangle2d {
            a: Point2d::new(0.0, 0.0),
            b: Point2d::new(1.0, 0.0),
            c: Point2d::new(0.0, 1.0),
        };
        let clone = triangle.clone();
        let correct = "Triangle2d { a: Point2d { x: 0.0, y: 0.0 }, b: Point2d { x: 1.0, y: 0.0 }, c: Point2d { x: 0.0, y: 1.0 } }";
        assert_eq!(format!("{:?}", triangle), correct);
        assert_eq!(format!("{:?}", clone), correct);

        let circle = Circle {
            center: [-1.0, -1.0],
            radius: 2.0,
        };
        let clone = circle.clone();
        let correct = "Circle { center: [-1.0, -1.0], radius: 2.0 }";
        assert_eq!(format!("{:?}", circle), correct);
        assert_eq!(format!("{:?}", clone), correct);
    }

    #[test]
    fn point2d_and_vector_2d_work() {
        let p1 = Point2d::new(1.0, 2.0);
        let p2 = Point2d::new(4.0, 3.0);
        let vec = Vector2d::from_points(&p1, &p2);
        assert_eq!(vec.ux, 3.0);
        assert_eq!(vec.uy, 1.0);
    }

    #[test]
    fn parallelogram_functions_work() {
        // Mathematica:
        // p = Parallelogram[{1, 2}, {{3, 1}, {1, 2}}];
        // Area[p]

        let u = Vector2d::new(3.0, 1.0);
        let v = Vector2d::new(1.0, 2.0);
        let area_vec = Parallelogram2d::from_vectors(&u, &v).signed_area();
        assert_eq!(area_vec, 5.0);

        // Sign reflects orientation (swap vectors -> negative area)
        let area_pos = Parallelogram2d::from_slices(&[3.0, 1.0], &[1.0, 2.0]).signed_area();
        let area_neg = Parallelogram2d::from_slices(&[1.0, 2.0], &[3.0, 1.0]).signed_area();
        assert_eq!(area_pos, -area_neg);

        // Colinear vectors -> zero area
        let zero_area = Parallelogram2d::from_slices(&[2.0, 4.0], &[1.0, 2.0]).signed_area();
        assert!((zero_area).abs() < 1e-12);
    }

    #[test]
    fn triangle_functions_work() {
        let p1 = &[0.0, 0.0];
        let p2 = &[3.0, 1.0];
        let p3 = &[1.0, 2.0];

        let triangle = Triangle2d::from_slices(p1, p2, p3);
        let tri_area = triangle.signed_area();
        assert_eq!(tri_area, 2.5);

        let ccw = Triangle2d::from_slices(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0]);
        assert!(ccw.signed_area() > 0.0);

        let cw = Triangle2d::from_slices(&[0.0, 0.0], &[0.0, 1.0], &[1.0, 0.0]);
        assert!(cw.signed_area() < 0.0);

        let colinear = Triangle2d::from_slices(&[0.0, 0.0], &[1.0, 1.0], &[2.0, 2.0]);
        assert_eq!(colinear.signed_area(), 0.0);
    }
}
