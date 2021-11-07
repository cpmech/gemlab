//! Contains basic geometric features and algorithms

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Circle {
    pub center: [f64; 2],
    pub radius: f64,
    pub tolerance: f64,
}

mod algorithms;
pub use crate::geometry::algorithms::*;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circle_traits_work() {
        let circle = Circle {
            center: [-1.0, -1.0],
            radius: 2.0,
            tolerance: 1e-3,
        };
        let clone = circle.clone();
        assert_eq!(circle, clone);
        assert_eq!(
            format!("{:?}", circle),
            "Circle { center: [-1.0, -1.0], radius: 2.0, tolerance: 0.001 }"
        );
    }
}
