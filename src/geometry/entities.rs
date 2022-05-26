/// Holds data defining a circle in 2D
#[derive(Clone, Debug)]
pub struct Circle {
    pub center: [f64; 2],
    pub radius: f64,
    pub tolerance: f64,
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Circle;

    #[test]
    fn circle_traits_work() {
        let circle = Circle {
            center: [-1.0, -1.0],
            radius: 2.0,
            tolerance: 1e-3,
        };
        let clone = circle.clone();
        let correct = "Circle { center: [-1.0, -1.0], radius: 2.0, tolerance: 0.001 }";
        assert_eq!(format!("{:?}", circle), correct);
        assert_eq!(format!("{:?}", clone), correct);
    }
}