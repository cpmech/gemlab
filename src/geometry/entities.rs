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
    use super::Circle;

    #[test]
    fn circle_traits_work() {
        let circle = Circle {
            center: [-1.0, -1.0],
            radius: 2.0,
        };
        let clone = circle.clone();
        let correct = "Circle { center: [-1.0, -1.0], radius: 2.0 }";
        assert_eq!(format!("{:?}", circle), correct);
        assert_eq!(format!("{:?}", clone), correct);
    }
}
