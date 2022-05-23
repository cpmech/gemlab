use crate::shapes::{Shape, StateOfShape};
use crate::util::SQRT_3;

/// Implements functions to assist in verifications
pub struct Verification;

impl Verification {
    /// Returns Shape and StateOfShape for a Lin2 line segment
    pub fn line_segment_lin2(l: f64) -> (Shape, StateOfShape) {
        let xa = 3.0;
        let xb = xa + l;
        let y = 4.0;
        let shape = Shape::new(2, 1, 2).unwrap();
        let state = StateOfShape::new(shape.geo_ndim, &[[xa, y], [xb, y]]).unwrap();
        (shape, state)
    }

    /// Returns (Shape,StateOfShape,area) for a Tri3 equilateral triangle
    ///
    /// ```text
    /// Equilateral triangle
    ///
    ///           /
    ///        2   \
    ///       / \   \
    ///      / ↑ \   l
    ///     /  h  \   \
    ///    /   ↓   \   \
    ///   /         \   /
    ///  0-----------1
    ///
    ///  |--s--|--s--|
    ///
    ///  |-----l-----|
    /// ```
    pub fn equilateral_triangle_tri3(l: f64) -> (Shape, StateOfShape, f64) {
        let s = l / 2.0;
        let h = l * SQRT_3 / 2.0;
        let (x0, y0) = (3.0, 4.0);
        let (x1, y1) = (x0 + l, y0);
        let (x2, y2) = (x0 + s, y0 + h);
        let shape = Shape::new(2, 2, 3).unwrap();
        let state = StateOfShape::new(shape.geo_ndim, &[[x0, y0], [x1, y1], [x2, y2]]).unwrap();
        let area = l * h / 2.0;
        (shape, state, area)
    }

    /// Returns (Shape,StateOfShape,area) for a Tri6 equilateral triangle
    ///
    /// ```text
    /// Equilateral triangle
    ///
    ///           /
    ///        2   \
    ///       / \   \
    ///      / ↑ \   l
    ///     5  h  4   \
    ///    /   ↓   \   \
    ///   /         \   /
    ///  0-----3-----1
    ///
    ///  |--s--|--s--|
    ///
    ///  |-----l-----|
    /// ```
    pub fn equilateral_triangle_tri6(l: f64) -> (Shape, StateOfShape, f64) {
        let s = l / 2.0;
        let h = l * SQRT_3 / 2.0;
        let (x0, y0) = (3.0, 4.0);
        let (x1, y1) = (x0 + l, y0);
        let (x2, y2) = (x0 + s, y0 + h);
        let (x3, y3) = (x0 + s, y0);
        let (x4, y4) = (x0 + 1.5 * s, y0 + 0.5 * h);
        let (x5, y5) = (x0 + 0.5 * s, y0 + 0.5 * h);
        let shape = Shape::new(2, 2, 6).unwrap();
        let state = StateOfShape::new(
            shape.geo_ndim,
            &[[x0, y0], [x1, y1], [x2, y2], [x3, y3], [x4, y4], [x5, y5]],
        )
        .unwrap();
        let area = l * h / 2.0;
        (shape, state, area)
    }
}
