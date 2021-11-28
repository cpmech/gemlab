use russell_lab::{Matrix, Vector};

/// Defines a quadrilateral with 8 points (quadratic edges)
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of points
///
/// ```text
/// 3-----6-----2
/// |     s     |           r     s            r     s
/// |     |     |   p:0 [-1.0, -1.0]   p:4 [ 0.0, -1.0]
/// 7     +--r  5   p:1 [ 1.0, -1.0]   p:5 [ 1.0,  0.0]
/// |           |   p:2 [ 1.0,  1.0]   p:6 [ 0.0,  1.0]
/// |           |   p:3 [-1.0,  1.0]   p:7 [-1.0,  0.0]
/// 0-----4-----1
/// ```
///
/// # Local IDs of edges
///
/// ```text
///        2
///  +-----------+         p0 p1 p2
///  |           |     e:0 [0, 1, 4]
///  |           |     e:1 [1, 2, 5]
/// 3|           |1    e:2 [2, 3, 6]
///  |           |     e:3 [3, 0, 7]
///  |           |
///  +-----------+
///        0
/// ```
pub struct Qua8 {}

impl Qua8 {
    pub const NDIM: usize = 2;
    pub const NPOINT: usize = 8;
    pub const NEDGE: usize = 4;
    pub const NFACE: usize = 0;
    pub const EDGE_NPOINT: usize = 3;
    pub const FACE_NPOINT: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const EDGE_POINT_IDS: [[usize; Qua8::EDGE_NPOINT]; Qua8::NEDGE] = [
        [0, 1, 4],
        [1, 2, 5],
        [2, 3, 6],
        [3, 0, 7],
    ];

    #[rustfmt::skip]
    pub const POINT_REFERENCE_COORDS: [[f64; Qua8::NDIM]; Qua8::NPOINT] = [
        [-1.0, -1.0],
        [ 1.0, -1.0],
        [ 1.0,  1.0],
        [-1.0,  1.0],
        [ 0.0, -1.0],
        [ 1.0,  0.0],
        [ 0.0,  1.0],
        [-1.0,  0.0],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        interp[0] = (1.0 - r) * (1.0 - s) * (-r - s - 1.0) / 4.0;
        interp[1] = (1.0 + r) * (1.0 - s) * (r - s - 1.0) / 4.0;
        interp[2] = (1.0 + r) * (1.0 + s) * (r + s - 1.0) / 4.0;
        interp[3] = (1.0 - r) * (1.0 + s) * (-r + s - 1.0) / 4.0;
        interp[4] = (1.0 - s) * (1.0 - r * r) / 2.0;
        interp[5] = (1.0 + r) * (1.0 - s * s) / 2.0;
        interp[6] = (1.0 + s) * (1.0 - r * r) / 2.0;
        interp[7] = (1.0 - r) * (1.0 - s * s) / 2.0;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        deriv[0][0] = -(1.0 - s) * (-r - r - s) / 4.0;
        deriv[1][0] = (1.0 - s) * (r + r - s) / 4.0;
        deriv[2][0] = (1.0 + s) * (r + r + s) / 4.0;
        deriv[3][0] = -(1.0 + s) * (-r - r + s) / 4.0;
        deriv[4][0] = -(1.0 - s) * r;
        deriv[5][0] = (1.0 - s * s) / 2.0;
        deriv[6][0] = -(1.0 + s) * r;
        deriv[7][0] = -(1.0 - s * s) / 2.0;

        deriv[0][1] = -(1.0 - r) * (-s - s - r) / 4.0;
        deriv[1][1] = -(1.0 + r) * (-s - s + r) / 4.0;
        deriv[2][1] = (1.0 + r) * (s + s + r) / 4.0;
        deriv[3][1] = (1.0 - r) * (s + s - r) / 4.0;
        deriv[4][1] = -(1.0 - r * r) / 2.0;
        deriv[5][1] = -(1.0 + r) * s;
        deriv[6][1] = (1.0 - r * r) / 2.0;
        deriv[7][1] = -(1.0 - r) * s;
    }
}
