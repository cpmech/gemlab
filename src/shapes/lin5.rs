use russell_lab::{Matrix, Vector};

/// Defines a line (segment) with 5 points (quartic functions)
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of points
///
/// ```text
/// -1  -0.5    0    0.5  +1
///  @----@-----@-----@----@  --> r
///  0    3     2     4    1
/// ```
pub struct Lin5 {}

impl Lin5 {
    pub const NDIM: usize = 1;
    pub const NPOINT: usize = 5;
    pub const NEDGE: usize = 0;
    pub const NFACE: usize = 0;
    pub const EDGE_NPOINT: usize = 0;
    pub const FACE_NPOINT: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const POINT_REFERENCE_COORDS: [[f64; Lin5::NDIM]; Lin5::NPOINT] = [
        [-1.0],
        [ 1.0],
        [ 0.0],
        [-0.5],
        [ 0.5],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &Vector) {
        let r = ksi[0];

        interp[0] = (r - 1.0) * (1.0 - 2.0 * r) * r * (-1.0 - 2.0 * r) / 6.0;
        interp[1] = (1.0 - 2.0 * r) * r * (-1.0 - 2.0 * r) * (1.0 + r) / 6.0;
        interp[2] = (1.0 - r) * (1.0 - 2.0 * r) * (-1.0 - 2.0 * r) * (-1.0 - r);
        interp[3] = 4.0 * (1.0 - r) * (1.0 - 2.0 * r) * r * (-1.0 - r) / 3.0;
        interp[4] = 4.0 * (1.0 - r) * r * (-1.0 - 2.0 * r) * (-1.0 - r) / 3.0;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &Vector) {
        let r = ksi[0];

        deriv[0][0] = -((1.0 - 2.0 * r) * (r - 1.0) * r) / 3.0 - ((-2.0 * r - 1.0) * (r - 1.0) * r) / 3.0
            + ((-2.0 * r - 1.0) * (1.0 - 2.0 * r) * r) / 6.0
            + ((-2.0 * r - 1.0) * (1.0 - 2.0 * r) * (r - 1.0)) / 6.0;

        deriv[1][0] = -((1.0 - 2.0 * r) * r * (r + 1.0)) / 3.0 - ((-2.0 * r - 1.0) * r * (r + 1.0)) / 3.0
            + ((-2.0 * r - 1.0) * (1.0 - 2.0 * r) * (r + 1.0)) / 6.0
            + ((-2.0 * r - 1.0) * (1.0 - 2.0 * r) * r) / 6.0;

        deriv[2][0] = -2.0 * (1.0 - 2.0 * r) * (-r - 1.0) * (1.0 - r)
            - 2.0 * (-2.0 * r - 1.0) * (-r - 1.0) * (1.0 - r)
            - (-2.0 * r - 1.0) * (1.0 - 2.0 * r) * (1.0 - r)
            - (-2.0 * r - 1.0) * (1.0 - 2.0 * r) * (-r - 1.0);

        deriv[3][0] = -(8.0 * (-r - 1.0) * (1.0 - r) * r) / 3.0
            - (4.0 * (1.0 - 2.0 * r) * (1.0 - r) * r) / 3.0
            - (4.0 * (1.0 - 2.0 * r) * (-r - 1.0) * r) / 3.0
            + (4.0 * (1.0 - 2.0 * r) * (-r - 1.0) * (1.0 - r)) / 3.0;

        deriv[4][0] = -(8.0 * (-r - 1.0) * (1.0 - r) * r) / 3.0
            - (4.0 * (-2.0 * r - 1.0) * (1.0 - r) * r) / 3.0
            - (4.0 * (-2.0 * r - 1.0) * (-r - 1.0) * r) / 3.0
            + (4.0 * (-2.0 * r - 1.0) * (-r - 1.0) * (1.0 - r)) / 3.0;
    }
}
