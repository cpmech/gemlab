use russell_lab::{Matrix, Vector};

/// Defines a line (segment) with 4 points (cubic functions)
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of points
///
/// ```text
/// -1                    +1
///  @------@-------@------@  --> r
///  0      2       3      1
/// ```
pub struct Lin4 {}

impl Lin4 {
    pub const NDIM: usize = 1;
    pub const NPOINT: usize = 4;
    pub const NEDGE: usize = 0;
    pub const NFACE: usize = 0;
    pub const EDGE_NPOINT: usize = 0;
    pub const FACE_NPOINT: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const POINT_REFERENCE_COORDS: [[f64; Lin4::NDIM]; Lin4::NPOINT] = [
        [-1.0],
        [ 1.0],
        [-1.0 / 3.0],
        [ 1.0 / 3.0],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let r = ksi[0];

        interp[0] = (-9.0 * r * r * r + 9.0 * r * r + r - 1.0) / 16.0;
        interp[1] = (9.0 * r * r * r + 9.0 * r * r - r - 1.0) / 16.0;
        interp[2] = (27.0 * r * r * r - 9.0 * r * r - 27.0 * r + 9.0) / 16.0;
        interp[3] = (-27.0 * r * r * r - 9.0 * r * r + 27.0 * r + 9.0) / 16.0;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &[f64]) {
        let r = ksi[0];

        deriv[0][0] = 1.0 / 16.0 * (-27.0 * r * r + 18.0 * r + 1.0);
        deriv[1][0] = 1.0 / 16.0 * (27.0 * r * r + 18.0 * r - 1.0);
        deriv[2][0] = 1.0 / 16.0 * (81.0 * r * r - 18.0 * r - 27.0);
        deriv[3][0] = 1.0 / 16.0 * (-81.0 * r * r - 18.0 * r + 27.0);
    }
}
