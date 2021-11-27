use russell_lab::{Matrix, Vector};

/// Defines a line (segment) with 2 points (linear functions)
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of points
///
/// ```text
/// -1                    +1
///  @---------------------@  --> r
///  0                     1
/// ```
pub struct Lin2 {}

impl Lin2 {
    pub const NDIM: usize = 1;
    pub const NPOINT: usize = 2;
    pub const NEDGE: usize = 0;
    pub const NFACE: usize = 0;
    pub const EDGE_NPOINT: usize = 0;
    pub const FACE_NPOINT: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const POINT_REFERENCE_COORDS: [[f64; Lin2::NDIM]; Lin2::NPOINT] = [
        [-1.0],
        [ 1.0],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &Vector) {
        let r = ksi[0];

        interp[0] = 0.5 * (1.0 - r);
        interp[1] = 0.5 * (1.0 + r);
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, _: &Vector) {
        deriv[0][0] = -0.5;
        deriv[1][0] = 0.5;
    }
}
