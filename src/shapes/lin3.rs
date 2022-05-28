use russell_lab::{Matrix, Vector};

/// Defines a line (segment) with 3 nodes (quadratic functions)
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of nodes
///
/// ```text
/// -1                    +1
///  @----------@----------@  --> r
///  0          2          1
/// ```
pub struct Lin3 {}

impl Lin3 {
    pub const GEO_NDIM: usize = 1;
    pub const NNODE: usize = 3;
    pub const NEDGE: usize = 0;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 0;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Lin3::GEO_NDIM]; Lin3::NNODE] = [
        [-1.0],
        [ 1.0],
        [ 0.0],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let r = ksi[0];

        interp[0] = 0.5 * (r * r - r);
        interp[1] = 0.5 * (r * r + r);
        interp[2] = 1.0 - r * r;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &[f64]) {
        let r = ksi[0];

        deriv[0][0] = r - 0.5;
        deriv[1][0] = r + 0.5;
        deriv[2][0] = -2.0 * r;
    }
}
