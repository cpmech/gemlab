use russell_lab::{Matrix, Vector};

/// Defines a line (segment) with 2 nodes (linear functions)
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of nodes
///
/// ```text
/// -1                    +1
///  @---------------------@  --> r
///  0                     1
/// ```
pub struct Lin2 {}

impl Lin2 {
    pub const NDIM: usize = 1;
    pub const NNODE: usize = 2;
    pub const NEDGE: usize = 0;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 0;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Lin2::NDIM]; Lin2::NNODE] = [
        [-1.0],
        [ 1.0],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let r = ksi[0];

        interp[0] = 0.5 * (1.0 - r);
        interp[1] = 0.5 * (1.0 + r);
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, _: &[f64]) {
        deriv[0][0] = -0.5;
        deriv[1][0] = 0.5;
    }
}
