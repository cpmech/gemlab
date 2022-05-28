use russell_lab::{Matrix, Vector};

/// Defines a quadrilateral with 4 nodes (linear edges)
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of nodes
///
/// ```text
/// 3-----------2
/// |     s     |              r     s
/// |     |     |      p:0 [-1.0, -1.0]
/// |     +--r  |      p:1 [ 1.0, -1.0]
/// |           |      p:2 [ 1.0,  1.0]
/// |           |      p:3 [-1.0,  1.0]
/// 0-----------1
/// ```
///
/// # Local IDs of edges
///
/// ```text
///        2
///  +-----------+         p0  p1
///  |           |     e:0 [1, 0]
///  |           |     e:1 [2, 1]
/// 3|           |1    e:2 [3, 2]
///  |           |     e:3 [0, 3]
///  |           |
///  +-----------+
///        0
/// ```
///
/// # Note about edges
///
/// * The order of edge nodes is such that the normals are outward
/// * The order of edge nodes corresponds to **Lin2** nodes
pub struct Qua4 {}

impl Qua4 {
    pub const GEO_NDIM: usize = 2;
    pub const NNODE: usize = 4;
    pub const NEDGE: usize = 4;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 2;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Qua4::EDGE_NNODE]; Qua4::NEDGE] = [
        [1, 0],
        [2, 1],
        [3, 2],
        [0, 3],
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Qua4::GEO_NDIM]; Qua4::NNODE] = [
        [-1.0, -1.0],
        [ 1.0, -1.0],
        [ 1.0,  1.0],
        [-1.0,  1.0],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        interp[0] = (1.0 - r - s + r * s) / 4.0;
        interp[1] = (1.0 + r - s - r * s) / 4.0;
        interp[2] = (1.0 + r + s + r * s) / 4.0;
        interp[3] = (1.0 - r + s - r * s) / 4.0;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        deriv[0][0] = (-1.0 + s) / 4.0;
        deriv[0][1] = (-1.0 + r) / 4.0;

        deriv[1][0] = (1.0 - s) / 4.0;
        deriv[1][1] = (-1.0 - r) / 4.0;

        deriv[2][0] = (1.0 + s) / 4.0;
        deriv[2][1] = (1.0 + r) / 4.0;

        deriv[3][0] = (-1.0 - s) / 4.0;
        deriv[3][1] = (1.0 - r) / 4.0;
    }
}
