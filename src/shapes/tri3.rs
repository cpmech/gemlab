use russell_lab::{Matrix, Vector};

/// Defines a triangle with 3 nodes (linear edges)
///
/// # Local IDs of nodes
///
/// ```text
/// s
/// |
/// 2, (0,1)
/// | ',
/// |   ',         node    r    s
/// |     ',          0  0.0  0.0
/// |       ',        1  1.0  0.0
/// |         ',      2  0.0  1.0
/// |           ',
/// |             ',
/// |               ',
/// | (0,0)           ', (1,0)
/// 0-------------------1 ---- r
/// ```
///
/// # Local IDs of edges
///
/// ```text
///  |\             p0  p1
///  | \        e:0 [1, 0]
///  |  \ 1     e:1 [2, 1]
/// 2|   \      e:2 [0, 2]
///  |    \
///  |_____\
///     0
/// ```
///
/// * The order of edge nodes is such that the normals are outward
/// * The order of edge nodes corresponds to **Lin2** nodes
pub struct Tri3 {}

impl Tri3 {
    pub const GEO_NDIM: usize = 2;
    pub const NNODE: usize = 3;
    pub const NEDGE: usize = 3;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 2;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Tri3::EDGE_NNODE]; Tri3::NEDGE] = [
        [1, 0],
        [2, 1],
        [0, 2],
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Tri3::GEO_NDIM]; Tri3::NNODE] = [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        interp[0] = 1.0 - r - s;
        interp[1] = r;
        interp[2] = s;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, _: &[f64]) {
        deriv[0][0] = -1.0;
        deriv[1][0] = 1.0;
        deriv[2][0] = 0.0;

        deriv[0][1] = -1.0;
        deriv[1][1] = 0.0;
        deriv[2][1] = 1.0;
    }
}
