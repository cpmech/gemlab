use russell_lab::{Matrix, Vector};

/// Defines a quadrilateral with 12 nodes (cubic edges)
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of nodes
///
/// ```text
///  3      10       6        2
///    @-----@-------@------@
///    |               (1,1)|
///    |       s ^          |
///  7 @         |          @ 9
///    |         |          |
///    |         +----> r   |
///    |       (0,0)        |
/// 11 @                    @ 5
///    |                    |
///    |(-1,-1)             |
///    @-----@-------@------@
///  0       4       8        1
/// ```
///
/// # Local IDs of edges
///
/// ```text
///        2
///  +-----------+         p0 p1  p2 p3
///  |           |     e:0 [1, 0,  8, 4]
///  |           |     e:1 [2, 1,  9, 5]
/// 3|           |1    e:2 [3, 2, 10, 6]
///  |           |     e:3 [0, 3, 11, 7]
///  |           |
///  +-----------+
///        0
/// ```
///
/// * The order of edge nodes is such that the normals are outward
/// * The order of edge nodes corresponds to **Lin4** nodes
pub struct Qua12 {}

impl Qua12 {
    pub const GEO_NDIM: usize = 2;
    pub const NNODE: usize = 12;
    pub const NEDGE: usize = 4;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 4;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Qua12::EDGE_NNODE]; Qua12::NEDGE] = [
        [1, 0,  8, 4],
        [2, 1,  9, 5],
        [3, 2, 10, 6],
        [0, 3, 11, 7]
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Qua12::GEO_NDIM]; Qua12::NNODE] = [
        [-1.0       , -1.0       ],
        [ 1.0       , -1.0       ],
        [ 1.0       ,  1.0       ],
        [-1.0       ,  1.0       ],
        [-1.0 / 3.0 , -1.0       ],
        [ 1.0       , -1.0 / 3.0 ],
        [ 1.0 / 3.0 ,  1.0       ],
        [-1.0       ,  1.0 / 3.0 ],
        [ 1.0 / 3.0 , -1.0       ],
        [ 1.0       ,  1.0 / 3.0 ],
        [-1.0 / 3.0 ,  1.0       ],
        [-1.0       , -1.0 / 3.0 ],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        let rm = 1.0 - r;
        let rp = 1.0 + r;
        let sm = 1.0 - s;
        let sp = 1.0 + s;

        interp[0] = rm * sm * (9.0 * (r * r + s * s) - 10.0) / 32.0;
        interp[1] = rp * sm * (9.0 * (r * r + s * s) - 10.0) / 32.0;
        interp[2] = rp * sp * (9.0 * (r * r + s * s) - 10.0) / 32.0;
        interp[3] = rm * sp * (9.0 * (r * r + s * s) - 10.0) / 32.0;
        interp[4] = 9.0 * (1.0 - r * r) * (1.0 - 3.0 * r) * sm / 32.0;
        interp[5] = 9.0 * (1.0 - s * s) * (1.0 - 3.0 * s) * rp / 32.0;
        interp[6] = 9.0 * (1.0 - r * r) * (1.0 + 3.0 * r) * sp / 32.0;
        interp[7] = 9.0 * (1.0 - s * s) * (1.0 + 3.0 * s) * rm / 32.0;
        interp[8] = 9.0 * (1.0 - r * r) * (1.0 + 3.0 * r) * sm / 32.0;
        interp[9] = 9.0 * (1.0 - s * s) * (1.0 + 3.0 * s) * rp / 32.0;
        interp[10] = 9.0 * (1.0 - r * r) * (1.0 - 3.0 * r) * sp / 32.0;
        interp[11] = 9.0 * (1.0 - s * s) * (1.0 - 3.0 * s) * rm / 32.0;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        let rm = 1.0 - r;
        let rp = 1.0 + r;
        let sm = 1.0 - s;
        let sp = 1.0 + s;

        deriv[0][0] = sm * (9.0 * (2.0 * r - 3.0 * r * r - s * s) + 10.0) / 32.0;
        deriv[1][0] = sm * (9.0 * (2.0 * r + 3.0 * r * r + s * s) - 10.0) / 32.0;
        deriv[2][0] = sp * (9.0 * (2.0 * r + 3.0 * r * r + s * s) - 10.0) / 32.0;
        deriv[3][0] = sp * (9.0 * (2.0 * r - 3.0 * r * r - s * s) + 10.0) / 32.0;
        deriv[4][0] = 9.0 * sm * (9.0 * r * r - 2.0 * r - 3.0) / 32.0;
        deriv[5][0] = 9.0 * (1.0 - s * s) * (1.0 - 3.0 * s) / 32.0;
        deriv[6][0] = 9.0 * sp * (-9.0 * r * r - 2.0 * r + 3.0) / 32.0;
        deriv[7][0] = -9.0 * (1.0 - s * s) * (1.0 + 3.0 * s) / 32.0;
        deriv[8][0] = 9.0 * sm * (-9.0 * r * r - 2.0 * r + 3.0) / 32.0;
        deriv[9][0] = 9.0 * (1.0 - s * s) * (1.0 + 3.0 * s) / 32.0;
        deriv[10][0] = 9.0 * sp * (9.0 * r * r - 2.0 * r - 3.0) / 32.0;
        deriv[11][0] = -9.0 * (1.0 - s * s) * (1.0 - 3.0 * s) / 32.0;

        deriv[0][1] = rm * (9.0 * (2.0 * s - 3.0 * s * s - r * r) + 10.0) / 32.0;
        deriv[1][1] = rp * (9.0 * (2.0 * s - 3.0 * s * s - r * r) + 10.0) / 32.0;
        deriv[2][1] = rp * (9.0 * (2.0 * s + 3.0 * s * s + r * r) - 10.0) / 32.0;
        deriv[3][1] = rm * (9.0 * (2.0 * s + 3.0 * s * s + r * r) - 10.0) / 32.0;
        deriv[4][1] = -9.0 * (1.0 - r * r) * (1.0 - 3.0 * r) / 32.0;
        deriv[5][1] = 9.0 * rp * (9.0 * s * s - 2.0 * s - 3.0) / 32.0;
        deriv[6][1] = 9.0 * (1.0 - r * r) * (1.0 + 3.0 * r) / 32.0;
        deriv[7][1] = 9.0 * rm * (-9.0 * s * s - 2.0 * s + 3.0) / 32.0;
        deriv[8][1] = -9.0 * (1.0 - r * r) * (1.0 + 3.0 * r) / 32.0;
        deriv[9][1] = 9.0 * rp * (-9.0 * s * s - 2.0 * s + 3.0) / 32.0;
        deriv[10][1] = 9.0 * (1.0 - r * r) * (1.0 - 3.0 * r) / 32.0;
        deriv[11][1] = 9.0 * rm * (9.0 * s * s - 2.0 * s - 3.0) / 32.0;
    }
}
