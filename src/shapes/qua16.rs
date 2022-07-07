use russell_lab::{Matrix, Vector};

/// Defines a quadrilateral with 16 nodes (cubic edges; interior nodes)
///
/// ![qua16](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_qua16.svg)
///
/// # Local IDs of nodes
///
/// ```text
///  3----10-------6------2
///  |               (1,1)|
///  |       s ^          |
///  7    15   |   14     9
///  |         |          |
///  |         +----> r   |
///  |       (0,0)        |
/// 11    12      13      5
///  |                    |
///  |(-1,-1)             |
///  0-----4-------8------1
/// ```
///
/// # Local IDs of edges
///
/// ```text
///              2
///    3----10-------6------2
///    |                    |
///    |                    |
///    7    15      14      9          p0 p1  p2 p3
///    |                    |      e:0 [1, 0,  8, 4]
/// 3  |                    | 1    e:1 [2, 1,  9, 5]
///    |                    |      e:2 [3, 2, 10, 6]
///   11    12      13      5      e:3 [0, 3, 11, 7]
///    |                    |
///    |                    |
///    0-----4-------8------1
///              0
/// ```
///
/// # Notes
///
/// * The reference coordinates range from -1 to +1 with the geometry centred @ 0
/// * The order of edge nodes is such that the normals are outward
/// * The order of edge nodes corresponds to [super::Lin4] nodes
pub struct Qua16 {}

impl Qua16 {
    pub const GEO_NDIM: usize = 2;
    pub const NNODE: usize = 16;
    pub const NEDGE: usize = 4;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 4;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;
    pub const N_INTERIOR_NODE: usize = 4;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Qua16::EDGE_NNODE]; Qua16::NEDGE] = [
        [1, 0,  8, 4],
        [2, 1,  9, 5],
        [3, 2, 10, 6],
        [0, 3, 11, 7]
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Qua16::GEO_NDIM]; Qua16::NNODE] = [
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
        [-1.0 / 3.0 , -1.0 / 3.0 ],
        [ 1.0 / 3.0 , -1.0 / 3.0 ],
        [ 1.0 / 3.0 ,  1.0 / 3.0 ],
        [-1.0 / 3.0 ,  1.0 / 3.0 ],
    ];

    pub const INTERIOR_NODES: [usize; Qua16::N_INTERIOR_NODE] = [12, 13, 14, 15];

    /// Computes the interpolation functions
    ///
    /// # Output
    ///
    /// * `interp` -- interpolation function evaluated at ksi (nnode)
    ///
    /// # Input
    ///
    /// * `ksi` -- reference coordinates with length ≥ geo_ndim
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        // interp fn of Lin4 along r
        let nn_r0 = (-9.0 * r * r * r + 9.0 * r * r + r - 1.0) / 16.0;
        let nn_r1 = (9.0 * r * r * r + 9.0 * r * r - r - 1.0) / 16.0;
        let nn_r2 = (27.0 * r * r * r - 9.0 * r * r - 27.0 * r + 9.0) / 16.0;
        let nn_r3 = (-27.0 * r * r * r - 9.0 * r * r + 27.0 * r + 9.0) / 16.0;

        // interp fn of Lin4 along s
        let nn_s0 = (-9.0 * s * s * s + 9.0 * s * s + s - 1.0) / 16.0;
        let nn_s1 = (9.0 * s * s * s + 9.0 * s * s - s - 1.0) / 16.0;
        let nn_s2 = (27.0 * s * s * s - 9.0 * s * s - 27.0 * s + 9.0) / 16.0;
        let nn_s3 = (-27.0 * s * s * s - 9.0 * s * s + 27.0 * s + 9.0) / 16.0;

        interp[0] = nn_r0 * nn_s0;
        interp[1] = nn_r1 * nn_s0;
        interp[2] = nn_r1 * nn_s1;
        interp[3] = nn_r0 * nn_s1;
        interp[4] = nn_r2 * nn_s0;
        interp[5] = nn_r1 * nn_s2;
        interp[6] = nn_r3 * nn_s1;
        interp[7] = nn_r0 * nn_s3;
        interp[8] = nn_r3 * nn_s0;
        interp[9] = nn_r1 * nn_s3;
        interp[10] = nn_r2 * nn_s1;
        interp[11] = nn_r0 * nn_s2;
        interp[12] = nn_r2 * nn_s2;
        interp[13] = nn_r3 * nn_s2;
        interp[14] = nn_r3 * nn_s3;
        interp[15] = nn_r2 * nn_s3;
    }

    /// Computes the derivatives of interpolation functions with respect to the reference coordinates
    ///
    /// # Output
    ///
    /// * `deriv` -- derivatives of the interpolation function with respect to
    ///   the reference coordinates ksi, evaluated at ksi (nnode,geo_ndim)
    ///
    /// # Input
    ///
    /// * `ksi` -- reference coordinates with length ≥ geo_ndim
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        // interp fn of Lin4 along r
        let nn_r0 = (-9.0 * r * r * r + 9.0 * r * r + r - 1.0) / 16.0;
        let nn_r1 = (9.0 * r * r * r + 9.0 * r * r - r - 1.0) / 16.0;
        let nn_r2 = (27.0 * r * r * r - 9.0 * r * r - 27.0 * r + 9.0) / 16.0;
        let nn_r3 = (-27.0 * r * r * r - 9.0 * r * r + 27.0 * r + 9.0) / 16.0;

        // interp fn of Lin4 along s
        let nn_s0 = (-9.0 * s * s * s + 9.0 * s * s + s - 1.0) / 16.0;
        let nn_s1 = (9.0 * s * s * s + 9.0 * s * s - s - 1.0) / 16.0;
        let nn_s2 = (27.0 * s * s * s - 9.0 * s * s - 27.0 * s + 9.0) / 16.0;
        let nn_s3 = (-27.0 * s * s * s - 9.0 * s * s + 27.0 * s + 9.0) / 16.0;

        // derivatives of Lin4 interp with respect to r
        let dnn_dr0 = (-27.0 * r * r + 18.0 * r + 1.0) / 16.0;
        let dnn_dr1 = (27.0 * r * r + 18.0 * r - 1.0) / 16.0;
        let dnn_dr2 = (81.0 * r * r - 18.0 * r - 27.0) / 16.0;
        let dnn_dr3 = (-81.0 * r * r - 18.0 * r + 27.0) / 16.0;

        // derivatives of Lin4 interp with respect to s
        let dnn_ds0 = (-27.0 * s * s + 18.0 * s + 1.0) / 16.0;
        let dnn_ds1 = (27.0 * s * s + 18.0 * s - 1.0) / 16.0;
        let dnn_ds2 = (81.0 * s * s - 18.0 * s - 27.0) / 16.0;
        let dnn_ds3 = (-81.0 * s * s - 18.0 * s + 27.0) / 16.0;

        deriv[0][0] = dnn_dr0 * nn_s0;
        deriv[1][0] = dnn_dr1 * nn_s0;
        deriv[2][0] = dnn_dr1 * nn_s1;
        deriv[3][0] = dnn_dr0 * nn_s1;
        deriv[4][0] = dnn_dr2 * nn_s0;
        deriv[5][0] = dnn_dr1 * nn_s2;
        deriv[6][0] = dnn_dr3 * nn_s1;
        deriv[7][0] = dnn_dr0 * nn_s3;
        deriv[8][0] = dnn_dr3 * nn_s0;
        deriv[9][0] = dnn_dr1 * nn_s3;
        deriv[10][0] = dnn_dr2 * nn_s1;
        deriv[11][0] = dnn_dr0 * nn_s2;
        deriv[12][0] = dnn_dr2 * nn_s2;
        deriv[13][0] = dnn_dr3 * nn_s2;
        deriv[14][0] = dnn_dr3 * nn_s3;
        deriv[15][0] = dnn_dr2 * nn_s3;

        deriv[0][1] = nn_r0 * dnn_ds0;
        deriv[1][1] = nn_r1 * dnn_ds0;
        deriv[2][1] = nn_r1 * dnn_ds1;
        deriv[3][1] = nn_r0 * dnn_ds1;
        deriv[4][1] = nn_r2 * dnn_ds0;
        deriv[5][1] = nn_r1 * dnn_ds2;
        deriv[6][1] = nn_r3 * dnn_ds1;
        deriv[7][1] = nn_r0 * dnn_ds3;
        deriv[8][1] = nn_r3 * dnn_ds0;
        deriv[9][1] = nn_r1 * dnn_ds3;
        deriv[10][1] = nn_r2 * dnn_ds1;
        deriv[11][1] = nn_r0 * dnn_ds2;
        deriv[12][1] = nn_r2 * dnn_ds2;
        deriv[13][1] = nn_r3 * dnn_ds2;
        deriv[14][1] = nn_r3 * dnn_ds3;
        deriv[15][1] = nn_r2 * dnn_ds3;
    }
}
