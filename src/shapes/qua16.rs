use russell_lab::{Matrix, Vector};

/// Defines a quadrilateral with 8 points
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of points
///
/// ```text
///  3      10       6        2
///    @-----@-------@------@
///    |               (1,1)|
///    |       s ^          |
///  7 @   15@   |    @14   @ 9  
///    |         |          |
///    |         +----> r   |
///    |       (0,0)        |
/// 11 @   12@       @13    @ 5
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
///  |           |     e:0 [0, 1, 4,  8]
///  |           |     e:1 [1, 2, 5,  9]
/// 3|           |1    e:2 [2, 3, 6, 10]
///  |           |     e:3 [3, 0, 7, 11]
///  |           |
///  +-----------+
///        0
/// ```
pub struct Qua16 {}

impl Qua16 {
    pub const NDIM: usize = 2;
    pub const NPOINT: usize = 16;
    pub const NEDGE: usize = 4;
    pub const NFACE: usize = 0;
    pub const EDGE_NPOINT: usize = 4;
    pub const FACE_NPOINT: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const EDGE_POINT_IDS: [[usize; Qua16::EDGE_NPOINT]; Qua16::NEDGE] = [
        [0, 1, 4,  8],
        [1, 2, 5,  9],
        [2, 3, 6, 10],
        [3, 0, 7, 11]
    ];

    #[rustfmt::skip]
    pub const POINT_REFERENCE_COORDS: [[f64; Qua16::NDIM]; Qua16::NPOINT] = [
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

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &Vector) {
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

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &Vector) {
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

        // derivs of Lin4 interp w.r.t r
        let dnn_dr0 = (-27.0 * r * r + 18.0 * r + 1.0) / 16.0;
        let dnn_dr1 = (27.0 * r * r + 18.0 * r - 1.0) / 16.0;
        let dnn_dr2 = (81.0 * r * r - 18.0 * r - 27.0) / 16.0;
        let dnn_dr3 = (-81.0 * r * r - 18.0 * r + 27.0) / 16.0;

        // derivs of Lin4 interp w.r.t s
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
