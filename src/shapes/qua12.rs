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
///  +-----------+         p0 p1 p2  p3
///  |           |     e:0 [0, 1, 4,  8]
///  |           |     e:1 [1, 2, 5,  9]
/// 3|           |1    e:2 [2, 3, 6, 10]
///  |           |     e:3 [3, 0, 7, 11]
///  |           |
///  +-----------+
///        0
/// ```
pub struct Qua12 {}

impl Qua12 {
    pub const NDIM: usize = 2;
    pub const NPOINT: usize = 12;
    pub const NEDGE: usize = 4;
    pub const NFACE: usize = 0;
    pub const EDGE_NPOINT: usize = 4;
    pub const FACE_NPOINT: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const EDGE_POINT_IDS: [[usize; Qua12::EDGE_NPOINT]; Qua12::NEDGE] = [
        [0, 1, 4,  8],
        [1, 2, 5,  9],
        [2, 3, 6, 10],
        [3, 0, 7, 11]
    ];

    #[rustfmt::skip]
    pub const POINT_REFERENCE_COORDS: [[f64; Qua12::NDIM]; Qua12::NPOINT] = [
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
    pub fn calc_interp(interp: &mut Vector, ksi: &Vector) {
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
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &Vector) {
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
