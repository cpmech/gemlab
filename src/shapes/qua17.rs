use russell_lab::{Matrix, Vector};

/// Defines a quadrilateral with 17 nodes (quartic edge; interior node)
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of nodes
///
/// ```text
///   3      14    10     6     2
///    @-----@-----@-----@-----@
///    |                  (1,1)|
///    |                       |
///  7 @                       @ 13
///    |         s ^           |
///    |           |           |
/// 11 @           |16         @ 9
///    |           @----> r    |
///    |         (0,0)         |
/// 15 @                       @ 5
///    |                       |
///    |(-1,-1)                |
///    @-----@-----@-----@-----@
///   0      4     8    12      1
/// ```
///
/// # Local IDs of edges
///
/// ```text
///        2
///  +-----------+         p0 p1  p2 p3  p4
///  |           |     e:0 [0, 1,  8, 4, 12]
///  |           |     e:1 [1, 2,  9, 5, 13]
/// 3|           |1    e:2 [2, 3, 10, 6, 14]
///  |           |     e:3 [3, 0, 11, 7, 15]
///  |           |
///  +-----------+
///        0
/// ```
pub struct Qua17 {}

impl Qua17 {
    pub const NDIM: usize = 2;
    pub const NNODE: usize = 17;
    pub const NEDGE: usize = 4;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 5;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Qua17::EDGE_NNODE]; Qua17::NEDGE] = [
        [0, 1,  8, 4, 12],
        [1, 2,  9, 5, 13],
        [2, 3, 10, 6, 14],
        [3, 0, 11, 7, 15]
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Qua17::NDIM]; Qua17::NNODE] = [
        [-1.0, -1.0], //  0
        [ 1.0, -1.0], //  1
        [ 1.0,  1.0], //  2
        [-1.0,  1.0], //  3
        [-0.5, -1.0], //  4
        [ 1.0, -0.5], //  5
        [ 0.5,  1.0], //  6
        [-1.0,  0.5], //  7
        [ 0.0, -1.0], //  8
        [ 1.0,  0.0], //  9
        [ 0.0,  1.0], // 10
        [-1.0,  0.0], // 11
        [ 0.5, -1.0], // 12
        [ 1.0,  0.5], // 13
        [-0.5,  1.0], // 14
        [-1.0, -0.5], // 15
        [ 0.0,  0.0], // 16
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        let a = 2.0 / 3.0;
        let rr = r * r;
        let ss = s * s;
        let rs = r * s;
        let rp = 1.0 + r;
        let rm = 1.0 - r;
        let sp = 1.0 + s;
        let sm = 1.0 - s;

        interp[0] = rm * sm * (-4.0 * r * (rr - 1.0) - 4.0 * s * (ss - 1.0) + 3.0 * rs) / 12.0;
        interp[1] = rp * sm * (4.0 * r * (rr - 1.0) - 4.0 * s * (ss - 1.0) - 3.0 * rs) / 12.0;
        interp[2] = rp * sp * (4.0 * r * (rr - 1.0) + 4.0 * s * (ss - 1.0) + 3.0 * rs) / 12.0;
        interp[3] = rm * sp * (-4.0 * r * (rr - 1.0) + 4.0 * s * (ss - 1.0) - 3.0 * rs) / 12.0;
        interp[4] = -a * r * sm * rm * rp * (1.0 - 2.0 * r);
        interp[5] = -a * s * rp * sm * sp * (1.0 - 2.0 * s);
        interp[6] = a * r * sp * rm * rp * (1.0 + 2.0 * r);
        interp[7] = a * s * sm * sp * (1.0 + 2.0 * s) * rm;
        interp[8] = 0.5 * rm * rp * (-s - 4.0 * rr) * sm;
        interp[9] = 0.5 * sm * sp * (r - 4.0 * ss) * rp;
        interp[10] = 0.5 * rm * rp * (s - 4.0 * rr) * sp;
        interp[11] = 0.5 * sm * sp * (-r - 4.0 * ss) * rm;
        interp[12] = a * r * sm * rm * rp * (1.0 + 2.0 * r);
        interp[13] = a * s * rp * sm * sp * (1.0 + 2.0 * s);
        interp[14] = -a * r * sp * rm * rp * (1.0 - 2.0 * r);
        interp[15] = -a * s * rm * sm * sp * (1.0 - 2.0 * s);
        interp[16] = rm * rp * sm * sp
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        let a = 2.0 / 3.0;
        let rr = r * r;
        let ss = s * s;
        let rs = r * s;
        let rp = 1.0 + r;
        let rm = 1.0 - r;
        let sp = 1.0 + s;
        let sm = 1.0 - s;

        let b = 1.0 / 12.0;
        let r1 = r - 1.0;
        let rrr = rr * r;
        let sss = ss * s;

        deriv[0][0] = b * sm * (16.0 * rrr - 12.0 * rr - 6.0 * rs - 8.0 * r + 4.0 * sss - s + 4.0);
        deriv[1][0] = b * sm * (16.0 * rrr + 12.0 * rr - 6.0 * rs - 8.0 * r - 4.0 * sss + s - 4.0);
        deriv[2][0] = b * sp * (16.0 * rrr + 12.0 * rr + 6.0 * rs - 8.0 * r + 4.0 * sss - s - 4.0);
        deriv[3][0] = b * sp * (16.0 * rrr - 12.0 * rr + 6.0 * rs - 8.0 * r - 4.0 * sss + s + 4.0);
        deriv[4][0] = -a * (1.0 - 4.0 * r - 3.0 * rr + 8.0 * rrr) * sm;
        deriv[5][0] = a * s * sm * sp * (-1.0 + 2.0 * s);
        deriv[6][0] = -a * (-1.0 - 4.0 * r + 3.0 * rr + 8.0 * rrr) * sp;
        deriv[7][0] = -a * s * sm * sp * (1.0 + 2.0 * s);
        deriv[8][0] = r * sm * (8.0 * rr + s - 4.0);
        deriv[9][0] = 0.5 * sm * sp * (2.0 * r - 4.0 * ss + 1.0);
        deriv[10][0] = r * sp * (8.0 * rr - s - 4.0);
        deriv[11][0] = 0.5 * sm * sp * (2.0 * r - 1.0 + 4.0 * ss);
        deriv[12][0] = a * (1.0 + 4.0 * r - 3.0 * rr - 8.0 * rrr) * sm;
        deriv[13][0] = a * s * sm * sp * (1.0 + 2.0 * s);
        deriv[14][0] = -a * (1.0 - 4.0 * r - 3.0 * rr + 8.0 * rrr) * sp;
        deriv[15][0] = a * s * sm * sp * (1.0 - 2.0 * s);
        deriv[16][0] = -2.0 * r * sm * sp;

        deriv[0][1] = b * rm * (16.0 * sss - 12.0 * ss - 6.0 * rs - 8.0 * s + 4.0 * rrr - r + 4.0);
        deriv[1][1] = -b * rp * (-16.0 * sss + 12.0 * ss - 6.0 * rs + 8.0 * s + 4.0 * rrr - r - 4.0);
        deriv[2][1] = b * rp * (16.0 * sss + 12.0 * ss + 6.0 * rs - 8.0 * s + 4.0 * rrr - r - 4.0);
        deriv[3][1] = b * r1 * (-16.0 * sss - 12.0 * ss + 6.0 * rs + 8.0 * s + 4.0 * rrr - r + 4.0);
        deriv[4][1] = a * r * r1 * rp * (2.0 * r - 1.0);
        deriv[5][1] = -a * (1.0 - 4.0 * s - 3.0 * ss + 8.0 * sss) * rp;
        deriv[6][1] = -a * r * r1 * rp * (1.0 + 2.0 * r);
        deriv[7][1] = a * (-1.0 - 4.0 * s + 3.0 * ss + 8.0 * sss) * r1;
        deriv[8][1] = -0.5 * r1 * rp * (2.0 * s - 1.0 + 4.0 * rr);
        deriv[9][1] = -s * rp * (-8.0 * ss + r + 4.0);
        deriv[10][1] = 0.5 * r1 * rp * (-2.0 * s + 4.0 * rr - 1.0);
        deriv[11][1] = -s * r1 * (8.0 * ss + r - 4.0);
        deriv[12][1] = a * r * r1 * rp * (1.0 + 2.0 * r);
        deriv[13][1] = -a * (-1.0 - 4.0 * s + 3.0 * ss + 8.0 * sss) * rp;
        deriv[14][1] = -a * r * r1 * rp * (2.0 * r - 1.0);
        deriv[15][1] = a * (1.0 - 4.0 * s - 3.0 * ss + 8.0 * sss) * r1;
        deriv[16][1] = 2.0 * s * r1 * rp;
    }
}
