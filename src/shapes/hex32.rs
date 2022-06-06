use russell_lab::{Matrix, Vector};

/// Defines a hexahedron with 32 nodes (cubic faces)
///
/// The reference coordinates (r,s,t) range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of nodes
///
/// ```text
///                               node     r     s     t   node     r     s    t
///                                  0 [-1.0, -1.0, -1.0]     4 [-1.0, -1.0, 1.0]
///                                  1 [ 1.0, -1.0, -1.0]     5 [ 1.0, -1.0, 1.0]
///                                  2 [ 1.0,  1.0, -1.0]     6 [ 1.0,  1.0, 1.0]
///                                  3 [-1.0,  1.0, -1.0]     7 [-1.0,  1.0, 1.0]
///            4----23----22----7                                 
///          ,'|              ,'| node                            node
///        16  |           21'  |    8 [-1.0/3.0,     -1.0, -1.0]   16 [-1.0/3.0,     -1.0, 1.0]
///     17'    25         ,'   31    9 [ 1.0/3.0,     -1.0, -1.0]   17 [ 1.0/3.0,     -1.0, 1.0]
///    ,'      |        ,20     |   10 [     1.0, -1.0/3.0, -1.0]   18 [     1.0, -1.0/3.0, 1.0]
///  5-----18----19---6'       30   11 [     1.0,  1.0/3.0, -1.0]   19 [     1.0,  1.0/3.0, 1.0]
///  |         24     |         |   12 [ 1.0/3.0,      1.0, -1.0]   20 [ 1.0/3.0,      1.0, 1.0]
///  |         |     29         |   13 [-1.0/3.0,      1.0, -1.0]   21 [-1.0/3.0,      1.0, 1.0]
/// 27         0--15- | --14----3   14 [    -1.0,  1.0/3.0, -1.0]   22 [    -1.0,  1.0/3.0, 1.0]
///  |       ,'       |       ,'    15 [    -1.0, -1.0/3.0, -1.0]   23 [    -1.0, -1.0/3.0, 1.0]
/// 26     ,8        28     ,13     
///  |   9'           |   12      node     r     s     t
///  | ,'             | ,'          24 [-1.0, -1.0, -1.0/3.0]
///  1----10----11----2'            25 [-1.0, -1.0,  1.0/3.0]
///                                 26 [ 1.0, -1.0, -1.0/3.0]
///            t^                   27 [ 1.0, -1.0,  1.0/3.0]
///             |                   28 [ 1.0,  1.0, -1.0/3.0]
///             o--->s              29 [ 1.0,  1.0,  1.0/3.0]
///           ,'                    30 [-1.0,  1.0, -1.0/3.0]
///          r                      31 [-1.0,  1.0,  1.0/3.0]
/// ```
///
/// # Local IDs of edges
///
/// ```text
///                     7
///             +----------------+    edge  points...
///           ,'|              ,'|       0 [0, 1,  8,  9]
///       4 ,'  |8          6,'  |       1 [1, 2, 10, 11]
///       ,'    |          ,'    |       2 [2, 3, 12, 13]
///     ,'      |   5    ,'      |11     3 [3, 0, 14, 15]
///   +'===============+'        |       4 [4, 5, 16, 17]
///   |         |      |         |       5 [5, 6, 18, 19]
///   |         |      |  3      |       6 [6, 7, 20, 21]
///   |         .- - - | -  - - -+       7 [7, 4, 22, 23]
///  9|       ,'       |       ,'        8 [0, 4, 24, 25]
///   |    0,'         |10   ,'          9 [1, 5, 26, 27]
///   |   ,'           |   ,' 2         10 [2, 6, 28, 29]
///   | ,'             | ,'             11 [3, 7, 30, 31]
///   +----------------+'
///           1
/// ```
///
/// * The order of edge nodes corresponds to **Lin4** nodes.
///
/// # Local IDs of faces
///
/// ```text
///           4----23----22----7
///         ,'|              ,'|
///      16' 25  ___       21  |
///    17'    |,'5,'  [0],'    |  face  points...
///   ,'      |~~~     20     31     0 [0, 4, 7, 3, 24, 23, 31, 14, 25, 22, 30, 15]
/// 5'======18===19==6'  ,'|   |     1 [1, 2, 6, 5, 10, 28, 19, 27, 11, 29, 18, 26]
/// |   ,'|   24     |   |3|  30     2 [0, 1, 5, 4,  8, 26, 17, 25,  9, 27, 16, 24]
/// |   |2|   |     29   |,'   |     3 [2, 3, 7, 6, 12, 30, 21, 29, 13, 31, 20, 28]
/// 27  |,'   0- 15- | -- -14 -3     4 [0, 3, 2, 1, 15, 13, 11,  9, 14, 12, 10,  8]
/// |       ,'      28       ,'      5 [4, 5, 6, 7, 16, 18, 20, 22, 17, 19, 21, 23]
/// |     8' [1]  ___|     ,13
/// 26  9'      ,'4,'|   12
/// | ,'        ~~~  | ,'
/// 1----10----11----2'
/// ```
///
/// * The order of face nodes is such that the normals are outward
/// * The order of face nodes corresponds to **Qua12** nodes
pub struct Hex32 {}

impl Hex32 {
    pub const GEO_NDIM: usize = 3;
    pub const NNODE: usize = 32;
    pub const NEDGE: usize = 12;
    pub const NFACE: usize = 6;
    pub const EDGE_NNODE: usize = 4;
    pub const FACE_NNODE: usize = 12;
    pub const FACE_NEDGE: usize = 4;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Hex32::EDGE_NNODE]; Hex32::NEDGE] = [
        [0, 1,  8,  9],
        [1, 2, 10, 11],
        [2, 3, 12, 13],
        [3, 0, 14, 15],
        [4, 5, 16, 17],
        [5, 6, 18, 19],
        [6, 7, 20, 21],
        [7, 4, 22, 23],
        [0, 4, 24, 25],
        [1, 5, 26, 27],
        [2, 6, 28, 29],
        [3, 7, 30, 31],
    ];

    #[rustfmt::skip]
    pub const FACE_NODE_IDS: [[usize; Hex32::FACE_NNODE]; Hex32::NFACE] = [
        [0, 4, 7, 3, 24, 23, 31, 14, 25, 22, 30, 15],
        [1, 2, 6, 5, 10, 28, 19, 27, 11, 29, 18, 26],
        [0, 1, 5, 4,  8, 26, 17, 25,  9, 27, 16, 24],
        [2, 3, 7, 6, 12, 30, 21, 29, 13, 31, 20, 28],
        [0, 3, 2, 1, 15, 13, 11,  9, 14, 12, 10,  8],
        [4, 5, 6, 7, 16, 18, 20, 22, 17, 19, 21, 23],
    ];

    #[rustfmt::skip]
    pub const FACE_EDGE_NODE_IDS: [[[usize; Hex32::EDGE_NNODE]; Hex32::FACE_NEDGE]; Hex32::NFACE] = [
        [[0, 4, 24, 25], [4, 7, 23, 22], [7, 3, 31, 30], [3, 0, 14, 15]],
        [[1, 2, 10, 11], [2, 6, 28, 29], [6, 5, 19, 18], [5, 1, 27, 26]],
        [[0, 1,  8,  9], [1, 5, 26, 27], [5, 4, 17, 16], [4, 0, 25, 24]],
        [[2, 3, 12, 13], [3, 7, 30, 31], [7, 6, 21, 20], [6, 2, 29, 28]],
        [[0, 3, 15, 14], [3, 2, 13, 12], [2, 1, 11, 10], [1, 0,  9,  8]],
        [[4, 5, 16, 17], [5, 6, 18, 19], [6, 7, 20, 21], [7, 4, 22, 23]],
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Hex32::GEO_NDIM]; Hex32::NNODE] = [
        [-1.0, -1.0, -1.0], // 0
        [ 1.0, -1.0, -1.0], // 1
        [ 1.0,  1.0, -1.0], // 2
        [-1.0,  1.0, -1.0], // 3

        [-1.0, -1.0, 1.0], // 4
        [ 1.0, -1.0, 1.0], // 5
        [ 1.0,  1.0, 1.0], // 6
        [-1.0,  1.0, 1.0], // 7

        [-1.0/3.0,     -1.0, -1.0], //  8
        [ 1.0/3.0,     -1.0, -1.0], //  9
        [     1.0, -1.0/3.0, -1.0], // 10
        [     1.0,  1.0/3.0, -1.0], // 11
        [ 1.0/3.0,      1.0, -1.0], // 12
        [-1.0/3.0,      1.0, -1.0], // 13
        [    -1.0,  1.0/3.0, -1.0], // 14
        [    -1.0, -1.0/3.0, -1.0], // 15

        [-1.0/3.0,     -1.0, 1.0], // 16
        [ 1.0/3.0,     -1.0, 1.0], // 17
        [     1.0, -1.0/3.0, 1.0], // 18
        [     1.0,  1.0/3.0, 1.0], // 19
        [ 1.0/3.0,      1.0, 1.0], // 20
        [-1.0/3.0,      1.0, 1.0], // 21
        [    -1.0,  1.0/3.0, 1.0], // 22
        [    -1.0, -1.0/3.0, 1.0], // 23

        [-1.0, -1.0, -1.0/3.0], // 24
        [-1.0, -1.0,  1.0/3.0], // 25
        [ 1.0, -1.0, -1.0/3.0], // 26
        [ 1.0, -1.0,  1.0/3.0], // 27
        [ 1.0,  1.0, -1.0/3.0], // 28
        [ 1.0,  1.0,  1.0/3.0], // 29
        [-1.0,  1.0, -1.0/3.0], // 30
        [-1.0,  1.0,  1.0/3.0], // 31
    ];

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
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        let rm = 1.0 - r;
        let sm = 1.0 - s;
        let tm = 1.0 - t;
        let rp = 1.0 + r;
        let sp = 1.0 + s;
        let tp = 1.0 + t;
        let r3m = 1.0 - 3.0 * r;
        let s3m = 1.0 - 3.0 * s;
        let t3m = 1.0 - 3.0 * t;
        let r3p = 1.0 + 3.0 * r;
        let s3p = 1.0 + 3.0 * s;
        let t3p = 1.0 + 3.0 * t;
        let rr = 1.0 - r * r;
        let ss = 1.0 - s * s;
        let tt = 1.0 - t * t;
        let rr9 = 9.0 * r * r;
        let ss9 = 9.0 * s * s;
        let tt9 = 9.0 * t * t;

        interp[0] = (rm * sm * tm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        interp[1] = (rp * sm * tm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        interp[2] = (rp * sp * tm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        interp[3] = (rm * sp * tm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        interp[4] = (rm * sm * tp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        interp[5] = (rp * sm * tp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        interp[6] = (rp * sp * tp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        interp[7] = (rm * sp * tp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        interp[8] = (9.0 * r3m * rr * sm * tm) / 64.0;
        interp[9] = (9.0 * r3p * rr * sm * tm) / 64.0;
        interp[10] = (9.0 * rp * s3m * ss * tm) / 64.0;
        interp[11] = (9.0 * rp * s3p * ss * tm) / 64.0;
        interp[12] = (9.0 * r3p * rr * sp * tm) / 64.0;
        interp[13] = (9.0 * r3m * rr * sp * tm) / 64.0;
        interp[14] = (9.0 * rm * s3p * ss * tm) / 64.0;
        interp[15] = (9.0 * rm * s3m * ss * tm) / 64.0;
        interp[16] = (9.0 * r3m * rr * sm * tp) / 64.0;
        interp[17] = (9.0 * r3p * rr * sm * tp) / 64.0;
        interp[18] = (9.0 * rp * s3m * ss * tp) / 64.0;
        interp[19] = (9.0 * rp * s3p * ss * tp) / 64.0;
        interp[20] = (9.0 * r3p * rr * sp * tp) / 64.0;
        interp[21] = (9.0 * r3m * rr * sp * tp) / 64.0;
        interp[22] = (9.0 * rm * s3p * ss * tp) / 64.0;
        interp[23] = (9.0 * rm * s3m * ss * tp) / 64.0;
        interp[24] = (9.0 * rm * sm * t3m * tt) / 64.0;
        interp[25] = (9.0 * rm * sm * t3p * tt) / 64.0;
        interp[26] = (9.0 * rp * sm * t3m * tt) / 64.0;
        interp[27] = (9.0 * rp * sm * t3p * tt) / 64.0;
        interp[28] = (9.0 * rp * sp * t3m * tt) / 64.0;
        interp[29] = (9.0 * rp * sp * t3p * tt) / 64.0;
        interp[30] = (9.0 * rm * sp * t3m * tt) / 64.0;
        interp[31] = (9.0 * rm * sp * t3p * tt) / 64.0;
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
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        let rm = 1.0 - r;
        let sm = 1.0 - s;
        let tm = 1.0 - t;
        let rp = 1.0 + r;
        let sp = 1.0 + s;
        let tp = 1.0 + t;
        let r3m = 1.0 - 3.0 * r;
        let s3m = 1.0 - 3.0 * s;
        let t3m = 1.0 - 3.0 * t;
        let r3p = 1.0 + 3.0 * r;
        let s3p = 1.0 + 3.0 * s;
        let t3p = 1.0 + 3.0 * t;
        let rr = 1.0 - r * r;
        let ss = 1.0 - s * s;
        let tt = 1.0 - t * t;
        let rr9 = 9.0 * r * r;
        let ss9 = 9.0 * s * s;
        let tt9 = 9.0 * t * t;

        // with respect to r
        deriv[0][0] = (9.0 * r * rm * sm * tm) / 32.0 - (sm * tm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[1][0] = (9.0 * r * rp * sm * tm) / 32.0 + (sm * tm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[2][0] = (9.0 * r * rp * sp * tm) / 32.0 + (sp * tm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[3][0] = (9.0 * r * rm * sp * tm) / 32.0 - (sp * tm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[4][0] = (9.0 * r * rm * sm * tp) / 32.0 - (sm * tp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[5][0] = (9.0 * r * rp * sm * tp) / 32.0 + (sm * tp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[6][0] = (9.0 * r * rp * sp * tp) / 32.0 + (sp * tp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[7][0] = (9.0 * r * rm * sp * tp) / 32.0 - (sp * tp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[8][0] = (-9.0 * r * r3m * sm * tm) / 32.0 - (27.0 * rr * sm * tm) / 64.0;
        deriv[9][0] = (-9.0 * r * r3p * sm * tm) / 32.0 + (27.0 * rr * sm * tm) / 64.0;
        deriv[10][0] = (9.0 * s3m * ss * tm) / 64.0;
        deriv[11][0] = (9.0 * s3p * ss * tm) / 64.0;
        deriv[12][0] = (-9.0 * r * r3p * sp * tm) / 32.0 + (27.0 * rr * sp * tm) / 64.0;
        deriv[13][0] = (-9.0 * r * r3m * sp * tm) / 32.0 - (27.0 * rr * sp * tm) / 64.0;
        deriv[14][0] = (-9.0 * s3p * ss * tm) / 64.0;
        deriv[15][0] = (-9.0 * s3m * ss * tm) / 64.0;
        deriv[16][0] = (-9.0 * r * r3m * sm * tp) / 32.0 - (27.0 * rr * sm * tp) / 64.0;
        deriv[17][0] = (-9.0 * r * r3p * sm * tp) / 32.0 + (27.0 * rr * sm * tp) / 64.0;
        deriv[18][0] = (9.0 * s3m * ss * tp) / 64.0;
        deriv[19][0] = (9.0 * s3p * ss * tp) / 64.0;
        deriv[20][0] = (-9.0 * r * r3p * sp * tp) / 32.0 + (27.0 * rr * sp * tp) / 64.0;
        deriv[21][0] = (-9.0 * r * r3m * sp * tp) / 32.0 - (27.0 * rr * sp * tp) / 64.0;
        deriv[22][0] = (-9.0 * s3p * ss * tp) / 64.0;
        deriv[23][0] = (-9.0 * s3m * ss * tp) / 64.0;
        deriv[24][0] = (-9.0 * sm * t3m * tt) / 64.0;
        deriv[25][0] = (-9.0 * sm * t3p * tt) / 64.0;
        deriv[26][0] = (9.0 * sm * t3m * tt) / 64.0;
        deriv[27][0] = (9.0 * sm * t3p * tt) / 64.0;
        deriv[28][0] = (9.0 * sp * t3m * tt) / 64.0;
        deriv[29][0] = (9.0 * sp * t3p * tt) / 64.0;
        deriv[30][0] = (-9.0 * sp * t3m * tt) / 64.0;
        deriv[31][0] = (-9.0 * sp * t3p * tt) / 64.0;

        // with respect to s
        deriv[0][1] = (9.0 * rm * s * sm * tm) / 32.0 - (rm * tm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[1][1] = (9.0 * rp * s * sm * tm) / 32.0 - (rp * tm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[2][1] = (9.0 * rp * s * sp * tm) / 32.0 + (rp * tm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[3][1] = (9.0 * rm * s * sp * tm) / 32.0 + (rm * tm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[4][1] = (9.0 * rm * s * sm * tp) / 32.0 - (rm * tp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[5][1] = (9.0 * rp * s * sm * tp) / 32.0 - (rp * tp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[6][1] = (9.0 * rp * s * sp * tp) / 32.0 + (rp * tp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[7][1] = (9.0 * rm * s * sp * tp) / 32.0 + (rm * tp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[8][1] = (-9.0 * r3m * rr * tm) / 64.0;
        deriv[9][1] = (-9.0 * r3p * rr * tm) / 64.0;
        deriv[10][1] = (-9.0 * rp * s * s3m * tm) / 32.0 - (27.0 * rp * ss * tm) / 64.0;
        deriv[11][1] = (-9.0 * rp * s * s3p * tm) / 32.0 + (27.0 * rp * ss * tm) / 64.0;
        deriv[12][1] = (9.0 * r3p * rr * tm) / 64.0;
        deriv[13][1] = (9.0 * r3m * rr * tm) / 64.0;
        deriv[14][1] = (-9.0 * rm * s * s3p * tm) / 32.0 + (27.0 * rm * ss * tm) / 64.0;
        deriv[15][1] = (-9.0 * rm * s * s3m * tm) / 32.0 - (27.0 * rm * ss * tm) / 64.0;
        deriv[16][1] = (-9.0 * r3m * rr * tp) / 64.0;
        deriv[17][1] = (-9.0 * r3p * rr * tp) / 64.0;
        deriv[18][1] = (-9.0 * rp * s * s3m * tp) / 32.0 - (27.0 * rp * ss * tp) / 64.0;
        deriv[19][1] = (-9.0 * rp * s * s3p * tp) / 32.0 + (27.0 * rp * ss * tp) / 64.0;
        deriv[20][1] = (9.0 * r3p * rr * tp) / 64.0;
        deriv[21][1] = (9.0 * r3m * rr * tp) / 64.0;
        deriv[22][1] = (-9.0 * rm * s * s3p * tp) / 32.0 + (27.0 * rm * ss * tp) / 64.0;
        deriv[23][1] = (-9.0 * rm * s * s3m * tp) / 32.0 - (27.0 * rm * ss * tp) / 64.0;
        deriv[24][1] = (-9.0 * rm * t3m * tt) / 64.0;
        deriv[25][1] = (-9.0 * rm * t3p * tt) / 64.0;
        deriv[26][1] = (-9.0 * rp * t3m * tt) / 64.0;
        deriv[27][1] = (-9.0 * rp * t3p * tt) / 64.0;
        deriv[28][1] = (9.0 * rp * t3m * tt) / 64.0;
        deriv[29][1] = (9.0 * rp * t3p * tt) / 64.0;
        deriv[30][1] = (9.0 * rm * t3m * tt) / 64.0;
        deriv[31][1] = (9.0 * rm * t3p * tt) / 64.0;

        // with respect to t
        deriv[0][2] = (9.0 * rm * sm * t * tm) / 32.0 - (rm * sm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[1][2] = (9.0 * rp * sm * t * tm) / 32.0 - (rp * sm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[2][2] = (9.0 * rp * sp * t * tm) / 32.0 - (rp * sp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[3][2] = (9.0 * rm * sp * t * tm) / 32.0 - (rm * sp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[4][2] = (9.0 * rm * sm * t * tp) / 32.0 + (rm * sm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[5][2] = (9.0 * rp * sm * t * tp) / 32.0 + (rp * sm * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[6][2] = (9.0 * rp * sp * t * tp) / 32.0 + (rp * sp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[7][2] = (9.0 * rm * sp * t * tp) / 32.0 + (rm * sp * (-19.0 + rr9 + ss9 + tt9)) / 64.0;
        deriv[8][2] = (-9.0 * r3m * rr * sm) / 64.0;
        deriv[9][2] = (-9.0 * r3p * rr * sm) / 64.0;
        deriv[10][2] = (-9.0 * rp * s3m * ss) / 64.0;
        deriv[11][2] = (-9.0 * rp * s3p * ss) / 64.0;
        deriv[12][2] = (-9.0 * r3p * rr * sp) / 64.0;
        deriv[13][2] = (-9.0 * r3m * rr * sp) / 64.0;
        deriv[14][2] = (-9.0 * rm * s3p * ss) / 64.0;
        deriv[15][2] = (-9.0 * rm * s3m * ss) / 64.0;
        deriv[16][2] = (9.0 * r3m * rr * sm) / 64.0;
        deriv[17][2] = (9.0 * r3p * rr * sm) / 64.0;
        deriv[18][2] = (9.0 * rp * s3m * ss) / 64.0;
        deriv[19][2] = (9.0 * rp * s3p * ss) / 64.0;
        deriv[20][2] = (9.0 * r3p * rr * sp) / 64.0;
        deriv[21][2] = (9.0 * r3m * rr * sp) / 64.0;
        deriv[22][2] = (9.0 * rm * s3p * ss) / 64.0;
        deriv[23][2] = (9.0 * rm * s3m * ss) / 64.0;
        deriv[24][2] = (-9.0 * rm * sm * t * t3m) / 32.0 - (27.0 * rm * sm * tt) / 64.0;
        deriv[25][2] = (-9.0 * rm * sm * t * t3p) / 32.0 + (27.0 * rm * sm * tt) / 64.0;
        deriv[26][2] = (-9.0 * rp * sm * t * t3m) / 32.0 - (27.0 * rp * sm * tt) / 64.0;
        deriv[27][2] = (-9.0 * rp * sm * t * t3p) / 32.0 + (27.0 * rp * sm * tt) / 64.0;
        deriv[28][2] = (-9.0 * rp * sp * t * t3m) / 32.0 - (27.0 * rp * sp * tt) / 64.0;
        deriv[29][2] = (-9.0 * rp * sp * t * t3p) / 32.0 + (27.0 * rp * sp * tt) / 64.0;
        deriv[30][2] = (-9.0 * rm * sp * t * t3m) / 32.0 - (27.0 * rm * sp * tt) / 64.0;
        deriv[31][2] = (-9.0 * rm * sp * t * t3p) / 32.0 + (27.0 * rm * sp * tt) / 64.0;
    }
}
