use russell_lab::{Matrix, Vector};

/// Defines a hexahedron with 20 nodes (quadratic faces)
///
/// ![hex20](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_hex20.svg)
///
/// # Local IDs of nodes
///
/// ```text
///            4_______15_______7            r     s     t             r     s     t
///          ,'|              ,'|    p:0 [-1.0, -1.0, -1.0]   p:8  [ 0.0, -1.0, -1.0]
///       12'  |            ,'  |    p:1 [ 1.0, -1.0, -1.0]   p:9  [ 1.0,  0.0, -1.0]
///      ,'    16         ,14   |    p:2 [ 1.0,  1.0, -1.0]   p:10 [ 0.0,  1.0, -1.0]
///    ,'      |        ,'      19   p:3 [-1.0,  1.0, -1.0]   p:11 [-1.0,  0.0, -1.0]
///  5'=====13========6'        |                             
///  |         |      |         |    p:4 [-1.0, -1.0,  1.0]   p:12 [ 0.0, -1.0,  1.0]
///  |         |      |         |    p:5 [ 1.0, -1.0,  1.0]   p:13 [ 1.0,  0.0,  1.0]
///  |         0_____ | _11_____3    p:6 [ 1.0,  1.0,  1.0]   p:14 [ 0.0,  1.0,  1.0]
/// 17       ,'       |       ,'     p:7 [-1.0,  1.0,  1.0]   p:15 [-1.0,  0.0,  1.0]
///  |     8'        18     ,'       
///  |   ,'           |   ,10       p:16 [-1.0, -1.0,  0.0]
///  | ,'             | ,'          p:17 [ 1.0, -1.0,  0.0]
///  1_______9________2'            p:18 [ 1.0,  1.0,  0.0]
///                                 p:19 [-1.0,  1.0,  0.0]
/// ```
///
/// # Local IDs of edges
///
/// ```text
///                     7
///             +----------------+          p0  p1 p2
///           ,'|              ,'|     e:0  [0, 1,  8],
///       4 ,'  |8          6,'  |     e:1  [1, 2,  9],
///       ,'    |          ,'    |     e:2  [2, 3, 10],
///     ,'      |   5    ,'      |11   e:3  [3, 0, 11],
///   +'===============+'        |     e:4  [4, 5, 12],
///   |         |      |         |     e:5  [5, 6, 13],
///   |         |      |  3      |     e:6  [6, 7, 14],
///   |         .- - - | -  - - -+     e:7  [7, 4, 15],
///  9|       ,'       |       ,'      e:8  [0, 4, 16],
///   |    0,'         |10   ,'        e:9  [1, 5, 17],
///   |   ,'           |   ,' 2        e:10 [2, 6, 18],
///   | ,'             | ,'            e:11 [3, 7, 19],
///   +----------------+'
///           1
/// ```
///
/// * The order of edge nodes corresponds to **Lin3** nodes.
///
/// # Local IDs of faces
///
/// ```text
///           4-------15-------7
///         ,'|              ,'|
///      12' 16  ___       14  |
///     ,'    |,'5,'  [0],'    |
///   ,'      |~~~     ,'     19   f:0 [0, 4, 7, 3, 16, 15, 19, 11]
/// 5'======13=======6'  ,'|   |   f:1 [1, 2, 6, 5,  9, 18, 13, 17]
/// |   ,'|   |      |   |3|   |   f:2 [0, 1, 5, 4,  8, 17, 12, 16]
/// |   |2|   |      |   |,'   |   f:3 [2, 3, 7, 6, 10, 19, 14, 18]
/// |   |,'   0- - - | 11 - - -3   f:4 [0, 3, 2, 1, 11, 10,  9,  8]
/// 17      ,'      18       ,'    f:5 [4, 5, 6, 7, 12, 13, 14, 15]
/// |     8' [1]  ___|     10
/// |   ,'      ,'4,'|   ,'
/// | ,'        ~~~  | ,'
/// 1-------9--------2'
/// ```
///
/// # Notes
///
/// * The reference coordinates range from -1 to +1 with the geometry centred @ 0
/// * The order of face nodes is such that the normals are outward
/// * The order of face nodes corresponds to [super::Qua8] nodes
/// * This shape is a higher-order version of [super::Hex8]
pub struct Hex20 {}

impl Hex20 {
    pub const GEO_NDIM: usize = 3;
    pub const NNODE: usize = 20;
    pub const NEDGE: usize = 12;
    pub const NFACE: usize = 6;
    pub const EDGE_NNODE: usize = 3;
    pub const FACE_NNODE: usize = 8;
    pub const FACE_NEDGE: usize = 4;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Hex20::EDGE_NNODE]; Hex20::NEDGE] = [
        [0, 1,  8],
        [1, 2,  9],
        [2, 3, 10],
        [3, 0, 11],
        [4, 5, 12],
        [5, 6, 13],
        [6, 7, 14],
        [7, 4, 15],
        [0, 4, 16],
        [1, 5, 17],
        [2, 6, 18],
        [3, 7, 19],
    ];

    #[rustfmt::skip]
    pub const FACE_NODE_IDS: [[usize; Hex20::FACE_NNODE]; Hex20::NFACE] = [
        [0, 4, 7, 3, 16, 15, 19, 11],
        [1, 2, 6, 5,  9, 18, 13, 17],
        [0, 1, 5, 4,  8, 17, 12, 16],
        [2, 3, 7, 6, 10, 19, 14, 18],
        [0, 3, 2, 1, 11, 10,  9,  8],
        [4, 5, 6, 7, 12, 13, 14, 15],
    ];

    #[rustfmt::skip]
    pub const FACE_EDGE_NODE_IDS: [[[usize; Hex20::EDGE_NNODE]; Hex20::FACE_NEDGE]; Hex20::NFACE] = [
        [[0, 4, 16], [4, 7, 15], [7, 3, 19], [3, 0, 11]],
        [[1, 2,  9], [2, 6, 18], [6, 5, 13], [5, 1, 17]],
        [[0, 1,  8], [1, 5, 17], [5, 4, 12], [4, 0, 16]],
        [[2, 3, 10], [3, 7, 19], [7, 6, 14], [6, 2, 18]],
        [[0, 3, 11], [3, 2, 10], [2, 1,  9], [1, 0,  8]],
        [[4, 5, 12], [5, 6, 13], [6, 7, 14], [7, 4, 15]],
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Hex20::GEO_NDIM]; Hex20::NNODE] = [
        [-1.0, -1.0, -1.0],
        [ 1.0, -1.0, -1.0],
        [ 1.0,  1.0, -1.0],
        [-1.0,  1.0, -1.0],

        [-1.0, -1.0,  1.0],
        [ 1.0, -1.0,  1.0],
        [ 1.0,  1.0,  1.0],
        [-1.0,  1.0,  1.0],

        [ 0.0, -1.0, -1.0],
        [ 1.0,  0.0, -1.0],
        [ 0.0,  1.0, -1.0],
        [-1.0,  0.0, -1.0],

        [ 0.0, -1.0,  1.0],
        [ 1.0,  0.0,  1.0],
        [ 0.0,  1.0,  1.0],
        [-1.0,  0.0,  1.0],

        [-1.0, -1.0,  0.0],
        [ 1.0, -1.0,  0.0],
        [ 1.0,  1.0,  0.0],
        [-1.0,  1.0,  0.0],
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
        let rr = 1.0 - r * r;
        let ss = 1.0 - s * s;
        let tt = 1.0 - t * t;

        interp[0] = rm * sm * tm * (-2.0 - r - s - t) / 8.0;
        interp[1] = rp * sm * tm * (-2.0 + r - s - t) / 8.0;
        interp[2] = rp * sp * tm * (-2.0 + r + s - t) / 8.0;
        interp[3] = rm * sp * tm * (-2.0 - r + s - t) / 8.0;
        interp[4] = rm * sm * tp * (-2.0 - r - s + t) / 8.0;
        interp[5] = rp * sm * tp * (-2.0 + r - s + t) / 8.0;
        interp[6] = rp * sp * tp * (-2.0 + r + s + t) / 8.0;
        interp[7] = rm * sp * tp * (-2.0 - r + s + t) / 8.0;
        interp[8] = rr * sm * tm / 4.0;
        interp[9] = rp * ss * tm / 4.0;
        interp[10] = rr * sp * tm / 4.0;
        interp[11] = rm * ss * tm / 4.0;
        interp[12] = rr * sm * tp / 4.0;
        interp[13] = rp * ss * tp / 4.0;
        interp[14] = rr * sp * tp / 4.0;
        interp[15] = rm * ss * tp / 4.0;
        interp[16] = rm * sm * tt / 4.0;
        interp[17] = rp * sm * tt / 4.0;
        interp[18] = rp * sp * tt / 4.0;
        interp[19] = rm * sp * tt / 4.0;
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
        let rr = 1.0 - r * r;
        let ss = 1.0 - s * s;
        let tt = 1.0 - t * t;

        // with respect to r
        deriv[0][0] = -rm * sm * tm / 8.0 - sm * (-2.0 - r - s - t) * tm / 8.0;
        deriv[1][0] = rp * sm * tm / 8.0 + sm * (-2.0 + r - s - t) * tm / 8.0;
        deriv[2][0] = rp * sp * tm / 8.0 + sp * (-2.0 + r + s - t) * tm / 8.0;
        deriv[3][0] = -rm * sp * tm / 8.0 - sp * (-2.0 - r + s - t) * tm / 8.0;
        deriv[4][0] = -rm * sm * tp / 8.0 - sm * (-2.0 - r - s + t) * tp / 8.0;
        deriv[5][0] = rp * sm * tp / 8.0 + sm * (-2.0 + r - s + t) * tp / 8.0;
        deriv[6][0] = rp * sp * tp / 8.0 + sp * (-2.0 + r + s + t) * tp / 8.0;
        deriv[7][0] = -rm * sp * tp / 8.0 - sp * (-2.0 - r + s + t) * tp / 8.0;
        deriv[8][0] = -r * sm * tm / 2.0;
        deriv[9][0] = ss * tm / 4.0;
        deriv[10][0] = -r * sp * tm / 2.0;
        deriv[11][0] = -ss * tm / 4.0;
        deriv[12][0] = -r * sm * tp / 2.0;
        deriv[13][0] = ss * tp / 4.0;
        deriv[14][0] = -r * sp * tp / 2.0;
        deriv[15][0] = -ss * tp / 4.0;
        deriv[16][0] = -sm * tt / 4.0;
        deriv[17][0] = sm * tt / 4.0;
        deriv[18][0] = sp * tt / 4.0;
        deriv[19][0] = -sp * tt / 4.0;

        // with respect to s
        deriv[0][1] = -rm * sm * tm / 8.0 - rm * (-2.0 - r - s - t) * tm / 8.0;
        deriv[1][1] = -rp * sm * tm / 8.0 - rp * (-2.0 + r - s - t) * tm / 8.0;
        deriv[2][1] = rp * sp * tm / 8.0 + rp * (-2.0 + r + s - t) * tm / 8.0;
        deriv[3][1] = rm * sp * tm / 8.0 + rm * (-2.0 - r + s - t) * tm / 8.0;
        deriv[4][1] = -rm * sm * tp / 8.0 - rm * (-2.0 - r - s + t) * tp / 8.0;
        deriv[5][1] = -rp * sm * tp / 8.0 - rp * (-2.0 + r - s + t) * tp / 8.0;
        deriv[6][1] = rp * sp * tp / 8.0 + rp * (-2.0 + r + s + t) * tp / 8.0;
        deriv[7][1] = rm * sp * tp / 8.0 + rm * (-2.0 - r + s + t) * tp / 8.0;
        deriv[8][1] = -rr * tm / 4.0;
        deriv[9][1] = -rp * s * tm / 2.0;
        deriv[10][1] = rr * tm / 4.0;
        deriv[11][1] = -rm * s * tm / 2.0;
        deriv[12][1] = -rr * tp / 4.0;
        deriv[13][1] = -rp * s * tp / 2.0;
        deriv[14][1] = rr * tp / 4.0;
        deriv[15][1] = -rm * s * tp / 2.0;
        deriv[16][1] = -rm * tt / 4.0;
        deriv[17][1] = -rp * tt / 4.0;
        deriv[18][1] = rp * tt / 4.0;
        deriv[19][1] = rm * tt / 4.0;

        // with respect to s
        deriv[0][2] = -rm * sm * (-2.0 - r - s - t) / 8.0 - rm * sm * tm / 8.0;
        deriv[1][2] = -rp * sm * (-2.0 + r - s - t) / 8.0 - rp * sm * tm / 8.0;
        deriv[2][2] = -rp * sp * (-2.0 + r + s - t) / 8.0 - rp * sp * tm / 8.0;
        deriv[3][2] = -rm * sp * (-2.0 - r + s - t) / 8.0 - rm * sp * tm / 8.0;
        deriv[4][2] = rm * sm * (-2.0 - r - s + t) / 8.0 + rm * sm * tp / 8.0;
        deriv[5][2] = rp * sm * (-2.0 + r - s + t) / 8.0 + rp * sm * tp / 8.0;
        deriv[6][2] = rp * sp * (-2.0 + r + s + t) / 8.0 + rp * sp * tp / 8.0;
        deriv[7][2] = rm * sp * (-2.0 - r + s + t) / 8.0 + rm * sp * tp / 8.0;
        deriv[8][2] = -rr * sm / 4.0;
        deriv[9][2] = -rp * ss / 4.0;
        deriv[10][2] = -rr * sp / 4.0;
        deriv[11][2] = -rm * ss / 4.0;
        deriv[12][2] = rr * sm / 4.0;
        deriv[13][2] = rp * ss / 4.0;
        deriv[14][2] = rr * sp / 4.0;
        deriv[15][2] = rm * ss / 4.0;
        deriv[16][2] = -rm * sm * t / 2.0;
        deriv[17][2] = -rp * sm * t / 2.0;
        deriv[18][2] = -rp * sp * t / 2.0;
        deriv[19][2] = -rm * sp * t / 2.0;
    }
}
