use russell_lab::{Matrix, Vector};

/// Defines a hexahedron with 20 points (quadratic faces)
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of points
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
/// # Local IDs of faces
///
/// Note: the order of points is such that the right-hand rule generates outward normals.
/// Also, the order of face points corresponds to Qua8 points.
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
pub struct Hex20 {}

impl Hex20 {
    pub const NDIM: usize = 3;
    pub const NPOINT: usize = 20;
    pub const NEDGE: usize = 12;
    pub const NFACE: usize = 6;
    pub const EDGE_NPOINT: usize = 3;
    pub const FACE_NPOINT: usize = 8;
    pub const FACE_NEDGE: usize = 4;

    #[rustfmt::skip]
    pub const EDGE_POINT_IDS: [[usize; Hex20::EDGE_NPOINT]; Hex20::NEDGE] = [
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
    pub const FACE_POINT_IDS: [[usize; Hex20::FACE_NPOINT]; Hex20::NFACE] = [
        [0, 4, 7, 3, 16, 15, 19, 11],
        [1, 2, 6, 5,  9, 18, 13, 17],
        [0, 1, 5, 4,  8, 17, 12, 16],
        [2, 3, 7, 6, 10, 19, 14, 18],
        [0, 3, 2, 1, 11, 10,  9,  8],
        [4, 5, 6, 7, 12, 13, 14, 15],
    ];

    #[rustfmt::skip]
    pub const FACE_EDGE_POINT_IDS: [[[usize; Hex20::EDGE_NPOINT]; Hex20::FACE_NEDGE]; Hex20::NFACE] = [
        [[0, 4, 16], [4, 7, 15], [7, 3, 19], [3, 0, 11]],
        [[1, 2,  9], [2, 6, 18], [6, 5, 13], [5, 1, 17]],
        [[0, 1,  8], [1, 5, 17], [5, 4, 12], [4, 0, 16]],
        [[2, 3, 10], [3, 7, 19], [7, 6, 14], [6, 2, 18]],
        [[0, 3, 11], [3, 2, 10], [2, 1,  9], [1, 0,  8]],
        [[4, 5, 12], [5, 6, 13], [6, 7, 14], [7, 4, 15]],
    ];

    #[rustfmt::skip]
    pub const POINT_REFERENCE_COORDS: [[f64; Hex20::NDIM]; Hex20::NPOINT] = [
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
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        let rp1 = 1.0 + r;
        let rm1 = 1.0 - r;
        let sp1 = 1.0 + s;
        let sm1 = 1.0 - s;
        let tp1 = 1.0 + t;
        let tm1 = 1.0 - t;

        interp[0] = rm1 * sm1 * tm1 * (-r - s - t - 2.0) / 8.0;
        interp[1] = rp1 * sm1 * tm1 * (r - s - t - 2.0) / 8.0;
        interp[2] = rp1 * sp1 * tm1 * (r + s - t - 2.0) / 8.0;
        interp[3] = rm1 * sp1 * tm1 * (-r + s - t - 2.0) / 8.0;
        interp[4] = rm1 * sm1 * tp1 * (-r - s + t - 2.0) / 8.0;
        interp[5] = rp1 * sm1 * tp1 * (r - s + t - 2.0) / 8.0;
        interp[6] = rp1 * sp1 * tp1 * (r + s + t - 2.0) / 8.0;
        interp[7] = rm1 * sp1 * tp1 * (-r + s + t - 2.0) / 8.0;
        interp[8] = (1.0 - r * r) * sm1 * tm1 / 4.0;
        interp[9] = rp1 * (1.0 - s * s) * tm1 / 4.0;
        interp[10] = (1.0 - r * r) * sp1 * tm1 / 4.0;
        interp[11] = rm1 * (1.0 - s * s) * tm1 / 4.0;
        interp[12] = (1.0 - r * r) * sm1 * tp1 / 4.0;
        interp[13] = rp1 * (1.0 - s * s) * tp1 / 4.0;
        interp[14] = (1.0 - r * r) * sp1 * tp1 / 4.0;
        interp[15] = rm1 * (1.0 - s * s) * tp1 / 4.0;
        interp[16] = rm1 * sm1 * (1.0 - t * t) / 4.0;
        interp[17] = rp1 * sm1 * (1.0 - t * t) / 4.0;
        interp[18] = rp1 * sp1 * (1.0 - t * t) / 4.0;
        interp[19] = rm1 * sp1 * (1.0 - t * t) / 4.0;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &[f64]) {
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        let rp1 = 1.0 + r;
        let rm1 = 1.0 - r;
        let sp1 = 1.0 + s;
        let sm1 = 1.0 - s;
        let tp1 = 1.0 + t;
        let tm1 = 1.0 - t;

        deriv[0][0] = -0.125 * sm1 * tm1 * (-r - s - t - 2.0) - 0.125 * rm1 * sm1 * tm1;
        deriv[1][0] = 0.125 * sm1 * tm1 * (r - s - t - 2.0) + 0.125 * rp1 * sm1 * tm1;
        deriv[2][0] = 0.125 * sp1 * tm1 * (r + s - t - 2.0) + 0.125 * rp1 * sp1 * tm1;
        deriv[3][0] = -0.125 * sp1 * tm1 * (-r + s - t - 2.0) - 0.125 * rm1 * sp1 * tm1;
        deriv[4][0] = -0.125 * sm1 * tp1 * (-r - s + t - 2.0) - 0.125 * rm1 * sm1 * tp1;
        deriv[5][0] = 0.125 * sm1 * tp1 * (r - s + t - 2.0) + 0.125 * rp1 * sm1 * tp1;
        deriv[6][0] = 0.125 * sp1 * tp1 * (r + s + t - 2.0) + 0.125 * rp1 * sp1 * tp1;
        deriv[7][0] = -0.125 * sp1 * tp1 * (-r + s + t - 2.0) - 0.125 * rm1 * sp1 * tp1;
        deriv[8][0] = -0.5 * r * sm1 * tm1;
        deriv[9][0] = 0.25 * (1.0 - s * s) * tm1;
        deriv[10][0] = -0.5 * r * sp1 * tm1;
        deriv[11][0] = -0.25 * (1.0 - s * s) * tm1;
        deriv[12][0] = -0.5 * r * sm1 * tp1;
        deriv[13][0] = 0.25 * (1.0 - s * s) * tp1;
        deriv[14][0] = -0.5 * r * sp1 * tp1;
        deriv[15][0] = -0.25 * (1.0 - s * s) * tp1;
        deriv[16][0] = -0.25 * sm1 * (1.0 - t * t);
        deriv[17][0] = 0.25 * sm1 * (1.0 - t * t);
        deriv[18][0] = 0.25 * sp1 * (1.0 - t * t);
        deriv[19][0] = -0.25 * sp1 * (1.0 - t * t);

        deriv[0][1] = -0.125 * rm1 * tm1 * (-r - s - t - 2.0) - 0.125 * rm1 * sm1 * tm1;
        deriv[1][1] = -0.125 * rp1 * tm1 * (r - s - t - 2.0) - 0.125 * rp1 * sm1 * tm1;
        deriv[2][1] = 0.125 * rp1 * tm1 * (r + s - t - 2.0) + 0.125 * rp1 * sp1 * tm1;
        deriv[3][1] = 0.125 * rm1 * tm1 * (-r + s - t - 2.0) + 0.125 * rm1 * sp1 * tm1;
        deriv[4][1] = -0.125 * rm1 * tp1 * (-r - s + t - 2.0) - 0.125 * rm1 * sm1 * tp1;
        deriv[5][1] = -0.125 * rp1 * tp1 * (r - s + t - 2.0) - 0.125 * rp1 * sm1 * tp1;
        deriv[6][1] = 0.125 * rp1 * tp1 * (r + s + t - 2.0) + 0.125 * rp1 * sp1 * tp1;
        deriv[7][1] = 0.125 * rm1 * tp1 * (-r + s + t - 2.0) + 0.125 * rm1 * sp1 * tp1;
        deriv[8][1] = -0.25 * (1.0 - r * r) * tm1;
        deriv[9][1] = -0.5 * s * rp1 * tm1;
        deriv[10][1] = 0.25 * (1.0 - r * r) * tm1;
        deriv[11][1] = -0.5 * s * rm1 * tm1;
        deriv[12][1] = -0.25 * (1.0 - r * r) * tp1;
        deriv[13][1] = -0.5 * s * rp1 * tp1;
        deriv[14][1] = 0.25 * (1.0 - r * r) * tp1;
        deriv[15][1] = -0.5 * s * rm1 * tp1;
        deriv[16][1] = -0.25 * rm1 * (1.0 - t * t);
        deriv[17][1] = -0.25 * rp1 * (1.0 - t * t);
        deriv[18][1] = 0.25 * rp1 * (1.0 - t * t);
        deriv[19][1] = 0.25 * rm1 * (1.0 - t * t);

        deriv[0][2] = -0.125 * rm1 * sm1 * (-r - s - t - 2.0) - 0.125 * rm1 * sm1 * tm1;
        deriv[1][2] = -0.125 * rp1 * sm1 * (r - s - t - 2.0) - 0.125 * rp1 * sm1 * tm1;
        deriv[2][2] = -0.125 * rp1 * sp1 * (r + s - t - 2.0) - 0.125 * rp1 * sp1 * tm1;
        deriv[3][2] = -0.125 * rm1 * sp1 * (-r + s - t - 2.0) - 0.125 * rm1 * sp1 * tm1;
        deriv[4][2] = 0.125 * rm1 * sm1 * (-r - s + t - 2.0) + 0.125 * rm1 * sm1 * tp1;
        deriv[5][2] = 0.125 * rp1 * sm1 * (r - s + t - 2.0) + 0.125 * rp1 * sm1 * tp1;
        deriv[6][2] = 0.125 * rp1 * sp1 * (r + s + t - 2.0) + 0.125 * rp1 * sp1 * tp1;
        deriv[7][2] = 0.125 * rm1 * sp1 * (-r + s + t - 2.0) + 0.125 * rm1 * sp1 * tp1;
        deriv[8][2] = -0.25 * (1.0 - r * r) * sm1;
        deriv[9][2] = -0.25 * rp1 * (1.0 - s * s);
        deriv[10][2] = -0.25 * (1.0 - r * r) * sp1;
        deriv[11][2] = -0.25 * rm1 * (1.0 - s * s);
        deriv[12][2] = 0.25 * (1.0 - r * r) * sm1;
        deriv[13][2] = 0.25 * rp1 * (1.0 - s * s);
        deriv[14][2] = 0.25 * (1.0 - r * r) * sp1;
        deriv[15][2] = 0.25 * rm1 * (1.0 - s * s);
        deriv[16][2] = -0.5 * t * rm1 * sm1;
        deriv[17][2] = -0.5 * t * rp1 * sm1;
        deriv[18][2] = -0.5 * t * rp1 * sp1;
        deriv[19][2] = -0.5 * t * rm1 * sp1;
    }
}
