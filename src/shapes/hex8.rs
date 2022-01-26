use russell_lab::{Matrix, Vector};

/// Defines a hexahedron with 8 nodes (bilinear faces)
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0.
///
/// # Local IDs of nodes
///
/// ```text
///           4________________7
///         ,'|              ,'|
///       ,'  |            ,'  |             r     s     t
///     ,'    |          ,'    |     p:0 [-1.0, -1.0, -1.0]
///   ,'      |        ,'      |     p:1 [ 1.0, -1.0, -1.0]
/// 5'===============6'        |     p:2 [ 1.0,  1.0, -1.0]
/// |         |      |         |     p:3 [-1.0,  1.0, -1.0]
/// |         |      |         |     p:4 [-1.0, -1.0,  1.0]
/// |         0_____ | ________3     p:5 [ 1.0, -1.0,  1.0]
/// |       ,'       |       ,'      p:6 [ 1.0,  1.0,  1.0]
/// |     ,'         |     ,'        p:7 [-1.0,  1.0,  1.0]
/// |   ,'           |   ,'
/// | ,'             | ,'
/// 1________________2'
/// ```
///
/// # Local IDs of edges
///
/// ```text
///                     7                   p0  p1
///             +----------------+     e:0  [0, 1]
///           ,'|              ,'|     e:1  [1, 2]
///       4 ,'  |8          6,'  |     e:2  [2, 3]
///       ,'    |          ,'    |     e:3  [3, 0]
///     ,'      |   5    ,'      |11   e:4  [4, 5]
///   +'===============+'        |     e:5  [5, 6]
///   |         |      |         |     e:6  [6, 7]
///   |         |      |  3      |     e:7  [7, 4]
///   |         .- - - | -  - - -+     e:8  [0, 4]
///  9|       ,'       |       ,'      e:9  [1, 5]
///   |    0,'         |10   ,'        e:10 [2, 6]
///   |   ,'           |   ,' 2        e:11 [3, 7]
///   | ,'             | ,'
///   +----------------+'
///           1
/// ```
///
/// # Local IDs of faces
///
/// Note: the order of nodes is such that the right-hand rule generates outward normals.
/// Also, the order of face nodes corresponds to Qua4 nodes.
///
/// ```text
///           4----------------7
///         ,'|              ,'|
///       ,'  |  ___       ,'  |
///     ,'    |,'5,'  [0],'    |        p0 p1 p2 p3
///   ,'      |~~~     ,'      |    f:0 [0, 4, 7, 3]
/// 5'===============6'  ,'|   |    f:1 [1, 2, 6, 5]
/// |   ,'|   |      |   |3|   |    f:2 [0, 1, 5, 4]
/// |   |2|   |      |   |,'   |    f:3 [2, 3, 7, 6]
/// |   |,'   0- - - | +- - - -3    f:4 [0, 3, 2, 1]
/// |       ,'       |       ,'     f:5 [4, 5, 6, 7]
/// |     ,' [1]  ___|     ,'
/// |   ,'      ,'4,'|   ,'
/// | ,'        ~~~  | ,'
/// 1----------------2'
/// ```
pub struct Hex8 {}

impl Hex8 {
    pub const NDIM: usize = 3;
    pub const NNODE: usize = 8;
    pub const NEDGE: usize = 12;
    pub const NFACE: usize = 6;
    pub const EDGE_NNODE: usize = 2;
    pub const FACE_NNODE: usize = 4;
    pub const FACE_NEDGE: usize = 4;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Hex8::EDGE_NNODE]; Hex8::NEDGE] = [
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 0],
        [4, 5],
        [5, 6],
        [6, 7],
        [7, 4],
        [0, 4],
        [1, 5],
        [2, 6],
        [3, 7],
    ];

    #[rustfmt::skip]
    pub const FACE_NODE_IDS: [[usize; Hex8::FACE_NNODE]; Hex8::NFACE] = [
        [0, 4, 7, 3],
        [1, 2, 6, 5],
        [0, 1, 5, 4],
        [2, 3, 7, 6],
        [0, 3, 2, 1],
        [4, 5, 6, 7],
    ];

    #[rustfmt::skip]
    pub const FACE_EDGE_NODE_IDS: [[[usize; Hex8::EDGE_NNODE]; Hex8::FACE_NEDGE]; Hex8::NFACE] = [
        [[0, 4], [4, 7], [7, 3], [3, 0]],
        [[1, 2], [2, 6], [6, 5], [5, 1]],
        [[0, 1], [1, 5], [5, 4], [4, 0]],
        [[2, 3], [3, 7], [7, 6], [6, 2]],
        [[0, 3], [3, 2], [2, 1], [1, 0]],
        [[4, 5], [5, 6], [6, 7], [7, 4]],
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Hex8::NDIM]; Hex8::NNODE] = [
        [-1.0, -1.0, -1.0],
        [ 1.0, -1.0, -1.0],
        [ 1.0,  1.0, -1.0],
        [-1.0,  1.0, -1.0],
        [-1.0, -1.0,  1.0],
        [ 1.0, -1.0,  1.0],
        [ 1.0,  1.0,  1.0],
        [-1.0,  1.0,  1.0],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        interp[0] = (1.0 - r - s + r * s - t + s * t + r * t - r * s * t) / 8.0;
        interp[1] = (1.0 + r - s - r * s - t + s * t - r * t + r * s * t) / 8.0;
        interp[2] = (1.0 + r + s + r * s - t - s * t - r * t - r * s * t) / 8.0;
        interp[3] = (1.0 - r + s - r * s - t - s * t + r * t + r * s * t) / 8.0;
        interp[4] = (1.0 - r - s + r * s + t - s * t - r * t + r * s * t) / 8.0;
        interp[5] = (1.0 + r - s - r * s + t - s * t + r * t - r * s * t) / 8.0;
        interp[6] = (1.0 + r + s + r * s + t + s * t + r * t + r * s * t) / 8.0;
        interp[7] = (1.0 - r + s - r * s + t + s * t - r * t - r * s * t) / 8.0;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &[f64]) {
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        deriv[0][0] = (-1.0 + s + t - s * t) / 8.0;
        deriv[0][1] = (-1.0 + r + t - r * t) / 8.0;
        deriv[0][2] = (-1.0 + r + s - r * s) / 8.0;

        deriv[1][0] = (1.0 - s - t + s * t) / 8.0;
        deriv[1][1] = (-1.0 - r + t + r * t) / 8.0;
        deriv[1][2] = (-1.0 - r + s + r * s) / 8.0;

        deriv[2][0] = (1.0 + s - t - s * t) / 8.0;
        deriv[2][1] = (1.0 + r - t - r * t) / 8.0;
        deriv[2][2] = (-1.0 - r - s - r * s) / 8.0;

        deriv[3][0] = (-1.0 - s + t + s * t) / 8.0;
        deriv[3][1] = (1.0 - r - t + r * t) / 8.0;
        deriv[3][2] = (-1.0 + r - s + r * s) / 8.0;

        deriv[4][0] = (-1.0 + s - t + s * t) / 8.0;
        deriv[4][1] = (-1.0 + r - t + r * t) / 8.0;
        deriv[4][2] = (1.0 - r - s + r * s) / 8.0;

        deriv[5][0] = (1.0 - s + t - s * t) / 8.0;
        deriv[5][1] = (-1.0 - r - t - r * t) / 8.0;
        deriv[5][2] = (1.0 + r - s - r * s) / 8.0;

        deriv[6][0] = (1.0 + s + t + s * t) / 8.0;
        deriv[6][1] = (1.0 + r + t + r * t) / 8.0;
        deriv[6][2] = (1.0 + r + s + r * s) / 8.0;

        deriv[7][0] = (-1.0 - s - t - s * t) / 8.0;
        deriv[7][1] = (1.0 - r + t - r * t) / 8.0;
        deriv[7][2] = (1.0 - r + s - r * s) / 8.0;
    }
}
