use russell_lab::{Matrix, Vector};

/// Defines a hexahedron with 8 nodes (bilinear faces)
///
/// ![hex8](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_hex8.svg)
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
/// * The order of edge nodes corresponds to **Lin2** nodes.
///
/// # Local IDs of faces
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
///
/// # Notes
///
/// * The reference coordinates range from -1 to +1 with the geometry centred @ 0.
/// * The order of face nodes is such that the normals are outward
/// * The order of face nodes corresponds to [super::Qua4] nodes
/// * This shape is a lower-order version of [super::Hex20]
pub struct Hex8 {}

impl Hex8 {
    pub const GEO_NDIM: usize = 3;
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
    pub const NODE_REFERENCE_COORDS: [[f64; Hex8::GEO_NDIM]; Hex8::NNODE] = [
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

        interp[0] = rm * sm * tm / 8.0;
        interp[1] = rp * sm * tm / 8.0;
        interp[2] = rp * sp * tm / 8.0;
        interp[3] = rm * sp * tm / 8.0;
        interp[4] = rm * sm * tp / 8.0;
        interp[5] = rp * sm * tp / 8.0;
        interp[6] = rp * sp * tp / 8.0;
        interp[7] = rm * sp * tp / 8.0;
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

        // with respect to r
        deriv.set(0, 0, -sm * tm / 8.0);
        deriv.set(1, 0, sm * tm / 8.0);
        deriv.set(2, 0, sp * tm / 8.0);
        deriv.set(3, 0, -sp * tm / 8.0);
        deriv.set(4, 0, -sm * tp / 8.0);
        deriv.set(5, 0, sm * tp / 8.0);
        deriv.set(6, 0, sp * tp / 8.0);
        deriv.set(7, 0, -sp * tp / 8.0);

        // with respect to s
        deriv.set(0, 1, -rm * tm / 8.0);
        deriv.set(1, 1, -rp * tm / 8.0);
        deriv.set(2, 1, rp * tm / 8.0);
        deriv.set(3, 1, rm * tm / 8.0);
        deriv.set(4, 1, -rm * tp / 8.0);
        deriv.set(5, 1, -rp * tp / 8.0);
        deriv.set(6, 1, rp * tp / 8.0);
        deriv.set(7, 1, rm * tp / 8.0);

        // with respect to t
        deriv.set(0, 2, -rm * sm / 8.0);
        deriv.set(1, 2, -rp * sm / 8.0);
        deriv.set(2, 2, -rp * sp / 8.0);
        deriv.set(3, 2, -rm * sp / 8.0);
        deriv.set(4, 2, rm * sm / 8.0);
        deriv.set(5, 2, rp * sm / 8.0);
        deriv.set(6, 2, rp * sp / 8.0);
        deriv.set(7, 2, rm * sp / 8.0);
    }
}
