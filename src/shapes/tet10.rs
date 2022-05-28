use russell_lab::{Matrix, Vector};

/// Defines a tetrahedron with 10 nodes (quadratic faces)
///
/// # Local IDs of nodes
///
/// ```text
///               t               r    s    t
///               |         p:0 [0.0, 0.0, 0.0]
///               @ 3       p:1 [1.0, 0.0, 0.0]
///              /|`.       p:2 [0.0, 1.0, 0.0]
///              ||  `,     p:3 [0.0, 0.0, 1.0]
///             / |    ',
///             | |      \
///            /  |       `.
///            |  |         `,  9
///           /   @ 7         @,
///           |   |             \
///          /    |              `.
///        8 @    |                ',
///         /     |                  \
///         |     @.,,_               `.
///        |     / 0   ``'-.,@__        `.
///        |    /           6  ``''-.,,_  ',
///       |    /                        ``-@-,_s
///       |  @'                       ,.-`` 2
///      |  , 4                  _,-'`
///      ' /            5    ,.'`
///     | /             _@-``      r    s    t
///     '/          ,-'`     p:4 [0.5, 0.0, 0.0]
///    |/      ,.-``         p:5 [0.5, 0.5, 0.0]
///    /  _,-``              p:6 [0.0, 0.5, 0.0]
///   @ '`                   p:7 [0.0, 0.0, 0.5]
///  / 1                     p:8 [0.5, 0.0, 0.5]
/// r                        p:9 [0.0, 0.5, 0.5]
/// ```
///
/// # Local IDs of edges
///
/// The order of edge nodes corresponds to Lin2 nodes.
///
/// ```text
///               t                    p0  p1 p2
///               |                e:0 [0, 1, 4]
///               3                e:1 [1, 2, 5]
///              /|`.              e:2 [2, 0, 6]
///              ||  `,            e:3 [0, 3, 7]
///             / |    ',          e:4 [1, 3, 8]
///             | |      \         e:5 [2, 3, 9]
///            /  |       `.
///            |  |         `,
///           /   7(3)        `9(5)
///           |   |             \
///          /    |              `.
///       (4)8    |                ',
///         /     |                  \
///         |     0.,,_               `.
///        |     /     ``'-.,6__        `.
///        |    /          (2) ``''-.,,_  ',
///       |    /                        ``'2-,_s
///       |  4'                       ,.-``
///      |  ,(0)                 _,-'`
///      ' /                 ,.'`
///     | /             _5-``
///     '/          ,-'` (1)
///    |/      ,.-``
///    /  _,-``
///   1,'`
///  /
/// r
/// ```
///
/// # Local IDs of faces
///
/// ```text
///               t           p0 p1 p2 p3 p4 p5
///               |       f:0 [0, 3, 2, 7, 9, 6]
///               3       f:1 [0, 1, 3, 4, 8, 7]
///              /|`.     f:2 [0, 2, 1, 6, 5, 4]
///              ||  `,   f:3 [1, 2, 3, 5, 9, 8]
///             / |    ',
///             | |      \
///            /  |       `.
///            |  |         `,
///           /   7           `9
///           |   |     [0]     \
///          /    |              `.
///          8    |                ',
///         / [1] |    ,             \
///         |     0. ,' \             `.
///        |     /   \ 3 \ .,6__        `.
///        |    /     \ ,'     ``''-.,,_  ',
///       |    /       '                ``'2-,_s
///       |  4'                       ,.-``
///      |  ,      [2]           _,-'`
///      ' /                 ,.'`
///     | /             _5-``
///     '/          ,-'`
///    |/      ,.-``
///    /  _,-``
///   1,'`
///  /
/// r
/// ```
///
/// # Note about face nodes
///
/// * The order of face nodes is such that the normals are outward
/// * The order of face nodes corresponds to **Tri6** nodes
pub struct Tet10 {}

impl Tet10 {
    pub const GEO_NDIM: usize = 3;
    pub const NNODE: usize = 10;
    pub const NEDGE: usize = 6;
    pub const NFACE: usize = 4;
    pub const EDGE_NNODE: usize = 3;
    pub const FACE_NNODE: usize = 6;
    pub const FACE_NEDGE: usize = 3;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Tet10::EDGE_NNODE]; Tet10::NEDGE] = [
        [0, 1, 4],
        [1, 2, 5],
        [2, 0, 6],
        [0, 3, 7],
        [1, 3, 8],
        [2, 3, 9],
    ];

    #[rustfmt::skip]
    pub const FACE_NODE_IDS: [[usize; Tet10::FACE_NNODE]; Tet10::NFACE] = [
        [0, 3, 2, 7, 9, 6],
        [0, 1, 3, 4, 8, 7],
        [0, 2, 1, 6, 5, 4],
        [1, 2, 3, 5, 9, 8],
    ];

    #[rustfmt::skip]
    pub const FACE_EDGE_NODE_IDS: [[[usize; Tet10::EDGE_NNODE]; Tet10::FACE_NEDGE]; Tet10::NFACE] = [
        [[0, 3, 7], [3, 2, 9], [2, 0, 6]],
        [[0, 1, 4], [1, 3, 8], [3, 0, 7]],
        [[0, 2, 6], [2, 1, 5], [1, 0, 4]],
        [[1, 2, 5], [2, 3, 9], [3, 1, 8]],
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Tet10::GEO_NDIM]; Tet10::NNODE] = [
        [0.0, 0.0, 0.0], // 0
        [1.0, 0.0, 0.0], // 1
        [0.0, 1.0, 0.0], // 2
        [0.0, 0.0, 1.0], // 3
        [0.5, 0.0, 0.0], // 4
        [0.5, 0.5, 0.0], // 5
        [0.0, 0.5, 0.0], // 6
        [0.0, 0.0, 0.5], // 7
        [0.5, 0.0, 0.5], // 8
        [0.0, 0.5, 0.5], // 9
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        let u = 1.0 - r - s - t;

        interp[0] = u * (2.0 * u - 1.0);
        interp[1] = r * (2.0 * r - 1.0);
        interp[2] = s * (2.0 * s - 1.0);
        interp[3] = t * (2.0 * t - 1.0);
        interp[4] = 4.0 * u * r;
        interp[5] = 4.0 * r * s;
        interp[6] = 4.0 * s * u;
        interp[7] = 4.0 * u * t;
        interp[8] = 4.0 * r * t;
        interp[9] = 4.0 * s * t;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &[f64]) {
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        deriv[0][0] = 4.0 * (r + s + t) - 3.0;
        deriv[1][0] = 4.0 * r - 1.0;
        deriv[2][0] = 0.0;
        deriv[3][0] = 0.0;
        deriv[4][0] = 4.0 - 8.0 * r - 4.0 * s - 4.0 * t;
        deriv[5][0] = 4.0 * s;
        deriv[6][0] = -4.0 * s;
        deriv[7][0] = -4.0 * t;
        deriv[8][0] = 4.0 * t;
        deriv[9][0] = 0.0;

        deriv[0][1] = 4.0 * (r + s + t) - 3.0;
        deriv[1][1] = 0.0;
        deriv[2][1] = 4.0 * s - 1.0;
        deriv[3][1] = 0.0;
        deriv[4][1] = -4.0 * r;
        deriv[5][1] = 4.0 * r;
        deriv[6][1] = 4.0 - 4.0 * r - 8.0 * s - 4.0 * t;
        deriv[7][1] = -4.0 * t;
        deriv[8][1] = 0.0;
        deriv[9][1] = 4.0 * t;

        deriv[0][2] = 4.0 * (r + s + t) - 3.0;
        deriv[1][2] = 0.0;
        deriv[2][2] = 0.0;
        deriv[3][2] = 4.0 * t - 1.0;
        deriv[4][2] = -4.0 * r;
        deriv[5][2] = 0.0;
        deriv[6][2] = -4.0 * s;
        deriv[7][2] = 4.0 - 4.0 * r - 4.0 * s - 8.0 * t;
        deriv[8][2] = 4.0 * r;
        deriv[9][2] = 4.0 * s;
    }
}
