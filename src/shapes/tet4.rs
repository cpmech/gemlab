use russell_lab::{Matrix, Vector};

/// Defines a tetrahedron with 4 points (linear faces)
///
/// # Local IDs of points
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
///            |  |         `,
///           /   |           `,
///           |   |             \
///          /    |              `.
///          |    |                ',
///         /     |                  \
///         |     @.,,_               `.
///        |     / 0   ``'-.,,__        `.
///        |    /              ``''-.,,_  ',
///       |    /                        ``-@-,_s
///       |  ,'                       ,.-`` 2
///      |  ,                    _,-'`
///      ' /                 ,.'`
///     | /             _.-``
///     '/          ,-'`
///    |/      ,.-``
///    /  _,-``
///   @ '`
///  / 1
/// r
/// ```
///
/// # Local IDs of edges
///
/// The order of edge points corresponds to Lin2 points.
///
/// ```text
///               t                    p0  p1
///               |                e:0 [0, 1]
///               3                e:1 [1, 2]
///              /|`.              e:2 [2, 0]
///              ||  `,            e:3 [0, 3]
///             / |    ',          e:4 [1, 3]
///             | |      \         e:5 [2, 3]
///            /  |       `.
///            |  |         `,
///           /   |(3)        `,(5)
///           |   |             \
///          /    |              `.
///       (4)|    |                ',
///         /     |                  \
///         |     0.,,_               `.
///        |     /     ``'-.,,__        `.
///        |    /          (2) ``''-.,,_  ',
///       |    /                        ``'2-,_s
///       |  ,'                       ,.-``
///      |  ,(0)                 _,-'`
///      ' /                 ,.'`
///     | /             _.-``
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
/// Note: the order of points is such that the right-hand rule generates outward normals.
/// Also, the order of face points corresponds to Tri3 points.
///
/// ```text
///               t
///               |
///               3                   p0 p1 p2
///              /|`.             f:0 [0, 3, 2]
///              ||  `,           f:1 [0, 1, 3]
///             / |    ',         f:2 [0, 2, 1]
///             | |      \        f:3 [1, 2, 3]
///            /  |       `.
///            |  |         `,
///           /   |           `,
///           |   |     [0]     \
///          /    |              `.
///          |    |                ',
///         / [1] |    ,             \
///         |     0. ,' \             `.
///        |     /   \ 3 \ .,,__        `.
///        |    /     \ ,'     ``''-.,,_  ',
///       |    /       '                ``'2-,_s
///       |  ,'                       ,.-``
///      |  ,      [2]           _,-'`
///      ' /                 ,.'`
///     | /             _.-``
///     '/          ,-'`
///    |/      ,.-``
///    /  _,-``
///   1,'`
///  /
/// r
/// ```
pub struct Tet4 {}

impl Tet4 {
    pub const NDIM: usize = 3;
    pub const NPOINT: usize = 4;
    pub const NEDGE: usize = 6;
    pub const NFACE: usize = 4;
    pub const EDGE_NPOINT: usize = 2;
    pub const FACE_NPOINT: usize = 3;
    pub const FACE_NEDGE: usize = 3;

    #[rustfmt::skip]
    pub const EDGE_POINT_IDS: [[usize; Tet4::EDGE_NPOINT]; Tet4::NEDGE] = [
        [0, 1],
        [1, 2],
        [2, 0],
        [0, 3],
        [1, 3],
        [2, 3],
    ];

    #[rustfmt::skip]
    pub const FACE_POINT_IDS: [[usize; Tet4::FACE_NPOINT]; Tet4::NFACE] = [
        [0, 3, 2],
        [0, 1, 3],
        [0, 2, 1],
        [1, 2, 3],
    ];

    #[rustfmt::skip]
    pub const FACE_EDGE_POINT_IDS: [[[usize; Tet4::EDGE_NPOINT]; Tet4::FACE_NEDGE]; Tet4::NFACE] = [
        [[0, 3], [3, 2], [2, 0]],
        [[0, 1], [1, 3], [3, 0]],
        [[0, 2], [2, 1], [1, 0]],
        [[1, 2], [2, 3], [3, 1]],
    ];

    #[rustfmt::skip]
    pub const POINT_REFERENCE_COORDS: [[f64; Tet4::NDIM]; Tet4::NPOINT] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &Vector) {
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        interp[0] = 1.0 - r - s - t;
        interp[1] = r;
        interp[2] = s;
        interp[3] = t;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, _: &Vector) {
        deriv[0][0] = -1.0;
        deriv[1][0] = 1.0;
        deriv[2][0] = 0.0;
        deriv[3][0] = 0.0;

        deriv[0][1] = -1.0;
        deriv[1][1] = 0.0;
        deriv[2][1] = 1.0;
        deriv[3][1] = 0.0;

        deriv[0][2] = -1.0;
        deriv[1][2] = 0.0;
        deriv[2][2] = 0.0;
        deriv[3][2] = 1.0;
    }
}