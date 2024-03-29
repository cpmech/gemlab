use russell_lab::{Matrix, Vector};

/// Defines a tetrahedron with 4 nodes (linear faces)
///
/// ![tet4](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_tet4.svg)
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
/// * The order of edge nodes corresponds to *Lin2* nodes.
///
/// # Local IDs of faces
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
///
/// # Notes
///
/// * The order of face nodes is such that the normals are outward
/// * The order of face nodes corresponds to [super::Tri3] nodes
pub struct Tet4 {}

impl Tet4 {
    pub const GEO_NDIM: usize = 3;
    pub const NNODE: usize = 4;
    pub const NEDGE: usize = 6;
    pub const NFACE: usize = 4;
    pub const EDGE_NNODE: usize = 2;
    pub const FACE_NNODE: usize = 3;
    pub const FACE_NEDGE: usize = 3;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Tet4::EDGE_NNODE]; Tet4::NEDGE] = [
        [0, 1],
        [1, 2],
        [2, 0],
        [0, 3],
        [1, 3],
        [2, 3],
    ];

    #[rustfmt::skip]
    pub const FACE_NODE_IDS: [[usize; Tet4::FACE_NNODE]; Tet4::NFACE] = [
        [0, 3, 2],
        [0, 1, 3],
        [0, 2, 1],
        [1, 2, 3],
    ];

    #[rustfmt::skip]
    pub const FACE_EDGE_NODE_IDS: [[[usize; Tet4::EDGE_NNODE]; Tet4::FACE_NEDGE]; Tet4::NFACE] = [
        [[0, 3], [3, 2], [2, 0]],
        [[0, 1], [1, 3], [3, 0]],
        [[0, 2], [2, 1], [1, 0]],
        [[1, 2], [2, 3], [3, 1]],
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Tet4::GEO_NDIM]; Tet4::NNODE] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
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

        interp[0] = 1.0 - r - s - t;
        interp[1] = r;
        interp[2] = s;
        interp[3] = t;
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
    pub fn calc_deriv(deriv: &mut Matrix, _: &[f64]) {
        deriv.set(0, 0, -1.0);
        deriv.set(1, 0, 1.0);
        deriv.set(2, 0, 0.0);
        deriv.set(3, 0, 0.0);

        deriv.set(0, 1, -1.0);
        deriv.set(1, 1, 0.0);
        deriv.set(2, 1, 1.0);
        deriv.set(3, 1, 0.0);

        deriv.set(0, 2, -1.0);
        deriv.set(1, 2, 0.0);
        deriv.set(2, 2, 0.0);
        deriv.set(3, 2, 1.0);
    }
}
