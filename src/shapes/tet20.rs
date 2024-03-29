use russell_lab::math::{ONE_BY_3, TWO_BY_3};
use russell_lab::{Matrix, Vector};

/// Defines a tetrahedron with 20 nodes (cubic faces)
///
/// ![tet20](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_tet20.svg)
///
/// # Local IDs of nodes
///
/// ```text
///               t
///               |                             node    r    s    t
///               3                                0 [0.0, 0.0, 0.0]
///              /|`.                              1 [1.0, 0.0, 0.0]
///              ||  `,                            2 [0.0, 1.0, 0.0]
///             / |    ',                          3 [0.0, 0.0, 1.0]
///             | |      \
///            /  9       11                       4 [1/3, 0.0, 0.0]
///            |  |         `,                     5 [2/3, 0.0, 0.0]
///          13   |           ',                   6 [0.0, 1/3, 0.0]
///           |   |             \                  7 [0.0, 2/3, 0.0]
///          /    8              `,                8 [0.0, 0.0, 1/3]
///          |    |  ,'\ 16        ',              9 [0.0, 0.0, 2/3]
///         / 17  |  \19\           10
///       12      0.,,\,'             `.          10 [0.0, 2/3, 1/3]
///        |     /     ``'6.,,__        `.        11 [0.0, 1/3, 2/3]
///        |    /              ``''7.,,_  ',      12 [2/3, 0.0, 1/3]
///       |    4      18                ``-2-,s   13 [1/3, 0.0, 2/3]
///       |   /                       ,.-``       14 [2/3, 1/3, 0.0]
///      |   /                   15-'`            15 [1/3, 2/3, 0.0]
///      ' 5                 ,.'`
///     | /             _,-``                     16 [0.0, 1/3, 1/3]
///     '/          14'`                          17 [1/3, 0.0, 1/3]
///    |/      ,.-``                              18 [1/3, 1/3, 0.0]
///    /  _,-``                                   19 [1/3, 1/3, 1/3]
///   1 '`
///  /
/// r
/// ```
///
/// # Local IDs of edges
///
/// ```text
///               t                    p0  p1 p2  p3
///               |                e:0 [0, 1,  4,  5]
///               3                e:1 [1, 2, 14, 15]
///              /|`.              e:2 [2, 0,  7,  6]
///              ||  `,            e:3 [0, 3,  8,  9]
///             / |    ',          e:4 [1, 3, 12, 13]
///             | |      \         e:5 [2, 3, 10, 11]
///            /  9       11
///            |  |         `,
///          13   |(3)        `.(5)
///           |   |             \
///          /    8              `.
///       (4)|    |                ',
///         /     |                  10
///       12      0.,,_               `.
///        |     /     ``'6.,,__        `.
///        |    /          (2) ``''7.,,_  ',
///       |    4                        ``'2-,_s
///       |   /                       ,.-``
///      |   (0)                 15-'`
///      ' 5                 ,.'`
///     | /             _,-``
///     '/          14'` (1)
///    |/      ,.-``
///    /  _,-``
///   1,'`
///  /
/// r
/// ```
///
/// * The order of edge nodes corresponds to **Lin4** nodes.
///
/// # Local IDs of faces
///
/// ```text
///               t           p0 p1 p2  p3  p4  p5  p6  p7  p8  p9
///               |       f:0 [0, 3, 2,  8, 11,  7,  9, 10,  6, 16]
///               3       f:1 [0, 1, 3,  4, 12,  9,  5, 13,  8, 17]
///              /|`.     f:2 [0, 2, 1,  6, 15,  5,  7, 14,  4, 18]
///              ||  `,   f:3 [1, 2, 3, 14, 10, 13, 15, 11, 12, 19]
///             / |    ',
///             | |      \
///            /  |       `.
///            |  |         `,
///           /   |           `.
///           |   |     [0]     \
///          /    |              `.
///          |    |                ',
///         / [1] |    ,             \
///         |     0. ,' \             `.
///        |     /   \ 3 \ .,.__        `.
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
/// * The order of face nodes corresponds to [super::Tri10] nodes
pub struct Tet20 {}

impl Tet20 {
    pub const GEO_NDIM: usize = 3;
    pub const NNODE: usize = 20;
    pub const NEDGE: usize = 6;
    pub const NFACE: usize = 4;
    pub const EDGE_NNODE: usize = 4;
    pub const FACE_NNODE: usize = 10;
    pub const FACE_NEDGE: usize = 3;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Tet20::EDGE_NNODE]; Tet20::NEDGE] = [
        [0, 1,  4,  5],
        [1, 2, 14, 15],
        [2, 0,  7,  6],
        [0, 3,  8,  9],
        [1, 3, 12, 13],
        [2, 3, 10, 11],
    ];

    #[rustfmt::skip]
    pub const FACE_NODE_IDS: [[usize; Tet20::FACE_NNODE]; Tet20::NFACE] = [
        [0, 3, 2,  8, 11,  7,  9, 10,  6, 16],
        [0, 1, 3,  4, 12,  9,  5, 13,  8, 17],
        [0, 2, 1,  6, 15,  5,  7, 14,  4, 18],
        [1, 2, 3, 14, 10, 13, 15, 11, 12, 19],
    ];

    #[rustfmt::skip]
    pub const FACE_EDGE_NODE_IDS: [[[usize; Tet20::EDGE_NNODE]; Tet20::FACE_NEDGE]; Tet20::NFACE] = [
        [[0, 3,  8,  9], [3, 2, 11, 10], [2, 0,  7,  6]],
        [[0, 1,  4,  5], [1, 3, 12, 13], [3, 0,  9,  8]],
        [[0, 2,  6,  7], [2, 1, 15, 14], [1, 0,  5,  4]],
        [[1, 2, 14, 15], [2, 3, 10, 11], [3, 1, 13, 12]],
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Tet20::GEO_NDIM]; Tet20::NNODE] = [
        [0.0, 0.0, 0.0], //  0
        [1.0, 0.0, 0.0], //  1
        [0.0, 1.0, 0.0], //  2
        [0.0, 0.0, 1.0], //  3

        [ONE_BY_3, 0.0, 0.0], //  4
        [TWO_BY_3, 0.0, 0.0], //  5
        [0.0, ONE_BY_3, 0.0], //  6
        [0.0, TWO_BY_3, 0.0], //  7
        [0.0, 0.0, ONE_BY_3], //  8
        [0.0, 0.0, TWO_BY_3], //  9

        [0.0, TWO_BY_3, ONE_BY_3], // 10
        [0.0, ONE_BY_3, TWO_BY_3], // 11
        [TWO_BY_3, 0.0, ONE_BY_3], // 12
        [ONE_BY_3, 0.0, TWO_BY_3], // 13
        [TWO_BY_3, ONE_BY_3, 0.0], // 14
        [ONE_BY_3, TWO_BY_3, 0.0], // 15

        [0.0, ONE_BY_3, ONE_BY_3], // 16
        [ONE_BY_3, 0.0, ONE_BY_3], // 17
        [ONE_BY_3, ONE_BY_3, 0.0], // 18

        [ONE_BY_3, ONE_BY_3, ONE_BY_3], // 19
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

        let u = 1.0 - r - s - t;

        interp[0] = (3.0 * u - 1.0) * (3.0 * u - 2.0) * u / 2.0; // corner r=0,s=0,t=0,u=1
        interp[1] = (3.0 * r - 1.0) * (3.0 * r - 2.0) * r / 2.0; // corner r=1,s=0,t=0,u=0
        interp[2] = (3.0 * s - 1.0) * (3.0 * s - 2.0) * s / 2.0; // corner r=0,s=1,t=0,u=0
        interp[3] = (3.0 * t - 1.0) * (3.0 * t - 2.0) * t / 2.0; // corner r=0,s=0,t=1,u=0
        interp[4] = 9.0 * (3.0 * u - 1.0) * u * r / 2.0; // along r r=1/3,u=2/3
        interp[5] = 9.0 * (3.0 * r - 1.0) * u * r / 2.0; // along r r=2/3,u=1/3
        interp[6] = 9.0 * (3.0 * u - 1.0) * u * s / 2.0; // along s s=1/3,u=2/3
        interp[7] = 9.0 * (3.0 * s - 1.0) * u * s / 2.0; // along s s=2/3,u=1/3
        interp[8] = 9.0 * (3.0 * u - 1.0) * u * t / 2.0; // along t t=1/3,u=2/3
        interp[9] = 9.0 * (3.0 * t - 1.0) * u * t / 2.0; // along t t=2/3,u=1/3
        interp[10] = 9.0 * (3.0 * s - 1.0) * s * t / 2.0; // diagonal of plane normal r s=2/3,t=1/3
        interp[11] = 9.0 * (3.0 * t - 1.0) * s * t / 2.0; // diagonal of plane normal r s=1/3,t=2/3
        interp[12] = 9.0 * (3.0 * r - 1.0) * r * t / 2.0; // diagonal of plane normal s r=2/3,t=1/3
        interp[13] = 9.0 * (3.0 * t - 1.0) * r * t / 2.0; // diagonal of plane normal s r=1/3,t=2/3
        interp[14] = 9.0 * (3.0 * r - 1.0) * r * s / 2.0; // diagonal of plane normal t r=2/3,s=1/3
        interp[15] = 9.0 * (3.0 * s - 1.0) * r * s / 2.0; // diagonal of plane normal t r=1/3,s=2/3
        interp[16] = 27.0 * s * t * u; // plane normal r s=1/3,t=1/3,u=1/3
        interp[17] = 27.0 * r * t * u; // plane normal s r=1/3,t=1/3,u=1/3
        interp[18] = 27.0 * r * s * u; // plane normal t r=1/3,s=1/3,u=1/3
        interp[19] = 27.0 * r * s * t; // plane normal (1,1,1) r=1/3,s=1/3,t=1/3,u=0
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

        let u = 1.0 - r - s - t;

        // with respect to r
        deriv.set(
            0,
            0,
            (-3.0 * u * (-2.0 + 3.0 * u)) / 2.0
                - (3.0 * u * (-1.0 + 3.0 * u)) / 2.0
                - ((-2.0 + 3.0 * u) * (-1.0 + 3.0 * u)) / 2.0,
        );
        deriv.set(
            1,
            0,
            (3.0 * r * (-2.0 + 3.0 * r)) / 2.0
                + (3.0 * r * (-1.0 + 3.0 * r)) / 2.0
                + ((-2.0 + 3.0 * r) * (-1.0 + 3.0 * r)) / 2.0,
        );
        deriv.set(2, 0, 0.0);
        deriv.set(3, 0, 0.0);
        deriv.set(
            4,
            0,
            (-27.0 * r * u) / 2.0 - (9.0 * r * (-1.0 + 3.0 * u)) / 2.0 + (9.0 * u * (-1.0 + 3.0 * u)) / 2.0,
        );
        deriv.set(
            5,
            0,
            (-9.0 * r * (-1.0 + 3.0 * r)) / 2.0 + (27.0 * r * u) / 2.0 + (9.0 * (-1.0 + 3.0 * r) * u) / 2.0,
        );
        deriv.set(6, 0, (-27.0 * s * u) / 2.0 - (9.0 * s * (-1.0 + 3.0 * u)) / 2.0);
        deriv.set(7, 0, (-9.0 * s * (-1.0 + 3.0 * s)) / 2.0);
        deriv.set(8, 0, (-27.0 * t * u) / 2.0 - (9.0 * t * (-1.0 + 3.0 * u)) / 2.0);
        deriv.set(9, 0, (-9.0 * t * (-1.0 + 3.0 * t)) / 2.0);
        deriv.set(10, 0, 0.0);
        deriv.set(11, 0, 0.0);
        deriv.set(12, 0, (27.0 * r * t) / 2.0 + (9.0 * (-1.0 + 3.0 * r) * t) / 2.0);
        deriv.set(13, 0, (9.0 * t * (-1.0 + 3.0 * t)) / 2.0);
        deriv.set(14, 0, (27.0 * r * s) / 2.0 + (9.0 * (-1.0 + 3.0 * r) * s) / 2.0);
        deriv.set(15, 0, (9.0 * s * (-1.0 + 3.0 * s)) / 2.0);
        deriv.set(16, 0, -27.0 * s * t);
        deriv.set(17, 0, -27.0 * r * t + 27.0 * t * u);
        deriv.set(18, 0, -27.0 * r * s + 27.0 * s * u);
        deriv.set(19, 0, 27.0 * s * t);

        // with respect to s
        deriv.set(
            0,
            1,
            (-3.0 * u * (-2.0 + 3.0 * u)) / 2.0
                - (3.0 * u * (-1.0 + 3.0 * u)) / 2.0
                - ((-2.0 + 3.0 * u) * (-1.0 + 3.0 * u)) / 2.0,
        );
        deriv.set(1, 1, 0.0);
        deriv.set(
            2,
            1,
            (3.0 * s * (-2.0 + 3.0 * s)) / 2.0
                + (3.0 * s * (-1.0 + 3.0 * s)) / 2.0
                + ((-2.0 + 3.0 * s) * (-1.0 + 3.0 * s)) / 2.0,
        );
        deriv.set(3, 1, 0.0);
        deriv.set(4, 1, (-27.0 * r * u) / 2.0 - (9.0 * r * (-1.0 + 3.0 * u)) / 2.0);
        deriv.set(5, 1, (-9.0 * r * (-1.0 + 3.0 * r)) / 2.0);
        deriv.set(
            6,
            1,
            (-27.0 * s * u) / 2.0 - (9.0 * s * (-1.0 + 3.0 * u)) / 2.0 + (9.0 * u * (-1.0 + 3.0 * u)) / 2.0,
        );
        deriv.set(
            7,
            1,
            (-9.0 * s * (-1.0 + 3.0 * s)) / 2.0 + (27.0 * s * u) / 2.0 + (9.0 * (-1.0 + 3.0 * s) * u) / 2.0,
        );
        deriv.set(8, 1, (-27.0 * t * u) / 2.0 - (9.0 * t * (-1.0 + 3.0 * u)) / 2.0);
        deriv.set(9, 1, (-9.0 * t * (-1.0 + 3.0 * t)) / 2.0);
        deriv.set(10, 1, (27.0 * s * t) / 2.0 + (9.0 * (-1.0 + 3.0 * s) * t) / 2.0);
        deriv.set(11, 1, (9.0 * t * (-1.0 + 3.0 * t)) / 2.0);
        deriv.set(12, 1, 0.0);
        deriv.set(13, 1, 0.0);
        deriv.set(14, 1, (9.0 * r * (-1.0 + 3.0 * r)) / 2.0);
        deriv.set(15, 1, (27.0 * r * s) / 2.0 + (9.0 * r * (-1.0 + 3.0 * s)) / 2.0);
        deriv.set(16, 1, -27.0 * s * t + 27.0 * t * u);
        deriv.set(17, 1, -27.0 * r * t);
        deriv.set(18, 1, -27.0 * r * s + 27.0 * r * u);
        deriv.set(19, 1, 27.0 * r * t);

        // with respect to t
        deriv.set(
            0,
            2,
            (-3.0 * u * (-2.0 + 3.0 * u)) / 2.0
                - (3.0 * u * (-1.0 + 3.0 * u)) / 2.0
                - ((-2.0 + 3.0 * u) * (-1.0 + 3.0 * u)) / 2.0,
        );
        deriv.set(1, 2, 0.0);
        deriv.set(2, 2, 0.0);
        deriv.set(
            3,
            2,
            (3.0 * t * (-2.0 + 3.0 * t)) / 2.0
                + (3.0 * t * (-1.0 + 3.0 * t)) / 2.0
                + ((-2.0 + 3.0 * t) * (-1.0 + 3.0 * t)) / 2.0,
        );
        deriv.set(4, 2, (-27.0 * r * u) / 2.0 - (9.0 * r * (-1.0 + 3.0 * u)) / 2.0);
        deriv.set(5, 2, (-9.0 * r * (-1.0 + 3.0 * r)) / 2.0);
        deriv.set(6, 2, (-27.0 * s * u) / 2.0 - (9.0 * s * (-1.0 + 3.0 * u)) / 2.0);
        deriv.set(7, 2, (-9.0 * s * (-1.0 + 3.0 * s)) / 2.0);
        deriv.set(
            8,
            2,
            (-27.0 * t * u) / 2.0 - (9.0 * t * (-1.0 + 3.0 * u)) / 2.0 + (9.0 * u * (-1.0 + 3.0 * u)) / 2.0,
        );
        deriv.set(
            9,
            2,
            (-9.0 * t * (-1.0 + 3.0 * t)) / 2.0 + (27.0 * t * u) / 2.0 + (9.0 * (-1.0 + 3.0 * t) * u) / 2.0,
        );
        deriv.set(10, 2, (9.0 * s * (-1.0 + 3.0 * s)) / 2.0);
        deriv.set(11, 2, (27.0 * s * t) / 2.0 + (9.0 * s * (-1.0 + 3.0 * t)) / 2.0);
        deriv.set(12, 2, (9.0 * r * (-1.0 + 3.0 * r)) / 2.0);
        deriv.set(13, 2, (27.0 * r * t) / 2.0 + (9.0 * r * (-1.0 + 3.0 * t)) / 2.0);
        deriv.set(14, 2, 0.0);
        deriv.set(15, 2, 0.0);
        deriv.set(16, 2, -27.0 * s * t + 27.0 * s * u);
        deriv.set(17, 2, -27.0 * r * t + 27.0 * r * u);
        deriv.set(18, 2, -27.0 * r * s);
        deriv.set(19, 2, 27.0 * r * s);
    }
}
