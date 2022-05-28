use russell_lab::{Matrix, Vector};

/// Defines a triangle with 15 nodes (quartic edges; interior nodes)
///
/// # Local IDs of nodes
///
/// ```text
///  s
///  |
///  2, (0,1)
///  | ',
///  |   ',
/// 10     9,
///  |       ',
///  |         ',
///  5   14      4,
///  |             ',
///  |               ',
/// 11    12    13     8,
///  |                   ',
///  |(0,0)                ', (1,0)
///  0-----6-----3-----7-----1 ---- r
/// ```
///
/// # Local IDs of edges
///
/// ```text
///  |\            p0  p1 p2 p3  p4
///  | \       e:0 [1, 0, 3,  7,  6]
///  |  \ 1    e:1 [2, 1, 4,  9,  8]
/// 2|   \     e:2 [0, 2, 5, 11, 10]
///  |    \
///  |_____\
///     0
/// ```
///
/// # Note about edges
///
/// * The order of edge nodes is such that the normals are outward
/// * The order of edge nodes corresponds to **Lin5** nodes
pub struct Tri15 {}

impl Tri15 {
    pub const GEO_NDIM: usize = 2;
    pub const NNODE: usize = 15;
    pub const NEDGE: usize = 3;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 5;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Tri15::EDGE_NNODE]; Tri15::NEDGE] = [
        [1, 0, 3,  7,  6],
        [2, 1, 4,  9,  8],
        [0, 2, 5, 11, 10],
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Tri15::GEO_NDIM]; Tri15::NNODE] = [
        [0.0  , 0.0 ], //  0
        [1.0  , 0.0 ], //  1
        [0.0  , 1.0 ], //  2
        [0.5  , 0.0 ], //  3
        [0.5  , 0.5 ], //  4
        [0.0  , 0.5 ], //  5
        [0.25 , 0.0 ], //  6
        [0.75 , 0.0 ], //  7
        [0.75 , 0.25], //  8
        [0.25 , 0.75], //  9
        [0.0  , 0.75], // 10
        [0.0  , 0.25], // 11
        [0.25 , 0.25], // 12
        [0.5  , 0.25], // 13
        [0.25 , 0.5 ], // 14
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        let pt1 = 128.0 / 3.0;
        let pt2 = 32.0 / 3.0;
        let cc = 1.0 - r - s;
        let t1 = r - 0.25;
        let t2 = r - 0.5;
        let t3 = r - 0.75;
        let t4 = s - 0.25;
        let t5 = s - 0.5;
        let t6 = s - 0.75;
        let t7 = cc - 0.25;
        let t8 = cc - 0.5;
        let t9 = cc - 0.75;

        interp[0] = pt2 * cc * t7 * t8 * t9;
        interp[1] = pt2 * r * t1 * t2 * t3;
        interp[2] = pt2 * s * t4 * t5 * t6;
        interp[3] = 64.0 * cc * r * t1 * t7;
        interp[4] = 64.0 * r * s * t1 * t4;
        interp[5] = 64.0 * s * cc * t4 * t7;
        interp[6] = pt1 * cc * r * t7 * t8;
        interp[7] = pt1 * cc * r * t1 * t2;
        interp[8] = pt1 * r * s * t1 * t2;
        interp[9] = pt1 * r * s * t4 * t5;
        interp[10] = pt1 * s * cc * t4 * t5;
        interp[11] = pt1 * s * cc * t7 * t8;
        interp[12] = 128.0 * r * s * cc * t7;
        interp[13] = 128.0 * r * s * t1 * cc;
        interp[14] = 128.0 * r * s * cc * t4;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        let pt1 = 128.0 / 3.0;
        let pt2 = 32.0 / 3.0;
        let cc = 1.0 - r - s;
        let t1 = r - 0.25;
        let t2 = r - 0.5;
        let t3 = r - 0.75;
        let t4 = s - 0.25;
        let t5 = s - 0.5;
        let t6 = s - 0.75;
        let t7 = cc - 0.25;
        let t8 = cc - 0.5;
        let t9 = cc - 0.75;

        deriv[0][0] = -pt2 * (t8 * t9 * (t7 + cc) + cc * t7 * (t8 + t9));
        deriv[1][0] = pt2 * (t2 * t3 * (t1 + r) + r * t1 * (t3 + t2));
        deriv[2][0] = 0.0;
        deriv[3][0] = 64.0 * (cc * t7 * (t1 + r) - r * t1 * (t7 + cc));
        deriv[4][0] = 64.0 * s * t4 * (t1 + r);
        deriv[5][0] = -64.0 * s * t4 * (t7 + cc);
        deriv[6][0] = pt1 * (cc * t7 * t8 - r * (t8 * (t7 + cc) + cc * t7));
        deriv[7][0] = pt1 * (cc * (t2 * (t1 + r) + r * t1) - r * t1 * t2);
        deriv[8][0] = pt1 * s * (t2 * (t1 + r) + r * t1);
        deriv[9][0] = pt1 * s * t4 * t5;
        deriv[10][0] = -pt1 * s * t4 * t5;
        deriv[11][0] = -pt1 * s * (t8 * (t7 + cc) + cc * t7);
        deriv[12][0] = 128.0 * s * (cc * t7 - r * (t7 + cc));
        deriv[13][0] = 128.0 * s * (cc * (t1 + r) - r * t1);
        deriv[14][0] = 128.0 * s * t4 * (cc - r);

        deriv[0][1] = -pt2 * (t8 * t9 * (t7 + cc) + cc * t7 * (t8 + t9));
        deriv[1][1] = 0.0;
        deriv[2][1] = pt2 * (t5 * t6 * (t4 + s) + s * t4 * (t6 + t5));
        deriv[3][1] = -64.0 * r * t1 * (t7 + cc);
        deriv[4][1] = 64.0 * r * t1 * (t4 + s);
        deriv[5][1] = 64.0 * (cc * t7 * (t4 + s) - s * t4 * (t7 + cc));
        deriv[6][1] = -pt1 * r * (t8 * (t7 + cc) + cc * t7);
        deriv[7][1] = -pt1 * r * t1 * t2;
        deriv[8][1] = pt1 * r * t1 * t2;
        deriv[9][1] = pt1 * r * (t5 * (t4 + s) + s * t4);
        deriv[10][1] = pt1 * ((cc * (t5 * (t4 + s) + s * t4)) - s * t4 * t5);
        deriv[11][1] = pt1 * (cc * t7 * t8 - s * (t8 * (t7 + cc) + cc * t7));
        deriv[12][1] = 128.0 * r * (cc * t7 - s * (cc + t7));
        deriv[13][1] = 128.0 * r * t1 * (cc - s);
        deriv[14][1] = 128.0 * r * (cc * (t4 + s) - s * t4);
    }
}
