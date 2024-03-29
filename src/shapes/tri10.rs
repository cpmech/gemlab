use russell_lab::{Matrix, Vector};

/// Defines a triangle with 10 nodes (cubic edges; interior node)
///
/// ![tri10](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_tri10.svg)
///
/// # Local IDs of nodes
///
/// ```text
/// s
/// |
/// 2, (0,1)
/// | ',
/// |   ',
/// 5     7,
/// |       ',
/// |         ',
/// 8     9     4,
/// |             ',
/// | (0,0)         ', (1,0)
/// 0-----3-----6-----1 ---- r
/// ```
///
/// # Local IDs of edges
///
/// ```text
///   2,
///   | ',
///   |   ',              p0  p1 p2 p3
///   5     7,        e:0 [1, 0, 6, 3]
///   |       ', 1    e:1 [2, 1, 7, 4]
/// 2 |         ',    e:2 [0, 2, 8, 5]
///   8     9     4,
///   |             ',
///   |               ',
///   0-----3-----6-----1
///            0
/// ```
///
/// # Notes
///
/// * The order of edge nodes is such that the normals are outward
/// * The order of edge nodes corresponds to [super::Lin4] nodes
pub struct Tri10 {}

impl Tri10 {
    pub const GEO_NDIM: usize = 2;
    pub const NNODE: usize = 10;
    pub const NEDGE: usize = 3;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 4;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;
    pub const N_INTERIOR_NODE: usize = 1;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Tri10::EDGE_NNODE]; Tri10::NEDGE] = [
        [1, 0, 6, 3],
        [2, 1, 7, 4],
        [0, 2, 8, 5],
    ];

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS_INWARD: [[usize; Tri10::EDGE_NNODE]; Tri10::NEDGE] = [
        [0, 1, 3, 6],
        [1, 2, 4, 7],
        [2, 0, 5, 8],
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Tri10::GEO_NDIM]; Tri10::NNODE] = [
        [0.0       , 0.0      ], // 0
        [1.0       , 0.0      ], // 1
        [0.0       , 1.0      ], // 2
        [1.0 / 3.0 , 0.0      ], // 3
        [2.0 / 3.0 , 1.0 / 3.0], // 4
        [0.0       , 2.0 / 3.0], // 5
        [2.0 / 3.0 , 0.0      ], // 6
        [1.0 / 3.0 , 2.0 / 3.0], // 7
        [0.0       , 1.0 / 3.0], // 8
        [1.0 / 3.0 , 1.0 / 3.0], // 9
    ];

    pub const INTERIOR_NODES: [usize; Tri10::N_INTERIOR_NODE] = [9];

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
        let (r, s) = (ksi[0], ksi[1]);

        let z = 1.0 - r - s;
        let t1 = s * (3.0 * s - 1.0);
        let t2 = z * (3.0 * z - 1.0);
        let t3 = r * (3.0 * r - 1.0);

        interp[0] = 0.5 * t2 * (3.0 * z - 2.0);
        interp[1] = 0.5 * t3 * (3.0 * r - 2.0);
        interp[2] = 0.5 * t1 * (3.0 * s - 2.0);
        interp[3] = 4.5 * r * t2;
        interp[4] = 4.5 * s * t3;
        interp[5] = 4.5 * z * t1;
        interp[6] = 4.5 * z * t3;
        interp[7] = 4.5 * r * t1;
        interp[8] = 4.5 * s * t2;
        interp[9] = 27.0 * s * z * r;
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
        let (r, s) = (ksi[0], ksi[1]);

        let z = 1.0 - r - s;

        let q0 = 4.5 * (6.0 * z - 1.0);
        let q1 = 4.5 * s * (3.0 * s - 1.0);
        let q2 = 4.5 * z * (3.0 * z - 1.0);
        let q3 = 4.5 * r * (3.0 * r - 1.0);
        let q4 = 4.5 * (6.0 * s - 1.0);
        let q5 = 4.5 * (6.0 * r - 1.0);
        let q6 = q0 * s;
        let q7 = q0 * r;
        let q8 = -0.5 * (27.0 * z * z - 18.0 * z + 2.0);
        let q9 = 0.5 * (27.0 * s * s - 18.0 * s + 2.0);
        let q10 = 0.5 * (27.0 * r * r - 18.0 * r + 2.0);

        deriv.set(0, 0, q8);
        deriv.set(1, 0, q10);
        deriv.set(2, 0, 0.0);
        deriv.set(3, 0, q2 - q7);
        deriv.set(4, 0, s * q5);
        deriv.set(5, 0, -q1);
        deriv.set(6, 0, z * q5 - q3);
        deriv.set(7, 0, q1);
        deriv.set(8, 0, -q6);
        deriv.set(9, 0, 27.0 * s * (z - r));

        deriv.set(0, 1, q8);
        deriv.set(1, 1, 0.0);
        deriv.set(2, 1, q9);
        deriv.set(3, 1, -q7);
        deriv.set(4, 1, q3);
        deriv.set(5, 1, z * q4 - q1);
        deriv.set(6, 1, -q3);
        deriv.set(7, 1, r * q4);
        deriv.set(8, 1, q2 - q6);
        deriv.set(9, 1, 27.0 * r * (z - s));
    }
}
