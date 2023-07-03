use russell_lab::{Matrix, Vector};

/// Defines a quadrilateral with 17 nodes (quartic edge; interior node)
///
/// ![qua17](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_qua17.svg)
///
/// The reference coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of nodes
///
/// ```text
///  3----14-----6----13-----2
///  |                  (1,1)|
///  |                       |
/// 15         s ^          12
///  |           |           |
///  |                       |
///  7           8  --> r    5
///  |                       |
///  |         (0,0)         |
/// 16                      11
///  |                       |
///  |(-1,-1)                |
///  0-----9-----4----10-----1
/// ```
///
/// # Local IDs of edges
///
/// ```text
///                2
///    3----14-----6----13-----2
///    |                       |
///    |                       |
///   15                      12         p0 p1  p2 p3  p4
///    |                       |     e:0 [1, 0, 4, 10,  9]
///    |                       |     e:1 [2, 1, 5, 12, 11]
/// 3  7           8           5 1   e:2 [3, 2, 6, 14, 13]
///    |                       |     e:3 [0, 3, 7, 16, 15]
///    |                       |
///   16                      11
///    |                       |
///    |                       |
///    0-----9-----4----10-----1
///                0
/// ```
///
/// # Notes
///
/// * The order of edge nodes is such that the normals are outward
/// * The order of edge nodes corresponds to [super::Lin5] nodes
/// * This shape is a higher-order version of [super::Qua9]
pub struct Qua17 {}

impl Qua17 {
    pub const GEO_NDIM: usize = 2;
    pub const NNODE: usize = 17;
    pub const NEDGE: usize = 4;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 5;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;
    pub const N_INTERIOR_NODE: usize = 1;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Qua17::EDGE_NNODE]; Qua17::NEDGE] = [
        [1, 0, 4, 10, 9],
        [2, 1, 5, 12, 11],
        [3, 2, 6, 14, 13],
        [0, 3, 7, 16, 15]
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Qua17::GEO_NDIM]; Qua17::NNODE] = [
        [-1.0, -1.0], //  0
        [ 1.0, -1.0], //  1
        [ 1.0,  1.0], //  2
        [-1.0,  1.0], //  3
        [ 0.0, -1.0], //  4
        [ 1.0,  0.0], //  5
        [ 0.0,  1.0], //  6
        [-1.0,  0.0], //  7
        [ 0.0,  0.0], //  8
        [-0.5, -1.0], //  9
        [ 0.5, -1.0], // 10
        [ 1.0, -0.5], // 11
        [ 1.0,  0.5], // 12
        [ 0.5,  1.0], // 13
        [-0.5,  1.0], // 14
        [-1.0,  0.5], // 15
        [-1.0, -0.5], // 16
    ];

    pub const INTERIOR_NODES: [usize; Qua17::N_INTERIOR_NODE] = [/*16*/ 8];

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

        let a = 2.0 / 3.0;
        let rr = r * r;
        let ss = s * s;
        let rs = r * s;
        let rp = 1.0 + r;
        let rm = 1.0 - r;
        let sp = 1.0 + s;
        let sm = 1.0 - s;

        interp[0] = rm * sm * (-4.0 * r * (rr - 1.0) - 4.0 * s * (ss - 1.0) + 3.0 * rs) / 12.0;
        interp[1] = rp * sm * (4.0 * r * (rr - 1.0) - 4.0 * s * (ss - 1.0) - 3.0 * rs) / 12.0;
        interp[2] = rp * sp * (4.0 * r * (rr - 1.0) + 4.0 * s * (ss - 1.0) + 3.0 * rs) / 12.0;
        interp[3] = rm * sp * (-4.0 * r * (rr - 1.0) + 4.0 * s * (ss - 1.0) - 3.0 * rs) / 12.0;

        interp[4] = 0.5 * rm * rp * (-s - 4.0 * rr) * sm;
        interp[5] = 0.5 * sm * sp * (r - 4.0 * ss) * rp;
        interp[6] = 0.5 * rm * rp * (s - 4.0 * rr) * sp;
        interp[7] = 0.5 * sm * sp * (-r - 4.0 * ss) * rm;

        interp[8] = rm * rp * sm * sp;

        interp[9] = -a * r * sm * rm * rp * (1.0 - 2.0 * r);
        interp[10] = a * r * sm * rm * rp * (1.0 + 2.0 * r);
        interp[11] = -a * s * rp * sm * sp * (1.0 - 2.0 * s);
        interp[12] = a * s * rp * sm * sp * (1.0 + 2.0 * s);
        interp[13] = a * r * sp * rm * rp * (1.0 + 2.0 * r);
        interp[14] = -a * r * sp * rm * rp * (1.0 - 2.0 * r);
        interp[15] = a * s * sm * sp * (1.0 + 2.0 * s) * rm;
        interp[16] = -a * s * rm * sm * sp * (1.0 - 2.0 * s);
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

        let a = 2.0 / 3.0;
        let rr = r * r;
        let ss = s * s;
        let rs = r * s;
        let rp = 1.0 + r;
        let rm = 1.0 - r;
        let sp = 1.0 + s;
        let sm = 1.0 - s;

        let b = 1.0 / 12.0;
        let r1 = r - 1.0;
        let rrr = rr * r;
        let sss = ss * s;

        deriv.set(
            0,
            0,
            b * sm * (16.0 * rrr - 12.0 * rr - 6.0 * rs - 8.0 * r + 4.0 * sss - s + 4.0),
        );
        deriv.set(
            1,
            0,
            b * sm * (16.0 * rrr + 12.0 * rr - 6.0 * rs - 8.0 * r - 4.0 * sss + s - 4.0),
        );
        deriv.set(
            2,
            0,
            b * sp * (16.0 * rrr + 12.0 * rr + 6.0 * rs - 8.0 * r + 4.0 * sss - s - 4.0),
        );
        deriv.set(
            3,
            0,
            b * sp * (16.0 * rrr - 12.0 * rr + 6.0 * rs - 8.0 * r - 4.0 * sss + s + 4.0),
        );

        deriv.set(9, 0, -a * (1.0 - 4.0 * r - 3.0 * rr + 8.0 * rrr) * sm);
        deriv.set(11, 0, a * s * sm * sp * (-1.0 + 2.0 * s));
        deriv.set(13, 0, -a * (-1.0 - 4.0 * r + 3.0 * rr + 8.0 * rrr) * sp);
        deriv.set(15, 0, -a * s * sm * sp * (1.0 + 2.0 * s));
        deriv.set(4, 0, r * sm * (8.0 * rr + s - 4.0));
        deriv.set(5, 0, 0.5 * sm * sp * (2.0 * r - 4.0 * ss + 1.0));
        deriv.set(6, 0, r * sp * (8.0 * rr - s - 4.0));
        deriv.set(7, 0, 0.5 * sm * sp * (2.0 * r - 1.0 + 4.0 * ss));
        deriv.set(10, 0, a * (1.0 + 4.0 * r - 3.0 * rr - 8.0 * rrr) * sm);
        deriv.set(12, 0, a * s * sm * sp * (1.0 + 2.0 * s));
        deriv.set(14, 0, -a * (1.0 - 4.0 * r - 3.0 * rr + 8.0 * rrr) * sp);
        deriv.set(16, 0, a * s * sm * sp * (1.0 - 2.0 * s));
        deriv.set(8, 0, -2.0 * r * sm * sp);

        deriv.set(
            0,
            1,
            b * rm * (16.0 * sss - 12.0 * ss - 6.0 * rs - 8.0 * s + 4.0 * rrr - r + 4.0),
        );
        deriv.set(
            1,
            1,
            -b * rp * (-16.0 * sss + 12.0 * ss - 6.0 * rs + 8.0 * s + 4.0 * rrr - r - 4.0),
        );
        deriv.set(
            2,
            1,
            b * rp * (16.0 * sss + 12.0 * ss + 6.0 * rs - 8.0 * s + 4.0 * rrr - r - 4.0),
        );
        deriv.set(
            3,
            1,
            b * r1 * (-16.0 * sss - 12.0 * ss + 6.0 * rs + 8.0 * s + 4.0 * rrr - r + 4.0),
        );

        deriv.set(9, 1, a * r * r1 * rp * (2.0 * r - 1.0));
        deriv.set(11, 1, -a * (1.0 - 4.0 * s - 3.0 * ss + 8.0 * sss) * rp);
        deriv.set(13, 1, -a * r * r1 * rp * (1.0 + 2.0 * r));
        deriv.set(15, 1, a * (-1.0 - 4.0 * s + 3.0 * ss + 8.0 * sss) * r1);
        deriv.set(4, 1, -0.5 * r1 * rp * (2.0 * s - 1.0 + 4.0 * rr));
        deriv.set(5, 1, -s * rp * (-8.0 * ss + r + 4.0));
        deriv.set(6, 1, 0.5 * r1 * rp * (-2.0 * s + 4.0 * rr - 1.0));
        deriv.set(7, 1, -s * r1 * (8.0 * ss + r - 4.0));
        deriv.set(10, 1, a * r * r1 * rp * (1.0 + 2.0 * r));
        deriv.set(12, 1, -a * (-1.0 - 4.0 * s + 3.0 * ss + 8.0 * sss) * rp);
        deriv.set(14, 1, -a * r * r1 * rp * (2.0 * r - 1.0));
        deriv.set(16, 1, a * (1.0 - 4.0 * s - 3.0 * ss + 8.0 * sss) * r1);
        deriv.set(8, 1, 2.0 * s * r1 * rp);
    }
}
