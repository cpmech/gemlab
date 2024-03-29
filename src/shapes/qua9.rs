use russell_lab::{Matrix, Vector};

/// Defines a quadrilateral with 9 nodes (quadratic edges; interior node)
///
/// ![qua9](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_qua9.svg)
///
/// # Local IDs of nodes
///
/// ```text
/// 3-----6-----2
/// |     s     |           r     s            r     s
/// |     |     |   p:0 [-1.0, -1.0]   p:4 [ 0.0, -1.0]
/// 7     8--r  5   p:1 [ 1.0, -1.0]   p:5 [ 1.0,  0.0]
/// |           |   p:2 [ 1.0,  1.0]   p:6 [ 0.0,  1.0]
/// |           |   p:3 [-1.0,  1.0]   p:7 [-1.0,  0.0]
/// 0-----4-----1   p:8 [ 0.0,  0.0]
/// ```
///
/// # Local IDs of edges
///
/// ```text
///         2
///   3-----6-----2          p0 p1 p2
///   |           |      e:0 [1, 0, 4]
///   |           |      e:1 [2, 1, 5]
/// 3 7     8     5 1    e:2 [3, 2, 6]
///   |           |      e:3 [0, 3, 7]
///   |           |
///   0-----4-----1
///         0
/// ```
///
/// # Notes
///
/// * The reference coordinates range from -1 to +1 with the geometry centred @ 0
/// * The order of edge nodes is such that the normals are outward
/// * The order of edge nodes corresponds to [super::Lin3] nodes
/// * This shape is a lower-order version of [super::Qua17]
/// * This shape is a higher-order version of [super::Qua4]
pub struct Qua9 {}

impl Qua9 {
    pub const GEO_NDIM: usize = 2;
    pub const NNODE: usize = 9;
    pub const NEDGE: usize = 4;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 3;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;
    pub const N_INTERIOR_NODE: usize = 1;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Qua9::EDGE_NNODE]; Qua9::NEDGE] = [
        [1, 0, 4],
        [2, 1, 5],
        [3, 2, 6],
        [0, 3, 7],
    ];

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS_INWARD: [[usize; Qua9::EDGE_NNODE]; Qua9::NEDGE] = [
        [0, 1, 4],
        [1, 2, 5],
        [2, 3, 6],
        [3, 0, 7],
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Qua9::GEO_NDIM]; Qua9::NNODE] = [
        [-1.0, -1.0],
        [ 1.0, -1.0],
        [ 1.0,  1.0],
        [-1.0,  1.0],
        [ 0.0, -1.0],
        [ 1.0,  0.0],
        [ 0.0,  1.0],
        [-1.0,  0.0],
        [ 0.0,  0.0],
    ];

    pub const INTERIOR_NODES: [usize; Qua9::N_INTERIOR_NODE] = [8];

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

        interp[0] = r * (r - 1.0) * s * (s - 1.0) / 4.0;
        interp[1] = r * (r + 1.0) * s * (s - 1.0) / 4.0;
        interp[2] = r * (r + 1.0) * s * (s + 1.0) / 4.0;
        interp[3] = r * (r - 1.0) * s * (s + 1.0) / 4.0;

        interp[4] = -(r * r - 1.0) * s * (s - 1.0) / 2.0;
        interp[5] = -r * (r + 1.0) * (s * s - 1.0) / 2.0;
        interp[6] = -(r * r - 1.0) * s * (s + 1.0) / 2.0;
        interp[7] = -r * (r - 1.0) * (s * s - 1.0) / 2.0;

        interp[8] = (r * r - 1.0) * (s * s - 1.0);
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

        deriv.set(0, 0, (r + r - 1.0) * s * (s - 1.0) / 4.0);
        deriv.set(1, 0, (r + r + 1.0) * s * (s - 1.0) / 4.0);
        deriv.set(2, 0, (r + r + 1.0) * s * (s + 1.0) / 4.0);
        deriv.set(3, 0, (r + r - 1.0) * s * (s + 1.0) / 4.0);

        deriv.set(0, 1, r * (r - 1.0) * (s + s - 1.0) / 4.0);
        deriv.set(1, 1, r * (r + 1.0) * (s + s - 1.0) / 4.0);
        deriv.set(2, 1, r * (r + 1.0) * (s + s + 1.0) / 4.0);
        deriv.set(3, 1, r * (r - 1.0) * (s + s + 1.0) / 4.0);

        deriv.set(4, 0, -(r + r) * s * (s - 1.0) / 2.0);
        deriv.set(5, 0, -(r + r + 1.0) * (s * s - 1.0) / 2.0);
        deriv.set(6, 0, -(r + r) * s * (s + 1.0) / 2.0);
        deriv.set(7, 0, -(r + r - 1.0) * (s * s - 1.0) / 2.0);

        deriv.set(4, 1, -(r * r - 1.0) * (s + s - 1.0) / 2.0);
        deriv.set(5, 1, -r * (r + 1.0) * (s + s) / 2.0);
        deriv.set(6, 1, -(r * r - 1.0) * (s + s + 1.0) / 2.0);
        deriv.set(7, 1, -r * (r - 1.0) * (s + s) / 2.0);

        deriv.set(8, 0, 2.0 * r * (s * s - 1.0));
        deriv.set(8, 1, 2.0 * s * (r * r - 1.0));
    }
}
