use russell_lab::{Matrix, Vector};

/// Defines a triangle with 6 nodes (quadratic edges)
///
/// ![tri6](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_tri6.svg)
///
/// # Local IDs of nodes
///
/// ```text
/// s
/// |
/// 2, (0,1)
/// | ',
/// |   ',
/// |     ',
/// |       ',
/// 5         4,
/// |           ',
/// |             ',
/// |               ',
/// | (0,0)           ', (1,0)
/// 0---------3---------1 ---- r
/// ```
///
/// # Local IDs of edges
///
/// ```text
///   2,
///   | ',                 p0  p1 p2
///   |   ',           e:0 [1, 0, 3]
///   |     ',         e:1 [2, 1, 4]
///   |       ', 1     e:2 [0, 2, 5]
/// 2 5         4,
///   |           ',
///   |             ',
///   |               ',
///   |                 ',
///   0---------3---------1
///             0
/// ```
///
/// # Notes
///
/// * The order of edge nodes is such that the normals are outward
/// * The order of edge nodes corresponds to [super::Lin3] nodes
/// * This shape is a lower-order version of [super::Tri15]
/// * This shape is a higher-order version of [super::Tri3]
pub struct Tri6 {}

impl Tri6 {
    pub const GEO_NDIM: usize = 2;
    pub const NNODE: usize = 6;
    pub const NEDGE: usize = 3;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 3;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const EDGE_NODE_IDS: [[usize; Tri6::EDGE_NNODE]; Tri6::NEDGE] = [
        [1, 0, 3],
        [2, 1, 4],
        [0, 2, 5],
    ];

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Tri6::GEO_NDIM]; Tri6::NNODE] = [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.5, 0.0],
        [0.5, 0.5],
        [0.0, 0.5],
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
        let (r, s) = (ksi[0], ksi[1]);

        interp[0] = 1.0 - (r + s) * (3.0 - 2.0 * (r + s));
        interp[1] = r * (2.0 * r - 1.0);
        interp[2] = s * (2.0 * s - 1.0);
        interp[3] = 4.0 * r * (1.0 - (r + s));
        interp[4] = 4.0 * r * s;
        interp[5] = 4.0 * s * (1.0 - (r + s));
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

        deriv.set(0, 0, -3.0 + 4.0 * (r + s));
        deriv.set(1, 0, 4.0 * r - 1.0);
        deriv.set(2, 0, 0.0);
        deriv.set(3, 0, 4.0 - 8.0 * r - 4.0 * s);
        deriv.set(4, 0, 4.0 * s);
        deriv.set(5, 0, -4.0 * s);

        deriv.set(0, 1, -3.0 + 4.0 * (r + s));
        deriv.set(1, 1, 0.0);
        deriv.set(2, 1, 4.0 * s - 1.0);
        deriv.set(3, 1, -4.0 * r);
        deriv.set(4, 1, 4.0 * r);
        deriv.set(5, 1, 4.0 - 4.0 * r - 8.0 * s);
    }
}
