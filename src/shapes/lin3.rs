use russell_lab::{Matrix, Vector};

/// Defines a line (segment) with 3 nodes (quadratic functions)
///
/// ![lin3](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_lin3.svg)
///
/// # Local IDs of nodes
///
/// ```text
///   0-----------2-----------1  --> r
/// -1.0         0.0         1.0
/// ```
///
/// # Notes
///
/// * The reference coordinates range from -1 to +1 with the geometry centred @ 0
/// * This shape is a lower-order version of [super::Lin5]
/// * This shape is a higher-order version of [super::Lin2]
pub struct Lin3 {}

impl Lin3 {
    pub const GEO_NDIM: usize = 1;
    pub const NNODE: usize = 3;
    pub const NEDGE: usize = 0;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 0;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;
    pub const N_INTERIOR_NODE: usize = 1;

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Lin3::GEO_NDIM]; Lin3::NNODE] = [
        [-1.0],
        [ 1.0],
        [ 0.0],
    ];

    pub const INTERIOR_NODES: [usize; Lin3::N_INTERIOR_NODE] = [2];

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
        let r = ksi[0];

        interp[0] = 0.5 * (r * r - r);
        interp[1] = 0.5 * (r * r + r);
        interp[2] = 1.0 - r * r;
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
        let r = ksi[0];

        deriv[0][0] = r - 0.5;
        deriv[1][0] = r + 0.5;
        deriv[2][0] = -2.0 * r;
    }
}
