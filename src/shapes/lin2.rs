use russell_lab::{Matrix, Vector};

/// Defines a line (segment) with 2 nodes (linear functions)
///
/// ![lin2](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_lin2.svg)
///
/// # Local IDs of nodes
///
/// ```text
///   0-----------------------1  --> r
/// -1.0                     1.0
/// ```
///
/// # Notes
///
/// * The reference coordinates range from -1 to +1 with the geometry centred @ 0
/// * This shape is a lower-order version of [super::Lin3]
pub struct Lin2 {}

impl Lin2 {
    pub const GEO_NDIM: usize = 1;
    pub const NNODE: usize = 2;
    pub const NEDGE: usize = 0;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 0;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Lin2::GEO_NDIM]; Lin2::NNODE] = [
        [-1.0],
        [ 1.0],
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
        let r = ksi[0];

        interp[0] = 0.5 * (1.0 - r);
        interp[1] = 0.5 * (1.0 + r);
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
        deriv.set(0, 0, -0.5);
        deriv.set(1, 0, 0.5);
    }
}
