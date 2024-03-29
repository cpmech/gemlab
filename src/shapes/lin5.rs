use russell_lab::{Matrix, Vector};

/// Defines a line (segment) with 5 nodes (quartic functions)
///
/// ![lin5](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_lin5.svg)
///
/// # Local IDs of nodes
///
/// ```text
///   0-----3-----2-----4-----1  --> r
/// -1.0  -0.5   0.0   0.5   1.0
/// ```
///
/// # Notes
///
/// * The reference coordinates range from -1 to +1 with the geometry centred @ 0
/// * This shape is a higher-order version of [super::Lin3]
pub struct Lin5 {}

impl Lin5 {
    pub const GEO_NDIM: usize = 1;
    pub const NNODE: usize = 5;
    pub const NEDGE: usize = 0;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 0;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;
    pub const N_INTERIOR_NODE: usize = 3;

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Lin5::GEO_NDIM]; Lin5::NNODE] = [
        [-1.0],
        [ 1.0],
        [ 0.0],
        [-0.5],
        [ 0.5],
    ];

    pub const INTERIOR_NODES: [usize; Lin5::N_INTERIOR_NODE] = [2, 3, 4];

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

        interp[0] = (r - 1.0) * (1.0 - 2.0 * r) * r * (-1.0 - 2.0 * r) / 6.0;
        interp[1] = (1.0 - 2.0 * r) * r * (-1.0 - 2.0 * r) * (1.0 + r) / 6.0;
        interp[2] = (1.0 - r) * (1.0 - 2.0 * r) * (-1.0 - 2.0 * r) * (-1.0 - r);
        interp[3] = 4.0 * (1.0 - r) * (1.0 - 2.0 * r) * r * (-1.0 - r) / 3.0;
        interp[4] = 4.0 * (1.0 - r) * r * (-1.0 - 2.0 * r) * (-1.0 - r) / 3.0;
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

        deriv.set(
            0,
            0,
            -((1.0 - 2.0 * r) * (r - 1.0) * r) / 3.0 - ((-2.0 * r - 1.0) * (r - 1.0) * r) / 3.0
                + ((-2.0 * r - 1.0) * (1.0 - 2.0 * r) * r) / 6.0
                + ((-2.0 * r - 1.0) * (1.0 - 2.0 * r) * (r - 1.0)) / 6.0,
        );

        deriv.set(
            1,
            0,
            -((1.0 - 2.0 * r) * r * (r + 1.0)) / 3.0 - ((-2.0 * r - 1.0) * r * (r + 1.0)) / 3.0
                + ((-2.0 * r - 1.0) * (1.0 - 2.0 * r) * (r + 1.0)) / 6.0
                + ((-2.0 * r - 1.0) * (1.0 - 2.0 * r) * r) / 6.0,
        );

        deriv.set(
            2,
            0,
            -2.0 * (1.0 - 2.0 * r) * (-r - 1.0) * (1.0 - r)
                - 2.0 * (-2.0 * r - 1.0) * (-r - 1.0) * (1.0 - r)
                - (-2.0 * r - 1.0) * (1.0 - 2.0 * r) * (1.0 - r)
                - (-2.0 * r - 1.0) * (1.0 - 2.0 * r) * (-r - 1.0),
        );

        deriv.set(
            3,
            0,
            -(8.0 * (-r - 1.0) * (1.0 - r) * r) / 3.0
                - (4.0 * (1.0 - 2.0 * r) * (1.0 - r) * r) / 3.0
                - (4.0 * (1.0 - 2.0 * r) * (-r - 1.0) * r) / 3.0
                + (4.0 * (1.0 - 2.0 * r) * (-r - 1.0) * (1.0 - r)) / 3.0,
        );

        deriv.set(
            4,
            0,
            -(8.0 * (-r - 1.0) * (1.0 - r) * r) / 3.0
                - (4.0 * (-2.0 * r - 1.0) * (1.0 - r) * r) / 3.0
                - (4.0 * (-2.0 * r - 1.0) * (-r - 1.0) * r) / 3.0
                + (4.0 * (-2.0 * r - 1.0) * (-r - 1.0) * (1.0 - r)) / 3.0,
        );
    }
}
