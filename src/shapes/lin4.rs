use russell_lab::{Matrix, Vector};

/// Defines a line (segment) with 4 nodes (cubic functions)
///
/// ![lin4](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/test_draw_shape_lin4.svg)
///
/// # Local IDs of nodes
///
/// ```text
///   0-------2-------3-------1  --> r
/// -1.0    -1/3     1/3     1.0
/// ```
///
/// # Notes
///
/// * The reference coordinates range from -1 to +1 with the geometry centred @ 0
pub struct Lin4 {}

impl Lin4 {
    pub const GEO_NDIM: usize = 1;
    pub const NNODE: usize = 4;
    pub const NEDGE: usize = 0;
    pub const NFACE: usize = 0;
    pub const EDGE_NNODE: usize = 0;
    pub const FACE_NNODE: usize = 0;
    pub const FACE_NEDGE: usize = 0;
    pub const N_INTERIOR_NODE: usize = 2;

    #[rustfmt::skip]
    pub const NODE_REFERENCE_COORDS: [[f64; Lin4::GEO_NDIM]; Lin4::NNODE] = [
        [-1.0],
        [ 1.0],
        [-1.0 / 3.0],
        [ 1.0 / 3.0],
    ];

    pub const INTERIOR_NODES: [usize; Lin4::N_INTERIOR_NODE] = [2, 3];

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

        interp[0] = (-9.0 * r * r * r + 9.0 * r * r + r - 1.0) / 16.0;
        interp[1] = (9.0 * r * r * r + 9.0 * r * r - r - 1.0) / 16.0;
        interp[2] = (27.0 * r * r * r - 9.0 * r * r - 27.0 * r + 9.0) / 16.0;
        interp[3] = (-27.0 * r * r * r - 9.0 * r * r + 27.0 * r + 9.0) / 16.0;
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

        deriv.set(0, 0, 1.0 / 16.0 * (-27.0 * r * r + 18.0 * r + 1.0));
        deriv.set(1, 0, 1.0 / 16.0 * (27.0 * r * r + 18.0 * r - 1.0));
        deriv.set(2, 0, 1.0 / 16.0 * (81.0 * r * r - 18.0 * r - 27.0));
        deriv.set(3, 0, 1.0 / 16.0 * (-81.0 * r * r - 18.0 * r + 27.0));
    }
}
