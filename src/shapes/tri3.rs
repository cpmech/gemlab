use russell_lab::{Matrix, Vector};

/// Defines a triangle with 3 points (linear edges)
///
/// # Local IDs of points
///
/// ```text
/// s
/// |
/// 2, (0,1)
/// | ',
/// |   ',
/// |     ',
/// |       ',
/// |         ',
/// |           ',
/// |             ',
/// |               ',
/// | (0,0)           ', (1,0)
/// 0-------------------1 ---- r
/// ```
///
/// # Local IDs of edges
///
/// ```text
///  |\
///  | \
///  |  \ 1
/// 2|   \
///  |    \
///  |_____\
///     0
/// ```
pub struct Tri3 {}

impl Tri3 {
    pub const NDIM: usize = 2;
    pub const NPOINT: usize = 3;
    pub const NEDGE: usize = 3;
    pub const NFACE: usize = 0;
    pub const EDGE_NPOINT: usize = 2;
    pub const FACE_NPOINT: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const EDGE_POINT_IDS: [[usize; Tri3::EDGE_NPOINT]; Tri3::NEDGE] = [
        [0, 1],
        [1, 2],
        [2, 0],
    ];

    #[rustfmt::skip]
    pub const POINT_REFERENCE_COORDS: [[f64; Tri3::NDIM]; Tri3::NPOINT] = [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        interp[0] = 1.0 - r - s;
        interp[1] = r;
        interp[2] = s;
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, _: &[f64]) {
        deriv[0][0] = -1.0;
        deriv[1][0] = 1.0;
        deriv[2][0] = 0.0;

        deriv[0][1] = -1.0;
        deriv[1][1] = 0.0;
        deriv[2][1] = 1.0;
    }
}
