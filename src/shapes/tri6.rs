use russell_lab::{Matrix, Vector};

/// Defines a triangle with 6 points (quadratic edges)
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
///  |\
///  | \
///  |  \ 1
/// 2|   \
///  |    \
///  |_____\
///     0
/// ```
pub struct Tri6 {}

impl Tri6 {
    pub const NDIM: usize = 2;
    pub const NPOINT: usize = 6;
    pub const NEDGE: usize = 3;
    pub const NFACE: usize = 0;
    pub const EDGE_NPOINT: usize = 3;
    pub const FACE_NPOINT: usize = 0;
    pub const FACE_NEDGE: usize = 0;

    #[rustfmt::skip]
    pub const EDGE_POINT_IDS: [[usize; Tri6::EDGE_NPOINT]; Tri6::NEDGE] = [
        [0, 1, 3],
        [1, 2, 4],
        [2, 0, 5],
    ];

    #[rustfmt::skip]
    pub const POINT_REFERENCE_COORDS: [[f64; Tri6::NDIM]; Tri6::NPOINT] = [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.5, 0.0],
        [0.5, 0.5],
        [0.0, 0.5],
    ];

    /// Computes the interpolation functions
    pub fn calc_interp(interp: &mut Vector, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        interp[0] = 1.0 - (r + s) * (3.0 - 2.0 * (r + s));
        interp[1] = r * (2.0 * r - 1.0);
        interp[2] = s * (2.0 * s - 1.0);
        interp[3] = 4.0 * r * (1.0 - (r + s));
        interp[4] = 4.0 * r * s;
        interp[5] = 4.0 * s * (1.0 - (r + s));
    }

    /// Computes the derivatives of interpolation functions
    pub fn calc_deriv(deriv: &mut Matrix, ksi: &[f64]) {
        let (r, s) = (ksi[0], ksi[1]);

        deriv[0][0] = -3.0 + 4.0 * (r + s);
        deriv[1][0] = 4.0 * r - 1.0;
        deriv[2][0] = 0.0;
        deriv[3][0] = 4.0 - 8.0 * r - 4.0 * s;
        deriv[4][0] = 4.0 * s;
        deriv[5][0] = -4.0 * s;

        deriv[0][1] = -3.0 + 4.0 * (r + s);
        deriv[1][1] = 0.0;
        deriv[2][1] = 4.0 * s - 1.0;
        deriv[3][1] = -4.0 * r;
        deriv[4][1] = 4.0 * r;
        deriv[5][1] = 4.0 - 4.0 * r - 8.0 * s;
    }
}
