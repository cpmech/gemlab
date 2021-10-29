use super::Shape;
use russell_lab::{Matrix, Vector};

const NDIM: usize = 3;
const NPOINT: usize = 8;
const NEDGE: usize = 12;
const NFACE: usize = 6;

/// Implements a hexahedron with 8 points
///
/// The natural coordinates range from -1 to +1 with the geometry centred @ 0
///
/// ```text
///           4________________7
///         ,'|              ,'|
///       ,'  |            ,'  |
///     ,'    |          ,'    |
///   ,'      |        ,'      |
/// 5'===============6'        |
/// |         |      |         |
/// |         |      |         |
/// |         0_____ | ________3
/// |       ,'       |       ,'
/// |     ,'         |     ,'
/// |   ,'           |   ,'
/// | ,'             | ,'
/// 1________________2'
/// ```
pub struct Hex8 {
    coords: Vec<Vector>, // natural coordinates (npoint, ndim)
    interp: Vector,      // interpolation functions @ natural coordinate (npoint)
    deriv: Matrix,       // derivatives of interpolation functions w.r.t natural coordinate (npoint, ndim)
}

impl Hex8 {
    /// Creates a new object with pre-calculated interpolation fn and derivatives
    pub fn new() -> Self {
        Hex8 {
            #[rustfmt::skip]
            coords: vec![
                Vector::from(&[-1.0, -1.0, -1.0]),
                Vector::from(&[ 1.0, -1.0, -1.0]),
                Vector::from(&[ 1.0,  1.0, -1.0]),
                Vector::from(&[-1.0,  1.0, -1.0]),
                Vector::from(&[-1.0, -1.0,  1.0]),
                Vector::from(&[ 1.0, -1.0,  1.0]),
                Vector::from(&[ 1.0,  1.0,  1.0]),
                Vector::from(&[-1.0,  1.0,  1.0]),
            ],
            interp: Vector::new(NPOINT),
            deriv: Matrix::new(NPOINT, NDIM),
        }
    }
}

impl Shape for Hex8 {
    fn calc_interp(&mut self, ksi: &Vector) {
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        self.interp[0] = (1.0 - r - s + r * s - t + s * t + r * t - r * s * t) / 8.0;
        self.interp[1] = (1.0 + r - s - r * s - t + s * t - r * t + r * s * t) / 8.0;
        self.interp[2] = (1.0 + r + s + r * s - t - s * t - r * t - r * s * t) / 8.0;
        self.interp[3] = (1.0 - r + s - r * s - t - s * t + r * t + r * s * t) / 8.0;
        self.interp[4] = (1.0 - r - s + r * s + t - s * t - r * t + r * s * t) / 8.0;
        self.interp[5] = (1.0 + r - s - r * s + t - s * t + r * t - r * s * t) / 8.0;
        self.interp[6] = (1.0 + r + s + r * s + t + s * t + r * t + r * s * t) / 8.0;
        self.interp[7] = (1.0 - r + s - r * s + t + s * t - r * t - r * s * t) / 8.0;
    }

    fn calc_deriv(&mut self, ksi: &Vector) {
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        self.deriv[0][0] = (-1.0 + s + t - s * t) / 8.0;
        self.deriv[0][1] = (-1.0 + r + t - r * t) / 8.0;
        self.deriv[0][2] = (-1.0 + r + s - r * s) / 8.0;

        self.deriv[1][0] = (1.0 - s - t + s * t) / 8.0;
        self.deriv[1][1] = (-1.0 - r + t + r * t) / 8.0;
        self.deriv[1][2] = (-1.0 - r + s + r * s) / 8.0;

        self.deriv[2][0] = (1.0 + s - t - s * t) / 8.0;
        self.deriv[2][1] = (1.0 + r - t - r * t) / 8.0;
        self.deriv[2][2] = (-1.0 - r - s - r * s) / 8.0;

        self.deriv[3][0] = (-1.0 - s + t + s * t) / 8.0;
        self.deriv[3][1] = (1.0 - r - t + r * t) / 8.0;
        self.deriv[3][2] = (-1.0 + r - s + r * s) / 8.0;

        self.deriv[4][0] = (-1.0 + s - t + s * t) / 8.0;
        self.deriv[4][1] = (-1.0 + r - t + r * t) / 8.0;
        self.deriv[4][2] = (1.0 - r - s + r * s) / 8.0;

        self.deriv[5][0] = (1.0 - s + t - s * t) / 8.0;
        self.deriv[5][1] = (-1.0 - r - t - r * t) / 8.0;
        self.deriv[5][2] = (1.0 + r - s - r * s) / 8.0;

        self.deriv[6][0] = (1.0 + s + t + s * t) / 8.0;
        self.deriv[6][1] = (1.0 + r + t + r * t) / 8.0;
        self.deriv[6][2] = (1.0 + r + s + r * s) / 8.0;

        self.deriv[7][0] = (-1.0 - s - t - s * t) / 8.0;
        self.deriv[7][1] = (1.0 - r + t - r * t) / 8.0;
        self.deriv[7][2] = (1.0 - r + s - r * s) / 8.0;
    }

    fn get_interp(&self, m: usize) -> f64 {
        self.interp[m]
    }

    fn get_deriv(&self, deriv: &mut Vector, m: usize) {
        deriv[0] = self.deriv[m][0];
        deriv[1] = self.deriv[m][1];
        deriv[2] = self.deriv[m][2];
    }

    fn get_ndim(&self) -> usize {
        NDIM
    }

    fn get_npoint(&self) -> usize {
        NPOINT
    }

    fn get_nedge(&self) -> usize {
        NEDGE
    }

    fn get_nface(&self) -> usize {
        NFACE
    }

    fn get_ksi(&self, ksi: &mut Vector, m: usize) {
        ksi[0] = self.coords[m][0];
        ksi[1] = self.coords[m][1];
        ksi[2] = self.coords[m][2];
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_works() {
        let geo = Hex8::new();
        assert_eq!(geo.coords.len(), NPOINT);
        assert_eq!(geo.interp.dim(), NPOINT);
        assert_eq!(geo.deriv.dims(), (NPOINT, NDIM));
    }
}
