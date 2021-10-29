use super::Shape;
use russell_lab::{Matrix, Vector};

const NDIM: usize = 2;
const NPOINT: usize = 4;
const NEDGE: usize = 4;
const NFACE: usize = 0;

/// Implements a quadrilateral with 4 points
///
/// The natural coordinates range from -1 to +1 with the geometry centred @ 0
///
/// ```text
/// 3-----------2
/// |     s     |
/// |     |     |
/// |     +--r  |
/// |           |
/// |           |
/// 0-----------1
/// ```
pub struct Qua4 {
    coords: Vec<Vector>, // natural coordinates (npoint, ndim)
    interp: Vector,      // interpolation functions @ natural coordinate (npoint)
    deriv: Matrix,       // derivatives of interpolation functions w.r.t natural coordinate (npoint, ndim)
}

impl Qua4 {
    pub fn new() -> Self {
        Qua4 {
            #[rustfmt::skip]
            coords: vec![
                Vector::from(&[-1.0, -1.0]),
                Vector::from(&[ 1.0, -1.0]),
                Vector::from(&[ 1.0,  1.0]),
                Vector::from(&[-1.0,  1.0]),
            ],
            interp: Vector::new(NPOINT),
            deriv: Matrix::new(NPOINT, NDIM),
        }
    }
}

impl Shape for Qua4 {
    fn calc_interp(&mut self, ksi: &Vector) {
        let (r, s) = (ksi[0], ksi[1]);

        self.interp[0] = (1.0 - r - s + r * s) / 4.0;
        self.interp[1] = (1.0 + r - s - r * s) / 4.0;
        self.interp[2] = (1.0 + r + s + r * s) / 4.0;
        self.interp[3] = (1.0 - r + s - r * s) / 4.0;
    }

    fn calc_deriv(&mut self, ksi: &Vector) {
        let (r, s) = (ksi[0], ksi[1]);

        self.deriv[0][0] = (-1.0 + s) / 4.0;
        self.deriv[0][1] = (-1.0 + r) / 4.0;

        self.deriv[1][0] = (1.0 - s) / 4.0;
        self.deriv[1][1] = (-1.0 - r) / 4.0;

        self.deriv[2][0] = (1.0 + s) / 4.0;
        self.deriv[2][1] = (1.0 + r) / 4.0;

        self.deriv[3][0] = (-1.0 - s) / 4.0;
        self.deriv[3][1] = (1.0 - r) / 4.0;
    }

    fn get_interp(&self, m: usize) -> f64 {
        self.interp[m]
    }

    fn get_deriv(&self, deriv: &mut Vector, m: usize) {
        deriv[0] = self.deriv[m][0];
        deriv[1] = self.deriv[m][1];
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
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use russell_chk::*;

    // Holds arguments for numerical differentiation
    struct Arguments {
        shape: Qua4,    // shape
        at_ksi: Vector, // at nat coord value
        ksi: Vector,    // temporary nat coord
        m: usize,       // point index
        i: usize,       // dimension index
    }

    // Computes the interpolation fn @ m as a function of x = ksi[i]
    fn sm(x: f64, args: &mut Arguments) -> f64 {
        args.ksi[0] = args.at_ksi[0];
        args.ksi[1] = args.at_ksi[1];
        args.ksi[args.i] = x;
        args.shape.calc_interp(&args.ksi);
        args.shape.get_interp(args.m)
    }

    #[test]
    fn new_works() {
        let geo = Qua4::new();
        assert_eq!(geo.coords.len(), NPOINT);
        assert_eq!(geo.interp.dim(), NPOINT);
        assert_eq!(geo.deriv.dims(), (NPOINT, NDIM));
    }

    #[test]
    fn calc_interp_works() {
        let mut geo = Qua4::new();
        let mut ksi = Vector::new(NDIM);
        for m in 0..NPOINT {
            geo.get_ksi(&mut ksi, m);
            geo.calc_interp(&ksi);
            for n in 0..NPOINT {
                let smn = geo.get_interp(n);
                if m == n {
                    assert_approx_eq!(smn, 1.0, 1e-15);
                } else {
                    assert_approx_eq!(smn, 0.0, 1e-15);
                }
            }
        }
    }

    #[test]
    fn calc_deriv_works() {
        let args = &mut Arguments {
            shape: Qua4::new(),
            at_ksi: Vector::from(&[0.25, 0.25]),
            ksi: Vector::new(NDIM),
            m: 0,
            i: 0,
        };
        let mut geo = Qua4::new();
        let mut dsm_dksi = Vector::new(NDIM);
        geo.calc_deriv(&args.at_ksi);
        for m in 0..NPOINT {
            geo.get_deriv(&mut dsm_dksi, m);
            args.m = m;
            for i in 0..NDIM {
                args.i = i;
                assert_deriv_approx_eq!(dsm_dksi[i], args.at_ksi[i], sm, args, 1e-13);
            }
        }
    }

    #[test]
    fn getters_work() {
        let geo = Qua4::new();
        assert_eq!(geo.get_ndim(), NDIM);
        assert_eq!(geo.get_npoint(), NPOINT);
        assert_eq!(geo.get_nedge(), NEDGE);
        assert_eq!(geo.get_nface(), NFACE);
    }
}
