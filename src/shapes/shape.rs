use super::*;
use russell_lab::Vector;

/// Defines the kind of shape
#[derive(Clone, Copy, Debug)]
pub enum Kind {
    Lin2,
    Lin3,
    Lin5,
    Tri3,
    Tri6,
    Tri15,
    Qua4,
    Qua8,
    Qua9,
    Tet4,
    Tet10,
    Hex8,
    Hex20,
    Joint,
    Lin4,
    Tri10,
    Qua12,
    Qua16,
    Beam,
}

/// Defines the functionality of shape
pub trait Shape {
    /// Calculates the interpolation functions at natural coordinate ξ
    ///
    /// ```text
    /// interp[m](ξ) = Sᵐ(ξ)
    /// ```
    fn calc_interp(&mut self, ksi: &Vector);

    /// Calculates the derivatives of interpolation fn at natural coordinate ξ
    ///
    /// ```text
    /// deriv[m][i](ξ) = ({dSᵐ(ξ)/dξ}_ξ)[i]
    /// ```
    fn calc_deriv(&mut self, ksi: &Vector);

    /// Returns the previously computed interpolation fn for point m
    ///
    /// ```text
    /// interp[m](ξ) = Sᵐ(ξ)
    /// ```
    fn get_interp(&self, m: usize) -> f64;

    /// Returns the previously computed derivative of interpolation fn for point m
    ///
    /// ```text
    /// deriv[m][i](ξ) = ({dSᵐ(ξ)/dξ}_ξ)[i]
    /// ```
    fn get_deriv(&self, deriv: &mut Vector, m: usize);

    /// Returns the number of dimensions
    fn get_ndim(&self) -> usize;

    /// Returns the number of points
    fn get_npoint(&self) -> usize;

    /// Returns the number of edges
    fn get_nedge(&self) -> usize;

    /// Returns the number of faces
    fn get_nface(&self) -> usize;

    /// Returns natural coordinates @ point m
    ///
    /// ```text
    /// coord = ξᵐ = vector{r, s, t} @ point m
    /// ```
    fn get_ksi(&self, ksi: &mut Vector, m: usize);
}

/// Returns new Shape
pub fn new(kind: Kind) -> Box<dyn Shape> {
    match kind {
        Kind::Qua4 => Box::new(Qua4::new()),
        Kind::Hex8 => Box::new(Hex8::new()),
        _ => panic!("Shape kind {:?} is not available yet", kind),
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use russell_chk::*;

    // Holds arguments for numerical differentiation
    struct Arguments {
        shape: Box<dyn Shape>, // shape
        at_ksi: Vector,        // at nat coord value
        ksi: Vector,           // temporary nat coord
        m: usize,              // point index
        i: usize,              // dimension index
    }

    // Computes the interpolation fn @ m as a function of x = ksi[i]
    fn sm(x: f64, args: &mut Arguments) -> f64 {
        for i in 0..args.shape.get_ndim() {
            args.ksi[i] = args.at_ksi[i];
        }
        args.ksi[args.i] = x;
        args.shape.calc_interp(&args.ksi);
        args.shape.get_interp(args.m)
    }

    #[test]
    fn calc_interp_works() {
        let kinds = vec![Kind::Hex8, Kind::Qua4];
        for kind in kinds {
            let mut shape = new(kind);
            let ndim = shape.get_ndim();
            let npoint = shape.get_npoint();
            let mut ksi = Vector::new(ndim);
            for m in 0..npoint {
                shape.get_ksi(&mut ksi, m);
                shape.calc_interp(&ksi);
                for n in 0..npoint {
                    let smn = shape.get_interp(n);
                    if m == n {
                        assert_approx_eq!(smn, 1.0, 1e-15);
                    } else {
                        assert_approx_eq!(smn, 0.0, 1e-15);
                    }
                }
            }
        }
    }

    #[test]
    fn calc_deriv_works() {
        let kinds = vec![Kind::Hex8, Kind::Qua4];
        for kind in kinds {
            let mut shape = new(kind);
            let ndim = shape.get_ndim();
            let npoint = shape.get_npoint();
            let mut at_ksi = Vector::new(ndim);
            for i in 0..ndim {
                at_ksi[i] = 0.25;
            }
            let args = &mut Arguments {
                shape: new(kind),
                at_ksi,
                ksi: Vector::new(ndim),
                m: 0,
                i: 0,
            };
            let mut dsm_dksi = Vector::new(ndim);
            shape.calc_deriv(&args.at_ksi);
            for m in 0..npoint {
                shape.get_deriv(&mut dsm_dksi, m);
                args.m = m;
                for i in 0..ndim {
                    args.i = i;
                    assert_deriv_approx_eq!(dsm_dksi[i], args.at_ksi[i], sm, args, 1e-13);
                }
            }
        }
    }

    #[test]
    fn getters_work() {
        let kinds = vec![Kind::Hex8, Kind::Qua4];
        for kind in kinds {
            let shape = new(kind);
            match kind {
                Kind::Qua4 => {
                    assert_eq!(shape.get_ndim(), 2);
                    assert_eq!(shape.get_npoint(), 4);
                    assert_eq!(shape.get_nedge(), 4);
                    assert_eq!(shape.get_nface(), 0);
                }
                Kind::Hex8 => {
                    assert_eq!(shape.get_ndim(), 3);
                    assert_eq!(shape.get_npoint(), 8);
                    assert_eq!(shape.get_nedge(), 12);
                    assert_eq!(shape.get_nface(), 6);
                }
                _ => panic!("Shape kind \"{:?}\" is not available yet", kind),
            }
        }
    }
}
