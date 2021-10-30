use super::*;
use russell_lab::{Matrix, Vector};

/// Defines the kind of shape
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Kind {
    Lin2,
    Lin3,
    Lin4,
    Lin5,
    Tri3,
    Tri6,
    Tri10,
    Tri15,
    Qua4,
    Qua8,
    Qua9,
    Qua12,
    Qua16,
    Tet4,
    Tet10,
    Hex8,
    Hex20,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum KindQua {
    Qua4,
    Qua8,
    Qua9,
    Qua12,
    Qua16,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum KindHex {
    Hex8,
    Hex20,
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

    /// Returns the number of points on edge
    fn get_edge_npoint(&self) -> usize;

    /// Returns the number of points on face
    fn get_face_npoint(&self) -> usize;

    /// Returns the local ids of vertices on selected edge
    ///
    /// # Input
    ///
    /// * `e` -- the index of edge in [0, nedge-1]
    /// * `local_vertex_ids` -- ids of vertices on edge. len = edge_npoint
    fn get_edge(&self, local_vertex_ids: &mut Vec<usize>, e: usize);

    /// Returns the local ids of vertices on selected face
    ///
    /// # Input
    ///
    /// * `f` -- index of face in [0, nface-1]
    /// * `local_vertex_ids` -- ids of vertices on face. len = face_npoint
    fn get_face(&self, local_vertex_ids: &mut Vec<usize>, f: usize);

    /// Returns natural coordinates @ point m
    ///
    /// ```text
    /// coord = ξᵐ = vector{r, s, t} @ point m
    /// ```
    fn get_ksi(&self, ksi: &mut Vector, m: usize);

    /// Multiplies the interpolation vector (S) by a matrix (a)
    ///
    /// ```text
    /// v[n] = Σ_m Sᵐ ⋅ a[m][n]
    /// ```
    ///
    /// or
    ///
    /// ```text
    ///  v  =  S  ⋅   a
    /// (n)   (m)   (m,n)
    /// ```
    ///
    /// or
    ///
    /// ```text
    ///  v  =   aᵀ  ⋅  S
    /// (n)   (n,m)   (m)   
    /// ```
    ///
    /// # Note
    ///
    /// The interpolation vector must be computed first by calling `calc_interp`.
    fn mul_interp_by_matrix(&self, v: &mut Vector, a: &Matrix) -> Result<(), &'static str>;
}

/// Returns new Shape
pub fn new(kind: Kind) -> Box<dyn Shape> {
    match kind {
        Kind::Qua4 => Box::new(Qua4::new()),
        Kind::Qua8 => Box::new(Qua8::new()),
        Kind::Hex8 => Box::new(Hex8::new()),
        Kind::Hex20 => Box::new(Hex20::new()),
        _ => panic!("Shape kind {:?} is not available yet", kind),
    }
}

/// Returns new Qua Shape
pub fn new_qua(kind: KindQua) -> Box<dyn Shape> {
    match kind {
        KindQua::Qua4 => Box::new(Qua4::new()),
        KindQua::Qua8 => Box::new(Qua8::new()),
        _ => panic!("Shape kind {:?} is not available yet", kind),
    }
}

/// Returns new Hex Shape
pub fn new_hex(kind: KindHex) -> Box<dyn Shape> {
    match kind {
        KindHex::Hex8 => Box::new(Hex8::new()),
        KindHex::Hex20 => Box::new(Hex20::new()),
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use russell_chk::*;

    fn gen_kinds() -> Vec<Kind> {
        vec![Kind::Qua4, Kind::Qua8, Kind::Hex8, Kind::Hex20]
    }

    #[test]
    fn kind_enums_work() {
        let lin = Kind::Lin2;
        let clone = lin.clone();
        assert_eq!(lin, clone);
        assert_eq!(format!("{:?}", lin), "Lin2");

        let qua = KindQua::Qua12;
        let clone = qua.clone();
        assert_eq!(qua, clone);
        assert_eq!(format!("{:?}", qua), "Qua12");

        let hex = KindHex::Hex20;
        let clone = hex.clone();
        assert_eq!(hex, clone);
        assert_eq!(format!("{:?}", hex), "Hex20");
    }

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
        let kinds = gen_kinds();
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
        let kinds = gen_kinds();
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
                    assert_deriv_approx_eq!(dsm_dksi[i], args.at_ksi[i], sm, args, 1e-12);
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

    #[test]
    fn mul_interp_by_matrix_works() -> Result<(), &'static str> {
        // iso-parametric elements: xᵢ = Σ_m Sᵐ ⋅ cᵐᵢ
        // where c is a matrix of coordinates
        let kinds = gen_kinds();
        for kind in kinds {
            let mut shape = new(kind);
            let ndim = shape.get_ndim();
            let npoint = shape.get_npoint();
            let mut ksi = Vector::new(ndim);
            // create matrix of coordinates with shifted nat-coords: edges = [0, 2]
            let mut c = Matrix::new(npoint, ndim);
            for m in 0..npoint {
                shape.get_ksi(&mut ksi, m);
                for i in 0..ndim {
                    c[m][i] = 1.0 + ksi[i];
                }
            }
            // set @ksi such that x[i] = 1.25
            let mut at_ksi = Vector::new(ndim);
            for i in 0..ndim {
                at_ksi[i] = 0.25;
            }
            shape.calc_interp(&at_ksi);
            // multiply S by C matrix
            let mut x = Vector::new(ndim);
            shape.mul_interp_by_matrix(&mut x, &c)?;
            let x_correct = vec![1.25; ndim];
            assert_vec_approx_eq!(x.as_data(), x_correct, 1e-15);
        }
        Ok(())
    }
}
