use super::{Hex20, Hex8, Qua4, Qua8};
use crate::StrError;
use russell_lab::{vec_mat_mul, Matrix, Vector};

/// Defines the kind of shape
#[allow(dead_code)]
enum ShapeKind {
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
    Qua17,
    Tet4,
    Tet10,
    Hex8,
    Hex20,
}

// alias for interpolation functions
type FnInterp = fn(&mut Vector, &Vector);

// alias for derivative of interpolation functions
type FnDeriv = fn(&mut Matrix, &Vector);

/// Defines a geometric shape for numerical integration and more
pub struct Shape {
    // shape kind
    kind: ShapeKind, // shape kind

    // constants (from each specific shape)
    ndim: usize,        // shape ndim (not space ndim)
    npoint: usize,      // number of points
    nedge: usize,       // number of edges
    nface: usize,       // number of faces
    edge_npoint: usize, // edge's number of points
    face_npoint: usize, // face's number of points
    face_nedge: usize,  // face's number of edges

    // functions (from each specific shape)
    fn_interp: FnInterp, // function to calculate interpolation functions
    fn_deriv: FnDeriv,   // function to calculate derivatives of interpolation functions

    // sandbox (temporary) variables
    interp: Vector, // interpolation functions @ natural coordinate (npoint)
    deriv: Matrix,  // derivatives of interpolation functions w.r.t natural coordinate (npoint, ndim)
    coords: Matrix, // (npoint, space_ndim) real coordinates matrix
}

impl Shape {
    /// Creates a new geometric shape
    pub fn new(space_ndim: usize, shape_ndim: usize, npoint: usize) -> Result<Self, StrError> {
        match (shape_ndim, npoint) {
            (1, 2) => return Err("Lin2 is not available yet"),
            (1, 3) => return Err("Lin3 is not available yet"),
            (1, 4) => return Err("Lin4 is not available yet"),
            (1, 5) => return Err("Lin5 is not available yet"),
            (2, 3) => return Err("Tri3 is not available yet"),
            (2, 6) => return Err("Tri6 is not available yet"),
            (2, 10) => return Err("Tri10 is not available yet"),
            (2, 15) => return Err("Tri15 is not available yet"),
            (2, 4) => {
                return Ok(Shape {
                    kind: ShapeKind::Qua4,
                    ndim: Qua4::NDIM,
                    npoint: Qua4::NPOINT,
                    nedge: Qua4::NEDGE,
                    nface: Qua4::NFACE,
                    edge_npoint: Qua4::EDGE_NPOINT,
                    face_npoint: Qua4::FACE_NPOINT,
                    face_nedge: Qua4::FACE_NEDGE,
                    fn_interp: Qua4::calc_interp,
                    fn_deriv: Qua4::calc_deriv,
                    interp: Vector::new(Qua4::NPOINT),
                    deriv: Matrix::new(Qua4::NPOINT, Qua4::NDIM),
                    coords: Matrix::new(Qua4::NPOINT, space_ndim),
                });
            }
            (2, 8) => {
                return Ok(Shape {
                    kind: ShapeKind::Qua8,
                    ndim: Qua8::NDIM,
                    npoint: Qua8::NPOINT,
                    nedge: Qua8::NEDGE,
                    nface: Qua8::NFACE,
                    edge_npoint: Qua8::EDGE_NPOINT,
                    face_npoint: Qua8::FACE_NPOINT,
                    face_nedge: Qua8::FACE_NEDGE,
                    fn_interp: Qua8::calc_interp,
                    fn_deriv: Qua8::calc_deriv,
                    interp: Vector::new(Qua8::NPOINT),
                    deriv: Matrix::new(Qua8::NPOINT, Qua8::NDIM),
                    coords: Matrix::new(Qua8::NPOINT, space_ndim),
                });
            }
            (2, 9) => return Err("Qua9 is not available yet"),
            (2, 12) => return Err("Qua12 is not available yet"),
            (2, 16) => return Err("Qua16 is not available yet"),
            (2, 17) => return Err("Qua17 is not available yet"),
            (3, 4) => return Err("Tet4 is not available yet"),
            (3, 10) => return Err("Tet10 is not available yet"),
            (3, 8) => {
                return Ok(Shape {
                    kind: ShapeKind::Hex8,
                    ndim: Hex8::NDIM,
                    npoint: Hex8::NPOINT,
                    nedge: Hex8::NEDGE,
                    nface: Hex8::NFACE,
                    edge_npoint: Hex8::EDGE_NPOINT,
                    face_npoint: Hex8::FACE_NPOINT,
                    face_nedge: Hex8::FACE_NEDGE,
                    fn_interp: Hex8::calc_interp,
                    fn_deriv: Hex8::calc_deriv,
                    interp: Vector::new(Hex8::NPOINT),
                    deriv: Matrix::new(Hex8::NPOINT, Hex8::NDIM),
                    coords: Matrix::new(Hex8::NPOINT, space_ndim),
                });
            }
            (3, 20) => {
                return Ok(Shape {
                    kind: ShapeKind::Hex20,
                    ndim: Hex20::NDIM,
                    npoint: Hex20::NPOINT,
                    nedge: Hex20::NEDGE,
                    nface: Hex20::NFACE,
                    edge_npoint: Hex20::EDGE_NPOINT,
                    face_npoint: Hex20::FACE_NPOINT,
                    face_nedge: Hex20::FACE_NEDGE,
                    fn_interp: Hex20::calc_interp,
                    fn_deriv: Hex20::calc_deriv,
                    interp: Vector::new(Hex20::NPOINT),
                    deriv: Matrix::new(Hex20::NPOINT, Hex20::NDIM),
                    coords: Matrix::new(Hex20::NPOINT, space_ndim),
                });
            }
            _ => return Err("(space_ndim, shape_ndim, npoint) combination is invalid"),
        };
    }

    /// Calculates the interpolation functions at natural coordinate ξ
    ///
    /// ```text
    /// interp[m](ξ) = Sᵐ(ξ)
    /// ```
    pub fn calc_interp(&mut self, ksi: &Vector) {
        (self.fn_interp)(&mut self.interp, ksi);
    }

    /// Calculates the derivatives of interpolation fn at natural coordinate ξ
    ///
    /// ```text
    /// deriv[m][i](ξ) = ({dSᵐ(ξ)/dξ}_ξ)[i]
    /// ```
    pub fn calc_deriv(&mut self, ksi: &Vector) {
        (self.fn_deriv)(&mut self.deriv, ksi);
    }

    /// Returns the previously computed interpolation fn for point m
    ///
    /// ```text
    /// interp[m](ξ) = Sᵐ(ξ)
    /// ```
    pub fn get_interp(&self, m: usize) -> f64 {
        self.interp[m]
    }

    /// Returns the previously computed derivative of interpolation fn for point m
    ///
    /// ```text
    /// deriv[m][i](ξ) = ({dSᵐ(ξ)/dξ}_ξ)[i]
    /// ```
    pub fn get_deriv(&self, deriv: &mut Vector, m: usize) {
        deriv[0] = self.deriv[m][0];
        deriv[1] = self.deriv[m][1];
        if self.ndim > 2 {
            deriv[2] = self.deriv[m][2];
        }
    }

    /// Returns the number of dimensions of the shape (not space)
    pub fn get_ndim(&self) -> usize {
        self.ndim
    }

    /// Returns the number of points
    pub fn get_npoint(&self) -> usize {
        self.npoint
    }

    /// Returns the number of edges
    pub fn get_nedge(&self) -> usize {
        self.nedge
    }

    /// Returns the number of faces
    pub fn get_nface(&self) -> usize {
        self.nface
    }

    /// Returns the number of points on edge
    pub fn get_edge_npoint(&self) -> usize {
        self.edge_npoint
    }

    /// Returns the number of points on face
    pub fn get_face_npoint(&self) -> usize {
        self.face_npoint
    }

    /// Returns the number of edges on face (not the cell's edges)
    pub fn get_face_nedge(&self) -> usize {
        self.face_nedge
    }

    /// Returns the local id of point on edge
    ///
    /// # Input
    ///
    /// * `e` -- index of edge in [0, nedge-1]
    /// * `i` -- index of local point [0, edge_npoint-1]
    pub fn get_edge_point_id(&self, e: usize, i: usize) -> usize {
        match self.kind {
            ShapeKind::Qua4 => Qua4::EDGE_POINT_IDS[e][i],
            ShapeKind::Qua8 => Qua8::EDGE_POINT_IDS[e][i],
            ShapeKind::Hex8 => Hex8::EDGE_POINT_IDS[e][i],
            ShapeKind::Hex20 => Hex20::EDGE_POINT_IDS[e][i],
            _ => panic!("ShapeKind is not available yet"),
        }
    }

    /// Returns the local id of point on face
    ///
    /// # Input
    ///
    /// * `f` -- index of face in [0, nface-1]
    /// * `i` -- index of local point [0, face_npoint-1]
    pub fn get_face_point_id(&self, f: usize, i: usize) -> usize {
        match self.kind {
            ShapeKind::Qua4 => 0,
            ShapeKind::Qua8 => 0,
            ShapeKind::Hex8 => Hex8::FACE_POINT_IDS[f][i],
            ShapeKind::Hex20 => Hex20::FACE_POINT_IDS[f][i],
            _ => panic!("ShapeKind is not available yet"),
        }
    }

    /// Returns the local point id on an edge on the face
    ///
    /// # Input
    ///
    /// * `f` -- index of face in [0, nface-1]
    /// * `k` -- index of face's edge (not the index of cell's edge) in [0, face_nedge-1]
    /// * `i` -- index of local point [0, edge_npoint-1]
    pub fn get_face_edge_point_id(&self, f: usize, k: usize, i: usize) -> usize {
        match self.kind {
            ShapeKind::Qua4 => 0,
            ShapeKind::Qua8 => 0,
            ShapeKind::Hex8 => Hex8::FACE_EDGE_POINT_IDS[f][k][i],
            ShapeKind::Hex20 => Hex20::FACE_EDGE_POINT_IDS[f][k][i],
            _ => panic!("ShapeKind is not available yet"),
        }
    }

    /// Returns natural coordinates @ point m
    ///
    /// ```text
    /// coord = ξᵐ = vector{r, s, t} @ point m
    /// ```
    pub fn get_ksi(&self, ksi: &mut Vector, m: usize) {
        match self.kind {
            ShapeKind::Qua4 => {
                ksi[0] = Qua4::POINT_NATURAL_COORDS[m][0];
                ksi[1] = Qua4::POINT_NATURAL_COORDS[m][1];
            }
            ShapeKind::Qua8 => {
                ksi[0] = Qua8::POINT_NATURAL_COORDS[m][0];
                ksi[1] = Qua8::POINT_NATURAL_COORDS[m][1];
            }
            ShapeKind::Hex8 => {
                ksi[0] = Hex8::POINT_NATURAL_COORDS[m][0];
                ksi[1] = Hex8::POINT_NATURAL_COORDS[m][1];
                ksi[2] = Hex8::POINT_NATURAL_COORDS[m][2];
            }
            ShapeKind::Hex20 => {
                ksi[0] = Hex20::POINT_NATURAL_COORDS[m][0];
                ksi[1] = Hex20::POINT_NATURAL_COORDS[m][1];
                ksi[2] = Hex20::POINT_NATURAL_COORDS[m][2];
            }
            _ => panic!("ShapeKind is not available yet"),
        }
    }

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
    pub fn mul_interp_by_matrix(&self, v: &mut Vector, a: &Matrix) -> Result<(), StrError> {
        vec_mat_mul(v, 1.0, &self.interp, a)
    }

    /// Sets the real coordinates i (0..space_ndim-1) of a point (0..npoint-1)
    pub fn set_coords(&mut self, m: usize, i: usize, val: f64) {
        self.coords[m][i] = val;
    }

    /// Returns a read-only access to the real coordinates matrix (npoint, space_ndim)
    pub fn get_coords_matrix(&self) -> &Matrix {
        &self.coords
    }

    /// Calculates the real coordinates of ksi
    ///
    /// ```text
    /// x[i] := Σ_m Sᵐ(ksi) ⋅ point_coords[m][i]
    /// ```
    pub fn calc_real_ksi_coords(&mut self, x: &mut Vector, ksi: &Vector) -> Result<(), StrError> {
        self.calc_interp(ksi);
        self.mul_interp_by_matrix(x, self.get_coords_matrix())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use russell_chk::*;

    fn gen_ndim_npoint() -> Vec<(usize, usize)> {
        vec![(2, 4), (2, 8), (3, 8), (3, 20)]
    }

    // Holds arguments for numerical differentiation
    struct Arguments {
        shape: Shape,   // shape
        at_ksi: Vector, // at nat coord value
        ksi: Vector,    // temporary nat coord
        m: usize,       // point index
        i: usize,       // dimension index
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
    fn calc_interp_works() -> Result<(), StrError> {
        let ndim_npoint = gen_ndim_npoint();
        for (ndim, npoint) in ndim_npoint {
            let mut shape = Shape::new(ndim, ndim, npoint)?;
            assert_eq!(shape.get_ndim(), ndim);
            assert_eq!(shape.get_npoint(), npoint);
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
        Ok(())
    }

    #[test]
    fn calc_deriv_works() -> Result<(), StrError> {
        let ndim_npoint = gen_ndim_npoint();
        for (ndim, npoint) in ndim_npoint {
            let mut shape = Shape::new(ndim, ndim, npoint)?;
            let mut at_ksi = Vector::new(ndim);
            for i in 0..ndim {
                at_ksi[i] = 0.25;
            }
            let args = &mut Arguments {
                shape: Shape::new(ndim, ndim, npoint)?,
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
        Ok(())
    }

    #[test]
    fn mul_interp_by_matrix_works() -> Result<(), StrError> {
        // iso-parametric elements: xᵢ = Σ_m Sᵐ ⋅ cᵐᵢ
        // where c is a matrix of coordinates
        let ndim_npoint = gen_ndim_npoint();
        for (ndim, npoint) in ndim_npoint {
            let mut shape = Shape::new(ndim, ndim, npoint)?;
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
