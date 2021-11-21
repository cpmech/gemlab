use super::{Hex20, Hex8, Qua4, Qua8};
use crate::StrError;
use russell_lab::{inverse, mat_mat_mul, mat_vec_mul, Matrix, Vector};

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

/// Defines an isoparametric geometric shape for numerical integration and more
///
/// The isoparametric formulation establishes that
///
/// ```text
/// → →         →  →
/// x(ξ) = Σ Sᵐ(ξ) xᵐ
///        m         
/// ```
///
/// where x is the mapped vector of real coordinates and ξ is the vector of natural (reference) coordinates.
/// Above, xm are the coordinates of each m-point of the geometric shape.
///
/// Given the matrix of coordinates X, we can calculate the mapped coordinates vector x by means of
///
/// ```text
/// x = Xᵀ ⋅ S
/// ```
///
/// The Jacobian tensor is
///
/// ```text
///         →
///   →    dx     →    →
/// J(ξ) = —— = Σ xᵐ ⊗ gᵐ
///         →   m
///        dξ
/// ```
///
/// where
///
/// ```text
///             →
/// →       dSᵐ(ξ)
/// gᵐ(ξ) = ——————
///            →
///           dξ
/// ```
///
/// We can write the Jacobian in matrix notation as follows
///
/// ```text
/// J = Xᵀ · g
/// ```
///
/// where X is the matrix of coordinates and g is a (npoint,shape_ndim) matrix.
///
/// The gradient of interpolation functions (derivatives w.r.t real coordinates) is given by
///
/// ```text
///             →
/// →       dSᵐ(ξ)
/// Gᵐ(ξ) = ——————
///            →
///           dx
/// ```
///
/// The inverse Jacobian allows us to determine the gradient vectors G as follows
///
/// ```text
/// →       →  →        →
/// Gᵐ(ξ) = gᵐ(ξ) · J⁻¹(ξ)
/// ```
///
/// Or, in matrix notation,
///
/// ```text
/// G = g · J⁻¹
/// ```
///
/// where G is a (npoint,space_ndim) matrix.
pub struct Shape {
    // space ndim and shape kind
    space_ndim: usize, // space ndim (not shape ndim)
    kind: ShapeKind,   // shape kind

    // constants (from each specific shape)
    shape_ndim: usize,  // shape ndim (not space ndim)
    npoint: usize,      // number of points
    nedge: usize,       // number of edges
    nface: usize,       // number of faces
    edge_npoint: usize, // edge's number of points
    face_npoint: usize, // face's number of points
    face_nedge: usize,  // face's number of edges

    // functions (from each specific shape)
    fn_interp: FnInterp, // function to calculate interpolation functions
    fn_deriv: FnDeriv,   // function to calculate derivatives of interpolation functions

    // input data
    xx_tra: Matrix, // (space_ndim, npoint) transposed coordinates matrix (real space)
    ok_xx: bool,    // User has provided X matrix

    // sandbox (temporary) variables
    ss: Vector,     // (npoint) interpolation functions @ natural coordinate ksi
    g: Matrix,      // (npoint, shape_ndim) derivatives of interpolation functions w.r.t natural coordinate ksi
    jj: Matrix,     // (space_ndim, space_ndim) Jacobian matrix
    jj_inv: Matrix, // (space_ndim, space_ndim) Inverse Jacobian matrix
}

impl Shape {
    /// Creates a new geometric shape
    ///
    /// # Input
    ///
    /// * `space_ndim` -- the space dimension (1, 2, or 3)
    /// * `shape_ndim` -- the dimension of the shape; e.g. a 1D line in a 3D space
    /// * `npoint` -- the number of points defining the shape
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
                    space_ndim,
                    kind: ShapeKind::Qua4,
                    shape_ndim: Qua4::NDIM,
                    npoint: Qua4::NPOINT,
                    nedge: Qua4::NEDGE,
                    nface: Qua4::NFACE,
                    edge_npoint: Qua4::EDGE_NPOINT,
                    face_npoint: Qua4::FACE_NPOINT,
                    face_nedge: Qua4::FACE_NEDGE,
                    fn_interp: Qua4::calc_interp,
                    fn_deriv: Qua4::calc_deriv,
                    xx_tra: Matrix::new(space_ndim, Qua4::NPOINT),
                    ok_xx: false,
                    ss: Vector::new(Qua4::NPOINT),
                    g: Matrix::new(Qua4::NPOINT, Qua4::NDIM),
                    jj: Matrix::new(space_ndim, space_ndim),
                    jj_inv: Matrix::new(space_ndim, space_ndim),
                });
            }
            (2, 8) => {
                return Ok(Shape {
                    space_ndim,
                    kind: ShapeKind::Qua8,
                    shape_ndim: Qua8::NDIM,
                    npoint: Qua8::NPOINT,
                    nedge: Qua8::NEDGE,
                    nface: Qua8::NFACE,
                    edge_npoint: Qua8::EDGE_NPOINT,
                    face_npoint: Qua8::FACE_NPOINT,
                    face_nedge: Qua8::FACE_NEDGE,
                    fn_interp: Qua8::calc_interp,
                    fn_deriv: Qua8::calc_deriv,
                    xx_tra: Matrix::new(space_ndim, Qua8::NPOINT),
                    ok_xx: false,
                    ss: Vector::new(Qua8::NPOINT),
                    g: Matrix::new(Qua8::NPOINT, Qua8::NDIM),
                    jj: Matrix::new(space_ndim, space_ndim),
                    jj_inv: Matrix::new(space_ndim, space_ndim),
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
                    space_ndim,
                    kind: ShapeKind::Hex8,
                    shape_ndim: Hex8::NDIM,
                    npoint: Hex8::NPOINT,
                    nedge: Hex8::NEDGE,
                    nface: Hex8::NFACE,
                    edge_npoint: Hex8::EDGE_NPOINT,
                    face_npoint: Hex8::FACE_NPOINT,
                    face_nedge: Hex8::FACE_NEDGE,
                    fn_interp: Hex8::calc_interp,
                    fn_deriv: Hex8::calc_deriv,
                    xx_tra: Matrix::new(space_ndim, Hex8::NPOINT),
                    ok_xx: false,
                    ss: Vector::new(Hex8::NPOINT),
                    g: Matrix::new(Hex8::NPOINT, Hex8::NDIM),
                    jj: Matrix::new(space_ndim, space_ndim),
                    jj_inv: Matrix::new(space_ndim, space_ndim),
                });
            }
            (3, 20) => {
                return Ok(Shape {
                    space_ndim,
                    kind: ShapeKind::Hex20,
                    shape_ndim: Hex20::NDIM,
                    npoint: Hex20::NPOINT,
                    nedge: Hex20::NEDGE,
                    nface: Hex20::NFACE,
                    edge_npoint: Hex20::EDGE_NPOINT,
                    face_npoint: Hex20::FACE_NPOINT,
                    face_nedge: Hex20::FACE_NEDGE,
                    fn_interp: Hex20::calc_interp,
                    fn_deriv: Hex20::calc_deriv,
                    xx_tra: Matrix::new(space_ndim, Hex20::NPOINT),
                    ok_xx: false,
                    ss: Vector::new(Hex20::NPOINT),
                    g: Matrix::new(Hex20::NPOINT, Hex20::NDIM),
                    jj: Matrix::new(space_ndim, space_ndim),
                    jj_inv: Matrix::new(space_ndim, space_ndim),
                });
            }
            _ => return Err("(space_ndim, shape_ndim, npoint) combination is invalid"),
        };
    }

    /// Sets the coordinates matrix
    ///
    /// # Input
    ///
    /// * `xx` -- (npoint, space_ndim) matrix of coordinates
    pub fn set_coords_matrix(&mut self, xx: &Matrix) -> Result<(), StrError> {
        let (m, n) = xx.dims();
        if m != self.npoint || n != self.space_ndim {
            return Err("matrix of coordinates has invalid dimensions");
        }
        for i in 0..self.npoint {
            for j in 0..self.space_ndim {
                self.xx_tra[j][i] = xx[i][j];
            }
        }
        self.ok_xx = true;
        Ok(())
    }

    /// Calculates the real coordinates x from natural coordinates ξ
    ///
    /// The isoparametric formulation establishes:
    ///
    /// ```text
    /// → →         →  →
    /// x(ξ) = Σ Sᵐ(ξ) xᵐ
    ///        m         
    /// ```
    ///
    /// or
    ///
    /// ```text
    /// x := Xᵀ ⋅ S
    /// ```
    pub fn calc_coords(&mut self, x: &mut Vector, ksi: &Vector) -> Result<(), StrError> {
        if !self.ok_xx {
            return Err("the coordinates matrix (xx_mat) must be set first");
        }
        (self.fn_interp)(&mut self.ss, ksi);
        mat_vec_mul(x, 1.0, &self.xx_tra, &self.ss)
    }

    /// Calculates the Jacobian matrix of the mapping from general to natural space
    ///
    /// In terms of Cartesian ortho-normal components, the Jacobian is
    ///
    /// ```text
    ///        ∂xᵢ
    /// Jᵢⱼ := ——— = Σ X[m][i] * g[m][j]
    ///        ∂ξⱼ   m
    /// ```
    ///
    /// Therefore, in matrix notation,
    ///
    /// ```text
    ///               ∂xi
    /// J_mat[i][j] = ———
    ///               ∂rj
    ///          _                           _
    ///         |  ∂x0/∂ξ0  ∂x0/∂ξ1  ∂x0/∂ξ2  |     
    /// J_mat = |  ∂x1/∂ξ0  ∂x1/∂ξ1  ∂x1/∂ξ2  |     
    ///         |_ ∂x2/∂ξ0  ∂x2/∂ξ1  ∂x2/∂ξ2 _|     
    ///
    /// J_mat = Xᵀ · g
    /// ```
    ///
    /// # Input
    ///
    /// * `ksi` -- ξ (natural) coordinate
    ///
    /// # Output
    ///
    /// * Returns the determinant of the Jacobian
    /// * stored: `jj` -- (space_ndim,space_ndim) Jacobian matrix
    /// * stored: `jj_inv` -- (space_ndim,space_ndim) inverse of the Jacobian matrix
    pub fn calc_jacobian(&mut self, ksi: &Vector) -> Result<f64, StrError> {
        if !self.ok_xx {
            return Err("the coordinates matrix (xx_mat) must be set first");
        }
        (self.fn_deriv)(&mut self.g, ksi);
        mat_mat_mul(&mut self.jj, 1.0, &self.xx_tra, &self.g)?;
        inverse(&mut self.jj_inv, &self.jj)
    }

    // --- getters ------------------------------------------------------------------------------

    /// Returns the number of dimensions of space (not shape)
    #[inline]
    pub fn get_space_ndim(&self) -> usize {
        self.space_ndim
    }

    /// Returns the number of dimensions of the shape (not space)
    #[inline]
    pub fn get_shape_ndim(&self) -> usize {
        self.shape_ndim
    }

    /// Returns the number of points
    #[inline]
    pub fn get_npoint(&self) -> usize {
        self.npoint
    }

    /// Returns the number of edges
    #[inline]
    pub fn get_nedge(&self) -> usize {
        self.nedge
    }

    /// Returns the number of faces
    #[inline]
    pub fn get_nface(&self) -> usize {
        self.nface
    }

    /// Returns the number of points on edge
    #[inline]
    pub fn get_edge_npoint(&self) -> usize {
        self.edge_npoint
    }

    /// Returns the number of points on face
    #[inline]
    pub fn get_face_npoint(&self) -> usize {
        self.face_npoint
    }

    /// Returns the number of edges on face (not the cell's edges)
    #[inline]
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

    // --- private ------------------------------------------------------------------------------

    ///
    ///
    ///
    ///
    ///
    fn calc_ss(&mut self) {}

    /// ```text
    /// →       dSᵐ(ξ)
    /// gᵐ(ξ) = ——————
    ///            →
    ///           dξ
    /// ```
    fn calc_g(&mut self) {}
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

    // Computes Sᵐ(ξ) with variable ξ
    fn sm(v: f64, args: &mut Arguments) -> f64 {
        for i in 0..args.shape.get_shape_ndim() {
            args.ksi[i] = args.at_ksi[i];
        }
        args.ksi[args.i] = v;
        (args.shape.fn_interp)(&mut args.shape.ss, &args.ksi);
        args.shape.ss[args.m]
    }

    #[test]
    fn fn_interp_works() -> Result<(), StrError> {
        let ndim_npoint = gen_ndim_npoint();
        for (ndim, npoint) in ndim_npoint {
            let mut shape = Shape::new(ndim, ndim, npoint)?;
            assert_eq!(shape.get_shape_ndim(), ndim);
            assert_eq!(shape.get_npoint(), npoint);
            let mut ksi = Vector::new(ndim);
            for m in 0..npoint {
                shape.get_ksi(&mut ksi, m);
                (shape.fn_interp)(&mut shape.ss, &ksi);
                for n in 0..npoint {
                    let smn = shape.ss[n];
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
    fn fn_deriv_works() -> Result<(), StrError> {
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
            (shape.fn_deriv)(&mut shape.g, &args.at_ksi);
            for m in 0..npoint {
                args.m = m;
                for i in 0..ndim {
                    args.i = i;
                    // gᵐᵢ := dSᵐ/dξᵢ
                    assert_deriv_approx_eq!(shape.g[m][i], args.at_ksi[i], sm, args, 1e-12);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn calc_coords_works() -> Result<(), StrError> {
        let ndim_npoint = gen_ndim_npoint();
        for (ndim, npoint) in ndim_npoint {
            let mut shape = Shape::new(ndim, ndim, npoint)?;
            let ndim = shape.get_shape_ndim();
            let npoint = shape.get_npoint();
            let mut ksi = Vector::new(ndim);
            // create matrix of coordinates with shifted nat-coords: edges = [0, 2]
            let mut xx = Matrix::new(npoint, ndim);
            for m in 0..npoint {
                shape.get_ksi(&mut ksi, m);
                for i in 0..ndim {
                    xx[m][i] = 1.0 + ksi[i];
                }
            }
            // set @ksi such that x[i] = 1.25
            let mut at_ksi = Vector::new(ndim);
            for i in 0..ndim {
                at_ksi[i] = 0.25;
            }
            // calculate coords
            let mut x = Vector::new(ndim);
            shape.calc_coords(&mut x, &at_ksi)?;
            // check
            let x_correct = vec![1.25; ndim];
            assert_vec_approx_eq!(x.as_data(), x_correct, 1e-15);
        }
        Ok(())
    }
}
