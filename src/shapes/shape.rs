use super::{Hex20, Hex8, Qua4, Qua8};
use crate::StrError;
use russell_lab::{Matrix, Vector};

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

    /// Returns the number of edges on face (not the cell's edges)
    fn get_face_nedge(&self) -> usize;

    /// Returns the local id of point on edge
    ///
    /// # Input
    ///
    /// * `e` -- index of edge in [0, nedge-1]
    /// * `i` -- index of local point [0, edge_npoint-1]
    fn get_edge_local_point_id(&self, e: usize, i: usize) -> usize;

    /// Returns the local id of point on face
    ///
    /// # Input
    ///
    /// * `f` -- index of face in [0, nface-1]
    /// * `i` -- index of local point [0, face_npoint-1]
    fn get_face_local_point_id(&self, f: usize, i: usize) -> usize;

    /// Returns the local point id on an edge on the face
    ///
    /// # Input
    ///
    /// * `f` -- index of face in [0, nface-1]
    /// * `k` -- index of face's edge (not the index of cell's edge) in [0, face_nedge-1]
    /// * `i` -- index of local point [0, edge_npoint-1]
    fn get_face_edge_local_point_id(&self, f: usize, k: usize, i: usize) -> usize;

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
    fn mul_interp_by_matrix(&self, v: &mut Vector, a: &Matrix) -> Result<(), StrError>;

    /// Sets the real coordinates i (0..space_ndim-1) of a point (0..npoint-1)
    fn set_coords(&mut self, m: usize, i: usize, val: f64);

    /// Returns a read-only access to the real coordinates matrix (npoint, space_ndim)
    fn get_coords_matrix(&self) -> &Matrix;

    /// Calculates the real coordinates of ksi
    ///
    /// ```text
    /// x[i] := Σ_m Sᵐ(ksi) ⋅ point_coords[m][i]
    /// ```
    fn calc_real_ksi_coords(&mut self, x: &mut Vector, ksi: &Vector) -> Result<(), StrError> {
        self.calc_interp(ksi);
        self.mul_interp_by_matrix(x, self.get_coords_matrix())
    }
}

/// Returns new Shape
pub fn new_shape(space_ndim: usize, shape_ndim: usize, npoint: usize) -> Result<Box<dyn Shape>, StrError> {
    match (shape_ndim, npoint) {
        (1, 2) => Err("Lin2 is not available yet"),
        (1, 3) => Err("Lin3 is not available yet"),
        (1, 4) => Err("Lin4 is not available yet"),
        (1, 5) => Err("Lin5 is not available yet"),
        (2, 3) => Err("Tri3 is not available yet"),
        (2, 6) => Err("Tri6 is not available yet"),
        (2, 10) => Err("Tri10 is not available yet"),
        (2, 15) => Err("Tri15 is not available yet"),
        (2, 4) => Ok(Box::new(Qua4::new(space_ndim)?)),
        (2, 8) => Ok(Box::new(Qua8::new(space_ndim)?)),
        (2, 9) => Err("Qua9 is not available yet"),
        (2, 12) => Err("Qua12 is not available yet"),
        (2, 16) => Err("Qua16 is not available yet"),
        (3, 4) => Err("Tet4 is not available yet"),
        (3, 10) => Err("Tet10 is not available yet"),
        (3, 8) => Ok(Box::new(Hex8::new(space_ndim)?)),
        (3, 20) => Ok(Box::new(Hex20::new(space_ndim)?)),
        _ => Err("(ndim, npoint) combination is invalid"),
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
    fn calc_interp_works() -> Result<(), StrError> {
        let ndim_npoint = gen_ndim_npoint();
        for (ndim, npoint) in ndim_npoint {
            let mut shape = new_shape(ndim, ndim, npoint)?;
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
            let mut shape = new_shape(ndim, ndim, npoint)?;
            let mut at_ksi = Vector::new(ndim);
            for i in 0..ndim {
                at_ksi[i] = 0.25;
            }
            let args = &mut Arguments {
                shape: new_shape(ndim, ndim, npoint)?,
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
            let mut shape = new_shape(ndim, ndim, npoint)?;
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
