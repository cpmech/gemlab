use crate::{Shape, StrError};
use russell_lab::{vec_mat_mul, Matrix, Vector};

const NDIM: usize = 2;
const NPOINT: usize = 4;
const NEDGE: usize = 4;
const NFACE: usize = 0;
const EDGE_NPOINT: usize = 2;
const FACE_NPOINT: usize = 0;
const FACE_NEDGE: usize = 0;

#[rustfmt::skip]
const EDGE_POINT_IDS: [[usize; EDGE_NPOINT]; NEDGE] = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 0],
];

#[rustfmt::skip]
const POINT_NATURAL_COORDS: [[f64; NDIM]; NPOINT] = [
    [-1.0, -1.0],
    [ 1.0, -1.0],
    [ 1.0,  1.0],
    [-1.0,  1.0],
];

/// Implements a quadrilateral with 4 points
///
/// The natural coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of points
///
/// ```text
/// 3-----------2
/// |     s     |              r     s
/// |     |     |      p:0 [-1.0, -1.0]
/// |     +--r  |      p:1 [ 1.0, -1.0]
/// |           |      p:2 [ 1.0,  1.0]
/// |           |      p:3 [-1.0,  1.0]
/// 0-----------1
/// ```
///
/// # Local IDs of edges
///
/// ```text
///        2
///  +-----------+         p0  p1
///  |           |     e:0 [0, 1],
///  |           |     e:1 [1, 2],
/// 3|           |1    e:2 [2, 3],
///  |           |     e:3 [3, 0],
///  |           |
///  +-----------+
///        0
/// ```
pub struct Qua4 {
    interp: Vector, // interpolation functions @ natural coordinate (npoint)
    deriv: Matrix,  // derivatives of interpolation functions w.r.t natural coordinate (npoint, ndim)
}

impl Qua4 {
    pub fn new() -> Self {
        Qua4 {
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

    fn get_edge_npoint(&self) -> usize {
        EDGE_NPOINT
    }

    fn get_face_npoint(&self) -> usize {
        FACE_NPOINT
    }

    fn get_face_nedge(&self) -> usize {
        FACE_NEDGE
    }

    fn get_edge_local_point_id(&self, e: usize, i: usize) -> usize {
        EDGE_POINT_IDS[e][i]
    }

    fn get_face_local_point_id(&self, _: usize, _: usize) -> usize {
        0 // none
    }

    fn get_face_edge_local_point_id(&self, _: usize, _: usize, _: usize) -> usize {
        0 // none
    }

    fn get_ksi(&self, ksi: &mut Vector, m: usize) {
        ksi[0] = POINT_NATURAL_COORDS[m][0];
        ksi[1] = POINT_NATURAL_COORDS[m][1];
    }

    fn mul_interp_by_matrix(&self, v: &mut Vector, a: &Matrix) -> Result<(), StrError> {
        vec_mat_mul(v, 1.0, &self.interp, a)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_works() {
        let shape = Qua4::new();
        assert_eq!(shape.interp.dim(), NPOINT);
        assert_eq!(shape.deriv.dims(), (NPOINT, NDIM));
    }
}
