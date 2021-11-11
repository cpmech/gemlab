use super::Shape;
use crate::StrError;
use russell_lab::{vec_mat_mul, Matrix, Vector};

const NDIM: usize = 3;
const NPOINT: usize = 20;
const NEDGE: usize = 12;
const NFACE: usize = 6;
const EDGE_NPOINT: usize = 3;
const FACE_NPOINT: usize = 8;
const FACE_NEDGE: usize = 4;

#[rustfmt::skip]
const EDGE_POINT_IDS: [[usize; EDGE_NPOINT]; NEDGE] = [
    [0, 1,  8],
    [1, 2,  9],
    [2, 3, 10],
    [3, 0, 11],
    [4, 5, 12],
    [5, 6, 13],
    [6, 7, 14],
    [7, 4, 15],
    [0, 4, 16],
    [1, 5, 17],
    [2, 6, 18],
    [3, 7, 19],
];

#[rustfmt::skip]
const FACE_POINT_IDS: [[usize; FACE_NPOINT]; NFACE] = [
    [0, 4, 7, 3, 16, 15, 19, 11],
    [1, 2, 6, 5,  9, 18, 13, 17],
    [0, 1, 5, 4,  8, 17, 12, 16],
    [2, 3, 7, 6, 10, 19, 14, 18],
    [0, 3, 2, 1, 11, 10,  9,  8],
    [4, 5, 6, 7, 12, 13, 14, 15],
];

#[rustfmt::skip]
const FACE_EDGE_POINT_IDS: [[[usize; EDGE_NPOINT]; FACE_NEDGE]; NFACE] = [
    [[0, 4, 16], [4, 7, 15], [7, 3, 19], [3, 0, 11]],
    [[1, 2,  9], [2, 6, 18], [6, 5, 13], [5, 1, 17]],
    [[0, 1,  8], [1, 5, 17], [5, 4, 12], [4, 0, 16]],
    [[2, 3, 10], [3, 7, 19], [7, 6, 14], [6, 2, 18]],
    [[0, 3, 11], [3, 2, 10], [2, 1,  9], [1, 0,  8]],
    [[4, 5, 12], [5, 6, 13], [6, 7, 14], [7, 4, 15]],
];

#[rustfmt::skip]
const POINT_NATURAL_COORDS: [[f64; NDIM]; NPOINT] = [
    [-1.0, -1.0, -1.0],
    [ 1.0, -1.0, -1.0],
    [ 1.0,  1.0, -1.0],
    [-1.0,  1.0, -1.0],

    [-1.0, -1.0,  1.0],
    [ 1.0, -1.0,  1.0],
    [ 1.0,  1.0,  1.0],
    [-1.0,  1.0,  1.0],

    [ 0.0, -1.0, -1.0],
    [ 1.0,  0.0, -1.0],
    [ 0.0,  1.0, -1.0],
    [-1.0,  0.0, -1.0],

    [ 0.0, -1.0,  1.0],
    [ 1.0,  0.0,  1.0],
    [ 0.0,  1.0,  1.0],
    [-1.0,  0.0,  1.0],

    [-1.0, -1.0,  0.0],
    [ 1.0, -1.0,  0.0],
    [ 1.0,  1.0,  0.0],
    [-1.0,  1.0,  0.0],
];

/// Implements a hexahedron with 20 points
///
/// The natural coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of points
///
/// ```text
///            4_______15_______7            r     s     t             r     s     t
///          ,'|              ,'|    p:0 [-1.0, -1.0, -1.0]   p:8  [ 0.0, -1.0, -1.0]
///       12'  |            ,'  |    p:1 [ 1.0, -1.0, -1.0]   p:9  [ 1.0,  0.0, -1.0]
///      ,'    16         ,14   |    p:2 [ 1.0,  1.0, -1.0]   p:10 [ 0.0,  1.0, -1.0]
///    ,'      |        ,'      19   p:3 [-1.0,  1.0, -1.0]   p:11 [-1.0,  0.0, -1.0]
///  5'=====13========6'        |                             
///  |         |      |         |    p:4 [-1.0, -1.0,  1.0]   p:12 [ 0.0, -1.0,  1.0]
///  |         |      |         |    p:5 [ 1.0, -1.0,  1.0]   p:13 [ 1.0,  0.0,  1.0]
///  |         0_____ | _11_____3    p:6 [ 1.0,  1.0,  1.0]   p:14 [ 0.0,  1.0,  1.0]
/// 17       ,'       |       ,'     p:7 [-1.0,  1.0,  1.0]   p:15 [-1.0,  0.0,  1.0]
///  |     8'        18     ,'       
///  |   ,'           |   ,10       p:16 [-1.0, -1.0,  0.0]
///  | ,'             | ,'          p:17 [ 1.0, -1.0,  0.0]
///  1_______9________2'            p:18 [ 1.0,  1.0,  0.0]
///                                 p:19 [-1.0,  1.0,  0.0]
/// ```
///
/// # Local IDs of edges
///
/// ```text
///                     7
///             +----------------+          p0  p1 p2
///           ,'|              ,'|     e:0  [0, 1,  8],
///       4 ,'  |8          6,'  |     e:1  [1, 2,  9],
///       ,'    |          ,'    |     e:2  [2, 3, 10],
///     ,'      |   5    ,'      |11   e:3  [3, 0, 11],
///   +'===============+'        |     e:4  [4, 5, 12],
///   |         |      |         |     e:5  [5, 6, 13],
///   |         |      |  3      |     e:6  [6, 7, 14],
///   |         .- - - | -  - - -+     e:7  [7, 4, 15],
///  9|       ,'       |       ,'      e:8  [0, 4, 16],
///   |    0,'         |10   ,'        e:9  [1, 5, 17],
///   |   ,'           |   ,' 2        e:10 [2, 6, 18],
///   | ,'             | ,'            e:11 [3, 7, 19],
///   +----------------+'
///           1
/// ```
///
/// # Local IDs of faces
///
/// Note: the order of points is such that the right-hand rule generates outward normals.
/// Also, the order of face points corresponds to Qua8 points.
///
/// ```text
///           4-------15-------7
///         ,'|              ,'|
///      12' 16  ___       14  |
///     ,'    |,'5,'  [0],'    |
///   ,'      |~~~     ,'     19   f:0 [0, 4, 7, 3, 16, 15, 19, 11]
/// 5'======13=======6'  ,'|   |   f:1 [1, 2, 6, 5,  9, 18, 13, 17]
/// |   ,'|   |      |   |3|   |   f:2 [0, 1, 5, 4,  8, 17, 12, 16]
/// |   |2|   |      |   |,'   |   f:3 [2, 3, 7, 6, 10, 19, 14, 18]
/// |   |,'   0- - - | 11 - - -3   f:4 [0, 3, 2, 1, 11, 10,  9,  8]
/// 17      ,'      18       ,'    f:5 [4, 5, 6, 7, 12, 13, 14, 15]
/// |     8' [1]  ___|     10
/// |   ,'      ,'4,'|   ,'
/// | ,'        ~~~  | ,'
/// 1-------9--------2'
/// ```
pub struct Hex20 {
    interp: Vector, // interpolation functions @ natural coordinate (npoint)
    deriv: Matrix,  // derivatives of interpolation functions w.r.t natural coordinate (npoint, ndim)
    coords: Matrix, // (npoint, space_ndim) real coordinates matrix
}

impl Hex20 {
    /// Creates a new Hex20
    ///
    /// **space_ndim** must be equal to 3
    pub fn new(space_ndim: usize) -> Result<Self, StrError> {
        if space_ndim != 3 {
            return Err("space_ndim must be 3 for Hex20");
        }
        Ok(Hex20 {
            interp: Vector::new(NPOINT),
            deriv: Matrix::new(NPOINT, NDIM),
            coords: Matrix::new(NPOINT, space_ndim),
        })
    }
}

impl Shape for Hex20 {
    fn calc_interp(&mut self, ksi: &Vector) {
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        let rp1 = 1.0 + r;
        let rm1 = 1.0 - r;
        let sp1 = 1.0 + s;
        let sm1 = 1.0 - s;
        let tp1 = 1.0 + t;
        let tm1 = 1.0 - t;

        self.interp[0] = rm1 * sm1 * tm1 * (-r - s - t - 2.0) / 8.0;
        self.interp[1] = rp1 * sm1 * tm1 * (r - s - t - 2.0) / 8.0;
        self.interp[2] = rp1 * sp1 * tm1 * (r + s - t - 2.0) / 8.0;
        self.interp[3] = rm1 * sp1 * tm1 * (-r + s - t - 2.0) / 8.0;
        self.interp[4] = rm1 * sm1 * tp1 * (-r - s + t - 2.0) / 8.0;
        self.interp[5] = rp1 * sm1 * tp1 * (r - s + t - 2.0) / 8.0;
        self.interp[6] = rp1 * sp1 * tp1 * (r + s + t - 2.0) / 8.0;
        self.interp[7] = rm1 * sp1 * tp1 * (-r + s + t - 2.0) / 8.0;
        self.interp[8] = (1.0 - r * r) * sm1 * tm1 / 4.0;
        self.interp[9] = rp1 * (1.0 - s * s) * tm1 / 4.0;
        self.interp[10] = (1.0 - r * r) * sp1 * tm1 / 4.0;
        self.interp[11] = rm1 * (1.0 - s * s) * tm1 / 4.0;
        self.interp[12] = (1.0 - r * r) * sm1 * tp1 / 4.0;
        self.interp[13] = rp1 * (1.0 - s * s) * tp1 / 4.0;
        self.interp[14] = (1.0 - r * r) * sp1 * tp1 / 4.0;
        self.interp[15] = rm1 * (1.0 - s * s) * tp1 / 4.0;
        self.interp[16] = rm1 * sm1 * (1.0 - t * t) / 4.0;
        self.interp[17] = rp1 * sm1 * (1.0 - t * t) / 4.0;
        self.interp[18] = rp1 * sp1 * (1.0 - t * t) / 4.0;
        self.interp[19] = rm1 * sp1 * (1.0 - t * t) / 4.0;
    }

    fn calc_deriv(&mut self, ksi: &Vector) {
        let (r, s, t) = (ksi[0], ksi[1], ksi[2]);

        let rp1 = 1.0 + r;
        let rm1 = 1.0 - r;
        let sp1 = 1.0 + s;
        let sm1 = 1.0 - s;
        let tp1 = 1.0 + t;
        let tm1 = 1.0 - t;

        self.deriv[0][0] = -0.125 * sm1 * tm1 * (-r - s - t - 2.0) - 0.125 * rm1 * sm1 * tm1;
        self.deriv[1][0] = 0.125 * sm1 * tm1 * (r - s - t - 2.0) + 0.125 * rp1 * sm1 * tm1;
        self.deriv[2][0] = 0.125 * sp1 * tm1 * (r + s - t - 2.0) + 0.125 * rp1 * sp1 * tm1;
        self.deriv[3][0] = -0.125 * sp1 * tm1 * (-r + s - t - 2.0) - 0.125 * rm1 * sp1 * tm1;
        self.deriv[4][0] = -0.125 * sm1 * tp1 * (-r - s + t - 2.0) - 0.125 * rm1 * sm1 * tp1;
        self.deriv[5][0] = 0.125 * sm1 * tp1 * (r - s + t - 2.0) + 0.125 * rp1 * sm1 * tp1;
        self.deriv[6][0] = 0.125 * sp1 * tp1 * (r + s + t - 2.0) + 0.125 * rp1 * sp1 * tp1;
        self.deriv[7][0] = -0.125 * sp1 * tp1 * (-r + s + t - 2.0) - 0.125 * rm1 * sp1 * tp1;
        self.deriv[8][0] = -0.5 * r * sm1 * tm1;
        self.deriv[9][0] = 0.25 * (1.0 - s * s) * tm1;
        self.deriv[10][0] = -0.5 * r * sp1 * tm1;
        self.deriv[11][0] = -0.25 * (1.0 - s * s) * tm1;
        self.deriv[12][0] = -0.5 * r * sm1 * tp1;
        self.deriv[13][0] = 0.25 * (1.0 - s * s) * tp1;
        self.deriv[14][0] = -0.5 * r * sp1 * tp1;
        self.deriv[15][0] = -0.25 * (1.0 - s * s) * tp1;
        self.deriv[16][0] = -0.25 * sm1 * (1.0 - t * t);
        self.deriv[17][0] = 0.25 * sm1 * (1.0 - t * t);
        self.deriv[18][0] = 0.25 * sp1 * (1.0 - t * t);
        self.deriv[19][0] = -0.25 * sp1 * (1.0 - t * t);

        self.deriv[0][1] = -0.125 * rm1 * tm1 * (-r - s - t - 2.0) - 0.125 * rm1 * sm1 * tm1;
        self.deriv[1][1] = -0.125 * rp1 * tm1 * (r - s - t - 2.0) - 0.125 * rp1 * sm1 * tm1;
        self.deriv[2][1] = 0.125 * rp1 * tm1 * (r + s - t - 2.0) + 0.125 * rp1 * sp1 * tm1;
        self.deriv[3][1] = 0.125 * rm1 * tm1 * (-r + s - t - 2.0) + 0.125 * rm1 * sp1 * tm1;
        self.deriv[4][1] = -0.125 * rm1 * tp1 * (-r - s + t - 2.0) - 0.125 * rm1 * sm1 * tp1;
        self.deriv[5][1] = -0.125 * rp1 * tp1 * (r - s + t - 2.0) - 0.125 * rp1 * sm1 * tp1;
        self.deriv[6][1] = 0.125 * rp1 * tp1 * (r + s + t - 2.0) + 0.125 * rp1 * sp1 * tp1;
        self.deriv[7][1] = 0.125 * rm1 * tp1 * (-r + s + t - 2.0) + 0.125 * rm1 * sp1 * tp1;
        self.deriv[8][1] = -0.25 * (1.0 - r * r) * tm1;
        self.deriv[9][1] = -0.5 * s * rp1 * tm1;
        self.deriv[10][1] = 0.25 * (1.0 - r * r) * tm1;
        self.deriv[11][1] = -0.5 * s * rm1 * tm1;
        self.deriv[12][1] = -0.25 * (1.0 - r * r) * tp1;
        self.deriv[13][1] = -0.5 * s * rp1 * tp1;
        self.deriv[14][1] = 0.25 * (1.0 - r * r) * tp1;
        self.deriv[15][1] = -0.5 * s * rm1 * tp1;
        self.deriv[16][1] = -0.25 * rm1 * (1.0 - t * t);
        self.deriv[17][1] = -0.25 * rp1 * (1.0 - t * t);
        self.deriv[18][1] = 0.25 * rp1 * (1.0 - t * t);
        self.deriv[19][1] = 0.25 * rm1 * (1.0 - t * t);

        self.deriv[0][2] = -0.125 * rm1 * sm1 * (-r - s - t - 2.0) - 0.125 * rm1 * sm1 * tm1;
        self.deriv[1][2] = -0.125 * rp1 * sm1 * (r - s - t - 2.0) - 0.125 * rp1 * sm1 * tm1;
        self.deriv[2][2] = -0.125 * rp1 * sp1 * (r + s - t - 2.0) - 0.125 * rp1 * sp1 * tm1;
        self.deriv[3][2] = -0.125 * rm1 * sp1 * (-r + s - t - 2.0) - 0.125 * rm1 * sp1 * tm1;
        self.deriv[4][2] = 0.125 * rm1 * sm1 * (-r - s + t - 2.0) + 0.125 * rm1 * sm1 * tp1;
        self.deriv[5][2] = 0.125 * rp1 * sm1 * (r - s + t - 2.0) + 0.125 * rp1 * sm1 * tp1;
        self.deriv[6][2] = 0.125 * rp1 * sp1 * (r + s + t - 2.0) + 0.125 * rp1 * sp1 * tp1;
        self.deriv[7][2] = 0.125 * rm1 * sp1 * (-r + s + t - 2.0) + 0.125 * rm1 * sp1 * tp1;
        self.deriv[8][2] = -0.25 * (1.0 - r * r) * sm1;
        self.deriv[9][2] = -0.25 * rp1 * (1.0 - s * s);
        self.deriv[10][2] = -0.25 * (1.0 - r * r) * sp1;
        self.deriv[11][2] = -0.25 * rm1 * (1.0 - s * s);
        self.deriv[12][2] = 0.25 * (1.0 - r * r) * sm1;
        self.deriv[13][2] = 0.25 * rp1 * (1.0 - s * s);
        self.deriv[14][2] = 0.25 * (1.0 - r * r) * sp1;
        self.deriv[15][2] = 0.25 * rm1 * (1.0 - s * s);
        self.deriv[16][2] = -0.5 * t * rm1 * sm1;
        self.deriv[17][2] = -0.5 * t * rp1 * sm1;
        self.deriv[18][2] = -0.5 * t * rp1 * sp1;
        self.deriv[19][2] = -0.5 * t * rm1 * sp1;
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

    fn get_face_local_point_id(&self, f: usize, i: usize) -> usize {
        FACE_POINT_IDS[f][i]
    }

    fn get_face_edge_local_point_id(&self, f: usize, k: usize, i: usize) -> usize {
        FACE_EDGE_POINT_IDS[f][k][i]
    }

    fn get_ksi(&self, ksi: &mut Vector, m: usize) {
        ksi[0] = POINT_NATURAL_COORDS[m][0];
        ksi[1] = POINT_NATURAL_COORDS[m][1];
        ksi[2] = POINT_NATURAL_COORDS[m][2];
    }

    fn mul_interp_by_matrix(&self, v: &mut Vector, a: &Matrix) -> Result<(), StrError> {
        vec_mat_mul(v, 1.0, &self.interp, a)
    }

    fn set_coords(&mut self, m: usize, i: usize, val: f64) {
        self.coords[m][i] = val;
    }

    fn get_coords_matrix(&self) -> &Matrix {
        &self.coords
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_works() -> Result<(), StrError> {
        let shape = Hex20::new(3)?;
        assert_eq!(shape.interp.dim(), NPOINT);
        assert_eq!(shape.deriv.dims(), (NPOINT, NDIM));
        Ok(())
    }
}
