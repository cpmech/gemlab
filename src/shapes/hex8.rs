use super::Shape;
use russell_lab::{vec_mat_mul, Matrix, Vector};

const NDIM: usize = 3;
const NPOINT: usize = 8;
const NEDGE: usize = 12;
const NFACE: usize = 6;
const EDGE_NPOINT: usize = 2;
const FACE_NPOINT: usize = 4;
const FACE_NEDGE: usize = 4;

#[rustfmt::skip]
const EDGE_POINT_IDS: [[usize; 2]; 12] = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 0],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 4],
    [0, 4],
    [1, 5],
    [2, 6],
    [3, 7],
];

#[rustfmt::skip]
const FACE_POINT_IDS: [[usize; 4]; 6] = [
    [0, 4, 7, 3],
    [1, 2, 6, 5],
    [0, 1, 5, 4],
    [2, 3, 7, 6],
    [0, 3, 2, 1],
    [4, 5, 6, 7],
];

#[rustfmt::skip]
const FACE_EDGE_POINT_IDS: [[[usize; 2]; 4]; 6] = [
    [[0, 4], [4, 7], [7, 3], [3, 0]],
    [[1, 2], [2, 6], [6, 5], [5, 1]],
    [[0, 1], [1, 5], [5, 4], [4, 0]],
    [[2, 3], [3, 7], [7, 6], [6, 2]],
    [[0, 3], [3, 2], [2, 1], [1, 0]],
    [[4, 5], [5, 6], [6, 7], [7, 4]],
];

#[rustfmt::skip]
const POINT_NATURAL_COORDS: [[f64; 3]; 8] = [
    [-1.0, -1.0, -1.0],
    [ 1.0, -1.0, -1.0],
    [ 1.0,  1.0, -1.0],
    [-1.0,  1.0, -1.0],
    [-1.0, -1.0,  1.0],
    [ 1.0, -1.0,  1.0],
    [ 1.0,  1.0,  1.0],
    [-1.0,  1.0,  1.0],
];

/// Implements a hexahedron with 8 points
///
/// The natural coordinates range from -1 to +1 with the geometry centred @ 0.
///
/// # Local IDs of points
///
/// ```text
///           4________________7
///         ,'|              ,'|
///       ,'  |            ,'  |             r     s     t
///     ,'    |          ,'    |     p:0 [-1.0, -1.0, -1.0]
///   ,'      |        ,'      |     p:1 [ 1.0, -1.0, -1.0]
/// 5'===============6'        |     p:2 [ 1.0,  1.0, -1.0]
/// |         |      |         |     p:3 [-1.0,  1.0, -1.0]
/// |         |      |         |     p:4 [-1.0, -1.0,  1.0]
/// |         0_____ | ________3     p:5 [ 1.0, -1.0,  1.0]
/// |       ,'       |       ,'      p:6 [ 1.0,  1.0,  1.0]
/// |     ,'         |     ,'        p:7 [-1.0,  1.0,  1.0]
/// |   ,'           |   ,'
/// | ,'             | ,'
/// 1________________2'
/// ```
///
/// # Local IDs of edges
///
/// ```text
///                     7                   p0  p1
///             +----------------+     e:0  [0, 1]
///           ,'|              ,'|     e:1  [1, 2]
///       4 ,'  |8          6,'  |     e:2  [2, 3]
///       ,'    |          ,'    |     e:3  [3, 0]
///     ,'      |   5    ,'      |11   e:4  [4, 5]
///   +'===============+'        |     e:5  [5, 6]
///   |         |      |         |     e:6  [6, 7]
///   |         |      |  3      |     e:7  [7, 4]
///   |         .- - - | -  - - -+     e:8  [0, 4]
///  9|       ,'       |       ,'      e:9  [1, 5]
///   |    0,'         |10   ,'        e:10 [2, 6]
///   |   ,'           |   ,' 2        e:11 [3, 7]
///   | ,'             | ,'
///   +----------------+'
///           1
/// ```
///
/// # Local IDs of faces
///
/// Note: the order of points is such that the right-hand rule generates outward normals.
/// Also, the order of face points corresponds to Qua4 points.
///
/// ```text
///           4----------------7
///         ,'|              ,'|
///       ,'  |  ___       ,'  |
///     ,'    |,'5,'  [0],'    |        p0 p1 p2 p3
///   ,'      |~~~     ,'      |    f:0 [0, 4, 7, 3]
/// 5'===============6'  ,'|   |    f:1 [1, 2, 6, 5]
/// |   ,'|   |      |   |3|   |    f:2 [0, 1, 5, 4]
/// |   |2|   |      |   |,'   |    f:3 [2, 3, 7, 6]
/// |   |,'   0- - - | +- - - -3    f:4 [0, 3, 2, 1]
/// |       ,'       |       ,'     f:5 [4, 5, 6, 7]
/// |     ,' [1]  ___|     ,'
/// |   ,'      ,'4,'|   ,'
/// | ,'        ~~~  | ,'
/// 1----------------2'
/// ```
pub struct Hex8 {
    interp: Vector, // interpolation functions @ natural coordinate (npoint)
    deriv: Matrix,  // derivatives of interpolation functions w.r.t natural coordinate (npoint, ndim)
}

impl Hex8 {
    /// Creates a new object with pre-calculated interpolation fn and derivatives
    pub fn new() -> Self {
        Hex8 {
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

    fn mul_interp_by_matrix(&self, v: &mut Vector, a: &Matrix) -> Result<(), &'static str> {
        vec_mat_mul(v, 1.0, &self.interp, a)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_works() {
        let shape = Hex8::new();
        assert_eq!(shape.interp.dim(), NPOINT);
        assert_eq!(shape.deriv.dims(), (NPOINT, NDIM));
    }
}
