use super::Shape;
use russell_lab::{vec_mat_mul, Matrix, Vector};

const NDIM: usize = 3;
const NPOINT: usize = 8;
const NEDGE: usize = 12;
const NFACE: usize = 6;
const EDGE_NPOINT: usize = 2;
const FACE_NPOINT: usize = 4;

/// Implements a hexahedron with 8 points
///
/// The natural coordinates range from -1 to +1 with the geometry centred @ 0.
///
/// # Local IDs of points
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
///
/// # Local IDs of edges
///
/// ```text
///                     7
///             +----------------+
///           ,'|              ,'|
///       4 ,'  |8          6,'  |
///       ,'    |          ,'    |   
///     ,'      |   5    ,'      |11
///   +'===============+'        |
///   |         |      |         |
///   |         |      |  3      |
///   |         .- - - | -  - - -+
///  9|       ,'       |       ,'
///   |    0,'         |10   ,'
///   |   ,'           |   ,' 2
///   | ,'             | ,'
///   +----------------+'
///           1
/// ```
///
/// # Local IDs of faces
///
/// ```text
///           +----------------+
///         ,'|              ,'|
///       ,'  |  ___       ,'  |
///     ,'    |,'5,'  [0],'    |
///   ,'      |~~~     ,'      |
/// +'===============+'  ,'|   |
/// |   ,'|   |      |   |3|   |
/// |   |2|   |      |   |,'   |
/// |   |,'   +- - - | +- - - -+
/// |       ,'       |       ,'
/// |     ,' [1]  ___|     ,'
/// |   ,'      ,'4,'|   ,'
/// | ,'        ~~~  | ,'
/// +----------------+'
/// ```
pub struct Hex8 {
    coords: Vec<Vector>,       // natural coordinates (npoint, ndim)
    interp: Vector,            // interpolation functions @ natural coordinate (npoint)
    deriv: Matrix,             // derivatives of interpolation functions w.r.t natural coordinate (npoint, ndim)
    edge_ids: Vec<Vec<usize>>, // ids of points on edges
    face_ids: Vec<Vec<usize>>, // ids of points on faces
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
            edge_ids: vec![
                vec![0, 1],
                vec![1, 2],
                vec![2, 3],
                vec![3, 0],
                vec![4, 5],
                vec![5, 6],
                vec![6, 7],
                vec![7, 4],
                vec![0, 4],
                vec![1, 5],
                vec![2, 6],
                vec![3, 7],
            ],
            face_ids: vec![
                vec![0, 4, 7, 3],
                vec![1, 2, 6, 5],
                vec![0, 1, 5, 4],
                vec![2, 3, 7, 6],
                vec![0, 3, 2, 1],
                vec![4, 5, 6, 7],
            ],
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

    fn get_edge(&self, local_point_ids: &mut Vec<usize>, e: usize) {
        for i in 0..EDGE_NPOINT {
            local_point_ids[i] = self.edge_ids[e][i];
        }
    }

    fn get_face(&self, local_point_ids: &mut Vec<usize>, f: usize) {
        for i in 0..FACE_NPOINT {
            local_point_ids[i] = self.face_ids[f][i];
        }
    }

    fn get_ksi(&self, ksi: &mut Vector, m: usize) {
        ksi[0] = self.coords[m][0];
        ksi[1] = self.coords[m][1];
        ksi[2] = self.coords[m][2];
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
        let geo = Hex8::new();
        assert_eq!(geo.coords.len(), NPOINT);
        assert_eq!(geo.interp.dim(), NPOINT);
        assert_eq!(geo.deriv.dims(), (NPOINT, NDIM));
    }
}
