use super::Shape;
use russell_lab::{vec_mat_mul, Matrix, Vector};

const NDIM: usize = 3;
const NPOINT: usize = 20;
const NEDGE: usize = 12;
const NFACE: usize = 6;
const EDGE_NPOINT: usize = 3;
const FACE_NPOINT: usize = 8;

/// Implements a hexahedron with 20 points
///
/// The natural coordinates range from -1 to +1 with the geometry centred @ 0
///
/// # Local IDs of points
///
/// ```text
///            4_______15_______7
///          ,'|              ,'|
///       12'  |            ,'  |
///      ,'    16         ,14   |
///    ,'      |        ,'      19
///  5'=====13========6'        |
///  |         |      |         |
///  |         |      |         |
///  |         0_____ | _11_____3
/// 17       ,'       |       ,'
///  |     8'        18     ,'
///  |   ,'           |   ,10
///  | ,'             | ,'
///  1_______9________2'
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
pub struct Hex20 {
    coords: Vec<Vector>,       // natural coordinates (npoint, ndim)
    interp: Vector,            // interpolation functions @ natural coordinate (npoint)
    deriv: Matrix,             // derivatives of interpolation functions w.r.t natural coordinate (npoint, ndim)
    edge_ids: Vec<Vec<usize>>, // ids of points on edges
    face_ids: Vec<Vec<usize>>, // ids of points on faces
}

impl Hex20 {
    /// Creates a new object with pre-calculated interpolation fn and derivatives
    pub fn new() -> Self {
        Hex20 {
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

                Vector::from(&[ 0.0, -1.0, -1.0]),
                Vector::from(&[ 1.0,  0.0, -1.0]),
                Vector::from(&[ 0.0,  1.0, -1.0]),
                Vector::from(&[-1.0,  0.0, -1.0]),

                Vector::from(&[ 0.0, -1.0,  1.0]),
                Vector::from(&[ 1.0,  0.0,  1.0]),
                Vector::from(&[ 0.0,  1.0,  1.0]),
                Vector::from(&[-1.0,  0.0,  1.0]),

                Vector::from(&[-1.0, -1.0,  0.0]),
                Vector::from(&[ 1.0, -1.0,  0.0]),
                Vector::from(&[ 1.0,  1.0,  0.0]),
                Vector::from(&[-1.0,  1.0,  0.0]),
            ],
            interp: Vector::new(NPOINT),
            deriv: Matrix::new(NPOINT, NDIM),
            edge_ids: vec![
                vec![0, 1, 8],
                vec![1, 2, 9],
                vec![2, 3, 10],
                vec![3, 0, 11],
                vec![4, 5, 12],
                vec![5, 6, 13],
                vec![6, 7, 14],
                vec![7, 4, 15],
                vec![0, 4, 16],
                vec![1, 5, 17],
                vec![2, 6, 18],
                vec![3, 7, 19],
            ],
            face_ids: vec![
                vec![0, 4, 7, 3, 16, 15, 19, 11],
                vec![1, 2, 6, 5, 9, 18, 13, 17],
                vec![0, 1, 5, 4, 8, 17, 12, 16],
                vec![2, 3, 7, 6, 10, 19, 14, 18],
                vec![0, 3, 2, 1, 11, 10, 9, 8],
                vec![4, 5, 6, 7, 12, 13, 14, 15],
            ],
        }
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

    fn get_edge_local_point_id(&self, e: usize, i: usize) -> usize {
        self.edge_ids[e][i]
    }

    fn get_face_local_point_id(&self, f: usize, i: usize) -> usize {
        self.face_ids[f][i]
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
        let geo = Hex20::new();
        assert_eq!(geo.coords.len(), NPOINT);
        assert_eq!(geo.interp.dim(), NPOINT);
        assert_eq!(geo.deriv.dims(), (NPOINT, NDIM));
    }
}
