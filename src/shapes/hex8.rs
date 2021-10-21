use russell_lab::{Matrix, Vector};

/// Defines a (bi-linear) Hexahedron with 8 nodes
///
/// ```text
///
///              4________________7
///            ,'|              ,'|
///          ,'  |            ,'  |
///        ,'    |          ,'    |
///      ,'      |        ,'      |
///    5'===============6'        |
///    |         |      |         |
///    |         |      |         |
///    |         0_____ | ________3
///    |       ,'       |       ,'
///    |     ,'         |     ,'
///    |   ,'           |   ,'
///    | ,'             | ,'
///    1________________2'
///
/// ```
pub struct Hex8 {
    coords: Matrix, // natural coordinates (nnode, ndim)
    shape: Vector,  // shape functions @ natural coordinate (nnode)
    deriv: Matrix,  // derivatives of shape functions w.r.t natural coordinate (nnode, ndim)
}

const NDIM: usize = 3;
const NNODE: usize = 8;

impl Hex8 {
    /// Creates a new object with pre-calculated shape fn and derivatives
    pub fn new() -> Self {
        Hex8 {
            #[rustfmt::skip]
            coords: Matrix::from(&[
                [-1.0, -1.0, -1.0],
                [ 1.0, -1.0, -1.0],
                [ 1.0,  1.0, -1.0],
                [-1.0,  1.0, -1.0],
                [-1.0, -1.0,  1.0],
                [ 1.0, -1.0,  1.0],
                [ 1.0,  1.0,  1.0],
                [-1.0,  1.0,  1.0],
            ]),
            shape: Vector::new(NNODE),
            deriv: Matrix::new(NNODE, NDIM),
        }
    }

    /// Calculates the shape functions at natural coordinate
    pub fn calc_shape(&mut self, coord: &Vector) {
        let (r, s, t) = (coord[0], coord[1], coord[2]);

        self.shape[0] = (1.0 - r - s + r * s - t + s * t + r * t - r * s * t) / 8.0;
        self.shape[1] = (1.0 + r - s - r * s - t + s * t - r * t + r * s * t) / 8.0;
        self.shape[2] = (1.0 + r + s + r * s - t - s * t - r * t - r * s * t) / 8.0;
        self.shape[3] = (1.0 - r + s - r * s - t - s * t + r * t + r * s * t) / 8.0;
        self.shape[4] = (1.0 - r - s + r * s + t - s * t - r * t + r * s * t) / 8.0;
        self.shape[5] = (1.0 + r - s - r * s + t - s * t + r * t - r * s * t) / 8.0;
        self.shape[6] = (1.0 + r + s + r * s + t + s * t + r * t + r * s * t) / 8.0;
        self.shape[7] = (1.0 - r + s - r * s + t + s * t - r * t - r * s * t) / 8.0;
    }

    /// Calculates the derivatives of shape functions at natural coordinate
    pub fn calc_deriv(&mut self, coord: &Vector) {
        let (r, s, t) = (coord[0], coord[1], coord[2]);

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

    /// Returns the previously computed shape function for node n
    pub fn get_shape(&self, n: usize) -> f64 {
        self.shape[n]
    }

    /// Returns the previously computed derivative of shape function for node n
    pub fn get_deriv(&self, deriv: &mut Vector, n: usize) {
        deriv[0] = self.deriv[n][0];
        deriv[1] = self.deriv[n][1];
        deriv[2] = self.deriv[n][2];
    }

    /// Returns the number of dimensions
    pub fn get_ndim(&self) -> usize {
        NDIM
    }

    /// Returns the number of nodes
    pub fn get_nnode(&self) -> usize {
        NNODE
    }

    /// Returns the natural coordinates @ node n
    pub fn get_coord(&self, coord: &mut Vector, n: usize) {
        coord[0] = self.coords[n][0];
        coord[1] = self.coords[n][1];
        coord[2] = self.coords[n][2];
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_works() {
        let geo = Hex8::new();
        assert_eq!(geo.coords.dims(), (NNODE, NDIM));
        assert_eq!(geo.shape.dim(), NNODE);
        assert_eq!(geo.deriv.dims(), (NNODE, NDIM));
    }
}
