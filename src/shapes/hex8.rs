use russell_lab::{Matrix, Vector};

/// Defines a (bi-linear) Hexahedron with 8 nodes
///
/// The natural coordinates range from -1 to +1 with the geometry centred @ 0
///
/// ```text
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
/// ```
///
/// # Example
///
/// ```
/// // import
/// use gemlab::Hex8;
/// use russell_lab::Vector;
///
/// // geometry object
/// let mut geo = Hex8::new();
/// let ndim = geo.get_ndim();
///
/// // compute shape fn and deriv @ ksi=(0,0,0) for all nodes
/// let ksi = Vector::new(ndim);
/// geo.calc_shape(&ksi);
/// geo.calc_deriv(&ksi);
///
/// // get shape fn and deriv for node m
/// let m = 0;
/// let sm = geo.get_shape(m);
/// let mut dsm_dksi = Vector::new(ndim);
/// geo.get_deriv(&mut dsm_dksi, m);
///
/// // check shape fn and deriv @ ksi for node m
/// assert_eq!(sm, 0.125);
/// assert_eq!(
///     format!("{}", dsm_dksi),
///     "┌        ┐\n\
///      │ -0.125 │\n\
///      │ -0.125 │\n\
///      │ -0.125 │\n\
///      └        ┘"
/// );
/// ```
pub struct Hex8 {
    // natural coordinates (nnode, ndim)
    //
    // ```text
    // ξᵐ = vector{r, s, t} @ node m
    // coords = [ξ⁰, ξ¹, ξ², ξ³, ξ⁴, ξ⁵, ξ⁶, ξ⁷, ξ⁸]
    // ```
    coords: Vec<Vector>,

    // shape functions @ natural coordinate (nnode)
    //
    // ```text
    // shape[m](ξ) = Sᵐ(ξ)
    // ```
    shape: Vector,

    // derivatives of shape functions w.r.t natural coordinate (nnode, ndim)
    //
    // ```text
    // deriv[m][i](ξ) = ({dSᵐ(ξ)/dξ}_ξ)[i]
    // ```
    deriv: Matrix,
}

const NDIM: usize = 3;
const NNODE: usize = 8;

impl Hex8 {
    /// Creates a new object with pre-calculated shape fn and derivatives
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
            shape: Vector::new(NNODE),
            deriv: Matrix::new(NNODE, NDIM),
        }
    }

    /// Calculates the shape functions at natural coordinate
    ///
    /// ```text
    /// shape[m](ξ) = Sᵐ(ξ)
    /// ```
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
    ///
    /// ```text
    /// deriv[m][i](ξ) = ({dSᵐ(ξ)/dξ}_ξ)[i]
    /// ```
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

    /// Returns the previously computed shape function for node m
    ///
    /// ```text
    /// shape[m](ξ) = Sᵐ(ξ)
    /// ```
    pub fn get_shape(&self, m: usize) -> f64 {
        self.shape[m]
    }

    /// Returns the previously computed derivative of shape function for node m
    //
    // ```text
    // deriv[m][i](ξ) = ({dSᵐ(ξ)/dξ}_ξ)[i]
    // ```
    pub fn get_deriv(&self, deriv: &mut Vector, m: usize) {
        deriv[0] = self.deriv[m][0];
        deriv[1] = self.deriv[m][1];
        deriv[2] = self.deriv[m][2];
    }

    /// Returns the number of dimensions
    pub fn get_ndim(&self) -> usize {
        NDIM
    }

    /// Returns the number of nodes
    pub fn get_nnode(&self) -> usize {
        NNODE
    }

    /// Returns the natural coordinates @ node m
    //
    // ```text
    // coord = ξᵐ = vector{r, s, t} @ node m
    // ```
    pub fn get_coord(&self, coord: &mut Vector, m: usize) {
        coord[0] = self.coords[m][0];
        coord[1] = self.coords[m][1];
        coord[2] = self.coords[m][2];
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use russell_chk::*;

    // Holds arguments for numerical differentiation
    struct Arguments {
        geo: Hex8,      // geometry
        at_ksi: Vector, // at nat coord value
        ksi: Vector,    // temporary nat coord
        m: usize,       // shape: node index
        i: usize,       // nat coord: dimension index
    }

    // Computes the shape function @ m as a function of x = ksi[i]
    fn sm(x: f64, args: &mut Arguments) -> f64 {
        args.ksi[0] = args.at_ksi[0];
        args.ksi[1] = args.at_ksi[1];
        args.ksi[2] = args.at_ksi[2];
        args.ksi[args.i] = x;
        args.geo.calc_shape(&args.ksi);
        args.geo.get_shape(args.m)
    }

    #[test]
    fn new_works() {
        let geo = Hex8::new();
        assert_eq!(geo.coords.len(), NNODE);
        assert_eq!(geo.shape.dim(), NNODE);
        assert_eq!(geo.deriv.dims(), (NNODE, NDIM));
    }

    #[test]
    fn calc_shape_works() {
        let mut geo = Hex8::new();
        let mut ksi = Vector::new(NDIM);
        for m in 0..NNODE {
            geo.get_coord(&mut ksi, m);
            geo.calc_shape(&ksi);
            for n in 0..NNODE {
                let smn = geo.get_shape(n);
                if m == n {
                    assert_approx_eq!(smn, 1.0, 1e-15);
                } else {
                    assert_approx_eq!(smn, 0.0, 1e-15);
                }
            }
        }
    }

    #[test]
    fn calc_deriv_works() {
        let args = &mut Arguments {
            geo: Hex8::new(),
            at_ksi: Vector::from(&[0.25, 0.25, 0.25]),
            ksi: Vector::new(NDIM),
            m: 0,
            i: 0,
        };
        let mut geo = Hex8::new();
        let mut dsm_dksi = Vector::new(NDIM);
        geo.calc_deriv(&args.at_ksi);
        for m in 0..NNODE {
            geo.get_deriv(&mut dsm_dksi, m);
            args.m = m;
            for i in 0..NDIM {
                args.i = i;
                assert_deriv_approx_eq!(dsm_dksi[i], args.at_ksi[i], sm, args, 1e-13);
            }
        }
    }

    #[test]
    fn getters_work() {
        let geo = Hex8::new();
        assert_eq!(geo.get_ndim(), NDIM);
        assert_eq!(geo.get_nnode(), NNODE);
    }
}
