use russell_lab::Matrix;

/// Performs analytical integrations on a Qua4
pub struct AnalyticalQua4 {
    a: f64,
    b: f64,
}

impl AnalyticalQua4 {
    /// Allocates a new instance
    ///
    /// `a` and `b` are half side lengths
    pub fn new(a: f64, b: f64) -> Self {
        AnalyticalQua4 { a, b }
    }

    /// Performs the n-s-n integration with constant s(x) field
    ///
    /// From @bhatti:05\page{332}
    /// @bhatti:05 Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
    pub fn mat_01_nsn(&self, s: f64, th: f64) -> Matrix {
        let c = th * s * self.a * self.b / 9.0;
        Matrix::from(&[
            [c * 4.0, c * 2.0, c * 1.0, c * 2.0],
            [c * 2.0, c * 4.0, c * 2.0, c * 1.0],
            [c * 1.0, c * 2.0, c * 4.0, c * 2.0],
            [c * 2.0, c * 1.0, c * 2.0, c * 4.0],
        ])
    }

    /// Performs the g-t-g integration with constant (and diagonal) tensor field
    ///
    /// From @bhatti:05\page{332}
    /// @bhatti:05 Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
    #[rustfmt::skip]
    pub fn mat_03_gtg(&self, kx: f64, ky: f64) -> Matrix {
        let a = self.a;
        let b = self.b;
        Matrix::from(&[
            [ (b*kx)/(3.0*a) + (a*ky)/(3.0*b), -(b*kx)/(3.0*a) + (a*ky)/(6.0*b), -(b*kx)/(6.0*a) - (a*ky)/(6.0*b),  (b*kx)/(6.0*a) - (a*ky)/(3.0*b)],
            [-(b*kx)/(3.0*a) + (a*ky)/(6.0*b),  (b*kx)/(3.0*a) + (a*ky)/(3.0*b),  (b*kx)/(6.0*a) - (a*ky)/(3.0*b), -(b*kx)/(6.0*a) - (a*ky)/(6.0*b)],
            [-(b*kx)/(6.0*a) - (a*ky)/(6.0*b),  (b*kx)/(6.0*a) - (a*ky)/(3.0*b),  (b*kx)/(3.0*a) + (a*ky)/(3.0*b), -(b*kx)/(3.0*a) + (a*ky)/(6.0*b)],
            [ (b*kx)/(6.0*a) - (a*ky)/(3.0*b), -(b*kx)/(6.0*a) - (a*ky)/(6.0*b), -(b*kx)/(3.0*a) + (a*ky)/(6.0*b),  (b*kx)/(3.0*a) + (a*ky)/(3.0*b)],
        ])
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::AnalyticalQua4;

    #[test]
    fn new_works() {
        let ana = AnalyticalQua4::new(2.0, 3.0);
        assert_eq!(ana.a, 2.0);
        assert_eq!(ana.b, 3.0);
    }
}
