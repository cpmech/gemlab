use russell_lab::Matrix;
use russell_tensor::Tensor2;

/// Performs analytical integrations on a Qua8
pub struct AnalyticalQua8 {
    a: f64,
    b: f64,
}

impl AnalyticalQua8 {
    /// Allocates a new instance
    ///
    /// `a` and `b` are half side lengths
    pub fn new(a: f64, b: f64) -> Self {
        AnalyticalQua8 { a, b }
    }

    /// Performs the n-s-n integration with constant s(x) field
    ///
    /// From @bhatti:05\page{348}
    /// @bhatti:05 Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
    #[rustfmt::skip]
    pub fn mat_01_nsn(&self, s: f64, th: f64) -> Matrix {
        let c = -th * s * self.a * self.b;
        Matrix::from(&[
            [(-(2.0/15.0))*c , (-(2.0/45.0))*c , (-(1.0/15.0))*c , (-(2.0/45.0))*c , (2.0*c)/15.0     , (8.0*c)/45.0     , (8.0*c)/45.0     , (2.0*c)/15.0]     , // bh1
            [(-(2.0/45.0))*c , (-(2.0/15.0))*c , (-(2.0/45.0))*c , (-(1.0/15.0))*c , (2.0*c)/15.0     , (2.0*c)/15.0     , (8.0*c)/45.0     , (8.0*c)/45.0]     , // bh3
            [(-(1.0/15.0))*c , (-(2.0/45.0))*c , (-(2.0/15.0))*c , (-(2.0/45.0))*c , (8.0*c)/45.0     , (2.0*c)/15.0     , (2.0*c)/15.0     , (8.0*c)/45.0]     , // bh5
            [(-(2.0/45.0))*c , (-(1.0/15.0))*c , (-(2.0/45.0))*c , (-(2.0/15.0))*c , (8.0*c)/45.0     , (8.0*c)/45.0     , (2.0*c)/15.0     , (2.0*c)/15.0]     , // bh7
            [(2.0*c)/15.0    , (2.0*c)/15.0    , (8.0*c)/45.0    , (8.0*c)/45.0    , (-(32.0/45.0))*c , (-(4.0/9.0))*c   , (-(16.0/45.0))*c , (-(4.0/9.0))*c]   , // bh2
            [(8.0*c)/45.0    , (2.0*c)/15.0    , (2.0*c)/15.0    , (8.0*c)/45.0    , (-(4.0/9.0))*c   , (-(32.0/45.0))*c , (-(4.0/9.0))*c   , (-(16.0/45.0))*c] , // bh4
            [(8.0*c)/45.0    , (8.0*c)/45.0    , (2.0*c)/15.0    , (2.0*c)/15.0    , (-(16.0/45.0))*c , (-(4.0/9.0))*c   , (-(32.0/45.0))*c , (-(4.0/9.0))*c]   , // bh6
            [(2.0*c)/15.0    , (8.0*c)/45.0    , (8.0*c)/45.0    , (2.0*c)/15.0    , (-(4.0/9.0))*c   , (-(16.0/45.0))*c , (-(4.0/9.0))*c   , (-(32.0/45.0))*c] , // bh8
            //            bh1               bh3               bh5               bh7                bh2                bh4                bh6                bh8
        ])
    }

    /// Performs the g-t-g integration with constant (and diagonal) stress field
    ///
    /// From @bhatti:05\page{348}
    /// @bhatti:05 Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
    #[rustfmt::skip]
    pub fn mat_03_btb(&self, kx: f64, ky: f64) -> Matrix {
        let a = self.a;
        let b = self.b;
        Matrix::from(&[
            [(26.0*b*kx)/(45.0*a) + (26.0*a*ky)/(45.0*b) , (14.0*b*kx)/(45.0*a) + (17.0*a*ky)/(90.0*b) , (23.0*b*kx)/(90.0*a) + (23.0*a*ky)/(90.0*b) , (17.0*b*kx)/(90.0*a) + (14.0*a*ky)/(45.0*b) , (-8.0*b*kx)/(9.0*a) + (a*ky)/(15.0*b)     , -(b*kx)/(15.0*a) - (4.0*a*ky)/(9.0*b)     , (-4.0*b*kx)/(9.0*a) - (a*ky)/(15.0*b)     , (b*kx)/(15.0*a) - (8.0*a*ky)/(9.0*b)]      , // bh1
            [(14.0*b*kx)/(45.0*a) + (17.0*a*ky)/(90.0*b) , (26.0*b*kx)/(45.0*a) + (26.0*a*ky)/(45.0*b) , (17.0*b*kx)/(90.0*a) + (14.0*a*ky)/(45.0*b) , (23.0*b*kx)/(90.0*a) + (23.0*a*ky)/(90.0*b) , (-8.0*b*kx)/(9.0*a) + (a*ky)/(15.0*b)     , (b*kx)/(15.0*a) - (8.0*a*ky)/(9.0*b)      , (-4.0*b*kx)/(9.0*a) - (a*ky)/(15.0*b)     , -(b*kx)/(15.0*a) - (4.0*a*ky)/(9.0*b)]     , // bh3
            [(23.0*b*kx)/(90.0*a) + (23.0*a*ky)/(90.0*b) , (17.0*b*kx)/(90.0*a) + (14.0*a*ky)/(45.0*b) , (26.0*b*kx)/(45.0*a) + (26.0*a*ky)/(45.0*b) , (14.0*b*kx)/(45.0*a) + (17.0*a*ky)/(90.0*b) , (-4.0*b*kx)/(9.0*a) - (a*ky)/(15.0*b)     , (b*kx)/(15.0*a) - (8.0*a*ky)/(9.0*b)      , (-8.0*b*kx)/(9.0*a) + (a*ky)/(15.0*b)     , -(b*kx)/(15.0*a) - (4.0*a*ky)/(9.0*b)]     , // bh5
            [(17.0*b*kx)/(90.0*a) + (14.0*a*ky)/(45.0*b) , (23.0*b*kx)/(90.0*a) + (23.0*a*ky)/(90.0*b) , (14.0*b*kx)/(45.0*a) + (17.0*a*ky)/(90.0*b) , (26.0*b*kx)/(45.0*a) + (26.0*a*ky)/(45.0*b) , (-4.0*b*kx)/(9.0*a) - (a*ky)/(15.0*b)     , -(b*kx)/(15.0*a) - (4.0*a*ky)/(9.0*b)     , (-8.0*b*kx)/(9.0*a) + (a*ky)/(15.0*b)     , (b*kx)/(15.0*a) - (8.0*a*ky)/(9.0*b)]      , // bh7
            [(-8.0*b*kx)/(9.0*a) + (a*ky)/(15.0*b)       , (-8.0*b*kx)/(9.0*a) + (a*ky)/(15.0*b)       , (-4.0*b*kx)/(9.0*a) - (a*ky)/(15.0*b)       , (-4.0*b*kx)/(9.0*a) - (a*ky)/(15.0*b)       , (16.0*b*kx)/(9.0*a) + (8.0*a*ky)/(15.0*b) , 0.0                                       , (8.0*b*kx)/(9.0*a) - (8.0*a*ky)/(15.0*b)  , 0.0]                                       , // bh2
            [-(b*kx)/(15.0*a) - (4.0*a*ky)/(9.0*b)       , (b*kx)/(15.0*a) - (8.0*a*ky)/(9.0*b)        , (b*kx)/(15.0*a) - (8.0*a*ky)/(9.0*b)        , -(b*kx)/(15.0*a) - (4.0*a*ky)/(9.0*b)       , 0.0                                       , (8.0*b*kx)/(15.0*a) + (16.0*a*ky)/(9.0*b) , 0.0                                       , (-8.0*b*kx)/(15.0*a) + (8.0*a*ky)/(9.0*b)] , // bh4
            [(-4.0*b*kx)/(9.0*a) - (a*ky)/(15.0*b)       , (-4.0*b*kx)/(9.0*a) - (a*ky)/(15.0*b)       , (-8.0*b*kx)/(9.0*a) + (a*ky)/(15.0*b)       , (-8.0*b*kx)/(9.0*a) + (a*ky)/(15.0*b)       , (8.0*b*kx)/(9.0*a) - (8.0*a*ky)/(15.0*b)  , 0.0                                       , (16.0*b*kx)/(9.0*a) + (8.0*a*ky)/(15.0*b) , 0.0]                                       , // bh6
            [(b*kx)/(15.0*a) - (8.0*a*ky)/(9.0*b)        , -(b*kx)/(15.0*a) - (4.0*a*ky)/(9.0*b)       , -(b*kx)/(15.0*a) - (4.0*a*ky)/(9.0*b)       , (b*kx)/(15.0*a) - (8.0*a*ky)/(9.0*b)        , 0.0                                       , (-8.0*b*kx)/(15.0*a) + (8.0*a*ky)/(9.0*b) , 0.0                                       , (8.0*b*kx)/(15.0*a) + (16.0*a*ky)/(9.0*b)] , // bh8
            //                                        bh1                                           bh3                                           bh5                                           bh7                                         bh2                                         bh4                                         bh6                                          bh8
        ])
    }

    /// Performs the coupled n-s-g integration with constant s(x) field
    #[rustfmt::skip]
    pub fn mat_04_nsb(&self, s: f64) -> Matrix {
        let a = self.a;
        let b = self.b;
        let c = s / 18.0;
        Matrix::from(&[
            [-7.0*b*c, -7.0*a*c,     -b*c, -2.0*a*c, -2.0*b*c, -2.0*a*c, -2.0*b*c,     -a*c,  8.0*b*c, -6.0*a*c, 6.0*b*c,  4.0*a*c,  4.0*b*c, 6.0*a*c, -6.0*b*c,  8.0*a*c],
            [     b*c, -2.0*a*c,  7.0*b*c, -7.0*a*c,  2.0*b*c,     -a*c,  2.0*b*c, -2.0*a*c, -8.0*b*c, -6.0*a*c, 6.0*b*c,  8.0*a*c, -4.0*b*c, 6.0*a*c, -6.0*b*c,  4.0*a*c],
            [ 2.0*b*c,  2.0*a*c,  2.0*b*c,      a*c,  7.0*b*c,  7.0*a*c,      b*c,  2.0*a*c, -4.0*b*c, -6.0*a*c, 6.0*b*c, -8.0*a*c, -8.0*b*c, 6.0*a*c, -6.0*b*c, -4.0*a*c],
            [-2.0*b*c,      a*c, -2.0*b*c,  2.0*a*c,     -b*c,  2.0*a*c, -7.0*b*c,  7.0*a*c,  4.0*b*c, -6.0*a*c, 6.0*b*c, -4.0*a*c,  8.0*b*c, 6.0*a*c, -6.0*b*c, -8.0*a*c],
        ])
    }

    /// Performs the coupled g-t-n integration with constant tensor field
    #[rustfmt::skip]
    pub fn mat_05_btn(&self, tt: &Tensor2) -> Matrix {
        let a = self.a;
        let b = self.b;
        let (t00,t01) = (tt.get(0,0), tt.get(0,1));
        let (t10,t11) = (tt.get(1,0), tt.get(1,1));
        Matrix::from(&[
            [(b*t00 + a*t10)/18.,(b*t01 + a*t11)/18.,(b*t00 + 2.*a*t10)/18.,(b*t01 + 2.*a*t11)/18.,(b*t00 + a*t10)/9.,(b*t01 + a*t11)/9.,(2.*b*t00 + a*t10)/18.,(2.*b*t01 + a*t11)/18., (-4.*b*t00 - 3.*a*t10)/9.,(-4.*b*t01 - 3.*a*t11)/9.,(-3.*b*t00 - 2.*a*t10)/9.,(-3.*b*t01 - 2.*a*t11)/9.,(-2.*b*t00 - 3.*a*t10)/9.,(-2.*b*t01 - 3.*a*t11)/9.,(-3.*b*t00 - 4.*a*t10)/9., (-3.*b*t01 - 4.*a*t11)/9.],
            [(-(b*t00) + 2.*a*t10)/18.,(-(b*t01) + 2.*a*t11)/18.,(-(b*t00) + a*t10)/18.,(-(b*t01) + a*t11)/18.,(-2.*b*t00 + a*t10)/18.,(-2.*b*t01 + a*t11)/18., (-(b*t00) + a*t10)/9.,(-(b*t01) + a*t11)/9.,(4.*b*t00 - 3.*a*t10)/9.,(4.*b*t01 - 3.*a*t11)/9.,(3.*b*t00 - 4.*a*t10)/9.,(3.*b*t01 - 4.*a*t11)/9.,(2.*b*t00 - 3.*a*t10)/9.,(2.*b*t01 - 3.*a*t11)/9., (3.*b*t00 - 2.*a*t10)/9.,(3.*b*t01 - 2.*a*t11)/9.],
            [(-(b*t00) - a*t10)/9.,(-(b*t01) - a*t11)/9.,(-2.*b*t00 - a*t10)/18.,(-2.*b*t01 - a*t11)/18.,(-(b*t00) - a*t10)/18.,(-(b*t01) - a*t11)/18., (-(b*t00) - 2.*a*t10)/18.,(-(b*t01) - 2.*a*t11)/18.,(2.*b*t00 + 3.*a*t10)/9.,(2.*b*t01 + 3.*a*t11)/9.,(3.*b*t00 + 4.*a*t10)/9.,(3.*b*t01 + 4.*a*t11)/9.,(4.*b*t00 + 3.*a*t10)/9.,(4.*b*t01 + 3.*a*t11)/9., (3.*b*t00 + 2.*a*t10)/9.,(3.*b*t01 + 2.*a*t11)/9.],
            [(2.*b*t00 - a*t10)/18.,(2.*b*t01 - a*t11)/18.,(b*t00 - a*t10)/9.,(b*t01 - a*t11)/9.,(b*t00 - 2.*a*t10)/18.,(b*t01 - 2.*a*t11)/18., (b*t00 - a*t10)/18.,(b*t01 - a*t11)/18.,(-2.*b*t00 + 3.*a*t10)/9.,(-2.*b*t01 + 3.*a*t11)/9.,(-3.*b*t00 + 2.*a*t10)/9.,(-3.*b*t01 + 2.*a*t11)/9.,(-4.*b*t00 + 3.*a*t10)/9.,(-4.*b*t01 + 3.*a*t11)/9., (-3.*b*t00 + 4.*a*t10)/9.,(-3.*b*t01 + 4.*a*t11)/9.],
        ])
    }

    /// Performs the coupled n-v-n integration with constant vector field
    #[rustfmt::skip]
    pub fn mat_06_nvn(&self, v0: f64, v1: f64) -> Matrix {
        let c = self.a * self.b / 9.0;
        Matrix::from(&[
            [     0.0,    -v0*c,    -v0*c,    -v0*c],
            [     0.0,    -v1*c,    -v1*c,    -v1*c],
            [   -v0*c,      0.0,    -v0*c,    -v0*c],
            [   -v1*c,      0.0,    -v1*c,    -v1*c],
            [   -v0*c,    -v0*c,      0.0,    -v0*c],
            [   -v1*c,    -v1*c,      0.0,    -v1*c],
            [   -v0*c,    -v0*c,    -v0*c,      0.0],
            [   -v1*c,    -v1*c,    -v1*c,      0.0],
            [4.0*v0*c, 4.0*v0*c, 2.0*v0*c, 2.0*v0*c],
            [4.0*v1*c, 4.0*v1*c, 2.0*v1*c, 2.0*v1*c],
            [2.0*v0*c, 4.0*v0*c, 4.0*v0*c, 2.0*v0*c],
            [2.0*v1*c, 4.0*v1*c, 4.0*v1*c, 2.0*v1*c],
            [2.0*v0*c, 2.0*v0*c, 4.0*v0*c, 4.0*v0*c],
            [2.0*v1*c, 2.0*v1*c, 4.0*v1*c, 4.0*v1*c],
            [4.0*v0*c, 2.0*v0*c, 2.0*v0*c, 4.0*v0*c],
            [4.0*v1*c, 2.0*v1*c, 2.0*v1*c, 4.0*v1*c],
        ])
    }

    /// Performs the coupled g-s-n integration with constant s(x) field
    #[rustfmt::skip]
    pub fn mat_07_bsn(&self, s: f64) -> Matrix {
        let a = self.a;
        let b = self.b;
        let c = s / 18.0;
        Matrix::from(&[
            [-7.0*b*c,      b*c,  2.0*b*c, -2.0*b*c],
            [-7.0*a*c, -2.0*a*c,  2.0*a*c,      a*c],
            [    -b*c,  7.0*b*c,  2.0*b*c, -2.0*b*c],
            [-2.0*a*c, -7.0*a*c,      a*c,  2.0*a*c],
            [-2.0*b*c,  2.0*b*c,  7.0*b*c,     -b*c],
            [-2.0*a*c,     -a*c,  7.0*a*c,  2.0*a*c],
            [-2.0*b*c,  2.0*b*c,      b*c, -7.0*b*c],
            [    -a*c, -2.0*a*c,  2.0*a*c,  7.0*a*c],
            [ 8.0*b*c, -8.0*b*c, -4.0*b*c,  4.0*b*c],
            [-6.0*a*c, -6.0*a*c, -6.0*a*c, -6.0*a*c],
            [ 6.0*b*c,  6.0*b*c,  6.0*b*c,  6.0*b*c],
            [ 4.0*a*c,  8.0*a*c, -8.0*a*c, -4.0*a*c],
            [ 4.0*b*c, -4.0*b*c, -8.0*b*c,  8.0*b*c],
            [ 6.0*a*c,  6.0*a*c,  6.0*a*c,  6.0*a*c],
            [-6.0*b*c, -6.0*b*c, -6.0*b*c, -6.0*b*c],
            [ 8.0*a*c,  4.0*a*c, -4.0*a*c, -8.0*a*c]
        ])
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::AnalyticalQua8;

    #[test]
    fn new_works() {
        let ana = AnalyticalQua8::new(2.0, 3.0);
        assert_eq!(ana.a, 2.0);
        assert_eq!(ana.b, 3.0);
    }
}
