use russell_lab::Matrix;

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

    /// Performs the nsn integration
    ///
    /// From [@bhatti:05]\page{348}
    /// [@bhatti:05] Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
    #[rustfmt::skip]
    pub fn integ_nsn(&self, s: f64, th: f64) -> Matrix {
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

    /// Performs the gtg integration
    ///
    /// From [@bhatti:05]\page{348}
    /// [@bhatti:05] Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
    #[rustfmt::skip]
    pub fn integ_gtg(&self, kx: f64, ky: f64) -> Matrix {
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
