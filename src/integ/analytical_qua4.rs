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

    /// Performs the nsn integration
    pub fn integ_nsn(&self, s: f64, th: f64) -> Matrix {
        let c = th * s * self.a * self.b / 9.0;
        Matrix::from(&[
            [c * 4.0, c * 2.0, c * 1.0, c * 2.0],
            [c * 2.0, c * 4.0, c * 2.0, c * 1.0],
            [c * 1.0, c * 2.0, c * 4.0, c * 2.0],
            [c * 2.0, c * 1.0, c * 2.0, c * 4.0],
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
