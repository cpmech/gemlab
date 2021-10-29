use super::Trait;
use russell_lab::Vector;

pub struct Qua4 {
    //
}

const NDIM: usize = 2;
const NPOINT: usize = 4;
const NEDGE: usize = 4;
const NFACE: usize = 0;

impl Qua4 {
    pub fn new() -> Self {
        Qua4 {}
    }
}

impl Trait for Qua4 {
    fn calc_interp(&mut self, coord: &Vector) {}

    fn calc_deriv(&mut self, coord: &Vector) {}

    fn get_interp(&self, m: usize) -> f64 {
        0.0
    }

    fn get_deriv(&self, deriv: &mut Vector, m: usize) {}

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

    fn get_ksi(&self, ksi: &mut Vector, m: usize) {}
}
