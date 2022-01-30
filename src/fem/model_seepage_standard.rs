use super::ModelSeepage;

pub struct ModelSeepageStandard {
    kx: f64,
    ky: f64,
    kz: f64,
}

impl ModelSeepageStandard {
    pub fn new(kx: f64, ky: f64, kz: f64) -> Self {
        ModelSeepageStandard { kx, ky, kz }
    }
}

impl ModelSeepage for ModelSeepageStandard {}
