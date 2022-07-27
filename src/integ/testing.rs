#[cfg(test)]
pub mod aux {
    use crate::shapes::{GeoKind, Scratchpad};

    // to test if variables are cleared before sum
    pub const NOISE: f64 = 1234.56;

    // generates pad Lin2 for tests
    pub fn gen_pad_lin2(l: f64) -> Scratchpad {
        let mut pad = Scratchpad::new(2, GeoKind::Lin2).unwrap();
        pad.set_xx(0, 0, 3.0);
        pad.set_xx(0, 1, 4.0);
        pad.set_xx(1, 0, 3.0 + l);
        pad.set_xx(1, 1, 4.0);
        pad
    }

    // generates pad Tri3 for tests
    pub fn gen_pad_tri3() -> Scratchpad {
        let mut pad = Scratchpad::new(2, GeoKind::Tri3).unwrap();
        pad.set_xx(0, 0, 3.0);
        pad.set_xx(0, 1, 4.0);
        pad.set_xx(1, 0, 8.0);
        pad.set_xx(1, 1, 4.0);
        pad.set_xx(2, 0, 5.0);
        pad.set_xx(2, 1, 9.0);
        pad
    }

    // generates pad Tet4 for tests
    pub fn gen_pad_tet4() -> Scratchpad {
        let mut pad = Scratchpad::new(3, GeoKind::Tet4).unwrap();
        pad.set_xx(0, 0, 2.0);
        pad.set_xx(0, 1, 3.0);
        pad.set_xx(0, 2, 4.0);
        pad.set_xx(1, 0, 6.0);
        pad.set_xx(1, 1, 3.0);
        pad.set_xx(1, 2, 2.0);
        pad.set_xx(2, 0, 2.0);
        pad.set_xx(2, 1, 5.0);
        pad.set_xx(2, 2, 1.0);
        pad.set_xx(3, 0, 4.0);
        pad.set_xx(3, 1, 3.0);
        pad.set_xx(3, 2, 6.0);
        pad
    }
}
