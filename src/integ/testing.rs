#[cfg(test)]
pub mod aux {
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::math::SQRT_2;

    // to test if variables are cleared before sum
    pub const NOISE: f64 = 1234.56;

    /// Generates pad Lin2 for tests
    pub fn gen_pad_lin2(l: f64) -> Scratchpad {
        let mut pad = Scratchpad::new(2, GeoKind::Lin2).unwrap();
        pad.set_xx(0, 0, 3.0);
        pad.set_xx(0, 1, 4.0);
        pad.set_xx(1, 0, 3.0 + l / SQRT_2);
        pad.set_xx(1, 1, 4.0 + l / SQRT_2);
        pad
    }

    /// Generates pad Tri3 for tests
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

    /// Generates pad Qua4 for tests
    pub fn gen_pad_qua4(xc: f64, yc: f64, a: f64, b: f64) -> Scratchpad {
        let mut pad = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        pad.set_xx(0, 0, xc - a);
        pad.set_xx(0, 1, yc - b);
        pad.set_xx(1, 0, xc + a);
        pad.set_xx(1, 1, yc - b);
        pad.set_xx(2, 0, xc + a);
        pad.set_xx(2, 1, yc + b);
        pad.set_xx(3, 0, xc - a);
        pad.set_xx(3, 1, yc + b);
        pad
    }

    /// Generates pad Qua8 for tests
    pub fn gen_pad_qua8(xc: f64, yc: f64, a: f64, b: f64) -> Scratchpad {
        let mut pad = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        pad.set_xx(0, 0, xc - a);
        pad.set_xx(0, 1, yc - b);
        pad.set_xx(1, 0, xc + a);
        pad.set_xx(1, 1, yc - b);
        pad.set_xx(2, 0, xc + a);
        pad.set_xx(2, 1, yc + b);
        pad.set_xx(3, 0, xc - a);
        pad.set_xx(3, 1, yc + b);
        pad.set_xx(4, 0, xc);
        pad.set_xx(4, 1, yc - b);
        pad.set_xx(5, 0, xc + a);
        pad.set_xx(5, 1, yc);
        pad.set_xx(6, 0, xc);
        pad.set_xx(6, 1, yc + b);
        pad.set_xx(7, 0, xc - a);
        pad.set_xx(7, 1, yc);
        pad
    }

    /// Generates pad Tet4 for tests
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

//
// NOTE: we don't need the functions below to be activated and making the coverage tool to go crazy.
// We only use these functions to create the figures that are already saved in data/figures.
//
#[cfg(test)]
mod tests {
    use super::aux;
    use plotpy::Plot;

    #[test]
    fn gen_pad_functions_figs() {
        let pad = aux::gen_pad_lin2(1.0);
        if false {
            let mut plot = Plot::new();
            pad.draw_shape(&mut plot, "blue", true, true).unwrap();
            plot.set_equal_axes(true)
                .grid_and_labels("x", "y")
                .save("/tmp/gemlab/test_gen_pad_functions_lin2.svg")
                .unwrap();
        }

        let pad = aux::gen_pad_tri3();
        if false {
            let mut plot = Plot::new();
            pad.draw_shape(&mut plot, "blue", true, true).unwrap();
            plot.set_equal_axes(true)
                .grid_and_labels("x", "y")
                .save("/tmp/gemlab/test_gen_pad_functions_tri3.svg")
                .unwrap();
        }

        let pad = aux::gen_pad_qua4(2.0, 1.0, 2.0, 1.5);
        if false {
            let mut plot = Plot::new();
            pad.draw_shape(&mut plot, "blue", true, true).unwrap();
            plot.set_equal_axes(true)
                .grid_and_labels("x", "y")
                .save("/tmp/gemlab/test_gen_pad_functions_qua4.svg")
                .unwrap();
        }

        let pad = aux::gen_pad_qua8(2.0, 1.0, 2.0, 1.5);
        if false {
            let mut plot = Plot::new();
            pad.draw_shape(&mut plot, "blue", true, true).unwrap();
            plot.set_equal_axes(true)
                .grid_and_labels("x", "y")
                .save("/tmp/gemlab/test_gen_pad_functions_qua8.svg")
                .unwrap();
        }

        let pad = aux::gen_pad_tet4();
        if false {
            let mut plot = Plot::new();
            pad.draw_shape(&mut plot, "blue", true, true).unwrap();
            plot.set_equal_axes(true)
                .save("/tmp/gemlab/test_gen_pad_functions_tet4.svg")
                .unwrap();
        }
    }
}
*/
