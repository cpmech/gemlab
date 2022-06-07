#[cfg(test)]
pub mod aux {
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::util::PI;
    use plotpy::{Canvas, Curve, Plot};
    use russell_lab::generate2d;
    use russell_lab::Vector;

    const RMIN: f64 = 5.0;
    const RMAX: f64 = 10.0;
    const AMIN: f64 = 30.0 * PI / 180.0;
    const AMAX: f64 = 60.0 * PI / 180.0;
    const ZMIN: f64 = 0.0;
    const ZMAX: f64 = 4.0;

    /// Maps point coordinates
    ///
    /// The shape is the area indicated with "?" or the edge with "%".
    /// If class == Tri, the shape is half of the highlighted wedge.
    /// In 3D, an extrusion is applied along the out-of-plane direction.
    ///
    /// ```text
    ///   |            /
    ///   |           / αmax
    ///   ***=---__  /
    ///   |         % _
    ///   |        % ? *._          ,
    ///   |       % ????? *.     ,-'
    ///   ***=-_ % ???????? *.,-' αmin
    ///   |     % - ?????? ,-'*
    ///   |    /    *.? ,-'    *
    ///   |   /      ,*'        *
    ///   |  /    ,-'  *         *
    ///   | /  ,-'      *         *
    ///   |/.-'         #         #
    ///   o ----------- # ------- # --> r
    ///               rmin       rmax
    /// ```
    ///
    /// Intermediary mapping:
    ///
    /// r(ξ₀,ξ₁,ξ₂) = rmin + (ξ₀ - ξ₀min) · Δr / Δξ₀
    /// α(ξ₀,ξ₁,ξ₂) = αmin + (ξ₁ - ξ₁min) · Δα / Δξ₁
    /// z(ξ₀,ξ₁,ξ₂) = ξ₂
    ///
    /// Cylindrical coordinates:
    ///
    /// x₀ := r · cos(α)
    /// x₁ := r · sin(α)
    /// x₂ := z
    pub fn map_point_coords(x: &mut Vector, ksi: &[f64], ksi_min: f64, ksi_del: f64) {
        assert_eq!(x.dim(), ksi.len());
        let r = RMIN + (ksi[0] - ksi_min) * (RMAX - RMIN) / ksi_del;
        let a = AMIN + (ksi[1] - ksi_min) * (AMAX - AMIN) / ksi_del;
        x[0] = r * f64::cos(a);
        x[1] = r * f64::sin(a);
        if x.dim() == 3 {
            x[2] = ZMIN + (ksi[2] - ksi_min) * (ZMAX - ZMIN) / ksi_del;
        }
    }

    /// Generates the Canvas for the mapping used in tests
    pub fn _gen_canvas_mapping() -> Canvas {
        let mut canvas = Canvas::new();
        let color = "#bfbfbf";
        canvas
            .set_face_color("None")
            .set_edge_color(color)
            .set_line_width(2.0)
            .draw_circle(0.0, 0.0, RMIN);
        canvas.draw_circle(0.0, 0.0, RMAX);
        canvas.set_edge_color(color).set_line_width(2.0);
        canvas.draw_polyline(&[[0.0, 0.0], [RMAX * f64::cos(AMIN), RMAX * f64::sin(AMIN)]], false);
        canvas.draw_polyline(&[[0.0, 0.0], [RMAX * f64::cos(AMAX), RMAX * f64::sin(AMAX)]], false);
        canvas
    }

    /// Draws the points in the natural (right) and real (left) spaces
    pub fn _draw_point_coords_2d(ksi_min: f64, ksi_del: f64) -> Plot {
        const N: usize = 11;
        let mut natural_space = Curve::new();
        natural_space.set_line_style("None").set_marker_style("o");
        let mut real_space = Curve::new();
        real_space.set_line_style("None").set_marker_style("o");
        let ksi_max = ksi_min + ksi_del;
        let (ksi_0, ksi_1) = generate2d(ksi_min, ksi_max, ksi_min, ksi_max, N, N);
        let mut x = Vector::new(2);
        natural_space.points_begin();
        real_space.points_begin();
        for i in 0..N {
            for j in 0..N {
                natural_space.points_add(ksi_0[i][j], ksi_1[i][j]);
                map_point_coords(&mut x, &[ksi_0[i][j], ksi_1[i][j]], ksi_min, ksi_del);
                real_space.points_add(x[0], x[1]);
            }
        }
        natural_space.points_end();
        real_space.points_end();
        let mut plot = Plot::new();
        plot.set_subplot(1, 2, 1)
            .add(&real_space)
            .add(&_gen_canvas_mapping())
            .set_equal_axes(true)
            .set_range(-0.1, RMAX + 0.1, -0.1, RMAX + 0.1)
            .set_ticks_x(1.0, 0.0, "")
            .set_ticks_y(1.0, 0.0, "")
            .grid_and_labels("x", "y");
        plot.set_subplot(1, 2, 2)
            .add(&natural_space)
            .set_equal_axes(true)
            .set_labels("ξ0", "ξ1");
        plot
    }

    /// Returns a new scratchpad with coordinates set
    pub fn gen_scratchpad_with_coords(space_ndim: usize, kind: GeoKind) -> Scratchpad {
        let geo_ndim = kind.ndim();
        let nnode = kind.nnode();
        let mut x = Vector::new(space_ndim);
        let mut ksi_aux = vec![0.0; space_ndim];
        let mut pad = Scratchpad::new(space_ndim, kind).unwrap();
        let (ksi_min, ksi_del) = kind.ksi_min_ksi_del();
        for m in 0..nnode {
            let ksi = kind.reference_coords(m);
            if geo_ndim == space_ndim {
                // 2D or 3D
                map_point_coords(&mut x, ksi, ksi_min, ksi_del);
            } else if geo_ndim == 1 && space_ndim == 2 {
                // line in 2D
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = 1.0;
                map_point_coords(&mut x, &ksi_aux, ksi_min, ksi_del);
            } else {
                // line or face in 3D
                ksi_aux[0] = ksi[0];
                ksi_aux[1] = ksi[1];
                ksi_aux[2] = 1.0;
                map_point_coords(&mut x, &ksi_aux, ksi_min, ksi_del);
            }
            for j in 0..space_ndim {
                pad.set_xx(m, j, x[j])
            }
        }
        pad
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::aux::_draw_point_coords_2d;
    use crate::StrError;

    // #[test]
    fn _draw_point_coords_2d_works() -> Result<(), StrError> {
        let mut plot = _draw_point_coords_2d(0.0, 1.0);
        plot.set_figure_size_points(1000.0, 600.0)
            .save("/tmp/gemlab/test_draw_point_coords_2d.svg")
    }
}
