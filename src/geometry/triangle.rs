/// Calculates the triangle (three) coordinates
///
/// # Output
///
/// * `zeta` -- the triangle coordinates (len = 3)
///
/// # Input
///
/// * `xa,xb,xc` -- corners (each has len = 2)
/// * `xp` -- the point to calculate the coordinates (len = 2)
///
/// # Panics
///
/// This function will panic if the array sizes are incorrect
pub fn triangle_coords(zeta: &mut [f64], xa: &[f64], xb: &[f64], xc: &[f64], xp: &[f64]) {
    let a2 = xa[0] * (xb[1] - xc[1]) + xb[0] * (xc[1] - xa[1]) + xc[0] * (xa[1] - xb[1]);
    zeta[0] = (xp[0] * (xb[1] - xc[1]) + xb[0] * (xc[1] - xp[1]) + xc[0] * (xp[1] - xb[1])) / a2;
    zeta[1] = (xa[0] * (xp[1] - xc[1]) + xp[0] * (xc[1] - xa[1]) + xc[0] * (xa[1] - xp[1])) / a2;
    zeta[2] = (xa[0] * (xb[1] - xp[1]) + xb[0] * (xp[1] - xa[1]) + xp[0] * (xa[1] - xb[1])) / a2;
}

/// Indicates if a point is inside a triangle by looking at its triangle coordinates (zeta)
///
/// Note: the point is inside (or on an edge) if all zeta are positive (or zero)
///
/// # Input
///
/// * `zeta` -- the triangle coordinates (len = 3)
///
/// # Panics
///
/// This function will panic if the array sizes are incorrect
#[inline]
pub fn in_triangle(zeta: &[f64]) -> bool {
    if zeta[0] < 0.0 || zeta[1] < 0.0 || zeta[2] < 0.0 {
        false
    } else {
        true
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{in_triangle, triangle_coords};
    use plotpy::{Canvas, Plot, PolyCode, Text};

    fn draw_triangles(plot: &mut Plot, triangles: &[[[f64; 2]; 3]]) {
        let mut ids = Text::new();
        let mut canvas = Canvas::new();
        canvas.set_face_color("#fefddc80").set_edge_color("#fcb827");
        for t in 0..triangles.len() {
            let (mut xc, mut yc) = (0.0, 0.0);
            canvas.polycurve_begin();
            for m in 0..3 {
                let code = if m == 0 { PolyCode::MoveTo } else { PolyCode::LineTo };
                canvas.polycurve_add(&triangles[t][m][0], &triangles[t][m][1], code);
                xc += triangles[t][m][0];
                yc += triangles[t][m][1];
            }
            canvas.polycurve_end(true);
            xc /= 3.0;
            yc /= 3.0;
            ids.draw(xc, yc, format!("{}", t).as_str());
        }
        plot.add(&canvas).add(&ids);
    }

    #[test]
    fn draw_triangles_works() {
        let mut plot = Plot::new();
        const T: [[[f64; 2]; 3]; 1] = [[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]];
        draw_triangles(&mut plot, &T);
        // if false {
        //     plot.set_equal_axes(true)
        //         .save("/tmp/gemlab/draw_triangles_works.svg")
        //         .unwrap();
        // }
    }

    #[test]
    fn triangle_coords_works() {
        let xa = &[0.0, 0.0];
        let xb = &[1.0, 0.0];
        let xc = &[0.0, 1.0];
        let mut zeta = vec![0.0; 3];
        triangle_coords(&mut zeta, xa, xb, xc, &[0.0, 0.0]);
        assert_eq!(zeta, &[1.0, 0.0, 0.0]);
        triangle_coords(&mut zeta, xa, xb, xc, &[1.0, 0.0]);
        assert_eq!(zeta, &[0.0, 1.0, 0.0]);
        triangle_coords(&mut zeta, xa, xb, xc, &[0.0, 1.0]);
        assert_eq!(zeta, &[0.0, 0.0, 1.0]);

        triangle_coords(&mut zeta, xa, xb, xc, &[0.5, 0.5]);
        assert_eq!(zeta, &[0.0, 0.5, 0.5]);
        triangle_coords(&mut zeta, xa, xb, xc, &[1.0, 1.0]);
        assert_eq!(zeta, &[-1.0, 1.0, 1.0]);
    }

    #[test]
    #[rustfmt::skip]
    fn in_triangle_works_1() {
        let xa = &[0.0, 0.0];
        let xb = &[1.0, 0.0];
        let xc = &[0.0, 1.0];
        let mut zeta = vec![0.0; 3];
        triangle_coords(&mut zeta, xa, xb, xc, &[0.0, 0.0]); assert!(in_triangle(&zeta));
        triangle_coords(&mut zeta, xa, xb, xc, &[1.0, 0.0]); assert!(in_triangle(&zeta));
        triangle_coords(&mut zeta, xa, xb, xc, &[0.0, 1.0]); assert!(in_triangle(&zeta));
        triangle_coords(&mut zeta, xa, xb, xc, &[0.5, 0.5]); assert!(in_triangle(&zeta));
        triangle_coords(&mut zeta, xa, xb, xc, &[0.5, 0.0]); assert!(in_triangle(&zeta));
        triangle_coords(&mut zeta, xa, xb, xc, &[0.0, 0.5]); assert!(in_triangle(&zeta));
        triangle_coords(&mut zeta, xa, xb, xc, &[0.3, 0.3]); assert!(in_triangle(&zeta));
        triangle_coords(&mut zeta, xa, xb, xc, &[1e-15, 1e-15]); assert!(in_triangle(&zeta));
        triangle_coords(&mut zeta, xa, xb, xc, &[1e-15, -1e-15]); assert!(!in_triangle(&zeta));
        triangle_coords(&mut zeta, xa, xb, xc, &[1.1, 0.0]); assert!(!in_triangle(&zeta));
        triangle_coords(&mut zeta, xa, xb, xc, &[0.0, 1.1]); assert!(!in_triangle(&zeta));
        triangle_coords(&mut zeta, xa, xb, xc, &[1.0, 1.0]); assert!(!in_triangle(&zeta));
        triangle_coords(&mut zeta, xa, xb, xc, &[-1e-15, -1e-15]); assert!(!in_triangle(&zeta));
    }

    #[test]
    fn in_triangle_works_2() {
        #[rustfmt::skip]
        const T: [[[f64; 2]; 3]; 12] = [
            [[0.230951,  0.558482], [0.133721,  0.348832],   [0.540745,  0.331184]],   //  0
            [[0.13928,   0.180603], [0.133721,  0.348832],   [0.0307942, 0.459123]],   //  1
            [[0.0307942, 0.459123], [0.230951,  0.558482],   [0.0980015, 0.981755]],   //  2
            [[0.230951,  0.558482], [0.0307942, 0.459123],   [0.133721,  0.348832]],   //  3
            [[0.0980015, 0.981755], [0.230951,  0.558482],   [0.578587,  0.760349]],   //  4
            [[0.133721,  0.348832], [0.13928,   0.180603],   [0.540745,  0.331184]],   //  5
            [[0.540745,  0.331184], [0.578587,  0.760349],   [0.230951,  0.558482]],   //  6
            [[0.540745,  0.331184], [0.478554,  0.00869692], [0.648071,  0.369534]],   //  7
            [[0.578587,  0.760349], [0.648071,  0.369534],   [0.903726,  0.975904]],   //  8
            [[0.648071,  0.369534], [0.578587,  0.760349],   [0.540745,  0.331184]],   //  9
            [[0.578587,  0.760349], [0.903726,  0.975904],   [0.0980015, 0.981755]],   // 10
            [[0.540745,  0.331184], [0.13928,   0.180603],   [0.478554,  0.00869692]], // 11
        ];
        let mut zeta = vec![0.0; 3];
        for tri in &T {
            for m in 0..3 {
                triangle_coords(&mut zeta, &tri[0], &tri[1], &tri[2], &tri[m]);
                assert!(in_triangle(&zeta));
            }
            triangle_coords(&mut zeta, &tri[0], &tri[1], &tri[2], &[0.1, 0.1]);
            assert!(!in_triangle(&zeta));
            triangle_coords(&mut zeta, &tri[0], &tri[1], &tri[2], &[0.6, 0.2]);
            assert!(!in_triangle(&zeta));
        }
        triangle_coords(&mut zeta, &T[0][0], &T[0][1], &T[0][2], &[0.3, 0.4]);
        assert!(in_triangle(&zeta));
        triangle_coords(&mut zeta, &T[8][0], &T[8][1], &T[8][2], &[0.7, 0.7]);
        assert!(in_triangle(&zeta));
        // if false {
        //     let mut plot = Plot::new();
        //     draw_triangles(&mut plot, &T);
        //     plot.grid_and_labels("x", "y")
        //         .set_range(0.0, 1.0, 0.0, 1.0)
        //         .set_ticks_x(0.1, 0.5, "")
        //         .set_ticks_y(0.1, 0.5, "")
        //         .set_equal_axes(true)
        //         .set_figure_size_points(600.0, 600.0)
        //         .save("/tmp/gemlab/test_is_point_inside_triangle_2.svg")
        //         .unwrap();
        // }
    }
}
