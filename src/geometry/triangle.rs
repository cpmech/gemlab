/// Returns true if the sign of the cross product vector is negative, including -0.0
///
/// Specifically, returns true if the sign of the component of the out-of-plane vector,
/// resulting from the cross product between u and v is negative
#[inline]
fn cross_is_negative(u0: f64, u1: f64, v0: f64, v1: f64) -> bool {
    f64::is_sign_negative(u0 * v1 - u1 * v0)
}

/// Checks whether a point is inside a triangle or not
///
/// ```text
/// →   →    →    →   →    →    →   →    →
/// a = AB × AP;  b = BC × BP;  c = CA × CP
///
/// The point is inside iff:
///     sign(a[2]) == sign(b[2]) == sign(c[2])
///
/// A,                      A, -------P   OUTSIDE
/// |\',                    | ',     / |
/// | \ ',                  |   ',  /    .  
/// |  \  ',                |     ',      .   
/// |   \   ',              |     / ',     .
/// |    P.   ',            |    /    ',    |
/// |   /  '.   ',          |   /       ',   .
/// |  /     `-.  ',        |  /          ',  |
/// | /  INSIDE `-. ',      | /             ', .
/// |/             `-.',    |/                ',|
/// C-------------------B   C-------------------B
/// ```
///
/// **Note:** This function returns true if the point is
/// on any boundary or coincides with any vertex
pub fn is_point_inside_triangle(xa: &[f64], xb: &[f64], xc: &[f64], xp: &[f64]) -> bool {
    if (xp[0] == xa[0] && xp[1] == xa[1]) || (xp[0] == xb[0] && xp[1] == xb[1]) || (xp[0] == xc[0] && xp[1] == xc[1]) {
        return true;
    }
    let na = cross_is_negative(xb[0] - xa[0], xb[1] - xa[1], xp[0] - xa[0], xp[1] - xa[1]);
    let nb = cross_is_negative(xc[0] - xb[0], xc[1] - xb[1], xp[0] - xb[0], xp[1] - xb[1]);
    let nc = cross_is_negative(xa[0] - xc[0], xa[1] - xc[1], xp[0] - xc[0], xp[1] - xc[1]);
    (na && nb && nc) || (!na && !nb && !nc)
}

/// Computes the signed area of a triangle given its vertices
/// 
/// The sign is positive if the vertices are given in counter-clockwise order.
/// Otherwise, the area is negative (clockwise order).
#[rustfmt::skip]
pub fn triangle_signed_area(xa: &[f64], xb: &[f64], xc: &[f64]) -> f64 {
    (   xa[0] * (xb[1] - xc[1])
      + xb[0] * (xc[1] - xa[1])
      + xc[0] * (xa[1] - xb[1])
    ) / 2.0
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{is_point_inside_triangle, triangle_signed_area};
    use crate::StrError;
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
    #[rustfmt::skip]
    fn is_point_inside_triangle_works_1() {
        assert!(is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[0.0, 0.0]));
        assert!(is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[1.0, 0.0]));
        assert!(is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[0.0, 1.0]));
        assert!(is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[0.5, 0.5]));
        assert!(is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[0.5, 0.0]));
        assert!(is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[0.0, 0.5]));
        assert!(is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[0.3, 0.3]));
        assert!(is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[1e-15, 1e-15]));
        // assert!(is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[1e-15, -1e-15])); // imprecision is not accepted
        assert!(!is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[1.1, 0.0]));
        assert!(!is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[0.0, 1.1]));
        assert!(!is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[1.0, 1.0]));
        assert!(!is_point_inside_triangle(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0], &[-1e-15, -1e-15]));
    }

    #[test]
    fn is_point_inside_triangle_works_2() -> Result<(), StrError> {
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
        for tri in &T {
            for m in 0..3 {
                assert!(is_point_inside_triangle(&tri[0], &tri[1], &tri[2], &tri[m]));
            }
            assert!(!is_point_inside_triangle(&tri[0], &tri[1], &tri[2], &[0.1, 0.1]));
            assert!(!is_point_inside_triangle(&tri[0], &tri[1], &tri[2], &[0.6, 0.2]));
        }
        assert!(is_point_inside_triangle(&T[0][0], &T[0][1], &T[0][2], &[0.3, 0.4]));
        assert!(is_point_inside_triangle(&T[8][0], &T[8][1], &T[8][2], &[0.7, 0.7]));
        if false {
            let mut plot = Plot::new();
            draw_triangles(&mut plot, &T);
            plot.grid_and_labels("x", "y")
                .set_range(0.0, 1.0, 0.0, 1.0)
                .set_ticks_x(0.1, 0.5, "")
                .set_ticks_y(0.1, 0.5, "")
                .set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .save("/tmp/gemlab/test_is_point_inside_triangle_2.svg")?;
        }
        Ok(())
    }

    #[test]
    fn triangle_signed_area_works() {
        assert_eq!(triangle_signed_area(&[0.0, 0.0], &[1.0, 0.0], &[0.0, 1.0]), 0.5);
        assert_eq!(triangle_signed_area(&[1.0, 0.0], &[0.0, 0.0], &[0.0, 1.0]), -0.5);
        assert_eq!(triangle_signed_area(&[-1.0, 2.0], &[4.0, -3.0], &[2.0, 3.0]), 10.0);
        assert_eq!(triangle_signed_area(&[-2.0, 3.0], &[-3.0, -1.0], &[3.0, -2.0]), 12.5);
    }
}
