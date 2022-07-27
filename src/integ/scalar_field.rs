use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;

/// Implements the integration of a scalar field over a geometric shape
///
/// ```text
///     ⌠   → →
/// I = │ s(x(ξ)) dΩ
///     ⌡
///     Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///     nip-1    →       →
/// I ≈   Σ    s(ιᵖ) |J|(ιᵖ) wᵖ
///      p=0
/// ```
///
/// # Input
///
/// * `pad` -- **modified** Scratchpad
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_s` -- Function `f(p)` corresponding to `s(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
///
/// # Output
///
/// * Returns `I`, the result of integration.
pub fn scalar_field<F>(pad: &mut Scratchpad, ips: IntegPointData, fn_s: F) -> Result<f64, StrError>
where
    F: Fn(usize) -> Result<f64, StrError>,
{
    // result from integration
    let mut ii = 0.0;

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate Jacobian
        let det_jac = pad.calc_jacobian(iota)?;

        // calculate s
        let s = fn_s(p)?;

        // loop over nodes and perform sum
        ii += s * det_jac * weight;
    }
    Ok(ii)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::scalar_field;
    use crate::integ;
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::StrError;
    use russell_chk::assert_approx_eq;

    #[allow(unused_imports)]
    use plotpy::Plot;

    #[test]
    fn scalar_fields_over_rotated_square() -> Result<(), StrError> {
        //       y         (1,1)
        //                   2
        //       ^        .'   `.
        //       |      .'   .   `.
        //       |    .'           `.
        //       |  .'       .       `.
        //        .'                   `.
        // (0,0) 0   .   .   .   .   .   1 (2,0)  ----> x
        //        '.                   .'
        //          '.       .       .'
        //            '.           .'
        //              '.   .   .'
        //                '.   .'
        //                  '1'
        //                 (1,-1)
        let mut pad = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, 1.0);
        pad.set_xx(1, 1, -1.0);
        pad.set_xx(2, 0, 2.0);
        pad.set_xx(2, 1, 0.0);
        pad.set_xx(3, 0, 1.0);
        pad.set_xx(3, 1, 1.0);

        // integration points
        let class = pad.kind.class();
        let selection: Vec<_> = [1, 4, 9, 16]
            .iter()
            .map(|n| integ::points(class, *n).unwrap())
            .collect();

        // s(x) is constant = 1; i.e., the integral will result in the area of the "diamond" shape
        let ii_correct = 2.0;
        let tolerances = [1e-15, 1e-15, 1e-15, 1e-14];
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let ii = scalar_field(&mut pad, ips, |_| Ok(1.0)).unwrap();
            assert_approx_eq!(ii, ii_correct, tol);
        });

        // ∫∫(x²+y²) dx dy
        let tolerances = [0.67, 1e-15, 1e-15, 1e-14];
        let ii_correct = 8.0 / 3.0;
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = integ::points_coords(&mut pad, ips).unwrap();
            let ii = scalar_field(&mut pad, ips, |p| {
                let x = x_ips[p][0];
                let y = x_ips[p][1];
                Ok(x * x + y * y)
            })
            .unwrap();
            assert_approx_eq!(ii, ii_correct, tol);
        });

        // ∫∫(x³+y³) dx dy
        let tolerances = [1.01, 1e-15, 1e-15, 1e-15];
        let ii_correct = 3.0;
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = integ::points_coords(&mut pad, ips).unwrap();
            let ii = scalar_field(&mut pad, ips, |p| {
                let x = x_ips[p][0];
                let y = x_ips[p][1];
                Ok(x * x * x + y * y * y)
            })
            .unwrap();
            assert_approx_eq!(ii, ii_correct, tol);
        });

        // Mathematica code:
        // region = Polygon[{{0, 0}, {1, -1}, {2, 0}, {1, 1}}];
        // Graphics[{Purple, region}]
        // Integrate[1, {x, y} \[Element] region]
        //     (* returns 2 *)
        // Integrate[x^2 + y^2, {x, y} \[Element] region]
        //     (* returns 8/3 *)
        // Integrate[x^3 + y^3, {x, y} \[Element] region]
        //     (* returns 3 *)
        Ok(())
    }

    #[test]
    fn scalar_fields_over_slanted_hex8() -> Result<(), StrError> {
        let mut pad = Scratchpad::new(3, GeoKind::Hex8).unwrap();
        // node 0
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(0, 2, 0.0);
        // node 1
        pad.set_xx(1, 0, 1.0);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(1, 2, 0.0);
        // node 2
        pad.set_xx(2, 0, 2.0);
        pad.set_xx(2, 1, 1.0);
        pad.set_xx(2, 2, 0.0);
        // node 3
        pad.set_xx(3, 0, 1.0);
        pad.set_xx(3, 1, 1.0);
        pad.set_xx(3, 2, 0.0);
        // node 4
        pad.set_xx(4, 0, 0.0);
        pad.set_xx(4, 1, 0.0);
        pad.set_xx(4, 2, 1.0);
        // node 5
        pad.set_xx(5, 0, 1.0);
        pad.set_xx(5, 1, 0.0);
        pad.set_xx(5, 2, 1.0);
        // node 6
        pad.set_xx(6, 0, 2.0);
        pad.set_xx(6, 1, 1.0);
        pad.set_xx(6, 2, 1.0);
        // node 7
        pad.set_xx(7, 0, 1.0);
        pad.set_xx(7, 1, 1.0);
        pad.set_xx(7, 2, 1.0);

        if false {
            let mut plot = Plot::new();
            pad.draw_shape(&mut plot, "", true, true)?;
            plot.set_equal_axes(true)
                .set_figure_size_points(400.0, 400.0)
                .save("/tmp/gemlab/test_scalar_fields_over_slanted_hex8.svg")?;
        }

        // integration points
        let class = pad.kind.class();
        let selection: Vec<_> = [6, 8, 14, 27, 64]
            .iter()
            .map(|n| integ::points(class, *n).unwrap())
            .collect();

        // s(x) is constant = 1; i.e., the integral will result in the volume
        let ii_correct = 1.0;
        let tolerances = [1e-15, 1e-15, 1e-15, 1e-15, 1e-15];
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let ii = scalar_field(&mut pad, ips, |_| Ok(1.0)).unwrap();
            assert_approx_eq!(ii, ii_correct, tol);
        });

        // ∫∫∫(x²+y²+z²) dx dy dz
        let ii_correct = 11.0 / 6.0;
        let tolerances = [1e-15, 1e-15, 1e-15, 1e-14, 1e-15];
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = integ::points_coords(&mut pad, ips).unwrap();
            let ii = scalar_field(&mut pad, ips, |p| {
                let x = x_ips[p][0];
                let y = x_ips[p][1];
                let z = x_ips[p][2];
                Ok(x * x + y * y + z * z)
            })
            .unwrap();
            assert_approx_eq!(ii, ii_correct, tol);
        });

        // ∫∫∫(x³+y³+z³) dx dy dz
        let ii_correct = 2.0;
        let tolerances = [1e-15, 1e-15, 1e-15, 1e-15, 1e-14];
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = integ::points_coords(&mut pad, ips).unwrap();
            let ii = scalar_field(&mut pad, ips, |p| {
                let x = x_ips[p][0];
                let y = x_ips[p][1];
                let z = x_ips[p][2];
                Ok(x * x * x + y * y * y + z * z * z)
            })
            .unwrap();
            assert_approx_eq!(ii, ii_correct, tol);
        });

        // mathematica code:
        // region = Hexahedron[{{0, 0, 0}, {1, 0, 0}, {2, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {2, 1, 1}, {1, 1, 1}}];
        // Graphics3D[{region}]
        // Integrate[1, {x, y, z} \[Element] region]
        //     (* returns 1 *)
        // Integrate[x^2 + y^2 + z^2, {x, y, z} \[Element] region]
        //     (* returns 11/6 *)
        // Integrate[x^3 + y^3 + z^3, {x, y, z} \[Element] region]
        //     (* returns 2 *)
        Ok(())
    }
}
