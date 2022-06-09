use crate::shapes::{op, Scratchpad};
use crate::StrError;

use super::IntegPointData;

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
/// * `pad` -- [modified] Scratchpad
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
        let det_jac = op::calc_jacobian(pad, iota)?;

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
    use crate::integ::{calc_ips_coords, select_integ_points};
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::StrError;
    use russell_chk::assert_approx_eq;

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
        let selection: Vec<_> = [1, 4, 5, 1_005, 8, 9, 16]
            .iter()
            .map(|n| select_integ_points(class, *n).unwrap())
            .collect();

        // s(x) is constant = 1; i.e., the integral will result in the area of the "diamond" shape
        let ii_correct = 2.0;
        let tolerances = [1e-15, 1e-15, 1e-15, 1e-15, 1e-15, 1e-15, 1e-14];
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let ii = scalar_field(&mut pad, ips, |_| Ok(1.0)).unwrap();
            assert_approx_eq!(ii, ii_correct, tol);
        });

        // ∫(x²+y²) dx dy
        //                   1      4      5 1_005      8      9     16
        let tolerances = [0.67, 1e-15, 1e-15, 1e-6, 1e-15, 1e-15, 1e-14];
        let ii_correct = 8.0 / 3.0;
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = calc_ips_coords(&mut pad, ips).unwrap();
            let ii = scalar_field(&mut pad, ips, |p| {
                let x = x_ips[p][0];
                let y = x_ips[p][1];
                Ok(x * x + y * y)
            })
            .unwrap();
            assert_approx_eq!(ii, ii_correct, tol);
        });

        // ∫(x³+y³) dx dy
        //                   1      4      5 1_005      8      9     16
        let tolerances = [1.01, 1e-15, 1e-15, 1e-6, 1e-15, 1e-15, 1e-15];
        let ii_correct = 3.0;
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = calc_ips_coords(&mut pad, ips).unwrap();
            let ii = scalar_field(&mut pad, ips, |p| {
                let x = x_ips[p][0];
                let y = x_ips[p][1];
                Ok(x * x * x + y * y * y)
            })
            .unwrap();
            assert_approx_eq!(ii, ii_correct, tol);
        });

        // ∫(y x² + x y²) dx dy
        //                   1      4      5 1_005      8      9     16
        let tolerances = [0.34, 1e-15, 1e-15, 1e-7, 1e-15, 1e-15, 1e-15];
        let ii_correct = 1.0 / 3.0;
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = calc_ips_coords(&mut pad, ips).unwrap();
            let ii = scalar_field(&mut pad, ips, |p| {
                let x = x_ips[p][0];
                let y = x_ips[p][1];
                Ok(y * x * x + x * y * y)
            })
            .unwrap();
            assert_approx_eq!(ii, ii_correct, tol);
        });

        // Mathematica codes:
        // Integrate[1, {x, y} \[Element] Polygon[{{0, 0}, {1, -1}, {2, 0}, {1, 1}}]]
        //     2
        // Integrate[x^2 + y^2, {x, y} \[Element] Polygon[{{0, 0}, {1, -1}, {2, 0}, {1, 1}}]]
        //     8/3
        // Integrate[x^3 + y^3, {x, y} \[Element] Polygon[{{0, 0}, {1, -1}, {2, 0}, {1, 1}}]]
        //     3
        // Integrate[y x^2 + x y^2, {x, y} \[Element] Polygon[{{0, 0}, {1, -1}, {2, 0}, {1, 1}}]]
        //     1/3
        Ok(())
    }
}
