use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::Vector;

/// Implements the vector(V) dot gradient(G) integration case 03
///
/// Vector dot gradient:
///
/// ```text
///      ⌠ → →    →  → →
/// cᵐ = │ w(x) · Gᵐ(x(ξ)) dΩ
///      ⌡
///      Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///      nip-1  → →     →  →       →
/// cᵐ ≈   Σ    w(ιᵖ) · Gᵐ(ιᵖ) |J|(ιᵖ) wᵖ
///       p=0
/// ```
///
/// # Output
///
/// ```text
///     ┌     ┐
///     |  c⁰ |  ⟸  ii0 = 0
///     |  c¹ |
/// c = |  c² |
///     | ··· |
///     |  cᵐ |  ⟸  ii
///     └     ┘
/// ```
///
/// * `c` -- A vector containing all `cᵐ` values, one after another, and sequentially placed as shown above.
///   `m` is the index of the node. The length must be `c.len() ≥ ii0 + nnode`.
/// * `pad` -- Some members of the scratchpad will be modified
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output vector where to add components
/// * `clear_c` -- Fills `c` vector with zeros, otherwise accumulate values into `c`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_w` -- Function `f(w,p)` that computes `w(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`.
///    The dim of `w` is equal to `space_ndim`.
///
/// # Examples
///
/// See also the `examples` directory.
///
/// ```
/// use gemlab::integ;
/// use gemlab::shapes::{GeoKind, Scratchpad};
/// use gemlab::StrError;
/// use russell_chk::assert_vec_approx_eq;
/// use russell_lab::Vector;
///
/// fn main() -> Result<(), StrError> {
///     let space_ndim = 2;
///     let mut pad = Scratchpad::new(space_ndim, GeoKind::Tri3)?;
///     pad.set_xx(0, 0, 2.0);
///     pad.set_xx(0, 1, 3.0);
///     pad.set_xx(1, 0, 6.0);
///     pad.set_xx(1, 1, 3.0);
///     pad.set_xx(2, 0, 2.0);
///     pad.set_xx(2, 1, 6.0);
///     let ips = integ::default_points(pad.kind);
///     let mut c = Vector::filled(pad.kind.nnode(), 0.0);
///     integ::vec_03_vg(&mut c, &mut pad, 0, true, ips, |w, _| {
///         w[0] = 1.0;
///         w[1] = 2.0;
///         Ok(())
///     })?;
///     // solution (A = 6):
///     // cᵐ = (w₀ Gᵐ₀ + w₁ Gᵐ₁) A
///     //     ┌       ┐
///     //     │ -¼ -⅓ │
///     // G = │  ¼  0 │
///     //     │  0  ⅓ │
///     //     └       ┘
///     assert_vec_approx_eq!(c.as_data(), &[-5.5, 1.5, 4.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_03_vg<F>(
    c: &mut Vector,
    pad: &mut Scratchpad,
    ii0: usize,
    clear_c: bool,
    ips: IntegPointData,
    fn_w: F,
) -> Result<(), StrError>
where
    F: Fn(&mut Vector, usize) -> Result<(), StrError>,
{
    // check
    let (space_ndim, nnode) = pad.xxt.dims();
    if c.dim() < ii0 + nnode {
        return Err("c.len() must be ≥ ii0 + nnode");
    }

    // allocate auxiliary vector
    let mut w = Vector::new(space_ndim);

    // clear output vector
    if clear_c {
        c.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate Jacobian and Gradient
        let det_jac = pad.calc_gradient(iota)?;

        // calculate w
        fn_w(&mut w, p)?;

        // add contribution to c vector
        let coef = det_jac * weight;
        let g = &pad.gradient;
        if space_ndim == 2 {
            for m in 0..nnode {
                c[ii0 + m] += coef * (w[0] * g[m][0] + w[1] * g[m][1]);
            }
        } else {
            for m in 0..nnode {
                c[ii0 + m] += coef * (w[0] * g[m][0] + w[1] * g[m][1] + w[2] * g[m][2]);
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, AnalyticalTet4, AnalyticalTri3};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut c = Vector::new(2);
        assert_eq!(
            integ::vec_03_vg(&mut c, &mut pad, 1, false, &[], |_, _| Ok(())).err(),
            Some("c.len() must be ≥ ii0 + nnode")
        );
    }

    #[test]
    fn vec_03_vg_works_tri3_constant() {
        // constant vector function: w(x) = {w₀, w₁}
        const W0: f64 = 2.0;
        const W1: f64 = 3.0;
        let mut pad = aux::gen_pad_tri3();

        // solution
        let ana = AnalyticalTri3::new(&pad);
        let c_correct = ana.vec_03_vg(W0, W1);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [1, 3].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let mut c = Vector::filled(pad.kind.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::vec_03_vg(&mut c, &mut pad, 0, true, ips, |w, _| {
                w[0] = W0;
                w[1] = W1;
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(c.as_data(), c_correct, tol);
        });
    }

    #[test]
    fn vec_03_vg_works_tri3_bilinear() {
        // bilinear vector function: w(x) = {x, y}
        let mut pad = aux::gen_pad_tri3();

        // solution
        let ana = AnalyticalTri3::new(&pad);
        let c_correct = ana.vec_03_vg_bilinear(&pad);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [1, 3].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let mut c = Vector::filled(pad.kind.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = integ::points_coords(&mut pad, ips).unwrap();
            integ::vec_03_vg(&mut c, &mut pad, 0, true, ips, |w, p| {
                w[0] = x_ips[p][0];
                w[1] = x_ips[p][1];
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(c.as_data(), c_correct, tol);
        });
    }

    #[test]
    fn vec_03_vg_works_tet4_constant() {
        // tet 4 with constant vector  w(x) = {w0, w1, w2}
        let mut pad = aux::gen_pad_tet4();

        // solution
        const W0: f64 = 2.0;
        const W1: f64 = 3.0;
        const W2: f64 = 4.0;
        let ana = AnalyticalTet4::new(&pad);
        let c_correct = ana.vec_03_vg(W0, W1, W2);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14];
        let selection: Vec<_> = [1, 4, 5, 8, 14, 15, 24]
            .iter()
            .map(|n| integ::points(class, *n).unwrap())
            .collect();

        // check
        let mut c = Vector::filled(pad.kind.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::vec_03_vg(&mut c, &mut pad, 0, true, ips, |w, _| {
                w[0] = W0;
                w[1] = W1;
                w[2] = W2;
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(c.as_data(), c_correct, tol);
        });
    }
}
