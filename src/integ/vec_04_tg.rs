use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the tensor(T) dot gradient(G) integration case 04
///
/// Tensor dot gradient:
///
/// ```text
/// →    ⌠   →    →  → →
/// dᵐ = │ σ(x) · Gᵐ(x(ξ)) dΩ
///      ⌡ ▔
///      Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
/// →    nip-1    →     →  →       →
/// dᵐ ≈   Σ    σ(ιᵖ) · Gᵐ(ιᵖ) |J|(ιᵖ) wᵖ
///       p=0   ▔
/// ```
///
/// # Output
///
/// ```text
///     ┌     ┐
///     | d⁰₀ |  ⟸  ii0 = 0
///     | d⁰₁ |
///     | d¹₀ |
/// d = | d¹₁ |
///     | d²₀ |
///     | d²₁ |
///     | ··· |
///     | dᵐᵢ |  ⟸  ii := i + m ⋅ space_ndim
///     └     ┘
///
/// m = ii / space_ndim
/// i = ii % space_ndim
/// ```
///
/// * `d` -- A vector containing all `dᵐᵢ` values, one after another, and sequentially placed
///   as shown above (in 2D). `m` is the index of the node and `i` corresponds to `space_ndim`.
///   The length must be `d.len() ≥ ii0 + nnode ⋅ space_ndim`
/// * `pad` -- Some members of the scratchpad will be modified
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output vector where to add components
/// * `clear_d` -- Fills `d` vector with zeros, otherwise accumulate values into `d`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_sig` -- Function `f(sig,p,G)` that computes `σ(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   and the gradients G(ιᵖ). `T` is set for `space_ndim`.
///
/// # Examples
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
///     let mut d = Vector::filled(pad.kind.nnode() * space_ndim, 0.0);
///     integ::vec_04_tg(&mut d, &mut pad, 0, true, ips, |sig, _, _| {
///         sig.sym_set(0, 0, 1.0);
///         sig.sym_set(1, 1, 2.0);
///         sig.sym_set(0, 1, 3.0);
///         Ok(())
///     })?;
///     // solution (A = 6):
///     // dᵐ₀ = (σ₀₀ Gᵐ₀ + σ₀₁ Gᵐ₁) A
///     // dᵐ₁ = (σ₁₀ Gᵐ₀ + σ₁₁ Gᵐ₁) A
///     //     ┌       ┐
///     //     │ -¼ -⅓ │
///     // G = │  ¼  0 │
///     //     │  0  ⅓ │
///     //     └       ┘
///     assert_vec_approx_eq!(d.as_data(), &[-7.5, -8.5, 1.5, 4.5, 6.0, 4.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_04_tg<F>(
    d: &mut Vector,
    pad: &mut Scratchpad,
    ii0: usize,
    clear_d: bool,
    ips: IntegPointData,
    mut fn_sig: F,
) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor2, usize, &Matrix) -> Result<(), StrError>,
{
    // check
    let (space_ndim, nnode) = pad.xxt.dims();
    if d.dim() < ii0 + nnode * space_ndim {
        return Err("d.len() must be ≥ ii0 + nnode ⋅ space_ndim");
    }

    // allocate auxiliary tensor
    let mut sig = Tensor2::new(true, space_ndim == 2);

    // clear output vector
    if clear_d {
        d.fill(0.0);
    }

    // loop over integration points
    let s = SQRT_2;
    for index in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[index];
        let weight = ips[index][3];

        // calculate Jacobian and Gradient
        let det_jac = pad.calc_gradient(iota)?; // G

        // calculate σ tensor
        let gg = &pad.gradient;
        fn_sig(&mut sig, index, gg)?;

        // add contribution to d vector
        let coef = det_jac * weight;
        let t = &sig.vec;
        if space_ndim == 2 {
            for m in 0..nnode {
                d[ii0 + 0 + m * 2] += coef * (t[0] * gg[m][0] + t[3] * gg[m][1] / s);
                d[ii0 + 1 + m * 2] += coef * (t[3] * gg[m][0] / s + t[1] * gg[m][1]);
            }
        } else {
            for m in 0..nnode {
                d[ii0 + 0 + m * 3] += coef * (t[0] * gg[m][0] + t[3] * gg[m][1] / s + t[5] * gg[m][2] / s);
                d[ii0 + 1 + m * 3] += coef * (t[3] * gg[m][0] / s + t[1] * gg[m][1] + t[4] * gg[m][2] / s);
                d[ii0 + 2 + m * 3] += coef * (t[5] * gg[m][0] / s + t[4] * gg[m][1] / s + t[2] * gg[m][2]);
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
    use russell_lab::{copy_vector, Vector};
    use russell_tensor::Tensor2;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut d = Vector::new(4);
        assert_eq!(
            integ::vec_04_tg(&mut d, &mut pad, 1, false, &[], |_, _, _| Ok(())).err(),
            Some("d.len() must be ≥ ii0 + nnode ⋅ space_ndim")
        );
    }

    #[test]
    fn vec_04_tri3_works_constant() {
        // constant tensor function: σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2}
        // solution:
        //    dᵐ₀ = ½ (σ₀₀ bₘ + σ₀₁ cₘ)
        //    dᵐ₁ = ½ (σ₁₀ bₘ + σ₁₁ cₘ)
        let mut pad = aux::gen_pad_tri3();

        // solution
        const S00: f64 = 2.0;
        const S11: f64 = 3.0;
        const S22: f64 = 4.0;
        const S01: f64 = 5.0;
        let ana = AnalyticalTri3::new(&pad);
        let d_correct = ana.vec_04_tg(S00, S11, S01);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14, 1e-14, 1e-13, 1e-14];
        let selection: Vec<_> = [1, 3, 4, 12, 16]
            .iter()
            .map(|n| integ::points(class, *n).unwrap())
            .collect();

        // check
        let (space_ndim, nnode) = pad.xxt.dims();
        let mut d = Vector::filled(nnode * space_ndim, aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::vec_04_tg(&mut d, &mut pad, 0, true, ips, |sig, _, _| {
                sig.sym_set(0, 0, S00);
                sig.sym_set(1, 1, S11);
                sig.sym_set(2, 2, S22);
                sig.sym_set(0, 1, S01);
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(d.as_data(), d_correct, tol);
        });
    }

    #[test]
    fn vec_04_tet4_works_constant() {
        // constant tensor function: σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2, σ₁₂√2, σ₀₂√2}
        let mut pad = aux::gen_pad_tet4();

        // solution
        #[rustfmt::skip]
        let tt = Tensor2::from_matrix(&[
            [2.0, 5.0, 7.0],
            [5.0, 3.0, 6.0],
            [7.0, 6.0, 4.0],
        ], true, false).unwrap();
        let ana = AnalyticalTet4::new(&pad);
        let d_correct = ana.vec_04_tg(&tt);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14, 1e-13, 1e-14, 1e-13, 1e-13, 1e-13];
        let selection: Vec<_> = [1, 4, 5, 8, 14, 15, 24]
            .iter()
            .map(|n| integ::points(class, *n).unwrap())
            .collect();

        // check
        let (space_ndim, nnode) = pad.xxt.dims();
        let mut d = Vector::filled(nnode * space_ndim, aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::vec_04_tg(&mut d, &mut pad, 0, true, ips, |sig, _, _| {
                copy_vector(&mut sig.vec, &tt.vec)
            })
            .unwrap();
            assert_vec_approx_eq!(d.as_data(), d_correct, tol);
        });
    }
}
