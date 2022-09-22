use super::CommonArgs;
use crate::StrError;
use russell_lab::math::SQRT_2;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the tensor(T) dot gradient(B) integration case 04
///
/// Callback function: `f(σ, p, N, B)`
///
/// Tensor dot gradient:
///
/// ```text
/// →    ⌠   →    →  → →
/// dᵐ = │ σ(x) · Bᵐ(x(ξ)) α dΩ
///      ⌡ ▔
///      Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
/// →    nip-1    →     →  →       →
/// dᵐ ≈   Σ    σ(ιᵖ) · Bᵐ(ιᵖ) |J|(ιᵖ) wᵖ α
///       p=0   ▔
/// ```
///
/// # Results
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
/// # Arguments
///
/// * `d` -- A vector containing all `dᵐᵢ` values, one after another, and sequentially placed
///   as shown above (in 2D). `m` is the index of the node and `i` corresponds to `space_ndim`.
///   The length must be `d.len() ≥ ii0 + nnode ⋅ space_ndim`
/// * `args` --- Common arguments
/// * `fn_sig` -- Function `f(σ,p,N,B)` that computes `σ(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions N(ιᵖ), and the gradients B(ιᵖ). `σ` is set for `space_ndim`.
///
/// # Examples
///
/// ```
/// use gemlab::integ;
/// use gemlab::shapes::{GeoKind, Scratchpad};
/// use gemlab::StrError;
/// use russell_chk::vec_approx_eq;
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
///     let mut args = integ::CommonArgs::new(&mut pad, ips);
///     integ::vec_04_tb(&mut d, &mut args, |sig, _, _, _| {
///         sig.sym_set(0, 0, 1.0);
///         sig.sym_set(1, 1, 2.0);
///         sig.sym_set(0, 1, 3.0);
///         Ok(())
///     })?;
///     // solution (A = 6):
///     // dᵐ₀ = (σ₀₀ Bᵐ₀ + σ₀₁ Bᵐ₁) A
///     // dᵐ₁ = (σ₁₀ Bᵐ₀ + σ₁₁ Bᵐ₁) A
///     //     ┌       ┐
///     //     │ -¼ -⅓ │
///     // B = │  ¼  0 │
///     //     │  0  ⅓ │
///     //     └       ┘
///     vec_approx_eq(d.as_data(), &[-7.5, -8.5, 1.5, 4.5, 6.0, 4.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_04_tb<F>(d: &mut Vector, args: &mut CommonArgs, mut fn_sig: F) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor2, usize, &Vector, &Matrix) -> Result<(), StrError>,
{
    // check
    let (space_ndim, nnode) = args.pad.xxt.dims();
    if d.dim() < args.ii0 + nnode * space_ndim {
        return Err("d.len() must be ≥ ii0 + nnode ⋅ space_ndim");
    }
    if args.axisymmetric && space_ndim != 2 {
        return Err("axisymmetric requires space_ndim = 2");
    }

    // allocate auxiliary tensor
    let mut sig = Tensor2::new(true, space_ndim == 2);

    // clear output vector
    if args.clear {
        d.fill(0.0);
    }

    // loop over integration points
    for index in 0..args.ips.len() {
        // ksi coordinates and weight
        let iota = &args.ips[index];
        let weight = args.ips[index][3];

        // calculate Jacobian and Gradient
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B

        // calculate σ tensor
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        fn_sig(&mut sig, index, nn, bb)?;

        // add contribution to d vector
        let c = det_jac * weight * args.alpha;
        if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt[0][m];
            }
            add_to_d_axisymmetric(d, nnode, c, r, &sig, args);
        } else {
            add_to_d(d, space_ndim, nnode, c, &sig, args);
        }
    }
    Ok(())
}

/// Adds contribution to the d-vector in vec_04_tb
#[inline]
fn add_to_d(d: &mut Vector, ndim: usize, nnode: usize, c: f64, sig: &Tensor2, args: &mut CommonArgs) {
    let t = &sig.vec;
    let s = SQRT_2;
    let b = &args.pad.gradient;
    let ii0 = args.ii0;
    if ndim == 2 {
        for m in 0..nnode {
            d[ii0 + 0 + m * 2] += c * (t[0] * b[m][0] + t[3] * b[m][1] / s);
            d[ii0 + 1 + m * 2] += c * (t[3] * b[m][0] / s + t[1] * b[m][1]);
        }
    } else {
        for m in 0..nnode {
            d[ii0 + 0 + m * 3] += c * (t[0] * b[m][0] + t[3] * b[m][1] / s + t[5] * b[m][2] / s);
            d[ii0 + 1 + m * 3] += c * (t[3] * b[m][0] / s + t[1] * b[m][1] + t[4] * b[m][2] / s);
            d[ii0 + 2 + m * 3] += c * (t[5] * b[m][0] / s + t[4] * b[m][1] / s + t[2] * b[m][2]);
        }
    }
}

/// Adds contribution to the d-vector in vec_04_tb (axisymmetric case)
#[inline]
fn add_to_d_axisymmetric(d: &mut Vector, nnode: usize, c: f64, r: f64, sig: &Tensor2, args: &mut CommonArgs) {
    let t = &sig.vec;
    let s = SQRT_2;
    let nn = &args.pad.interp;
    let b = &args.pad.gradient;
    let ii0 = args.ii0;
    for m in 0..nnode {
        d[ii0 + 0 + m * 2] += c * r * (t[0] * b[m][0] + t[3] * b[m][1] / s) + nn[m] * t[2];
        d[ii0 + 1 + m * 2] += c * r * (t[3] * b[m][0] / s + t[1] * b[m][1]);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, AnalyticalTet4, AnalyticalTri3, CommonArgs};
    use russell_chk::vec_approx_eq;
    use russell_lab::{Matrix, Vector};
    use russell_tensor::{copy_tensor2, Tensor2};

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut d = Vector::new(4);
        let mut sig = Tensor2::new(true, true);
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let f = |_: &mut Tensor2, _, _: &Vector, _: &Matrix| Ok(());
        f(&mut sig, 0, &nn, &bb).unwrap();
        let mut args = CommonArgs::new(&mut pad, &[]);
        args.ii0 = 1;
        assert_eq!(
            integ::vec_04_tb(&mut d, &mut args, f).err(),
            Some("d.len() must be ≥ ii0 + nnode ⋅ space_ndim")
        );
    }

    #[test]
    fn tri3_constant_works() {
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
        let sig = Tensor2::from_matrix(&[[S00, S01, 0.0], [S01, S11, 0.0], [0.0, 0.0, S22]], true, true).unwrap();
        let d_correct = ana.vec_04_tb(&sig, false);

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
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_04_tb(&mut d, &mut args, |sig, _, _, _| {
                sig.sym_set(0, 0, S00);
                sig.sym_set(1, 1, S11);
                sig.sym_set(2, 2, S22);
                sig.sym_set(0, 1, S01);
                Ok(())
            })
            .unwrap();
            vec_approx_eq(d.as_data(), d_correct.as_data(), tol);
        });
    }

    #[test]
    fn tet4_constant_works() {
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
        let d_correct = ana.vec_04_tb(&tt);

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
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_04_tb(&mut d, &mut args, |sig, _, _, _| {
                copy_tensor2(sig, &tt).unwrap();
                Ok(())
            })
            .unwrap();
            vec_approx_eq(d.as_data(), &d_correct, tol);
        });
    }
}
