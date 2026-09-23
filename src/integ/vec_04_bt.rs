use super::CommonArgs;
use crate::StrError;
use russell_lab::math::SQRT_2;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the gradient(B) dot transpose tensor(T) integration case 04
///
/// Callback function: `f(σ, p, N, B)`
///
/// Tensor dot gradient:
///
/// ```text
/// →    ⌠  →  → →        →         ⌠   →    →  → →
/// dᵐ = │  Bᵐ(x(ξ)) · σᵀ(x) α dΩ = │ σ(x) · Bᵐ(x(ξ)) α dΩ
///      ⌡             ▔            ⌡ ▔
///      Ωₑ                         Ωₑ
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
/// * `fn_sig` -- Function `f(σ,p,N,B)` that computes `σ(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`,
///   shape functions N(ιᵖ), and the gradients B(ιᵖ). `σ` is set for `space_ndim`.
///
/// # Examples
///
/// ```
/// use gemlab::integ::{self, CommonArgs, Gauss};
/// use gemlab::shapes::{GeoKind, Scratchpad};
/// use gemlab::StrError;
/// use russell_lab::{Vector, vec_approx_eq};
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
///     let gauss = Gauss::new(pad.kind);
///     let mut d = Vector::filled(pad.kind.nnode() * space_ndim, 0.0);
///     let mut args = CommonArgs::new(&mut pad, &gauss);
///     integ::vec_04_bt::<4, _>(&mut d, &mut args, |sig, _, _, _| {
///         sig.sym_set_std(0, 0, 1.0);
///         sig.sym_set_std(1, 1, 2.0);
///         sig.sym_set_std(0, 1, 3.0);
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
///     vec_approx_eq(&d, &[-7.5, -8.5, 1.5, 4.5, 6.0, 4.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_04_bt<const N: usize, F>(d: &mut Vector, args: &mut CommonArgs, mut fn_sig: F) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor2<N>, usize, &Vector, &Matrix) -> Result<(), StrError>,
{
    // check
    let (space_ndim, nnode) = args.pad.xxt.dims();
    if d.dim() < args.ii0 + nnode * space_ndim {
        return Err("d.len() must be ≥ ii0 + nnode ⋅ space_ndim");
    }
    if args.axisymmetric && space_ndim != 2 {
        return Err("axisymmetric requires space_ndim = 2");
    }
    if N != 2 * space_ndim {
        return Err("tensor dimension N must equal 2 * space_ndim");
    }

    // allocate auxiliary tensor
    let mut sig = Tensor2::<N>::new();

    // clear output vector
    if args.clear {
        d.fill(0.0);
    }

    // loop over integration points
    for p in 0..args.gauss.npoint() {
        // ksi coordinates and weight
        let iota = args.gauss.coords(p);
        let weight = args.gauss.weight(p);

        // calculate Jacobian and Gradient
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B

        // calculate σ tensor
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        fn_sig(&mut sig, p, nn, bb)?;

        // add contribution to d vector
        let c = det_jac * weight * args.alpha;
        if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt.get(0, m);
            }
            add_to_d_axisymmetric(d, nnode, c, r, &sig, args);
        } else {
            add_to_d(d, space_ndim, nnode, c, &sig, args);
        }
    }
    Ok(())
}

/// Adds contribution to the d-vector in vec_04_bt
#[inline]
fn add_to_d<const N: usize>(
    d: &mut Vector,
    ndim: usize,
    nnode: usize,
    c: f64,
    sig: &Tensor2<N>,
    args: &mut CommonArgs,
) {
    let s = SQRT_2;
    let b = &args.pad.gradient;
    let ii0 = args.ii0;
    if ndim == 2 {
        for m in 0..nnode {
            d[ii0 + 0 + m * 2] += c * (sig.get(0) * b.get(m, 0) + sig.get(3) * b.get(m, 1) / s);
            d[ii0 + 1 + m * 2] += c * (sig.get(3) * b.get(m, 0) / s + sig.get(1) * b.get(m, 1));
        }
    } else {
        for m in 0..nnode {
            d[ii0 + 0 + m * 3] +=
                c * (sig.get(0) * b.get(m, 0) + sig.get(3) * b.get(m, 1) / s + sig.get(5) * b.get(m, 2) / s);
            d[ii0 + 1 + m * 3] +=
                c * (sig.get(3) * b.get(m, 0) / s + sig.get(1) * b.get(m, 1) + sig.get(4) * b.get(m, 2) / s);
            d[ii0 + 2 + m * 3] +=
                c * (sig.get(5) * b.get(m, 0) / s + sig.get(4) * b.get(m, 1) / s + sig.get(2) * b.get(m, 2));
        }
    }
}

/// Adds contribution to the d-vector in vec_04_bt (axisymmetric case)
#[inline]
fn add_to_d_axisymmetric<const N: usize>(
    d: &mut Vector,
    nnode: usize,
    c: f64,
    r: f64,
    sig: &Tensor2<N>,
    args: &mut CommonArgs,
) {
    let s = SQRT_2;
    let nn = &args.pad.interp;
    let b = &args.pad.gradient;
    let ii0 = args.ii0;
    for m in 0..nnode {
        d[ii0 + 0 + m * 2] +=
            c * r * (sig.get(0) * b.get(m, 0) + sig.get(3) * b.get(m, 1) / s) + c * nn[m] * sig.get(2);
        d[ii0 + 1 + m * 2] += c * r * (sig.get(3) * b.get(m, 0) / s + sig.get(1) * b.get(m, 1));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, AnalyticalTet4, AnalyticalTri3, CommonArgs, Gauss};
    use russell_lab::{vec_approx_eq, Matrix, Vector};
    use russell_tensor::Tensor2;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut d = Vector::new(4);
        let mut sig = Tensor2::<4>::new();
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let f = |_: &mut Tensor2<4>, _, _: &Vector, _: &Matrix| Ok(());
        f(&mut sig, 0, &nn, &bb).unwrap();
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::vec_04_bt::<4, _>(&mut d, &mut args, f).err(),
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
        let sig = Tensor2::<4>::from_std_matrix(&[[S00, S01, 0.0], [S01, S11, 0.0], [0.0, 0.0, S22]]).unwrap();
        let d_correct = ana.vec_04_bt(&sig, false);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14, 1e-14, 1e-13, 1e-14];
        let selection: Vec<_> = [1, 3, 4, 12, 16]
            .iter()
            .map(|n| Gauss::new_sized(class, *n).unwrap())
            .collect();

        // check
        let (space_ndim, nnode) = pad.xxt.dims();
        let mut d = Vector::filled(nnode * space_ndim, aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_04_bt::<4, _>(&mut d, &mut args, |sig, _, _, _| {
                sig.sym_set_std(0, 0, S00);
                sig.sym_set_std(1, 1, S11);
                sig.sym_set_std(2, 2, S22);
                sig.sym_set_std(0, 1, S01);
                Ok(())
            })
            .unwrap();
            vec_approx_eq(&d, &d_correct, tol);
        });
    }

    #[test]
    fn tri3_constant_axisymmetric_works() {
        // constant tensor function: σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2}
        let mut pad = aux::gen_pad_tri3();

        // solution
        const S00: f64 = 2.0;
        const S11: f64 = 3.0;
        const S22: f64 = 4.0;
        const S01: f64 = 5.0;
        let ana = AnalyticalTri3::new(&pad);
        let sig = Tensor2::<4>::from_std_matrix(&[[S00, S01, 0.0], [S01, S11, 0.0], [0.0, 0.0, S22]]).unwrap();
        let d_correct = ana.vec_04_bt(&sig, true);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-13, 1e-13];
        let selection: Vec<_> = [1, 3].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();

        // check
        let (space_ndim, nnode) = pad.xxt.dims();
        let mut d = Vector::filled(nnode * space_ndim, aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.data.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            args.axisymmetric = true;
            integ::vec_04_bt::<4, _>(&mut d, &mut args, |sig, _, _, _| {
                sig.sym_set_std(0, 0, S00);
                sig.sym_set_std(1, 1, S11);
                sig.sym_set_std(2, 2, S22);
                sig.sym_set_std(0, 1, S01);
                Ok(())
            })
            .unwrap();
            vec_approx_eq(&d, &d_correct, tol);
        });
    }

    #[test]
    fn tet4_constant_works() {
        // constant tensor function: σ(x) = {σ₀₀, σ₁₁, σ₂₂, σ₀₁√2, σ₁₂√2, σ₀₂√2}
        let mut pad = aux::gen_pad_tet4();

        // solution
        #[rustfmt::skip]
        let tt = Tensor2::<6>::from_std_matrix(&[
            [2.0, 5.0, 7.0],
            [5.0, 3.0, 6.0],
            [7.0, 6.0, 4.0],
        ]).unwrap();
        let ana = AnalyticalTet4::new(&pad);
        let d_correct = ana.vec_04_bt(&tt);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14, 1e-13, 1e-14, 1e-13, 1e-13, 1e-13];
        let selection: Vec<_> = [1, 4, 5, 8, 14, 15, 24]
            .iter()
            .map(|n| Gauss::new_sized(class, *n).unwrap())
            .collect();

        // check
        let (space_ndim, nnode) = pad.xxt.dims();
        let mut d = Vector::filled(nnode * space_ndim, aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_04_bt::<6, _>(&mut d, &mut args, |sig, _, _, _| {
                sig.set_tensor(1.0, &tt);
                Ok(())
            })
            .unwrap();
            vec_approx_eq(&d, &d_correct, tol);
        });
    }
}
