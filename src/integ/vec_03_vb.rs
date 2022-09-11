use super::CommonArgs;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the vector(V) dot gradient(B) integration case 03
///
/// Callback function: `f(w, p, N, B)`
///
/// Vector dot gradient:
///
/// ```text
///      ⌠ → →    →  → →
/// cᵐ = │ w(x) · Bᵐ(x(ξ)) α dΩ
///      ⌡
///      Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///      nip-1  → →     →  →       →
/// cᵐ ≈   Σ    w(ιᵖ) · Bᵐ(ιᵖ) |J|(ιᵖ) wᵖ α
///       p=0
/// ```
///
/// # Results
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
/// # Arguments
///
/// * `c` -- A vector containing all `cᵐ` values, one after another, and sequentially placed as shown above.
///   `m` is the index of the node. The length must be `c.len() ≥ ii0 + nnode`.
/// * `args` --- Common arguments
/// * `fn_w` -- Function `f(w,p,N,B)` that computes `w(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions N(ιᵖ), and gradients B(ιᵖ). `w.dim() = space_ndim`.
///
/// # Examples
///
/// See also the `examples` directory.
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
///     let mut c = Vector::filled(pad.kind.nnode(), 0.0);
///     let mut args = integ::CommonArgs::new(&mut pad, ips);
///     integ::vec_03_vb(&mut c, &mut args, |w, _, _, _| {
///         w[0] = 1.0;
///         w[1] = 2.0;
///         Ok(())
///     })?;
///     // solution (A = 6):
///     // cᵐ = (w₀ Bᵐ₀ + w₁ Bᵐ₁) A
///     //     ┌       ┐
///     //     │ -¼ -⅓ │
///     // B = │  ¼  0 │
///     //     │  0  ⅓ │
///     //     └       ┘
///     vec_approx_eq(c.as_data(), &[-5.5, 1.5, 4.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_03_vb<F>(c: &mut Vector, args: &mut CommonArgs, mut fn_w: F) -> Result<(), StrError>
where
    F: FnMut(&mut Vector, usize, &Vector, &Matrix) -> Result<(), StrError>,
{
    // check
    let (space_ndim, nnode) = args.pad.xxt.dims();
    let ii0 = args.ii0;
    if c.dim() < ii0 + nnode {
        return Err("c.len() must be ≥ ii0 + nnode");
    }

    // allocate auxiliary vector
    let mut w = Vector::new(space_ndim);

    // clear output vector
    if args.clear {
        c.fill(0.0);
    }

    // loop over integration points
    for p in 0..args.ips.len() {
        // ksi coordinates and weight
        let iota = &args.ips[p];
        let weight = args.ips[p][3];

        // calculate Jacobian and Gradient
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B

        // calculate w
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        fn_w(&mut w, p, nn, bb)?;

        // calculate coefficient
        let coef = if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt[0][m];
            }
            det_jac * weight * args.alpha * r
        } else {
            det_jac * weight * args.alpha
        };

        // add contribution to c vector
        if space_ndim == 2 {
            for m in 0..nnode {
                c[ii0 + m] += coef * (w[0] * bb[m][0] + w[1] * bb[m][1]);
            }
        } else {
            for m in 0..nnode {
                c[ii0 + m] += coef * (w[0] * bb[m][0] + w[1] * bb[m][1] + w[2] * bb[m][2]);
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, AnalyticalTet4, AnalyticalTri3, CommonArgs};
    use russell_chk::vec_approx_eq;
    use russell_lab::{Matrix, Vector};

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut c = Vector::new(2);
        let mut w = Vector::new(0);
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let f = |_: &mut Vector, _: usize, _: &Vector, _: &Matrix| Ok(());
        f(&mut w, 0, &nn, &bb).unwrap();
        let mut args = CommonArgs::new(&mut pad, &[]);
        args.ii0 = 1;
        assert_eq!(
            integ::vec_03_vb(&mut c, &mut args, f).err(),
            Some("c.len() must be ≥ ii0 + nnode")
        );
    }

    #[test]
    fn tri3_constant_works() {
        // constant vector function: w(x) = {w₀, w₁}
        const W0: f64 = 2.0;
        const W1: f64 = 3.0;
        let mut pad = aux::gen_pad_tri3();

        // solution
        let ana = AnalyticalTri3::new(&pad);
        let c_correct = ana.vec_03_vb(W0, W1);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [1, 3].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let mut c = Vector::filled(pad.kind.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_03_vb(&mut c, &mut args, |w, _, _, _| {
                w[0] = W0;
                w[1] = W1;
                Ok(())
            })
            .unwrap();
            vec_approx_eq(c.as_data(), c_correct.as_data(), tol);
        });
    }

    #[test]
    fn tri3_bilinear_works() {
        // bilinear vector function: w(x) = {x, y}
        let mut pad = aux::gen_pad_tri3();

        // solution
        let ana = AnalyticalTri3::new(&pad);
        let c_correct = ana.vec_03_vb_bilinear(&pad);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [1, 3].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let mut c = Vector::filled(pad.kind.nnode(), aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            let x_ips = integ::points_coords(&mut args.pad, ips).unwrap();
            integ::vec_03_vb(&mut c, &mut args, |w, p, _, _| {
                w[0] = x_ips[p][0];
                w[1] = x_ips[p][1];
                Ok(())
            })
            .unwrap();
            vec_approx_eq(c.as_data(), c_correct.as_data(), tol);
        });
    }

    #[test]
    fn tet4_constant_works() {
        // tet 4 with constant vector  w(x) = {w0, w1, w2}
        let mut pad = aux::gen_pad_tet4();

        // solution
        const W0: f64 = 2.0;
        const W1: f64 = 3.0;
        const W2: f64 = 4.0;
        let ana = AnalyticalTet4::new(&pad);
        let c_correct = ana.vec_03_vb(W0, W1, W2);

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
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_03_vb(&mut c, &mut args, |w, _, _, _| {
                w[0] = W0;
                w[1] = W1;
                w[2] = W2;
                Ok(())
            })
            .unwrap();
            vec_approx_eq(c.as_data(), &c_correct, tol);
        });
    }
}
