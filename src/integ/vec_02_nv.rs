use super::CommonArgs;
use crate::StrError;
use russell_lab::Vector;

/// Implements the the shape(N) times vector(V) integration case 02
///
/// Callback function: `f(v, p, N)`
///
/// Interpolation functions times vector field:
///
/// ```text
/// →    ⌠    → →   → →
/// bᵐ = │ Nᵐ(x(ξ)) v(x) α dΩ
///      ⌡
///      Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
/// →    nip-1     →   → →       →
/// bᵐ ≈   Σ    Nᵐ(ιᵖ) v(ιᵖ) |J|(ιᵖ) wᵖ α
///       p=0
/// ```
///
/// # Results
///
/// ```text
///     ┌     ┐
///     | b⁰₀ |  ⟸  ii0 = 0
///     | b⁰₁ |
///     | b¹₀ |
/// b = | b¹₁ |
///     | b²₀ |
///     | b²₁ |
///     | ··· |
///     | bᵐᵢ |  ⟸  ii := i + m ⋅ space_ndim
///     └     ┘       
///
/// m = ii / space_ndim
/// i = ii % space_ndim
/// ```
///
/// # Arguments
///
/// * `b` -- A vector containing all `bᵐᵢ` values, one after another, and sequentially placed
///   as shown above (in 2D). `m` is the index of the node and `i` corresponds to `space_ndim`.
///   The length must be `b.len() ≥ ii0 + nnode ⋅ space_ndim`
/// * `args` --- Common arguments
/// * `fn_v` -- Function `f(v,p,N)` that computes `v(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   and shape functions N(ιᵖ). `v.dim() = space_ndim`.
///
/// # Examples
///
/// ```
/// use gemlab::integ;
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
///     let ips = integ::default_points(pad.kind);
///     let mut b = Vector::filled(pad.kind.nnode() * space_ndim, 0.0);
///     let mut args = integ::CommonArgs::new(&mut pad, ips);
///     integ::vec_02_nv(&mut b, &mut args, |v, _, _| {
///         v[0] = 1.0;
///         v[1] = 2.0;
///         Ok(())
///     })?;
///     // solution (A = 6):
///     // bᵐ₀ = v₀ A / 3
///     // bᵐ₁ = v₁ A / 3
///     vec_approx_eq(&b, &[2.0, 4.0, 2.0, 4.0, 2.0, 4.0], 1e-14);
///     Ok(())
/// }
/// ```
pub fn vec_02_nv<F>(b: &mut Vector, args: &mut CommonArgs, mut fn_v: F) -> Result<(), StrError>
where
    F: FnMut(&mut Vector, usize, &Vector) -> Result<(), StrError>,
{
    // check
    let (space_ndim, nnode) = args.pad.xxt.dims();
    let ii0 = args.ii0;
    if b.dim() < ii0 + nnode * space_ndim {
        return Err("b.len() must be ≥ ii0 + nnode ⋅ space_ndim");
    }

    // allocate auxiliary vector
    let mut v = Vector::new(space_ndim);

    // clear output vector
    if args.clear {
        b.fill(0.0);
    }

    // loop over integration points
    for p in 0..args.ips.len() {
        // ksi coordinates and weight
        let iota = &args.ips[p];
        let weight = args.ips[p][3];

        // calculate interpolation functions and Jacobian
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_jacobian(iota)?;

        // calculate v
        let nn = &args.pad.interp;
        fn_v(&mut v, p, nn)?;

        // calculate coefficient
        let coef = if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt.get(0, m);
            }
            det_jac * weight * args.alpha * r
        } else {
            det_jac * weight * args.alpha
        };

        // add contribution to b vector
        if space_ndim == 2 {
            for m in 0..nnode {
                b[ii0 + 0 + m * 2] += coef * nn[m] * v[0];
                b[ii0 + 1 + m * 2] += coef * nn[m] * v[1];
            }
        } else {
            for m in 0..nnode {
                b[ii0 + 0 + m * 3] += coef * nn[m] * v[0];
                b[ii0 + 1 + m * 3] += coef * nn[m] * v[1];
                b[ii0 + 2 + m * 3] += coef * nn[m] * v[2];
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
    use russell_lab::{vec_approx_eq, Vector};

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut b = Vector::new(4);
        let mut v = Vector::new(0);
        let nn = Vector::new(0);
        let f = |_v: &mut Vector, _p: usize, _nn: &Vector| Ok(());
        f(&mut v, 0, &nn).unwrap();
        let mut args = CommonArgs::new(&mut pad, &[]);
        args.ii0 = 1;
        assert_eq!(
            integ::vec_02_nv(&mut b, &mut args, f).err(),
            Some("b.len() must be ≥ ii0 + nnode ⋅ space_ndim")
        );
    }

    #[test]
    fn lin2_linear_works() {
        // This test is similar to the shape_times_scalar with lin2
        const L: f64 = 6.0;
        let mut pad = aux::gen_pad_lin2(L);

        // solution
        let cf = L / 6.0;
        let (xa, xb) = (pad.xxt.get(0, 0), pad.xxt.get(0, 1));
        let b_correct = &[
            cf * (2.0 * xa + xb),
            cf * (2.0 * xa + xb),
            cf * (xa + 2.0 * xb),
            cf * (xa + 2.0 * xb),
        ];

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [2, 3].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let (space_ndim, nnode) = pad.xxt.dims();
        let mut b = Vector::filled(nnode * space_ndim, aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let x_ips = integ::points_coords(&mut pad, ips).unwrap();
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_02_nv(&mut b, &mut args, |v, p, _| {
                v[0] = x_ips[p][0];
                v[1] = x_ips[p][0]; // << note use of x component here too
                Ok(())
            })
            .unwrap();
            vec_approx_eq(&b, b_correct, tol);
        });
    }

    #[test]
    fn tri3_constant_works() {
        // This test is similar to the shape_times_scalar with tri3, however using a vector
        // So, each component of `b` equals `Fₛ`
        let mut pad = aux::gen_pad_tri3();

        // solution
        let ana = AnalyticalTri3::new(&pad);
        const V0: f64 = -3.0;
        const V1: f64 = 8.0;
        let b_correct = ana.vec_02_nv(V0, V1);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [1, 3].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let (space_ndim, nnode) = pad.xxt.dims();
        let mut b = Vector::filled(nnode * space_ndim, aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_02_nv(&mut b, &mut args, |v, _, _| {
                v[0] = V0;
                v[1] = V1;
                Ok(())
            })
            .unwrap();
            vec_approx_eq(&b, &b_correct, tol);
        });
    }

    #[test]
    fn tet4_constant_works() {
        // tet 4 with constant vector
        const V0: f64 = 2.0;
        const V1: f64 = 3.0;
        const V2: f64 = 4.0;
        let mut pad = aux::gen_pad_tet4();

        // solution
        let ana = AnalyticalTet4::new(&pad);
        let b_correct = ana.vec_02_nv(V0, V1, V2);

        // integration points
        let class = pad.kind.class();
        let tolerances = [1e-15, 1e-15];
        let selection: Vec<_> = [1, 4].iter().map(|n| integ::points(class, *n).unwrap()).collect();

        // check
        let (space_ndim, nnode) = pad.xxt.dims();
        let mut b = Vector::filled(nnode * space_ndim, aux::NOISE);
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::vec_02_nv(&mut b, &mut args, |v, _, _| {
                v[0] = V0;
                v[1] = V1;
                v[2] = V2;
                Ok(())
            })
            .unwrap();
            vec_approx_eq(&b, &b_correct, tol);
        });
    }
}
