use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the shape(N) times vector(V) times shape(Nb) integration case 06 (e.g., coupling matrix)
///
/// Callback function: `α ← f(v, p, N, G, Nb)`
///
/// **Notes:**
///
/// * `m` ranges over the number of nodes of the *driver* shape specified by `pad` (for `Nᵐ`)
/// * `n` ranges over the number of nodes of the *lower-order* shape specified by `pad_b` (for `Nbⁿ`)
/// * For example, `1 ≤ m ≤ 8` for a `pad→Qua8` and `1 ≤ n ≤ 4` for `pad_b→Qua4`
/// * The determinant of the Jacobian is calculated for `pad` (`pad` is the driver of the calculations)
/// * The number of integration points must consider the nodes of `pad` and the expected order of the whole integrand
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠    →
/// Kᵐⁿ = │ Nᵐ v Nbⁿ α dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///        nip-1     →     →       →       →
/// Kᵐⁿᵢ ≈   Σ   Nᵐ(ιᵖ) vᵢ(ιᵖ) Nbⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
///         p=0
/// ```
///
/// # Output
///
/// ```text
///     ┌                         ┐
///     | K⁰⁰₀ K⁰¹₀ K⁰²₀ ··· K⁰ⁿ₀ |  ⟸  ii0
///     | K⁰⁰₁ K⁰¹₁ K⁰²₁ ··· K⁰ⁿ₁ |
///     | K¹⁰₀ K¹¹₀ K¹²₀ ··· K¹ⁿ₀ |
/// K = | K¹⁰₁ K¹¹₁ K¹²₁ ··· K¹ⁿ₁ |
///     | K²⁰₀ K²¹₀ K²²₀ ··· K²ⁿ₀ |
///     | K²⁰₁ K²¹₁ K²²₁ ··· K²ⁿ₁ |
///     |  ···  ···  ··· ···  ··· |     (pad)
///     | Kᵐ⁰ᵢ Kᵐ¹ᵢ Kᵐ²ᵢ ··· Kᵐⁿᵢ |  ⟸  ii := i + m ⋅ space_ndim
///     └                         ┘
///        ⇑                  ⇑
///       jj0         (pad_b) jj
///
/// m = ii / space_ndim
/// i = ii % space_ndim
/// ```
///
/// * `kk` -- A matrix containing all `Kᵐⁿᵢ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `i` corresponds to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode ⋅ space_ndim` and `ncol(K) ≥ jj0 + pad_b.nnode`.
/// * `pad` -- Driver scratchpad (modified) to compute G
/// * `pad_b` -- Lower-order scratchpad (modified) to compute Nb
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output matrix where to add components.
/// * `jj0` -- Stride marking the first column in the output matrix where to add components.
/// * `clear_kk` -- Fills `kk` matrix with zeros, otherwise accumulate values into `kk`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_v` -- Function `f(v,p,N,G,Nb)→α` that computes `v(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   shape functions N(ιᵖ), gradients G(ιᵖ), and shape functions Nb(ιᵖ). `v.dim() = space_ndim`.
///   `fn_v` returns α that can accommodate plane-strain simulations.
///   **NOTE:** the value α is ignored if axisymmetric = true, because the radius is calculated and used instead.
///
/// # Warning
///
/// The two [crate::shapes::Scratchpad]s mut be compatible, otherwise **calculation errors may occur**.
/// Therefore, `pad_b` must be either the lower-version of `pad` or have the same shape as `pad`.
pub fn mat_06_nvn<F>(
    kk: &mut Matrix,
    pad: &mut Scratchpad,
    pad_b: &mut Scratchpad,
    ii0: usize,
    jj0: usize,
    clear_kk: bool,
    axisymmetric: bool,
    ips: IntegPointData,
    mut fn_v: F,
) -> Result<(), StrError>
where
    F: FnMut(&mut Vector, usize, &Vector, &Matrix, &Vector) -> Result<f64, StrError>,
{
    // check
    let nnode_b = pad_b.interp.dim();
    let (space_ndim, nnode) = pad.xxt.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    if nrow_kk < ii0 + nnode * space_ndim {
        return Err("nrow(K) must be ≥ ii0 + pad.nnode ⋅ space_ndim");
    }
    if ncol_kk < jj0 + nnode_b {
        return Err("ncol(K) must be ≥ jj0 + pad_b.nnode");
    }

    // allocate auxiliary vector
    let mut v = Vector::new(space_ndim);

    // clear output matrix
    if clear_kk {
        kk.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate interpolation functions and Jacobian
        (pad.fn_interp)(&mut pad.interp, iota); // N
        let det_jac = pad.calc_gradient(iota)?; // G
        (pad_b.fn_interp)(&mut pad_b.interp, iota); // Nb

        // calculate v
        let nn = &pad.interp;
        let gg = &pad.gradient;
        let nnb = &pad_b.interp;
        let alpha = fn_v(&mut v, p, nn, gg, nnb)?;

        // calculate coefficient
        let c = if axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * pad.xxt[0][m];
            }
            r * det_jac * weight
        } else {
            alpha * det_jac * weight
        };

        // add contribution to K matrix
        if space_ndim == 2 {
            for m in 0..nnode {
                for n in 0..nnode_b {
                    kk[ii0 + 0 + m * 2][jj0 + n] += c * nn[m] * v[0] * nnb[n];
                    kk[ii0 + 1 + m * 2][jj0 + n] += c * nn[m] * v[1] * nnb[n];
                }
            }
        } else {
            for m in 0..nnode {
                for n in 0..nnode_b {
                    kk[ii0 + 0 + m * 3][jj0 + n] += c * nn[m] * v[0] * nnb[n];
                    kk[ii0 + 1 + m * 3][jj0 + n] += c * nn[m] * v[1] * nnb[n];
                    kk[ii0 + 2 + m * 3][jj0 + n] += c * nn[m] * v[2] * nnb[n];
                }
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, AnalyticalQua8, AnalyticalTet4};
    use russell_chk::vec_approx_eq;
    use russell_lab::{Matrix, Vector};

    #[test]
    fn capture_some_errors() {
        let (a, b) = (2.0, 3.0);
        let mut pad_b = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut kk = Matrix::new(8 * 2, 4);
        let mut v = Vector::new(0);
        let nn = Vector::new(0);
        let gg = Matrix::new(0, 0);
        let nnb = Vector::new(0);
        let f = |_v: &mut Vector, _p: usize, _nn: &Vector, _gg: &Matrix, _nnb: &Vector| Ok(1.0);
        assert_eq!(f(&mut v, 0, &nn, &gg, &nnb).unwrap(), 1.0);
        let (clear, axis) = (true, false);
        assert_eq!(
            integ::mat_06_nvn(&mut kk, &mut pad, &mut pad_b, 1, 0, clear, axis, &[], f).err(),
            Some("nrow(K) must be ≥ ii0 + pad.nnode ⋅ space_ndim")
        );
        assert_eq!(
            integ::mat_06_nvn(&mut kk, &mut pad, &mut pad_b, 0, 1, clear, axis, &[], f).err(),
            Some("ncol(K) must be ≥ jj0 + pad_b.nnode")
        );
    }

    #[test]
    fn qua4_qua8_works() {
        let (a, b) = (2.0, 3.0);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_b = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut kk = Matrix::new(8 * 2, 4);
        let (clear, axis) = (true, false);
        let ana = AnalyticalQua8::new(a, b);
        let (v0, v1) = (4.0, 5.0);
        let kk_correct = ana.mat_06_nvn(v0, v1);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [4, 9].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_06_nvn(
                &mut kk,
                &mut pad,
                &mut pad_b,
                0,
                0,
                clear,
                axis,
                ips,
                |v, _, _, _, _| {
                    v[0] = v0;
                    v[1] = v1;
                    Ok(1.0)
                },
            )
            .unwrap();
            // println!("{}", kk);
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }

    #[test]
    fn tet4_tet4_works() {
        let mut pad_b = aux::gen_pad_tet4();
        let mut pad = pad_b.clone();
        let mut kk = Matrix::new(4 * 3, 4);
        let (clear, axis) = (true, false);
        let ana = AnalyticalTet4::new(&pad);
        let (v0, v1, v2) = (4.0, 5.0, 6.0);
        let kk_correct = ana.mat_06_nvn(v0, v1, v2);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-15];
        let selection: Vec<_> = [4].iter().map(|n| integ::points(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_06_nvn(
                &mut kk,
                &mut pad,
                &mut pad_b,
                0,
                0,
                clear,
                axis,
                ips,
                |v, _, _, _, _| {
                    v[0] = v0;
                    v[1] = v1;
                    v[2] = v2;
                    Ok(1.0)
                },
            )
            .unwrap();
            // println!("{}", kk);
            vec_approx_eq(kk.as_data(), kk_correct.as_data(), tol);
        });
    }
}
