use super::CommonArgs;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::math::SQRT_2;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the gradient(B) dot tensor(T) times shape(Na) integration case 12 (e.g., coupling matrix)
///
/// Callback function: `T ← f(T, p, N, B, Na)`
///
/// **Notes:**
///
/// * `pad` (in `args`) is the **main** pad: it carries the geometry, the integration measure `|J|`,
///   and its own shape functions `N` and gradients `B`
/// * `pad_alt` is the **alternative** pad: only its shape functions `Na` are used, never its
///   coordinates, Jacobian, or gradients, and it may be of any order — lower, equal to, or higher
///   than that of `pad`
/// * `m` ranges over the nodes of `pad` (for `Bᵐ`) and `n` over the nodes of `pad_alt` (for `Naⁿ`),
///   as in the result layout below
/// * For example, `0 ≤ m < 8` with a `pad→Qua8` and `pad_alt→Qua4`, or the other way around
/// * The number of integration points must consider the nodes of both pads and the expected order of
///   the whole integrand
///
/// Coupling matrices:
///
/// ```text
/// →     ⌠ →
/// Kᵐⁿ = │ Bᵐ · T Naⁿ α dΩ
///       ⌡       ▔
///       Ωₑ
/// ```
///
/// The numerical integration is (note the sum over repeated lower indices):
///
/// ```text
///        nip-1       →       →        →       →
/// Kᵐⁿⱼ ≈   Σ   Bᵐᵢ(ιᵖ) Tᵢⱼ(ιᵖ) Naⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
///         p=0
/// ```
///
/// # Results
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
///       jj0         (pad_alt) jj
///
/// m = ii / space_ndim
/// i = ii % space_ndim
/// ```
///
/// # Arguments
///
/// * `kk` -- A matrix containing all `Kᵐⁿᵢ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `i` corresponds to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + pad.nnode ⋅ space_ndim` and `ncol(K) ≥ jj0 + pad_alt.nnode`
/// * `args` --- Common arguments (`pad` is the **main** scratchpad (modified), which carries the
///   geometry)
/// * `pad_alt` -- Scratchpad of the alternative shape (modified) to compute `Na`; its
///   coordinates are not used
/// * `fn_tt` -- Function `f(T,p,N,B,Na)` that computes `T(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`,
///   shape functions N(ιᵖ), gradients B(ιᵖ), and shape functions Na(ιᵖ). `T` is a full symmetric
///   tensor set for `space_ndim` and stored in the symmetric Mandel format `Tensor2<N>` with
///   `N = 2 ⋅ space_ndim` components (order `00, 11, 22, 01√2[, 12√2, 20√2]`)
///
/// # Axisymmetric models
///
/// When `space_ndim = 2` and `args.axisymmetric` is true, the integration is taken over one radian,
/// the in-plane terms are scaled by the radius `r = Σₖ Nᵏ x₀ᵏ`, and the radial displacement
/// component receives the hoop contribution `+ c Nᵐ T₂₂ Naⁿ` (because `ε_θθ = u_r / r`, so the
/// `1/r` cancels — exactly as in [crate::integ::vec_04_bt()]).
///
/// # Warning
///
/// Both [crate::shapes::Scratchpad]s must be functions on the same reference element, otherwise
/// **calculation errors will occur**: one `ιᵖ` must mean the same point in both. Their classes and
/// geometric dimensions are checked here; the orders are free, and `pad_alt` may be of lower, equal,
/// or higher order than `pad`.
///
/// # Examples
///
/// ```
/// use gemlab::integ::{self, CommonArgs, Gauss};
/// use gemlab::shapes::{GeoKind, Scratchpad};
/// use gemlab::StrError;
/// use russell_lab::{mat_approx_eq, Matrix};
///
/// fn main() -> Result<(), StrError> {
///     // triangle with area A = 6
///     let space_ndim = 2;
///     let mut pad = Scratchpad::new(space_ndim, GeoKind::Tri3)?;
///     pad.set_xx(0, 0, 2.0);
///     pad.set_xx(0, 1, 3.0);
///     pad.set_xx(1, 0, 6.0);
///     pad.set_xx(1, 1, 3.0);
///     pad.set_xx(2, 0, 2.0);
///     pad.set_xx(2, 1, 6.0);
///     let mut pad_alt = pad.clone();
///
///     // constant tensor T = diag(1, 2, 3)
///     let gauss = Gauss::new(pad.kind);
///     let mut kk = Matrix::new(pad.nnode() * space_ndim, pad.nnode());
///     let mut args = CommonArgs::new(&mut pad, &gauss);
///     integ::mat_12_btna::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, _, _, _, _| {
///         tt.sym_set_std(0, 0, 1.0);
///         tt.sym_set_std(1, 1, 2.0);
///         tt.sym_set_std(2, 2, 3.0);
///         Ok(())
///     })?;
///
///     // solution: since Σₙ Naⁿ = 1 and Bᵐ · T is constant, K = (A/3) (Bᵐ · T) in every column
///     //     ┌       ┐
///     //     │ -¼ -⅓ │
///     // B = │  ¼  0 │
///     //     │  0  ⅓ │
///     //     └       ┘
///     #[rustfmt::skip]
///     let kk_correct = Matrix::from(&[
///         [-0.5, -0.5, -0.5],
///         [-4.0 / 3.0, -4.0 / 3.0, -4.0 / 3.0],
///         [0.5, 0.5, 0.5],
///         [0.0, 0.0, 0.0],
///         [0.0, 0.0, 0.0],
///         [4.0 / 3.0, 4.0 / 3.0, 4.0 / 3.0],
///     ]);
///     mat_approx_eq(&kk, &kk_correct, 1e-14);
///     Ok(())
/// }
/// ```
pub fn mat_12_btna<const N: usize, F>(
    kk: &mut Matrix,
    args: &mut CommonArgs,
    pad_alt: &mut Scratchpad,
    mut fn_tt: F,
) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor2<N>, usize, &Vector, &Matrix, &Vector) -> Result<(), StrError>,
{
    // check
    let nnode_alt = pad_alt.interp.dim();
    let (space_ndim, _, nnode) = args.pad.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if pad_alt.kind.class() != args.pad.kind.class() || pad_alt.kind.ndim() != args.pad.kind.ndim() {
        return Err("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)");
    }
    if nrow_kk < ii0 + nnode * space_ndim {
        return Err("nrow(K) must be ≥ ii0 + pad.nnode ⋅ space_ndim");
    }
    if ncol_kk < jj0 + nnode_alt {
        return Err("ncol(K) must be ≥ jj0 + pad_alt.nnode");
    }
    if args.axisymmetric && space_ndim != 2 {
        return Err("axisymmetric requires space_ndim = 2");
    }
    if N != 2 * space_ndim {
        return Err("tensor dimension N must equal 2 * space_ndim");
    }

    // allocate auxiliary tensor
    let mut tt = Tensor2::<N>::new();

    // clear output matrix
    if args.clear {
        kk.fill(0.0);
    }

    // loop over integration points
    let s = SQRT_2;
    for p in 0..args.gauss.npoint() {
        // ksi coordinates and weight
        let iota = args.gauss.coords(p);
        let weight = args.gauss.weight(p);

        // interpolation functions of both pads and the Jacobian of the main pad
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let det_jac = args.pad.calc_gradient(iota)?; // B of the main pad
        (pad_alt.fn_interp)(&mut pad_alt.interp, iota); // Na

        // calculate T tensor
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        let na = &pad_alt.interp;
        fn_tt(&mut tt, p, nn, bb, na)?;

        // calculate coefficient
        let c = det_jac * weight * args.alpha;

        // add contribution to K matrix
        if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for k in 0..nnode {
                r += nn[k] * args.pad.xxt.get(0, k);
            }
            for m in 0..nnode {
                let q0 = c * r * (tt.get(0) * bb.get(m, 0) + tt.get(3) * bb.get(m, 1) / s) + c * nn[m] * tt.get(2);
                let q1 = c * r * (tt.get(3) * bb.get(m, 0) / s + tt.get(1) * bb.get(m, 1));
                for n in 0..nnode_alt {
                    kk.add(ii0 + 0 + m * 2, jj0 + n, q0 * na[n]);
                    kk.add(ii0 + 1 + m * 2, jj0 + n, q1 * na[n]);
                }
            }
        } else if space_ndim == 2 {
            for m in 0..nnode {
                for n in 0..nnode_alt {
                    kk.add(
                        ii0 + 0 + m * 2,
                        jj0 + n,
                        c * (tt.get(0) * bb.get(m, 0) + tt.get(3) * bb.get(m, 1) / s) * na[n],
                    );
                    kk.add(
                        ii0 + 1 + m * 2,
                        jj0 + n,
                        c * (tt.get(3) * bb.get(m, 0) / s + tt.get(1) * bb.get(m, 1)) * na[n],
                    );
                }
            }
        } else {
            // 3D: component 0 -> tt0,tt3,tt5 ; component 1 -> tt3,tt1,tt4 ; component 2 -> tt5,tt4,tt2
            // (same Mandel contraction as vec_04_bt::add_to_d, multiplied by na[n])
            for m in 0..nnode {
                for n in 0..nnode_alt {
                    kk.add(
                        ii0 + 0 + m * 3,
                        jj0 + n,
                        c * (tt.get(0) * bb.get(m, 0) + tt.get(3) * bb.get(m, 1) / s + tt.get(5) * bb.get(m, 2) / s)
                            * na[n],
                    );
                    kk.add(
                        ii0 + 1 + m * 3,
                        jj0 + n,
                        c * (tt.get(3) * bb.get(m, 0) / s + tt.get(1) * bb.get(m, 1) + tt.get(4) * bb.get(m, 2) / s)
                            * na[n],
                    );
                    kk.add(
                        ii0 + 2 + m * 3,
                        jj0 + n,
                        c * (tt.get(5) * bb.get(m, 0) / s + tt.get(4) * bb.get(m, 1) / s + tt.get(2) * bb.get(m, 2))
                            * na[n],
                    );
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
    use crate::integ::{self, CommonArgs, Gauss};
    use crate::recovery::get_points_coords;
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::StrError;
    use russell_lab::{approx_eq, mat_approx_eq, Matrix, Vector};
    use russell_tensor::{Tensor1, Tensor2};

    /// Callback that does nothing (used in the error tests)
    fn ignore_callback<const N: usize>(
        _tt: &mut Tensor2<N>,
        _p: usize,
        _nn: &Vector,
        _bb: &Matrix,
        _na: &Vector,
    ) -> Result<(), StrError> {
        Ok(())
    }

    /// Callback that fails (used in the error propagation test)
    fn fail_callback<const N: usize>(
        _tt: &mut Tensor2<N>,
        _p: usize,
        _nn: &Vector,
        _bb: &Matrix,
        _na: &Vector,
    ) -> Result<(), StrError> {
        Err("callback error")
    }

    /// Generates a Tet10 pad for tests (straight edges; midside nodes at the edge midpoints)
    fn gen_pad_tet10() -> Scratchpad {
        // corner coordinates as in aux::gen_pad_tet4
        let xx = [[2.0, 3.0, 4.0], [6.0, 3.0, 2.0], [2.0, 5.0, 1.0], [4.0, 3.0, 6.0]];
        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
        let mut pad = Scratchpad::new(3, GeoKind::Tet10).unwrap();
        for m in 0..4 {
            for j in 0..3 {
                pad.set_xx(m, j, xx[m][j]);
            }
        }
        for (i, (a, b)) in edges.iter().enumerate() {
            let m = 4 + i;
            for j in 0..3 {
                pad.set_xx(m, j, 0.5 * (xx[*a][j] + xx[*b][j]));
            }
        }
        pad
    }

    /// Returns the maximum absolute difference between two matrices
    fn max_abs_diff(a: &Matrix, b: &Matrix) -> f64 {
        let (nrow, ncol) = a.dims();
        assert_eq!(b.dims(), (nrow, ncol));
        let mut max_diff = 0.0;
        for i in 0..nrow {
            for j in 0..ncol {
                max_diff = f64::max(max_diff, f64::abs(a.get(i, j) - b.get(i, j)));
            }
        }
        max_diff
    }

    /// Returns the maximum absolute value of a matrix
    fn max_abs_value(a: &Matrix) -> f64 {
        let (nrow, ncol) = a.dims();
        let mut max_abs = 0.0;
        for i in 0..nrow {
            for j in 0..ncol {
                max_abs = f64::max(max_abs, f64::abs(a.get(i, j)));
            }
        }
        max_abs
    }

    /// Returns the spatially varying scalar s(x) = 2 + x₀ - x₁ used in the scalar-identity test (2D)
    fn s_2d(x: &Tensor1) -> f64 {
        2.0 + x.get(0) - x.get(1)
    }

    /// Returns the spatially varying scalar s(x) = 2 + x₀ - x₁ + 0.5 x₂ used in the scalar-identity test (3D)
    fn s_3d(x: &Tensor1) -> f64 {
        2.0 + x.get(0) - x.get(1) + 0.5 * x.get(2)
    }

    /// Sets a spatially varying full symmetric tensor in 2D (T₂₂ is the hoop component)
    fn tensor_2d(tt: &mut Tensor2<4>, x: &Tensor1) {
        let (x0, x1) = (x.get(0), x.get(1));
        tt.sym_set_std(0, 0, 2.0 + x0 - 0.5 * x1);
        tt.sym_set_std(1, 1, 1.0 + 0.3 * x0 + 0.7 * x1);
        tt.sym_set_std(2, 2, 0.5 - 0.2 * x0 + 0.4 * x1);
        tt.sym_set_std(0, 1, 0.4 + 0.1 * x0 + 0.25 * x1);
    }

    /// Sets a spatially varying full symmetric tensor in 3D
    fn tensor_3d(tt: &mut Tensor2<6>, x: &Tensor1) {
        let (x0, x1, x2) = (x.get(0), x.get(1), x.get(2));
        tt.sym_set_std(0, 0, 2.0 + x0 - 0.5 * x1 + 0.25 * x2);
        tt.sym_set_std(1, 1, 1.0 + 0.3 * x0 + 0.7 * x1 - 0.2 * x2);
        tt.sym_set_std(2, 2, 0.5 - 0.4 * x0 + 0.2 * x1 + 0.6 * x2);
        tt.sym_set_std(0, 1, 0.4 + 0.1 * x0 + 0.25 * x1 - 0.15 * x2);
        tt.sym_set_std(1, 2, -0.3 + 0.2 * x0 - 0.1 * x1 + 0.35 * x2);
        tt.sym_set_std(2, 0, 0.6 - 0.25 * x0 + 0.3 * x1 + 0.05 * x2);
    }

    /// Sets a deviatoric tensor in 2D: tr(T) = 0 and non-zero shear (T₂₂ is the hoop component)
    fn deviatoric_2d(tt: &mut Tensor2<4>, _x: &Tensor1) {
        tt.sym_set_std(0, 0, 2.0);
        tt.sym_set_std(1, 1, -1.0);
        tt.sym_set_std(2, 2, -1.0);
        tt.sym_set_std(0, 1, 0.75);
    }

    /// Sets a deviatoric tensor in 3D: tr(T) = 0 and non-zero shear
    fn deviatoric_3d(tt: &mut Tensor2<6>, _x: &Tensor1) {
        tt.sym_set_std(0, 0, 2.0);
        tt.sym_set_std(1, 1, -0.5);
        tt.sym_set_std(2, 2, -1.5);
        tt.sym_set_std(0, 1, 0.75);
        tt.sym_set_std(1, 2, -0.4);
        tt.sym_set_std(2, 0, 0.6);
    }

    /// Asserts that summing the pad_alt block of mat_12_btna gives vec_04_bt exactly (Σₙ Naⁿ = 1)
    ///
    /// Both functions use the same Gauss rule and the same integration point coordinates. The matrix
    /// is returned so that the caller can also inspect its magnitude.
    fn check_partition_of_unity<const N: usize>(
        pad: &Scratchpad,
        pad_alt: &Scratchpad,
        axisymmetric: bool,
        tensor_fn: fn(&mut Tensor2<N>, &Tensor1),
    ) -> Matrix {
        let gauss = Gauss::new(pad.kind);
        let mut pad_main = pad.clone();
        let mut pad_alt_main = pad_alt.clone();
        let x_ips = get_points_coords(&mut pad_main, &gauss).unwrap();
        let (space_ndim, _, nnode) = pad.dims();
        let nnode_alt = pad_alt.nnode();

        // coupled matrix: nnode ⋅ space_ndim rows for Bᵐ times nnode_alt columns for Naⁿ
        let mut kk = Matrix::new(nnode * space_ndim, nnode_alt);
        let mut args = CommonArgs::new(&mut pad_main, &gauss);
        args.axisymmetric = axisymmetric;
        integ::mat_12_btna::<N, _>(&mut kk, &mut args, &mut pad_alt_main, |tt, p, _, _, _| {
            tensor_fn(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();

        // reference: the same tensor through vec_04_bt on the main pad
        let mut d = Vector::new(nnode * space_ndim);
        let mut args = CommonArgs::new(&mut pad_main, &gauss);
        args.axisymmetric = axisymmetric;
        integ::vec_04_bt::<N, _>(&mut d, &mut args, |sig, p, _, _| {
            tensor_fn(sig, &x_ips[p]);
            Ok(())
        })
        .unwrap();

        // since Σₙ Naⁿ = 1 at every integration point, the sum over the alternative nodes must match.
        // The 1e-13 tolerance (instead of 1e-14) accounts for the axisymmetric entries reaching O(100),
        // where one ulp is already ~3e-14; it is still ~1e-15 in relative terms.
        for m in 0..nnode {
            for i in 0..space_ndim {
                let mut sum = 0.0;
                for n in 0..nnode_alt {
                    sum += kk.get(i + m * space_ndim, n);
                }
                approx_eq(sum, d[i + m * space_ndim], 1e-13);
            }
        }
        kk
    }

    #[test]
    fn capture_some_errors() {
        // exercise the ignore_callback stub itself
        let mut tt = Tensor2::<4>::new();
        let empty_v = Vector::new(0);
        let empty_m = Matrix::new(0, 0);
        ignore_callback::<4>(&mut tt, 0, &empty_v, &empty_m, &empty_v).unwrap();

        let (a, b) = (2.0, 3.0);

        // class mismatch: Qua8 (main) + Tri3 (alt)
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_alt = aux::gen_pad_tri3();
        let mut kk = Matrix::new(16, 4);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::mat_12_btna::<4, _>(&mut kk, &mut args, &mut pad_alt, ignore_callback::<4>).err(),
            Some("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)")
        );

        // nrow(K) is too small
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_12_btna::<4, _>(&mut kk, &mut args, &mut pad_alt, ignore_callback::<4>).err(),
            Some("nrow(K) must be ≥ ii0 + pad.nnode ⋅ space_ndim")
        );

        // ncol(K) is too small
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.jj0 = 1;
        assert_eq!(
            integ::mat_12_btna::<4, _>(&mut kk, &mut args, &mut pad_alt, ignore_callback::<4>).err(),
            Some("ncol(K) must be ≥ jj0 + pad_alt.nnode")
        );

        // axisymmetric requires space_ndim = 2
        let mut pad_3d = gen_pad_tet10();
        let mut pad_alt_3d = aux::gen_pad_tet4();
        let gauss_3d = Gauss::new(pad_3d.kind);
        let mut kk_3d = Matrix::new(30, 4);
        let mut args = CommonArgs::new(&mut pad_3d, &gauss_3d);
        args.axisymmetric = true;
        assert_eq!(
            integ::mat_12_btna::<6, _>(&mut kk_3d, &mut args, &mut pad_alt_3d, ignore_callback::<6>).err(),
            Some("axisymmetric requires space_ndim = 2")
        );

        // the tensor dimension must match the space dimension
        let mut args = CommonArgs::new(&mut pad_3d, &gauss_3d);
        assert_eq!(
            integ::mat_12_btna::<4, _>(&mut kk_3d, &mut args, &mut pad_alt_3d, ignore_callback::<4>).err(),
            Some("tensor dimension N must equal 2 * space_ndim")
        );

        // the callback error is propagated
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let mut kk = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::mat_12_btna::<4, _>(&mut kk, &mut args, &mut pad_alt, fail_callback::<4>).err(),
            Some("callback error")
        );
    }

    #[test]
    fn equals_mat_07_bsn_when_t_is_scalar_identity() {
        // 2D: Qua8 (main) + Qua4 (alt)
        let (a, b) = (2.0, 3.0);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();
        let mut kk = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_12_btna::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            let s = s_2d(&x_ips[p]);
            tt.sym_set_std(0, 0, s);
            tt.sym_set_std(1, 1, s);
            tt.sym_set_std(2, 2, s);
            Ok(())
        })
        .unwrap();
        let mut kk_ref = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_07_bsn(&mut kk_ref, &mut args, &mut pad_alt, |p, _, _, _| Ok(s_2d(&x_ips[p]))).unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-14);

        // 3D: Tet10 (main) + Tet4 (alt)
        let mut pad = gen_pad_tet10();
        let mut pad_alt = aux::gen_pad_tet4();
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();
        let mut kk = Matrix::new(30, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_12_btna::<6, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            let s = s_3d(&x_ips[p]);
            tt.sym_set_std(0, 0, s);
            tt.sym_set_std(1, 1, s);
            tt.sym_set_std(2, 2, s);
            Ok(())
        })
        .unwrap();
        let mut kk_ref = Matrix::new(30, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_07_bsn(&mut kk_ref, &mut args, &mut pad_alt, |p, _, _, _| Ok(s_3d(&x_ips[p]))).unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-14);
    }

    #[test]
    fn curved_element_uses_the_main_jacobian() {
        // Curved 8-node element: the midside nodes are pushed off the straight edges, so the
        // element's map and the bilinear (4-node) map are different domains. The reference below
        // is assembled here from the main pad's interp, gradient, and det(J), independently of
        // mat_12_btna.
        #[rustfmt::skip]
        let xx = [
            [0.00, 0.00], [2.00, 0.10], [1.70, 1.30], [-0.20, 1.00],
            [1.00, 0.35], [2.05, 0.72], [0.75, 1.45], [-0.10, 0.55],
        ];
        let mut pad = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for m in 0..8 {
            pad.set_xx(m, 0, xx[m][0]);
            pad.set_xx(m, 1, xx[m][1]);
        }
        let mut pad_alt = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad_alt.set_xx(m, 0, xx[m][0]);
            pad_alt.set_xx(m, 1, xx[m][1]);
        }
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();

        // kk_ref uses the main pad's det(J); kk_wrong uses the alt pad's bilinear det(J)
        let mut kk_ref = Matrix::new(16, 4);
        let mut kk_wrong = Matrix::new(16, 4);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad_alt_calc.set_xx(m, 0, xx[m][0]);
            pad_alt_calc.set_xx(m, 1, xx[m][1]);
        }
        let mut tt = Tensor2::<4>::new();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            let det_jac_wrong = pad_alt_calc.calc_gradient(iota).unwrap();
            tensor_2d(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            let c_wrong = det_jac_wrong * gauss.weight(p);
            for m in 0..8 {
                for n in 0..4 {
                    // (Bᵐ · T)ᵢ = Σₖ Bᵐₖ Tₖᵢ, with Tₖᵢ the standard (non-Mandel) component
                    let mut q = [0.0; 2];
                    for i in 0..2 {
                        for k in 0..2 {
                            q[i] += pad_calc.gradient.get(m, k) * tt.get_std(k, i);
                        }
                    }
                    kk_ref.add(m * 2, n, pad_alt_calc.interp[n] * c * q[0]);
                    kk_ref.add(m * 2 + 1, n, pad_alt_calc.interp[n] * c * q[1]);
                    kk_wrong.add(m * 2, n, pad_alt_calc.interp[n] * c_wrong * q[0]);
                    kk_wrong.add(m * 2 + 1, n, pad_alt_calc.interp[n] * c_wrong * q[1]);
                }
            }
        }

        let mut kk = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_12_btna::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the coordinates of pad_alt are not used: the same matrix comes out of a pad without coordinates
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let mut kk_nocoords = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_12_btna::<4, _>(&mut kk_nocoords, &mut args, &mut pad_alt_nocoords, |tt, p, _, _, _| {
            tensor_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk_nocoords, &kk, 1e-15);

        // on this element the bilinear measure differs by an O(0.1) amount, so the check discriminates
        let max_diff = max_abs_diff(&kk, &kk_wrong);
        println!("max |K - K_wrong| = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn reversed_order_qua4_qua8() {
        // pad_alt = Qua8 (shapes only) and main = Qua4 (geometry): any order combination is allowed.
        // The reference below is assembled here from the main pad's interp, gradient, and det(J).
        #[rustfmt::skip]
        let xx = [
            [0.00, 0.00], [2.00, 0.10], [1.70, 1.30], [-0.20, 1.00],
            [1.00, 0.35], [2.05, 0.72], [0.75, 1.45], [-0.10, 0.55],
        ];
        let mut pad = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad.set_xx(m, 0, xx[m][0]);
            pad.set_xx(m, 1, xx[m][1]);
        }
        let mut pad_alt = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for m in 0..8 {
            pad_alt.set_xx(m, 0, xx[m][0]);
            pad_alt.set_xx(m, 1, xx[m][1]);
        }
        let gauss = Gauss::new_or_sized(GeoKind::Qua4, Some(9)).unwrap();
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();

        let mut kk_ref = Matrix::new(8, 8);
        let mut kk_wrong = Matrix::new(8, 8);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for m in 0..8 {
            pad_alt_calc.set_xx(m, 0, xx[m][0]);
            pad_alt_calc.set_xx(m, 1, xx[m][1]);
        }
        let mut tt = Tensor2::<4>::new();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            let det_jac_wrong = pad_alt_calc.calc_gradient(iota).unwrap();
            tensor_2d(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            let c_wrong = det_jac_wrong * gauss.weight(p);
            for m in 0..4 {
                for n in 0..8 {
                    // (Bᵐ · T)ᵢ = Σₖ Bᵐₖ Tₖᵢ, with Tₖᵢ the standard (non-Mandel) component
                    let mut q = [0.0; 2];
                    for i in 0..2 {
                        for k in 0..2 {
                            q[i] += pad_calc.gradient.get(m, k) * tt.get_std(k, i);
                        }
                    }
                    kk_ref.add(m * 2, n, pad_alt_calc.interp[n] * c * q[0]);
                    kk_ref.add(m * 2 + 1, n, pad_alt_calc.interp[n] * c * q[1]);
                    kk_wrong.add(m * 2, n, pad_alt_calc.interp[n] * c_wrong * q[0]);
                    kk_wrong.add(m * 2 + 1, n, pad_alt_calc.interp[n] * c_wrong * q[1]);
                }
            }
        }

        let mut kk = Matrix::new(8, 8);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_12_btna::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the coordinates of pad_alt are not used: the same matrix comes out of a pad without coordinates
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        let mut kk_nocoords = Matrix::new(8, 8);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_12_btna::<4, _>(&mut kk_nocoords, &mut args, &mut pad_alt_nocoords, |tt, p, _, _, _| {
            tensor_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk_nocoords, &kk, 1e-15);

        // on this element the quadratic measure differs by an O(0.1) amount, so the check discriminates
        let max_diff = max_abs_diff(&kk, &kk_wrong);
        println!("max |K - K_wrong| = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn axisymmetric_includes_the_hoop_term() {
        // translated so that the radius is positive at every integration point
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua8(5.0, 2.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(5.0, 2.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();
        let mut kk = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.axisymmetric = true;
        integ::mat_12_btna::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();

        // hand quadrature: the hoop term + c Nᵐ T₂₂ Naⁿ is added to the radial rows
        let mut kk_ref = Matrix::new(16, 4);
        let mut kk_no_hoop = Matrix::new(16, 4);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let mut tt = Tensor2::<4>::new();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_alt_calc.fn_interp)(&mut pad_alt_calc.interp, iota);
            tensor_2d(&mut tt, &x_ips[p]);
            let mut r = 0.0;
            for k in 0..8 {
                r += pad_calc.interp[k] * pad_calc.xxt.get(0, k);
            }
            assert!(r > 0.0);
            let c = det_jac * gauss.weight(p);
            for m in 0..8 {
                // (Bᵐ · T)ᵢ = Σₖ Bᵐₖ Tₖᵢ, with Tₖᵢ the standard (non-Mandel) component
                let mut q = [0.0; 2];
                for i in 0..2 {
                    for k in 0..2 {
                        q[i] += pad_calc.gradient.get(m, k) * tt.get_std(k, i);
                    }
                }
                let q0 = c * r * q[0];
                let q0_hoop = q0 + c * pad_calc.interp[m] * tt.get_std(2, 2);
                let q1 = c * r * q[1];
                for n in 0..4 {
                    kk_ref.add(m * 2, n, pad_alt_calc.interp[n] * q0_hoop);
                    kk_ref.add(m * 2 + 1, n, pad_alt_calc.interp[n] * q1);
                    kk_no_hoop.add(m * 2, n, pad_alt_calc.interp[n] * q0);
                    kk_no_hoop.add(m * 2 + 1, n, pad_alt_calc.interp[n] * q1);
                }
            }
        }
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // omitting the T₂₂ term differs by an O(0.1) amount, so the check discriminates
        let max_diff = max_abs_diff(&kk, &kk_no_hoop);
        println!("max |K - K_no_hoop| = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn anisotropic_3d_tet10_tet4() {
        // Tet10 main pad + Tet4 alternative pad with a full anisotropic tensor. The whole 30 x 4
        // matrix is compared against a reference assembled with the standard (non-Mandel)
        // components, so a swap of T₁₂/T₂₀ or a missing √2 in 3D cannot pass.
        let mut pad = gen_pad_tet10();
        let mut pad_alt = aux::gen_pad_tet4();
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();

        let mut kk_ref = Matrix::new(30, 4);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(3, GeoKind::Tet4).unwrap();
        let mut tt = Tensor2::<6>::new();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_alt_calc.fn_interp)(&mut pad_alt_calc.interp, iota);
            tensor_3d(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            for m in 0..10 {
                for n in 0..4 {
                    // (Bᵐ · T)ᵢ = Σₖ Bᵐₖ Tₖᵢ, with Tₖᵢ the standard (non-Mandel) component
                    let mut q = [0.0; 3];
                    for i in 0..3 {
                        for k in 0..3 {
                            q[i] += pad_calc.gradient.get(m, k) * tt.get_std(k, i);
                        }
                    }
                    for i in 0..3 {
                        kk_ref.add(i + m * 3, n, c * q[i] * pad_alt_calc.interp[n]);
                    }
                }
            }
        }

        let mut kk = Matrix::new(30, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_12_btna::<6, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_3d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);
    }

    #[test]
    fn partition_of_unity_matches_vec_04_bt() {
        // 2D: Qua8 (main) + Qua4 (alt)
        check_partition_of_unity(
            &aux::gen_pad_qua8(0.0, 0.0, 2.0, 3.0),
            &aux::gen_pad_qua4(0.0, 0.0, 2.0, 3.0),
            false,
            tensor_2d,
        );

        // axisymmetric 2D with r > 0 at every integration point
        check_partition_of_unity(
            &aux::gen_pad_qua8(5.0, 2.0, 2.0, 3.0),
            &aux::gen_pad_qua4(5.0, 2.0, 2.0, 3.0),
            true,
            tensor_2d,
        );

        // 3D: Tet10 (main) + Tet4 (alt)
        check_partition_of_unity(&gen_pad_tet10(), &aux::gen_pad_tet4(), false, tensor_3d);

        // Deviatoric 2D tensor: tr(T) = 0 with non-zero shear
        let kk_2d = check_partition_of_unity(
            &aux::gen_pad_qua8(0.0, 0.0, 2.0, 3.0),
            &aux::gen_pad_qua4(0.0, 0.0, 2.0, 3.0),
            false,
            deviatoric_2d,
        );
        // the coupled block is not the zero matrix even though tr(T) = 0
        let max_abs = max_abs_value(&kk_2d);
        println!("max |K| = {:.3e}", max_abs);
        assert!(max_abs > 1e-1);

        // the deviatoric coupling is invisible to the scalar case: with s = tr(T)/3 = 0,
        // mat_07_bsn sees nothing on the same fixture
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, 2.0, 3.0);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, 2.0, 3.0);
        let gauss = Gauss::new(pad.kind);
        let mut kk_scalar = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_07_bsn(&mut kk_scalar, &mut args, &mut pad_alt, |_, _, _, _| Ok(0.0)).unwrap();
        for i in 0..16 {
            for j in 0..4 {
                assert!(f64::abs(kk_scalar.get(i, j)) < 1e-15);
            }
        }

        // Deviatoric 3D tensor: tr(T) = 0 with non-zero shear
        let kk_3d = check_partition_of_unity(&gen_pad_tet10(), &aux::gen_pad_tet4(), false, deviatoric_3d);
        let max_abs = max_abs_value(&kk_3d);
        println!("max |K| = {:.3e}", max_abs);
        assert!(max_abs > 1e-1);

        // the deviatoric coupling is invisible to the scalar case: with s = tr(T)/3 = 0,
        // mat_07_bsn sees nothing on the same fixture
        let mut pad = gen_pad_tet10();
        let mut pad_alt = aux::gen_pad_tet4();
        let gauss = Gauss::new(pad.kind);
        let mut kk_scalar = Matrix::new(30, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_07_bsn(&mut kk_scalar, &mut args, &mut pad_alt, |_, _, _, _| Ok(0.0)).unwrap();
        for i in 0..30 {
            for j in 0..4 {
                assert!(f64::abs(kk_scalar.get(i, j)) < 1e-15);
            }
        }
    }

    #[test]
    fn offset_and_accumulation_work() {
        let (a, b) = (2.0, 3.0);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();

        // single run: the 16 x 4 block
        let mut kk_single = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_12_btna::<4, _>(&mut kk_single, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();

        // accumulate twice, with an offset, on a prefilled matrix
        let mut kk = Matrix::filled(17, 6, 0.5);
        for i in 1..17 {
            for j in 2..6 {
                kk.set(i, j, 0.0); // start the accumulated block from zero
            }
        }
        for _ in 0..2 {
            let mut args = CommonArgs::new(&mut pad, &gauss);
            args.clear = false;
            args.ii0 = 1;
            args.jj0 = 2;
            integ::mat_12_btna::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
                tensor_2d(tt, &x_ips[p]);
                Ok(())
            })
            .unwrap();
        }

        // the border is untouched
        for i in 0..17 {
            for j in 0..6 {
                let in_block = (1..17).contains(&i) && (2..6).contains(&j);
                if !in_block {
                    assert_eq!(kk.get(i, j), 0.5);
                }
            }
        }

        // the 16 x 4 block is twice the single run
        for i in 0..16 {
            for j in 0..4 {
                approx_eq(kk.get(1 + i, 2 + j), 2.0 * kk_single.get(i, j), 1e-14);
            }
        }
    }
}
