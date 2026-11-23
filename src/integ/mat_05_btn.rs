use super::CommonArgs;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::math::SQRT_2;
use russell_lab::{mat_mat_mul, Matrix, Vector};
use russell_tensor::Tensor2;

/// Implements the gradient(Ba) times tensor(T) times shape(N) integration case 05 (e.g., coupling matrix)
///
/// Callback function: `f(T, p, Ba, N, B)`
///
/// **Notes:**
///
/// * `pad` (in `args`) is the **main** pad: it carries the geometry, the integration measure `|J|`,
///   and its own shape functions `N` and gradients `B`
/// * `pad_alt` is the **alternative** pad: only its shape functions are used, never its coordinates,
///   and it may be of any order — lower, equal to, or higher than that of `pad`
/// * `m` ranges over the nodes of `pad_alt` (for `Baᵐ`), whose gradient is taken in the **main** pad's
///   mapping, `Ba = dNa/dξ · J⁻¹`; `n` ranges over the nodes of `pad` (for `Nⁿ`)
/// * For example, `1 ≤ m ≤ 4` with a `pad_alt→Qua4` and `pad→Qua8`, or the reverse, which
///   differentiates the quadratic shape in the 4-node map (a subparametric element)
/// * The number of integration points must consider the nodes of both pads and the expected order of
///   the whole integrand
///
/// Coupling vectors:
///
/// ```text
/// →     ⌠ →
/// Kᵐⁿ = │ Bᵃᵐ ⋅ T Nⁿ α dΩ
///       ⌡       ▔
///       Ωₑ
/// ```
///
/// The numerical integration is (assuming an implicit sum over repeated lower indices):
///
/// ```text
///        nip-1      →       →      →       →
/// Kᵐⁿⱼ ≈   Σ   Bᵃᵐᵢ(ιᵖ) Tᵢⱼ(ιᵖ) Nⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
///         p=0
/// ```
///
/// # Results
///
/// ```text
///     ┌                                        ┐
///     | K⁰⁰₀ K⁰⁰₁ K⁰¹₀ K⁰¹₁ K⁰²₀ K⁰²₁ ··· K⁰ⁿⱼ |  ⟸  ii0
///     | K¹⁰₀ K¹⁰₁ K¹¹₀ K¹¹₁ K¹²₀ K¹²₁ ··· K¹ⁿⱼ |
/// K = | K²⁰₀ K²⁰₁ K²¹₀ K²¹₁ K²²₀ K²²₁ ··· K²ⁿⱼ |
///     |  ··   ··   ··   ··   ··   ··  ···  ··  |    (pad_alt)
///     | Kᵐ⁰₀ Kᵐ⁰₁ Kᵐ¹₀ Kᵐ¹₁ Kᵐ²₀ Kᵐ²₁ ··· Kᵐⁿⱼ |  ⟸  ii
///     └                                        ┘
///        ⇑                                  ⇑
///       jj0                          (pad)  jj := j + n ⋅ space_ndim
///
/// n = jj / space_ndim
/// j = jj % space_ndim
/// ```
///
/// # Arguments
///
/// * `kk` -- A matrix containing all `Kᵐⁿⱼ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `j` corresponds to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + pad_alt.nnode` and `ncol(K) ≥ jj0 + pad.nnode ⋅ space_ndim`
/// * `args` --- Common arguments (`pad` is the **main** scratchpad (modified), which carries the
///   geometry)
/// * `pad_alt` -- Scratchpad of the alternative shape (modified) to compute `Na` and `Ba`; its
///   coordinates are not used
/// * `fn_tt` -- Function `f(T,p,Ba,N,B)` that computes `T(x(ιᵖ))`, given `0 ≤ p ≤ ngauss`,
///   the gradients Ba(ιᵖ) (in the main pad's mapping), shape functions N(ιᵖ), and gradients B(ιᵖ).
///   `T` is set for `space_ndim`.
///
/// # Warning
///
/// Both [crate::shapes::Scratchpad]s must be functions on the same reference element, otherwise
/// **calculation errors will occur**: one `ιᵖ` must mean the same point in both. Their classes and
/// geometric dimensions are checked here; the orders are free, and `pad_alt` may be of lower, equal,
/// or higher order than `pad`.
pub fn mat_05_btn<const N: usize, F>(
    kk: &mut Matrix,
    args: &mut CommonArgs,
    pad_alt: &mut Scratchpad,
    mut fn_tt: F,
) -> Result<(), StrError>
where
    F: FnMut(&mut Tensor2<N>, usize, &Matrix, &Vector, &Matrix) -> Result<(), StrError>,
{
    // check
    let nnode_alt = pad_alt.interp.dim();
    let (space_ndim, _, nnode) = args.pad.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if pad_alt.kind.class() != args.pad.kind.class() || pad_alt.kind.ndim() != args.pad.kind.ndim() {
        return Err("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)");
    }
    if pad_alt.dims().1 != space_ndim {
        return Err("mat_05_btn requires pad_alt with geo_ndim = space_ndim");
    }
    if nrow_kk < ii0 + nnode_alt {
        return Err("nrow(K) must be ≥ ii0 + pad_alt.nnode");
    }
    if ncol_kk < jj0 + nnode * space_ndim {
        return Err("ncol(K) must be ≥ jj0 + pad.nnode ⋅ space_ndim");
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
        pad_alt.calc_interp(iota); // Nb
        let det_jac = args.pad.calc_gradient(iota)?; // N, B, and J of the main pad

        // gradient of the alternative shape in the main pad's mapping: Ba = dNa/dξ · J⁻¹
        (pad_alt.fn_deriv)(&mut pad_alt.deriv, iota);
        mat_mat_mul(&mut pad_alt.gradient, 1.0, &pad_alt.deriv, &args.pad.inv_jacobian, 0.0).unwrap(); // cannot fail because the dimensions are checked above

        // radius @ x(ιᵖ) from the geometry
        let radius = if args.axisymmetric {
            args.pad.calc_radius(iota)?
        } else {
            0.0
        };

        // calculate T tensor
        let ggb = &pad_alt.gradient;
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        fn_tt(&mut tt, p, ggb, nn, bb)?;

        // calculate coefficient
        let c = if args.axisymmetric {
            det_jac * weight * args.alpha * radius
        } else {
            det_jac * weight * args.alpha
        };

        // add contribution to K matrix
        if space_ndim == 2 {
            for m in 0..nnode_alt {
                for n in 0..nnode {
                    kk.add(
                        ii0 + m,
                        jj0 + 0 + n * 2,
                        c * (ggb.get(m, 0) * tt.get(0) + ggb.get(m, 1) * tt.get(3) / s) * nn[n],
                    );
                    kk.add(
                        ii0 + m,
                        jj0 + 1 + n * 2,
                        c * (ggb.get(m, 0) * tt.get(3) / s + ggb.get(m, 1) * tt.get(1)) * nn[n],
                    );
                }
            }
        } else {
            for m in 0..nnode_alt {
                for n in 0..nnode {
                    kk.add(
                        ii0 + m,
                        jj0 + 0 + n * 3,
                        c * (ggb.get(m, 0) * tt.get(0) + ggb.get(m, 1) * tt.get(3) / s + ggb.get(m, 2) * tt.get(5) / s)
                            * nn[n],
                    );
                    kk.add(
                        ii0 + m,
                        jj0 + 1 + n * 3,
                        c * (ggb.get(m, 0) * tt.get(3) / s + ggb.get(m, 1) * tt.get(1) + ggb.get(m, 2) * tt.get(4) / s)
                            * nn[n],
                    );
                    kk.add(
                        ii0 + m,
                        jj0 + 2 + n * 3,
                        c * (ggb.get(m, 0) * tt.get(5) / s + ggb.get(m, 1) * tt.get(4) / s + ggb.get(m, 2) * tt.get(2))
                            * nn[n],
                    );
                }
            }
        }
    }
    Ok(())
}

/// Implements the gradient(B) dot tensor(T) times shape(Na) integration case 05 mapped into `Qᵉ`
///
/// Callback function: `f(T, p, N, B, Na)`
///
/// **Notes:**
///
/// * `pad` (in `args`) is the **main** pad: it carries the vector field addressing the rows of `K`,
///   the geometry, the integration measure `|J|`, and its own shape functions `N` and gradients `B`
/// * `pad_alt` is the **alternative** pad: it provides the scalar shape functions `Na` addressing the
///   columns of `K`; only `Na` is evaluated, never its coordinates, Jacobian, or gradients, and it
///   may be of any order — lower, equal to, or higher than that of `pad`
/// * `m` ranges over the nodes of `pad` (for `Bᵐ`) and `n` over the nodes of `pad_alt` (for `Naⁿ`),
///   as in the result layout below
/// * For example, `0 ≤ m < 8` with a `pad→Qua8` and `pad_alt→Qua4`, or the other way around
/// * The number of integration points must consider the nodes of both pads and the expected order of
///   the whole integrand
///
/// Coupling matrices:
///
/// ```text
/// →     ⌠ →        →
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
///     integ::mat_05_btn_qq::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, _, _, _, _| {
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
pub fn mat_05_btn_qq<const N: usize, F>(
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
        let det_jac = args.pad.calc_gradient(iota)?; // N, B, and J of the main pad
        pad_alt.calc_interp(iota); // Na

        // radius @ x(ιᵖ) from the geometry
        let radius = if args.axisymmetric {
            args.pad.calc_radius(iota)?
        } else {
            0.0
        };

        // calculate T tensor
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        let na = &pad_alt.interp;
        fn_tt(&mut tt, p, nn, bb, na)?;

        // calculate coefficient
        let c = det_jac * weight * args.alpha;

        // add contribution to K matrix
        if args.axisymmetric {
            for m in 0..nnode {
                let q0 = c * radius * (tt.get(0) * bb.get(m, 0) + tt.get(3) * bb.get(m, 1) / s) + c * nn[m] * tt.get(2);
                let q1 = c * radius * (tt.get(3) * bb.get(m, 0) / s + tt.get(1) * bb.get(m, 1));
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

/// Implements the shape(Na) times tensor(T) dot gradient(B) integration case 05 mapped into `Pᵉ`
///
/// Callback function: `f(T, p, N, B, Na)`
///
/// **Notes:**
///
/// * `pad` (in `args`) is the **main** pad: it carries the vector field addressing the columns of
///   `K`, the geometry, the integration measure `|J|`, and its own shape functions `N` and gradients
///   `B`
/// * `pad_alt` is the **alternative** pad: it provides the scalar shape functions `Na` addressing the
///   rows of `K`; only `Na` is evaluated, never its coordinates, Jacobian, or gradients, and it may
///   be of any order — lower, equal to, or higher than that of `pad`
/// * `m` ranges over the nodes of `pad_alt` (for `Naᵐ`) and `n` over the nodes of `pad` (for `Bⁿ`),
///   as in the result layout below
/// * For example, `0 ≤ m < 4` with a `pad_alt→Qua4` and `pad→Qua8`, or `0 ≤ m < 8` the other way around
/// * The number of integration points must consider the nodes of both pads and the expected order of
///   the whole integrand
///
/// Coupling matrices:
///
/// ```text
/// →     ⌠           →
/// Kᵐⁿ = │ Naᵐ T · Bⁿ α dΩ
///       ⌡     ▔
///       Ωₑ
/// ```
///
/// The numerical integration is (note the sum over repeated lower indices):
///
/// ```text
///        nip-1       →        →       →       →
/// Kᵐⁿⱼ ≈   Σ   Naᵐ(ιᵖ) Tⱼₖ(ιᵖ) Bⁿₖ(ιᵖ) |J|(ιᵖ) wᵖ α
///         p=0
/// ```
///
/// # Results
///
/// ```text
///     ┌                                        ┐
///     | K⁰⁰₀ K⁰⁰₁ K⁰¹₀ K⁰¹₁ K⁰²₀ K⁰²₁ ··· K⁰ⁿⱼ |  ⟸  ii0
///     | K¹⁰₀ K¹⁰₁ K¹¹₀ K¹¹₁ K¹²₀ K¹²₁ ··· K¹ⁿⱼ |
/// K = | K²⁰₀ K²⁰₁ K²¹₀ K²¹₁ K²²₀ K²²₁ ··· K²ⁿⱼ |
///     |  ··   ··   ··   ··   ··   ··  ···  ··  |    (pad_alt)
///     | Kᵐ⁰₀ Kᵐ⁰₁ Kᵐ¹₀ Kᵐ¹₁ Kᵐ²₀ Kᵐ²₁ ··· Kᵐⁿⱼ |  ⟸  ii
///     └                                        ┘
///        ⇑                                  ⇑
///       jj0                          (pad)  jj := j + n ⋅ space_ndim
///
/// n = jj / space_ndim
/// j = jj % space_ndim
/// ```
///
/// # Arguments
///
/// * `kk` -- A matrix containing all `Kᵐⁿⱼ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `j` corresponds to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + pad_alt.nnode` and `ncol(K) ≥ jj0 + pad.nnode ⋅ space_ndim`
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
/// component receives the hoop contribution `+ c Naᵐ T₂₂ Nⁿ` (because `ε_θθ = u_r / r`, so the
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
///     let mut kk = Matrix::new(pad.nnode(), pad.nnode() * space_ndim);
///     let mut args = CommonArgs::new(&mut pad, &gauss);
///     integ::mat_05_btn_nm::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, _, _, _, _| {
///         tt.sym_set_std(0, 0, 1.0);
///         tt.sym_set_std(1, 1, 2.0);
///         tt.sym_set_std(2, 2, 3.0);
///         Ok(())
///     })?;
///
///     // solution: since Σₘ Naᵐ = 1 and T · Bⁿ is constant, K = (A/3) (T · Bⁿ) in every row
///     //     ┌       ┐
///     //     │ -¼ -⅓ │
///     // B = │  ¼  0 │
///     //     │  0  ⅓ │
///     //     └       ┘
///     #[rustfmt::skip]
///     let kk_correct = Matrix::from(&[
///         [-0.5, -4.0 / 3.0, 0.5, 0.0, 0.0, 4.0 / 3.0],
///         [-0.5, -4.0 / 3.0, 0.5, 0.0, 0.0, 4.0 / 3.0],
///         [-0.5, -4.0 / 3.0, 0.5, 0.0, 0.0, 4.0 / 3.0],
///     ]);
///     mat_approx_eq(&kk, &kk_correct, 1e-14);
///     Ok(())
/// }
/// ```
pub fn mat_05_btn_nm<const N: usize, F>(
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
    if nrow_kk < ii0 + nnode_alt {
        return Err("nrow(K) must be ≥ ii0 + pad_alt.nnode");
    }
    if ncol_kk < jj0 + nnode * space_ndim {
        return Err("ncol(K) must be ≥ jj0 + pad.nnode ⋅ space_ndim");
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
        let det_jac = args.pad.calc_gradient(iota)?; // N, B, and J of the main pad
        pad_alt.calc_interp(iota); // Na

        // radius @ x(ιᵖ) from the geometry
        let radius = if args.axisymmetric {
            args.pad.calc_radius(iota)?
        } else {
            0.0
        };

        // calculate T tensor
        let nn = &args.pad.interp;
        let bb = &args.pad.gradient;
        let na = &pad_alt.interp;
        fn_tt(&mut tt, p, nn, bb, na)?;

        // calculate coefficient
        let c = det_jac * weight * args.alpha;

        // add contribution to K matrix
        if args.axisymmetric {
            for m in 0..nnode_alt {
                for n in 0..nnode {
                    let q0 =
                        c * radius * (tt.get(0) * bb.get(n, 0) + tt.get(3) * bb.get(n, 1) / s) + c * nn[n] * tt.get(2);
                    let q1 = c * radius * (tt.get(3) * bb.get(n, 0) / s + tt.get(1) * bb.get(n, 1));
                    kk.add(ii0 + m, jj0 + 0 + n * 2, na[m] * q0);
                    kk.add(ii0 + m, jj0 + 1 + n * 2, na[m] * q1);
                }
            }
        } else if space_ndim == 2 {
            for m in 0..nnode_alt {
                for n in 0..nnode {
                    kk.add(
                        ii0 + m,
                        jj0 + 0 + n * 2,
                        na[m] * c * (tt.get(0) * bb.get(n, 0) + tt.get(3) * bb.get(n, 1) / s),
                    );
                    kk.add(
                        ii0 + m,
                        jj0 + 1 + n * 2,
                        na[m] * c * (tt.get(3) * bb.get(n, 0) / s + tt.get(1) * bb.get(n, 1)),
                    );
                }
            }
        } else {
            // 3D: component 0 -> tt0,tt3,tt5 ; component 1 -> tt3,tt1,tt4 ; component 2 -> tt5,tt4,tt2
            // (same Mandel contraction as vec_04_bt::add_to_d, multiplied by na[m])
            for m in 0..nnode_alt {
                for n in 0..nnode {
                    kk.add(
                        ii0 + m,
                        jj0 + 0 + n * 3,
                        na[m]
                            * c
                            * (tt.get(0) * bb.get(n, 0) + tt.get(3) * bb.get(n, 1) / s + tt.get(5) * bb.get(n, 2) / s),
                    );
                    kk.add(
                        ii0 + m,
                        jj0 + 1 + n * 3,
                        na[m]
                            * c
                            * (tt.get(3) * bb.get(n, 0) / s + tt.get(1) * bb.get(n, 1) + tt.get(4) * bb.get(n, 2) / s),
                    );
                    kk.add(
                        ii0 + m,
                        jj0 + 2 + n * 3,
                        na[m]
                            * c
                            * (tt.get(5) * bb.get(n, 0) / s + tt.get(4) * bb.get(n, 1) / s + tt.get(2) * bb.get(n, 2)),
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
    use crate::integ::{self, AnalyticalQua8, AnalyticalTet4, CommonArgs, Gauss};
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

    /// Asserts that mat_05_btn_nm(T) is the transpose of mat_05_btn_qq(T)
    fn check_transpose<const N: usize>(
        pad: &mut Scratchpad,
        pad_alt: &mut Scratchpad,
        gauss: &Gauss,
        x_ips: &[Tensor1],
        tensor_fn: fn(&mut Tensor2<N>, &Tensor1),
        axisymmetric: bool,
        nrow: usize,
        ncol: usize,
    ) {
        let mut kk_nm = Matrix::new(nrow, ncol);
        let mut args = CommonArgs::new(&mut *pad, gauss);
        args.axisymmetric = axisymmetric;
        integ::mat_05_btn_nm::<N, _>(&mut kk_nm, &mut args, pad_alt, |tt, p, _, _, _| {
            tensor_fn(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        let mut kk_qq = Matrix::new(ncol, nrow);
        let mut args = CommonArgs::new(&mut *pad, gauss);
        args.axisymmetric = axisymmetric;
        integ::mat_05_btn_qq::<N, _>(&mut kk_qq, &mut args, pad_alt, |tt, p, _, _, _| {
            tensor_fn(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        let mut max_error = 0.0;
        for i in 0..nrow {
            for j in 0..ncol {
                let diff = f64::abs(kk_qq.get(j, i) - kk_nm.get(i, j));
                max_error = f64::max(max_error, diff);
                approx_eq(kk_qq.get(j, i), kk_nm.get(i, j), 1e-14);
            }
        }
        println!("max |K_nm - K_qqᵀ| = {:.3e}", max_error);
    }

    /// Asserts that summing the pad_alt block of mat_05_btn_nm gives vec_04_bt exactly (Σₘ Naᵐ = 1)
    ///
    /// Both functions use the same Gauss rule and the same integration point coordinates. The matrix
    /// is returned so that the caller can also inspect its magnitude.
    fn check_partition_of_unity_nm<const N: usize>(
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

        // coupled matrix: nnode_alt rows for Naᵐ times nnode ⋅ space_ndim columns for Bⁿ
        let mut kk = Matrix::new(nnode_alt, nnode * space_ndim);
        let mut args = CommonArgs::new(&mut pad_main, &gauss);
        args.axisymmetric = axisymmetric;
        integ::mat_05_btn_nm::<N, _>(&mut kk, &mut args, &mut pad_alt_main, |tt, p, _, _, _| {
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

        // since Σₘ Naᵐ = 1 at every integration point, the sum over the alternative nodes must match.
        // The 1e-13 tolerance (instead of 1e-14) accounts for the axisymmetric entries reaching O(100),
        // where one ulp is already ~3e-14; it is still ~1e-15 in relative terms.
        for n in 0..nnode {
            for j in 0..space_ndim {
                let mut sum = 0.0;
                for m in 0..nnode_alt {
                    sum += kk.get(m, j + n * space_ndim);
                }
                approx_eq(sum, d[j + n * space_ndim], 1e-13);
            }
        }
        kk
    }

    /// Asserts that summing the pad_alt block of mat_05_btn_qq gives vec_04_bt exactly (Σₙ Naⁿ = 1)
    ///
    /// Both functions use the same Gauss rule and the same integration point coordinates. The matrix
    /// is returned so that the caller can also inspect its magnitude.
    fn check_partition_of_unity_qq<const N: usize>(
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
        integ::mat_05_btn_qq::<N, _>(&mut kk, &mut args, &mut pad_alt_main, |tt, p, _, _, _| {
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
        let (a, b) = (2.0, 3.0);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut kk = Matrix::new(4, 8 * 2);
        let mut tt = Tensor2::<4>::new();
        let ggb = Matrix::new(0, 0);
        let nn = Vector::new(0);
        let bb = Matrix::new(0, 0);
        let f = |_tt: &mut Tensor2<4>, _p: usize, _ggb: &Matrix, _nn: &Vector, _bb: &Matrix| Ok(());
        f(&mut tt, 0, &ggb, &nn, &bb).unwrap();
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_05_btn::<4, _>(&mut kk, &mut args, &mut pad_alt, f).err(),
            Some("nrow(K) must be ≥ ii0 + pad_alt.nnode")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_05_btn::<4, _>(&mut kk, &mut args, &mut pad_alt, f).err(),
            Some("ncol(K) must be ≥ jj0 + pad.nnode ⋅ space_ndim")
        );
    }

    #[test]
    fn qua4_qua8_works() {
        let (a, b) = (2.0, 3.0);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut kk = Matrix::new(4, 8 * 2);
        let ana = AnalyticalQua8::new(a, b);
        #[rustfmt::skip]
        let tt = Tensor2::<4>::from_std_matrix(&[
            [2.0, 5.0, 0.0],
            [5.0, 3.0, 0.0],
            [0.0, 0.0, 4.0]]).unwrap();
        let kk_correct = ana.mat_05_btn(&tt);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14, 1e-14];
        let selection: Vec<_> = [4, 9].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_05_btn::<4, _>(&mut kk, &mut args, &mut pad_alt, |ten, _, _, _, _| {
                ten.set_tensor(1.0, &tt);
                Ok(())
            })
            .unwrap();
            // println!("{}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn tet4_tet4_works() {
        let mut pad_alt = aux::gen_pad_tet4();
        let mut pad = pad_alt.clone();
        let mut kk = Matrix::new(4, 4 * 3);
        let ana = AnalyticalTet4::new(&pad);
        #[rustfmt::skip]
        let tt = Tensor2::<6>::from_std_matrix(&[
            [1.0, 4.0, 6.0],
            [4.0, 2.0, 5.0],
            [6.0, 5.0, 3.0]]).unwrap();
        let kk_correct = ana.mat_05_btn(&tt);
        // println!("{}", kk_correct);
        let class = pad.kind.class();
        let tolerances = [1e-14];
        let selection: Vec<_> = [4].iter().map(|n| Gauss::new_sized(class, *n).unwrap()).collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            let mut args = CommonArgs::new(&mut pad, ips);
            integ::mat_05_btn::<6, _>(&mut kk, &mut args, &mut pad_alt, |ten, _, _, _, _| {
                ten.set_tensor(1.0, &tt);
                Ok(())
            })
            .unwrap();
            // println!("{}", kk);
            mat_approx_eq(&kk, &kk_correct, tol);
        });
    }

    #[test]
    fn curved_element_uses_the_main_jacobian() {
        // Curved 8-node element: the midside nodes are pushed off the straight edges, so the
        // element's map and the bilinear (4-node) map are different domains. The reference below
        // is assembled here from dNa/dξ and the main pad's J⁻¹, independently of mat_05_btn.
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
        #[rustfmt::skip]
        let tt = Tensor2::<4>::from_std_matrix(&[
            [2.0, 0.3, 0.0],
            [0.3, 1.0, 0.0],
            [0.0, 0.0, 0.0]]).unwrap();

        // kk_ref: Ba = dNa/dξ · J⁻¹ (main pad's J). kk_own: Ba from the bilinear map (old behavior).
        let mut kk_ref = Matrix::new(4, 16);
        let mut kk_own = Matrix::new(4, 16);
        let mut pad_calc = pad.clone();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            pad_alt.calc_gradient(iota).unwrap(); // bilinear map: its own gradient
            let c = det_jac * gauss.weight(p);
            for m in 0..4 {
                for n in 0..8 {
                    for j in 0..2 {
                        let mut g_ref = 0.0;
                        let mut g_own = 0.0;
                        for k in 0..2 {
                            g_ref += (pad_alt.deriv.get(m, 0) * pad_calc.inv_jacobian.get(0, k)
                                + pad_alt.deriv.get(m, 1) * pad_calc.inv_jacobian.get(1, k))
                                * tt.get_std(k, j);
                            g_own += pad_alt.gradient.get(m, k) * tt.get_std(k, j);
                        }
                        kk_ref.add(m, n * 2 + j, g_ref * pad_calc.interp[n] * c);
                        kk_own.add(m, n * 2 + j, g_own * pad_calc.interp[n] * c);
                    }
                }
            }
        }

        let mut kk = Matrix::new(4, 16);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn::<4, _>(&mut kk, &mut args, &mut pad_alt, |ten, _, _, _, _| {
            ten.set_tensor(1.0, &tt);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the coordinates of pad_alt are not used: the same matrix comes out of a pad without coordinates
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let mut kk_nocoords = Matrix::new(4, 16);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn::<4, _>(&mut kk_nocoords, &mut args, &mut pad_alt_nocoords, |ten, _, _, _, _| {
            ten.set_tensor(1.0, &tt);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk_nocoords, &kk_ref, 1e-13);

        // on this element the old convention differs by an O(0.1) amount, so the check discriminates
        let mut max_diff = 0.0;
        for i in 0..4 {
            for j in 0..16 {
                max_diff = f64::max(max_diff, f64::abs(kk.get(i, j) - kk_own.get(i, j)));
            }
        }
        println!("max |K - K_own| = {:.3e}", max_diff);
        assert!(
            max_diff > 1e-2,
            "test does not discriminate: max_diff = {:.3e}",
            max_diff
        );
    }

    #[test]
    fn reversed_order_gives_the_subparametric_element() {
        // pad_alt = Qua8 (differentiated), main = Qua4 (geometry): the quadratic shape functions are
        // differentiated in the bilinear map. The reference is assembled here from dNb/dξ · J⁻¹.
        #[rustfmt::skip]
        let xx = [
            [0.00, 0.00], [2.00, 0.10], [1.70, 1.30], [-0.20, 1.00],
            [1.00, 0.35], [2.05, 0.72], [0.75, 1.45], [-0.10, 0.55],
        ];
        let mut pad_alt = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for m in 0..8 {
            pad_alt.set_xx(m, 0, xx[m][0]);
            pad_alt.set_xx(m, 1, xx[m][1]);
        }
        let mut pad = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad.set_xx(m, 0, xx[m][0]);
            pad.set_xx(m, 1, xx[m][1]);
        }
        let gauss = Gauss::new_or_sized(GeoKind::Qua4, Some(9)).unwrap();
        #[rustfmt::skip]
        let tt = Tensor2::<4>::from_std_matrix(&[
            [2.0, 0.3, 0.0],
            [0.3, 1.0, 0.0],
            [0.0, 0.0, 0.0]]).unwrap();

        let mut kk_ref = Matrix::new(8, 4 * 2);
        let mut pad_calc = pad.clone();
        let mut pad_deriv = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_deriv.fn_deriv)(&mut pad_deriv.deriv, iota);
            let c = det_jac * gauss.weight(p);
            for m in 0..8 {
                for n in 0..4 {
                    for j in 0..2 {
                        let mut g = 0.0;
                        for k in 0..2 {
                            g += (pad_deriv.deriv.get(m, 0) * pad_calc.inv_jacobian.get(0, k)
                                + pad_deriv.deriv.get(m, 1) * pad_calc.inv_jacobian.get(1, k))
                                * tt.get_std(k, j);
                        }
                        kk_ref.add(m, n * 2 + j, g * pad_calc.interp[n] * c);
                    }
                }
            }
        }

        let mut kk = Matrix::new(8, 4 * 2);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn::<4, _>(&mut kk, &mut args, &mut pad_alt, |ten, _, _, _, _| {
            ten.set_tensor(1.0, &tt);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);
    }

    #[test]
    fn qq_capture_some_errors() {
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
            integ::mat_05_btn_qq::<4, _>(&mut kk, &mut args, &mut pad_alt, ignore_callback::<4>).err(),
            Some("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)")
        );

        // nrow(K) is too small
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_05_btn_qq::<4, _>(&mut kk, &mut args, &mut pad_alt, ignore_callback::<4>).err(),
            Some("nrow(K) must be ≥ ii0 + pad.nnode ⋅ space_ndim")
        );

        // ncol(K) is too small
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.jj0 = 1;
        assert_eq!(
            integ::mat_05_btn_qq::<4, _>(&mut kk, &mut args, &mut pad_alt, ignore_callback::<4>).err(),
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
            integ::mat_05_btn_qq::<6, _>(&mut kk_3d, &mut args, &mut pad_alt_3d, ignore_callback::<6>).err(),
            Some("axisymmetric requires space_ndim = 2")
        );

        // the tensor dimension must match the space dimension
        let mut args = CommonArgs::new(&mut pad_3d, &gauss_3d);
        assert_eq!(
            integ::mat_05_btn_qq::<4, _>(&mut kk_3d, &mut args, &mut pad_alt_3d, ignore_callback::<4>).err(),
            Some("tensor dimension N must equal 2 * space_ndim")
        );

        // the callback error is propagated
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let mut kk = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::mat_05_btn_qq::<4, _>(&mut kk, &mut args, &mut pad_alt, fail_callback::<4>).err(),
            Some("callback error")
        );
    }

    #[test]
    fn nm_capture_some_errors() {
        // exercise the ignore_callback stub itself
        let mut tt = Tensor2::<4>::new();
        let empty_v = Vector::new(0);
        let empty_m = Matrix::new(0, 0);
        ignore_callback::<4>(&mut tt, 0, &empty_v, &empty_m, &empty_v).unwrap();

        let (a, b) = (2.0, 3.0);

        // class mismatch: Qua8 (main) + Tri3 (alt)
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_alt = aux::gen_pad_tri3();
        let mut kk = Matrix::new(4, 16);
        let gauss = Gauss::new(pad.kind);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::mat_05_btn_nm::<4, _>(&mut kk, &mut args, &mut pad_alt, ignore_callback::<4>).err(),
            Some("pad_alt and pad must be of the same class and geo_ndim (e.g. Qua4 with Qua8)")
        );

        // nrow(K) is too small
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_05_btn_nm::<4, _>(&mut kk, &mut args, &mut pad_alt, ignore_callback::<4>).err(),
            Some("nrow(K) must be ≥ ii0 + pad_alt.nnode")
        );

        // ncol(K) is too small
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.jj0 = 1;
        assert_eq!(
            integ::mat_05_btn_nm::<4, _>(&mut kk, &mut args, &mut pad_alt, ignore_callback::<4>).err(),
            Some("ncol(K) must be ≥ jj0 + pad.nnode ⋅ space_ndim")
        );

        // axisymmetric requires space_ndim = 2
        let mut pad_3d = gen_pad_tet10();
        let mut pad_alt_3d = aux::gen_pad_tet4();
        let gauss_3d = Gauss::new(pad_3d.kind);
        let mut kk_3d = Matrix::new(4, 30);
        let mut args = CommonArgs::new(&mut pad_3d, &gauss_3d);
        args.axisymmetric = true;
        assert_eq!(
            integ::mat_05_btn_nm::<6, _>(&mut kk_3d, &mut args, &mut pad_alt_3d, ignore_callback::<6>).err(),
            Some("axisymmetric requires space_ndim = 2")
        );

        // the tensor dimension must match the space dimension
        let mut args = CommonArgs::new(&mut pad_3d, &gauss_3d);
        assert_eq!(
            integ::mat_05_btn_nm::<4, _>(&mut kk_3d, &mut args, &mut pad_alt_3d, ignore_callback::<4>).err(),
            Some("tensor dimension N must equal 2 * space_ndim")
        );

        // the callback error is propagated
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let mut kk = Matrix::new(4, 16);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        assert_eq!(
            integ::mat_05_btn_nm::<4, _>(&mut kk, &mut args, &mut pad_alt, fail_callback::<4>).err(),
            Some("callback error")
        );
    }

    #[test]
    fn qq_equals_mat_07_bsn_when_t_is_scalar_identity() {
        // 2D: Qua8 (main) + Qua4 (alt)
        let (a, b) = (2.0, 3.0);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();
        let mut kk = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn_qq::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
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
        integ::mat_05_btn_qq::<6, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
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
    fn nm_equals_mat_04_nsb_when_t_is_scalar_identity() {
        // 2D: Qua8 (main) + Qua4 (alt)
        let (a, b) = (2.0, 3.0);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();
        let mut kk = Matrix::new(4, 16);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn_nm::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            let s = s_2d(&x_ips[p]);
            tt.sym_set_std(0, 0, s);
            tt.sym_set_std(1, 1, s);
            tt.sym_set_std(2, 2, s);
            Ok(())
        })
        .unwrap();
        let mut kk_ref = Matrix::new(4, 16);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_04_nsb(&mut kk_ref, &mut args, &mut pad_alt, |p, _, _, _| Ok(s_2d(&x_ips[p]))).unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-14);

        // 3D: Tet10 (main) + Tet4 (alt)
        let mut pad = gen_pad_tet10();
        let mut pad_alt = aux::gen_pad_tet4();
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();
        let mut kk = Matrix::new(4, 30);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn_nm::<6, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            let s = s_3d(&x_ips[p]);
            tt.sym_set_std(0, 0, s);
            tt.sym_set_std(1, 1, s);
            tt.sym_set_std(2, 2, s);
            Ok(())
        })
        .unwrap();
        let mut kk_ref = Matrix::new(4, 30);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_04_nsb(&mut kk_ref, &mut args, &mut pad_alt, |p, _, _, _| Ok(s_3d(&x_ips[p]))).unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-14);
    }

    #[test]
    fn qq_curved_element_uses_the_main_jacobian() {
        // Curved 8-node element: the midside nodes are pushed off the straight edges, so the
        // element's map and the bilinear (4-node) map are different domains. The reference below
        // is assembled here from the main pad's interp, gradient, and det(J), independently of
        // mat_05_btn_qq.
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
        integ::mat_05_btn_qq::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the coordinates of pad_alt are not used: the same matrix comes out of a pad without coordinates
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let mut kk_nocoords = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn_qq::<4, _>(&mut kk_nocoords, &mut args, &mut pad_alt_nocoords, |tt, p, _, _, _| {
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
    fn nm_curved_element_uses_the_main_jacobian() {
        // Curved 8-node element: the midside nodes are pushed off the straight edges, so the
        // element's map and the bilinear (4-node) map are different domains. The reference below
        // is assembled here from the main pad's interp, gradient, and det(J), independently of
        // mat_05_btn_nm.
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
        let mut kk_ref = Matrix::new(4, 16);
        let mut kk_wrong = Matrix::new(4, 16);
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
            for m in 0..4 {
                for n in 0..8 {
                    // (T · Bⁿ)ⱼ = Σᵢ Tⱼᵢ Bⁿᵢ, with Tⱼᵢ the standard (non-Mandel) component
                    let mut q = [0.0; 2];
                    for j in 0..2 {
                        for i in 0..2 {
                            q[j] += tt.get_std(j, i) * pad_calc.gradient.get(n, i);
                        }
                    }
                    kk_ref.add(m, n * 2, pad_alt_calc.interp[m] * c * q[0]);
                    kk_ref.add(m, n * 2 + 1, pad_alt_calc.interp[m] * c * q[1]);
                    kk_wrong.add(m, n * 2, pad_alt_calc.interp[m] * c_wrong * q[0]);
                    kk_wrong.add(m, n * 2 + 1, pad_alt_calc.interp[m] * c_wrong * q[1]);
                }
            }
        }

        let mut kk = Matrix::new(4, 16);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn_nm::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the coordinates of pad_alt are not used: the same matrix comes out of a pad without coordinates
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let mut kk_nocoords = Matrix::new(4, 16);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn_nm::<4, _>(&mut kk_nocoords, &mut args, &mut pad_alt_nocoords, |tt, p, _, _, _| {
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
    fn qq_reversed_order_qua4_qua8() {
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
        integ::mat_05_btn_qq::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the coordinates of pad_alt are not used: the same matrix comes out of a pad without coordinates
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        let mut kk_nocoords = Matrix::new(8, 8);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn_qq::<4, _>(&mut kk_nocoords, &mut args, &mut pad_alt_nocoords, |tt, p, _, _, _| {
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
    fn nm_reversed_order_qua4_qua8() {
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
            for m in 0..8 {
                for n in 0..4 {
                    // (T · Bⁿ)ⱼ = Σᵢ Tⱼᵢ Bⁿᵢ, with Tⱼᵢ the standard (non-Mandel) component
                    let mut q = [0.0; 2];
                    for j in 0..2 {
                        for i in 0..2 {
                            q[j] += tt.get_std(j, i) * pad_calc.gradient.get(n, i);
                        }
                    }
                    kk_ref.add(m, n * 2, pad_alt_calc.interp[m] * c * q[0]);
                    kk_ref.add(m, n * 2 + 1, pad_alt_calc.interp[m] * c * q[1]);
                    kk_wrong.add(m, n * 2, pad_alt_calc.interp[m] * c_wrong * q[0]);
                    kk_wrong.add(m, n * 2 + 1, pad_alt_calc.interp[m] * c_wrong * q[1]);
                }
            }
        }

        let mut kk = Matrix::new(8, 8);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn_nm::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);

        // the coordinates of pad_alt are not used: the same matrix comes out of a pad without coordinates
        let mut pad_alt_nocoords = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        let mut kk_nocoords = Matrix::new(8, 8);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn_nm::<4, _>(&mut kk_nocoords, &mut args, &mut pad_alt_nocoords, |tt, p, _, _, _| {
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
    fn qq_axisymmetric_includes_the_hoop_term() {
        // translated so that the radius is positive at every integration point
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua8(5.0, 2.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(5.0, 2.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();
        let mut kk = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.axisymmetric = true;
        integ::mat_05_btn_qq::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
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
            let r = pad_calc.calc_radius(iota).unwrap();
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
    fn nm_axisymmetric_includes_the_hoop_term() {
        // translated so that the radius is positive at every integration point
        let (a, b) = (2.0, 1.5);
        let mut pad = aux::gen_pad_qua8(5.0, 2.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(5.0, 2.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();
        let mut kk = Matrix::new(4, 16);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        args.axisymmetric = true;
        integ::mat_05_btn_nm::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();

        // hand quadrature: the hoop term + c Nⁿ T₂₂ Naᵐ is added to the radial columns
        let mut kk_ref = Matrix::new(4, 16);
        let mut kk_no_hoop = Matrix::new(4, 16);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        let mut tt = Tensor2::<4>::new();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_alt_calc.fn_interp)(&mut pad_alt_calc.interp, iota);
            tensor_2d(&mut tt, &x_ips[p]);
            let r = pad_calc.calc_radius(iota).unwrap();
            assert!(r > 0.0);
            let c = det_jac * gauss.weight(p);
            for m in 0..4 {
                for n in 0..8 {
                    // (T · Bⁿ)ⱼ = Σᵢ Tⱼᵢ Bⁿᵢ, with Tⱼᵢ the standard (non-Mandel) component
                    let mut q = [0.0; 2];
                    for j in 0..2 {
                        for i in 0..2 {
                            q[j] += tt.get_std(j, i) * pad_calc.gradient.get(n, i);
                        }
                    }
                    let q0 = c * r * q[0];
                    let q0_hoop = q0 + c * pad_calc.interp[n] * tt.get_std(2, 2);
                    let q1 = c * r * q[1];
                    kk_ref.add(m, n * 2, pad_alt_calc.interp[m] * q0_hoop);
                    kk_ref.add(m, n * 2 + 1, pad_alt_calc.interp[m] * q1);
                    kk_no_hoop.add(m, n * 2, pad_alt_calc.interp[m] * q0);
                    kk_no_hoop.add(m, n * 2 + 1, pad_alt_calc.interp[m] * q1);
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
    fn transpose_nm_equals_qq() {
        let (a, b) = (2.0, 3.0);

        // 2D, non-axisymmetric
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();
        check_transpose(&mut pad, &mut pad_alt, &gauss, &x_ips, tensor_2d, false, 4, 16);

        // 2D, axisymmetric (r > 0)
        let mut pad = aux::gen_pad_qua8(5.0, 2.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(5.0, 2.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();
        check_transpose(&mut pad, &mut pad_alt, &gauss, &x_ips, tensor_2d, true, 4, 16);

        // 3D
        let mut pad = gen_pad_tet10();
        let mut pad_alt = aux::gen_pad_tet4();
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();
        check_transpose(&mut pad, &mut pad_alt, &gauss, &x_ips, tensor_3d, false, 4, 30);
    }

    #[test]
    fn qq_anisotropic_3d_tet10_tet4() {
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
        integ::mat_05_btn_qq::<6, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_3d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);
    }

    #[test]
    fn nm_anisotropic_3d_tet10_tet4() {
        // Tet10 main pad + Tet4 alternative pad with a full anisotropic tensor. The whole 4 x 30
        // matrix is compared against a reference assembled with the standard (non-Mandel)
        // components, so a swap of T₁₂/T₂₀ or a missing √2 in 3D cannot pass.
        let mut pad = gen_pad_tet10();
        let mut pad_alt = aux::gen_pad_tet4();
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();

        let mut kk_ref = Matrix::new(4, 30);
        let mut pad_calc = pad.clone();
        let mut pad_alt_calc = Scratchpad::new(3, GeoKind::Tet4).unwrap();
        let mut tt = Tensor2::<6>::new();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = pad_calc.calc_gradient(iota).unwrap();
            (pad_alt_calc.fn_interp)(&mut pad_alt_calc.interp, iota);
            tensor_3d(&mut tt, &x_ips[p]);
            let c = det_jac * gauss.weight(p);
            for m in 0..4 {
                for n in 0..10 {
                    // (T · Bⁿ)ⱼ = Σᵢ Tⱼᵢ Bⁿᵢ, with Tⱼᵢ the standard (non-Mandel) component
                    let mut q = [0.0; 3];
                    for j in 0..3 {
                        for i in 0..3 {
                            q[j] += tt.get_std(j, i) * pad_calc.gradient.get(n, i);
                        }
                    }
                    for j in 0..3 {
                        kk_ref.add(m, n * 3 + j, pad_alt_calc.interp[m] * c * q[j]);
                    }
                }
            }
        }

        let mut kk = Matrix::new(4, 30);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn_nm::<6, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_3d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();
        mat_approx_eq(&kk, &kk_ref, 1e-13);
    }

    #[test]
    fn qq_partition_of_unity_matches_vec_04_bt() {
        // 2D: Qua8 (main) + Qua4 (alt)
        check_partition_of_unity_qq(
            &aux::gen_pad_qua8(0.0, 0.0, 2.0, 3.0),
            &aux::gen_pad_qua4(0.0, 0.0, 2.0, 3.0),
            false,
            tensor_2d,
        );

        // axisymmetric 2D with r > 0 at every integration point
        check_partition_of_unity_qq(
            &aux::gen_pad_qua8(5.0, 2.0, 2.0, 3.0),
            &aux::gen_pad_qua4(5.0, 2.0, 2.0, 3.0),
            true,
            tensor_2d,
        );

        // 3D: Tet10 (main) + Tet4 (alt)
        check_partition_of_unity_qq(&gen_pad_tet10(), &aux::gen_pad_tet4(), false, tensor_3d);

        // Deviatoric 2D tensor: tr(T) = 0 with non-zero shear
        let kk_2d = check_partition_of_unity_qq(
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
        let kk_3d = check_partition_of_unity_qq(&gen_pad_tet10(), &aux::gen_pad_tet4(), false, deviatoric_3d);
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
    fn nm_partition_of_unity_matches_vec_04_bt() {
        // 2D: Qua8 (main) + Qua4 (alt)
        check_partition_of_unity_nm(
            &aux::gen_pad_qua8(0.0, 0.0, 2.0, 3.0),
            &aux::gen_pad_qua4(0.0, 0.0, 2.0, 3.0),
            false,
            tensor_2d,
        );

        // axisymmetric 2D with r > 0 at every integration point
        check_partition_of_unity_nm(
            &aux::gen_pad_qua8(5.0, 2.0, 2.0, 3.0),
            &aux::gen_pad_qua4(5.0, 2.0, 2.0, 3.0),
            true,
            tensor_2d,
        );

        // 3D: Tet10 (main) + Tet4 (alt)
        check_partition_of_unity_nm(&gen_pad_tet10(), &aux::gen_pad_tet4(), false, tensor_3d);

        // Deviatoric 2D tensor: tr(T) = 0 with non-zero shear
        let kk_2d = check_partition_of_unity_nm(
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
        // mat_04_nsb sees nothing on the same fixture
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, 2.0, 3.0);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, 2.0, 3.0);
        let gauss = Gauss::new(pad.kind);
        let mut kk_scalar = Matrix::new(4, 16);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_04_nsb(&mut kk_scalar, &mut args, &mut pad_alt, |_, _, _, _| Ok(0.0)).unwrap();
        for i in 0..4 {
            for j in 0..16 {
                assert!(f64::abs(kk_scalar.get(i, j)) < 1e-15);
            }
        }

        // Deviatoric 3D tensor: tr(T) = 0 with non-zero shear
        let kk_3d = check_partition_of_unity_nm(&gen_pad_tet10(), &aux::gen_pad_tet4(), false, deviatoric_3d);
        let max_abs = max_abs_value(&kk_3d);
        println!("max |K| = {:.3e}", max_abs);
        assert!(max_abs > 1e-1);

        // the deviatoric coupling is invisible to the scalar case: with s = tr(T)/3 = 0,
        // mat_04_nsb sees nothing on the same fixture
        let mut pad = gen_pad_tet10();
        let mut pad_alt = aux::gen_pad_tet4();
        let gauss = Gauss::new(pad.kind);
        let mut kk_scalar = Matrix::new(4, 30);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_04_nsb(&mut kk_scalar, &mut args, &mut pad_alt, |_, _, _, _| Ok(0.0)).unwrap();
        for i in 0..4 {
            for j in 0..30 {
                assert!(f64::abs(kk_scalar.get(i, j)) < 1e-15);
            }
        }
    }

    #[test]
    fn qq_offset_and_accumulation_work() {
        let (a, b) = (2.0, 3.0);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();

        // single run: the 16 x 4 block
        let mut kk_single = Matrix::new(16, 4);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn_qq::<4, _>(&mut kk_single, &mut args, &mut pad_alt, |tt, p, _, _, _| {
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
            integ::mat_05_btn_qq::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
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

    #[test]
    fn nm_offset_and_accumulation_work() {
        let (a, b) = (2.0, 3.0);
        let mut pad = aux::gen_pad_qua8(0.0, 0.0, a, b);
        let mut pad_alt = aux::gen_pad_qua4(0.0, 0.0, a, b);
        let gauss = Gauss::new(pad.kind);
        let x_ips = get_points_coords(&mut pad, &gauss).unwrap();

        // single run: the 4 x 16 block
        let mut kk_single = Matrix::new(4, 16);
        let mut args = CommonArgs::new(&mut pad, &gauss);
        integ::mat_05_btn_nm::<4, _>(&mut kk_single, &mut args, &mut pad_alt, |tt, p, _, _, _| {
            tensor_2d(tt, &x_ips[p]);
            Ok(())
        })
        .unwrap();

        // accumulate twice, with an offset, on a prefilled matrix
        let mut kk = Matrix::filled(6, 17, 0.5);
        for i in 2..6 {
            for j in 1..17 {
                kk.set(i, j, 0.0); // start the accumulated block from zero
            }
        }
        for _ in 0..2 {
            let mut args = CommonArgs::new(&mut pad, &gauss);
            args.clear = false;
            args.ii0 = 2;
            args.jj0 = 1;
            integ::mat_05_btn_nm::<4, _>(&mut kk, &mut args, &mut pad_alt, |tt, p, _, _, _| {
                tensor_2d(tt, &x_ips[p]);
                Ok(())
            })
            .unwrap();
        }

        // the border is untouched
        for i in 0..6 {
            for j in 0..17 {
                let in_block = (2..6).contains(&i) && (1..17).contains(&j);
                if !in_block {
                    assert_eq!(kk.get(i, j), 0.5);
                }
            }
        }

        // the 4 x 16 block is twice the single run
        for i in 0..4 {
            for j in 0..16 {
                approx_eq(kk.get(2 + i, 1 + j), 2.0 * kk_single.get(i, j), 1e-14);
            }
        }
    }
}
