use super::{GeoKind, Scratchpad};
use crate::StrError;
use russell_lab::mat_mat_mul;

/// Holds a **main** pad and an **alternative** pad evaluated at the same reference coordinates
///
/// This struct supports two-field elements whose fields are interpolated by different kinds of shape
/// functions over the same cell. At a given reference coordinate ξ, such an element needs:
///
/// * the **main** shape functions `N_main` and gradients `B_main = dN_main/dξ · J⁻¹`, together with
///   `det(J)` for the integration weight;
/// * the **alternative** shape functions `N_alt` and the gradient taken in the main pad's mapping
///
/// ```text
/// B_alt = dN_alt/dξ · J⁻¹
/// ```
///
/// where `J` is the **main** pad's Jacobian. Differentiating the alternative shape in its own map
/// would use two different maps for the two fields: on a straight-edged cell the two coincide, but on
/// a curved one they describe different domains.
///
/// # Kinds and orders
///
/// `main_kind` and `alt_kind` may be given in any combination of orders — the alternative may be of
/// lower, equal, or higher order than the main one — provided both are of the same class and
/// geometric dimension (for example `Qua8` with `Qua4`, or `Qua4` with `Qua8`). The main kind carries
/// the geometry and is usually the mesh cell's kind.
///
/// # Coordinates
///
/// [MixedPad::new] allocates both pads but leaves their coordinates unset; the caller must fill the
/// **main** pad's coordinates with [Scratchpad::set_xx], or use [crate::mesh::Mesh::get_mixed_pad],
/// which does it for a cell. The alternative pad's coordinates are never read.
pub struct MixedPad {
    /// The main pad: it carries the geometry, the integration measure and its own shape functions
    pub main: Scratchpad,

    /// The alternative pad: only its shape functions are used, differentiated in the main mapping
    pub alt: Scratchpad,
}

impl MixedPad {
    /// Allocates the pair of pads
    ///
    /// **Note:** The coordinates are left unset. Use [Scratchpad::set_xx] on `main` afterwards, or
    /// [crate::mesh::Mesh::get_mixed_pad], which does it for a cell. The coordinates of `alt` are not
    /// used by [MixedPad::calc].
    ///
    /// # Input
    ///
    /// * `space_ndim` -- space dimension (2 or 3); must equal `main_kind`'s geometric dimension
    /// * `main_kind` -- the kind of the main shape (e.g. [GeoKind::Qua8]); it carries the geometry
    /// * `alt_kind` -- the kind of the alternative shape (e.g. [GeoKind::Qua4]); any order, but of the
    ///   same class and geometric dimension as `main_kind`
    ///
    /// # Returns
    ///
    /// Returns an error if `main_kind` does not match `space_ndim`
    /// ("MixedPad requires main_kind with geo_ndim = space_ndim"), if `alt_kind` is not of the same
    /// class and geometric dimension as `main_kind`
    /// ("MixedPad requires main_kind and alt_kind of the same class and geo_ndim (e.g. Qua8 and
    /// Qua4)"), or if the scratchpads cannot be allocated (see [Scratchpad::new]).
    pub fn new(space_ndim: usize, main_kind: GeoKind, alt_kind: GeoKind) -> Result<Self, StrError> {
        if main_kind.ndim() != space_ndim {
            return Err("MixedPad requires main_kind with geo_ndim = space_ndim");
        }
        if alt_kind.class() != main_kind.class() || alt_kind.ndim() != main_kind.ndim() {
            return Err("MixedPad requires main_kind and alt_kind of the same class and geo_ndim (e.g. Qua8 and Qua4)");
        }
        Ok(MixedPad {
            main: Scratchpad::new(space_ndim, main_kind)?,
            alt: Scratchpad::new(space_ndim, alt_kind)?,
        })
    }

    /// Fills both pads at ξ and returns the main pad's `det(J)`
    ///
    /// After this call:
    ///
    /// ```text
    /// main.interp, main.gradient        -- the main shape
    /// main.jacobian, main.inv_jacobian  -- the main (geometric) map
    /// alt.interp                        -- the alternative shape values
    /// alt.deriv                         -- dN_alt/dξ
    /// alt.gradient                      -- dN_alt/dξ · J_main⁻¹
    /// ```
    ///
    /// The alternative pad's `jacobian` and `inv_jacobian` are neither computed nor used.
    ///
    /// # Input
    ///
    /// * `ksi` -- reference coordinates ξ with len ≥ `geo_ndim`, as in [Scratchpad::calc_gradient]
    ///
    /// # Returns
    ///
    /// Returns the determinant of the main pad's Jacobian (the integration weight is `det(J) · w`,
    /// with `w` from [crate::integ::Gauss]), or an error if the main coordinates have not been set or
    /// the Jacobian cannot be inverted.
    pub fn calc(&mut self, ksi: &[f64]) -> Result<f64, StrError> {
        // main pad: its own map, N, B and det(J)
        let det_jac = self.main.calc_gradient(ksi)?;

        // alternative shape: values and ξ-derivatives only (its coordinates are not used)
        (self.alt.fn_interp)(&mut self.alt.interp, ksi);
        (self.alt.fn_deriv)(&mut self.alt.deriv, ksi);

        // gradient in the main pad's mapping: B_alt = dN_alt/dξ · J⁻¹
        mat_mat_mul(
            &mut self.alt.gradient,
            1.0,
            &self.alt.deriv,
            &self.main.inv_jacobian,
            0.0,
        )
        .unwrap(); // cannot fail because the dimensions are fixed by the two scratchpads
        Ok(det_jac)
    }

    /// Returns the number of nodes of the main shape
    #[inline]
    pub fn nnode_main(&self) -> usize {
        self.main.nnode()
    }

    /// Returns the number of nodes of the alternative shape
    #[inline]
    pub fn nnode_alt(&self) -> usize {
        self.alt.nnode()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::MixedPad;
    use crate::integ::Gauss;
    use crate::mesh::Samples;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::{approx_eq, mat_approx_eq};

    /// Coefficients of the linear field used in the patch tests: a₀ + a · x
    const A: [f64; 3] = [2.0, -3.0, 0.5];

    /// Returns the linear field a₀ + a · x at the point with coordinates x
    fn linear_field(x: &[f64]) -> f64 {
        let mut value = 1.0;
        for j in 0..x.len() {
            value += A[j] * x[j];
        }
        value
    }

    /// Computes the area of a straight-edged 2D element with the shoelace formula on its `ncorners`
    fn shoelace(pad: &Scratchpad, ncorners: usize) -> f64 {
        let mut sum = 0.0;
        for m in 0..ncorners {
            let next = if m + 1 == ncorners { 0 } else { m + 1 };
            sum += pad.xxt.get(0, m) * pad.xxt.get(1, next) - pad.xxt.get(0, next) * pad.xxt.get(1, m);
        }
        0.5 * sum
    }

    /// Checks the invariants of the pair at every Gauss point of the main kind
    ///
    /// Requires the alternative nodes to be the first nodes of the main pad (the usual mixed element).
    ///
    /// 1. `Σ N = 1` for both shapes (machine precision)
    /// 2. `Σ B_alt = 0` per column (rigid-body translation)
    /// 3. a linear field is reproduced by `N_alt` and its gradient by `B_alt`
    /// 4. `Σ_i Σ_j ∫ N_i N_j dΩ = measure` and `Σ_i Σ_j ∫ N_i B_j dΩ = 0`
    fn check_invariants(mixed: &mut MixedPad, ngauss: usize, measure: f64) {
        let space_ndim = mixed.main.space_ndim();
        let n_main = mixed.nnode_main();
        let k = mixed.nnode_alt();
        assert!(k <= n_main);
        let gauss = Gauss::new_or_sized(mixed.main.kind, Some(ngauss)).unwrap();

        // nodal coordinates of the alternative shape (its nodes are the first k nodes of the main one)
        let mut xr = vec![vec![0.0; space_ndim]; k];
        for m in 0..k {
            for j in 0..space_ndim {
                xr[m][j] = mixed.main.xxt.get(j, m);
            }
        }

        // values of the linear field at the alternative nodes
        let phi: Vec<f64> = (0..k).map(|m| linear_field(&xr[m])).collect();

        // integral accumulators: Σ_ij ∫ N_i N_j dΩ and Σ_ij ∫ N_i B_j dΩ
        let mut norm_int = 0.0;
        let mut mixed_int = vec![0.0; space_ndim];

        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = mixed.calc(iota).unwrap();
            assert!(det_jac > 0.0);
            let w = det_jac * gauss.weight(p);

            // 1. partition of unity of both shapes
            let mut sum_main = 0.0;
            for m in 0..n_main {
                sum_main += mixed.main.interp[m];
            }
            approx_eq(sum_main, 1.0, 1e-14);
            let mut sum_alt = 0.0;
            for m in 0..k {
                sum_alt += mixed.alt.interp[m];
            }
            approx_eq(sum_alt, 1.0, 1e-14);

            // 2. the alternative gradient sums to zero over the nodes
            for j in 0..space_ndim {
                let mut sum = 0.0;
                for m in 0..k {
                    sum += mixed.alt.gradient.get(m, j);
                }
                approx_eq(sum, 0.0, 1e-14);
            }

            // coordinates of the current point, from the main shape's own interpolation
            let mut x = vec![0.0; space_ndim];
            for m in 0..n_main {
                for j in 0..space_ndim {
                    x[j] += mixed.main.interp[m] * mixed.main.xxt.get(j, m);
                }
            }

            // 3. the linear field is reproduced at x(ξ), and its gradient is exact
            let mut value = 0.0;
            for m in 0..k {
                value += mixed.alt.interp[m] * phi[m];
            }
            approx_eq(value, linear_field(&x), 1e-13);
            for j in 0..space_ndim {
                let mut g = 0.0;
                for m in 0..k {
                    g += mixed.alt.gradient.get(m, j) * phi[m];
                }
                approx_eq(g, A[j], 1e-13);
            }

            // 4. accumulate the two integrals
            for i in 0..k {
                for j in 0..k {
                    norm_int += mixed.alt.interp[i] * mixed.alt.interp[j] * w;
                }
                for j in 0..n_main {
                    for d in 0..space_ndim {
                        mixed_int[d] += mixed.alt.interp[i] * mixed.main.gradient.get(j, d) * w;
                    }
                }
            }
        }

        // Σ_ij ∫ N_i N_j dΩ = ∫ (Σ N)² dΩ = measure
        approx_eq(norm_int, measure, 1e-12);

        // Σ_ij ∫ N_i B_j dΩ = ∫ (Σ N) (Σ B) dΩ = 0
        for d in 0..space_ndim {
            approx_eq(mixed_int[d], 0.0, 1e-12);
        }
    }

    #[test]
    fn get_mixed_pad_qua8() {
        // unit square: Qua8 (main) + Qua4 (alt)
        let mesh = Samples::one_qua8();
        let mut mixed = mesh.get_mixed_pad(0, GeoKind::Qua4);
        assert_eq!(mixed.nnode_main(), 8);
        assert_eq!(mixed.nnode_alt(), 4);
        check_invariants(&mut mixed, 9, 1.0);
    }

    #[test]
    fn get_mixed_pad_tri6() {
        // triangle with base 1 and height 0.85: Tri6 (main) + Tri3 (alt), area = 0.425
        let mesh = Samples::one_tri6();
        let mut mixed = mesh.get_mixed_pad(0, GeoKind::Tri3);
        assert_eq!(mixed.nnode_main(), 6);
        assert_eq!(mixed.nnode_alt(), 3);
        check_invariants(&mut mixed, 7, 0.425);
    }

    #[test]
    fn get_mixed_pad_hex20() {
        // unit cube: Hex20 (main) + Hex8 (alt), volume = 1
        let mesh = Samples::one_hex20(0.0, 0.0, 0.0);
        let mut mixed = mesh.get_mixed_pad(0, GeoKind::Hex8);
        assert_eq!(mixed.nnode_main(), 20);
        assert_eq!(mixed.nnode_alt(), 8);
        check_invariants(&mut mixed, 27, 1.0);
    }

    #[test]
    fn distorted_qua8_matches_the_reference_values() {
        // Distorted quadrilateral with straight edges: the midside nodes are the edge midpoints, so
        // the area is exactly the shoelace area of the corners (2.195), as measured in the
        // investigation recorded in the plan note.
        let corners = [(0.0, 0.0), (2.0, 0.1), (1.7, 1.3), (-0.2, 1.0)];
        let mut mixed = MixedPad::new(2, GeoKind::Qua8, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            mixed.main.set_xx(m, 0, corners[m].0);
            mixed.main.set_xx(m, 1, corners[m].1);
        }
        // midside nodes: 4 = middle(0,1), 5 = middle(1,2), 6 = middle(2,3), 7 = middle(3,0)
        let edges = [(0, 1), (1, 2), (2, 3), (3, 0)];
        for (i, (a, b)) in edges.iter().enumerate() {
            let m = 4 + i;
            mixed.main.set_xx(m, 0, 0.5 * (corners[*a].0 + corners[*b].0));
            mixed.main.set_xx(m, 1, 0.5 * (corners[*a].1 + corners[*b].1));
        }

        // the area comes from the corners only (straight edges)
        approx_eq(shoelace(&mixed.main, 4), 2.195, 1e-14);

        // 9 Gauss points, as in the investigation
        check_invariants(&mut mixed, 9, 2.195);
    }

    #[test]
    fn equal_order_works() {
        // Qua4 with Qua4: the alternative shape is the main one and the mixed gradient is the plain
        // gradient
        let mut mixed = MixedPad::new(2, GeoKind::Qua4, GeoKind::Qua4).unwrap();
        assert_eq!(mixed.nnode_main(), 4);
        assert_eq!(mixed.nnode_alt(), 4);

        // skewed quadrilateral
        let xx = [(0.0, 0.0), (2.0, 0.1), (1.7, 1.3), (-0.2, 1.0)];
        for m in 0..4 {
            mixed.main.set_xx(m, 0, xx[m].0);
            mixed.main.set_xx(m, 1, xx[m].1);
        }

        // compare with a plain scratchpad at the same ξ
        let at_ksi = [-0.5, 0.3];
        let det_jac = mixed.calc(&at_ksi).unwrap();
        let mut pad = Scratchpad::new(2, GeoKind::Qua4).unwrap();
        for m in 0..4 {
            pad.set_xx(m, 0, xx[m].0);
            pad.set_xx(m, 1, xx[m].1);
        }
        let det_jac_ref = pad.calc_gradient(&at_ksi).unwrap();
        approx_eq(det_jac, det_jac_ref, 1e-15);
        mat_approx_eq(&mixed.main.gradient, &pad.gradient, 1e-15);
        mat_approx_eq(&mixed.alt.gradient, &pad.gradient, 1e-15);
        for m in 0..4 {
            approx_eq(mixed.alt.interp[m], pad.interp[m], 1e-15);
        }

        // the invariants also hold for the equal-order case
        check_invariants(&mut mixed, 9, 2.195);
    }

    #[test]
    fn subparametric_main_qua4_alt_qua8() {
        // The main pad is the geometry (Qua4); the alternative is quadratic. Its gradient must be
        // dN_alt/dξ · J_main⁻¹, and its coordinates are irrelevant (they are never set here).
        let corners = [(0.0, 0.0), (2.0, 0.1), (1.7, 1.3), (-0.2, 1.0)];
        let mut mixed = MixedPad::new(2, GeoKind::Qua4, GeoKind::Qua8).unwrap();
        assert_eq!(mixed.nnode_main(), 4);
        assert_eq!(mixed.nnode_alt(), 8);
        for m in 0..4 {
            mixed.main.set_xx(m, 0, corners[m].0);
            mixed.main.set_xx(m, 1, corners[m].1);
        }

        // the alternative derivatives are dNa/dξ; the reference multiplies them by the main's J⁻¹
        let mut pad_alt = Scratchpad::new(2, GeoKind::Qua8).unwrap();
        let gauss = Gauss::new_or_sized(GeoKind::Qua4, Some(9)).unwrap();
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            mixed.calc(iota).unwrap();
            (pad_alt.fn_deriv)(&mut pad_alt.deriv, iota);

            let mut sum_alt = 0.0;
            for m in 0..8 {
                sum_alt += mixed.alt.interp[m];
            }
            approx_eq(sum_alt, 1.0, 1e-14);

            for m in 0..8 {
                for d in 0..2 {
                    let mut b = 0.0;
                    for j in 0..2 {
                        b += pad_alt.deriv.get(m, j) * mixed.main.inv_jacobian.get(j, d);
                    }
                    approx_eq(mixed.alt.gradient.get(m, d), b, 1e-15);
                }
            }
            for d in 0..2 {
                let mut sum = 0.0;
                for m in 0..8 {
                    sum += mixed.alt.gradient.get(m, d);
                }
                approx_eq(sum, 0.0, 1e-13);
            }
        }
    }

    #[test]
    fn mixed_pad_handles_errors() {
        // the main kind must match the space dimension
        assert_eq!(
            MixedPad::new(2, GeoKind::Lin2, GeoKind::Lin2).err(),
            Some("MixedPad requires main_kind with geo_ndim = space_ndim")
        );
        assert_eq!(
            MixedPad::new(3, GeoKind::Tri3, GeoKind::Tri3).err(),
            Some("MixedPad requires main_kind with geo_ndim = space_ndim")
        );
        assert_eq!(
            MixedPad::new(4, GeoKind::Hex8, GeoKind::Hex8).err(),
            Some("MixedPad requires main_kind with geo_ndim = space_ndim")
        );

        // the alternative must be of the same class and geometric dimension
        assert_eq!(
            MixedPad::new(2, GeoKind::Qua4, GeoKind::Tri3).err(),
            Some("MixedPad requires main_kind and alt_kind of the same class and geo_ndim (e.g. Qua8 and Qua4)")
        );

        // error from Scratchpad::new is propagated (a line in 1D)
        assert_eq!(
            MixedPad::new(1, GeoKind::Lin2, GeoKind::Lin2).err(),
            Some("space_ndim must be 2 or 3")
        );

        // coordinates not set
        let mut mixed = MixedPad::new(2, GeoKind::Qua8, GeoKind::Qua4).unwrap();
        assert_eq!(
            mixed.calc(&[0.0, 0.0]).err(),
            Some("all components of the coordinates matrix must be set first")
        );
    }
}
