use super::{GeoKind, Scratchpad};
use crate::StrError;
use russell_lab::mat_mat_mul;

/// Holds a full-order Scratchpad and its reduced (LBB) counterpart evaluated at the same point
///
/// This struct supports mixed-order (two-field) elements whose lower-order field lives on the
/// corner nodes of a quadratic cell. At a given reference coordinate ξ, such an element needs:
///
/// * the **full** shape functions `N_u` and the full gradient `B_u = L_u · J⁻¹`, together with
///   `det(J)` for the integration weight;
/// * the **reduced** shape functions `N_φ` and the **mixed** gradient
///
/// ```text
/// B_φ = L_φ · J⁻¹
/// ```
///
/// where `L_φ = dN_φ/dξ` comes from the reduced shape but `J⁻¹` comes from the **full** shape.
/// Using the reduced shape's own Jacobian for `B_φ` is wrong: the two fields would then be
/// differentiated with respect to two different mappings.
///
/// # Reduced kind
///
/// The reduced kind is `full_kind.lower_order().unwrap_or(full_kind)`, and it takes the first
/// `k = nnode_reduced()` nodes of the full shape. gemlab numbers the corner nodes first in every
/// quadratic kind, thus the first `k` nodes are the corners (e.g. Qua8 → Qua4 with nodes 0..3,
/// Tri6 → Tri3 with nodes 0..2, Hex20 → Hex8 with nodes 0..7). When `lower_order()` returns
/// `None` (a linear kind such as Qua4), the reduced kind **is** the full kind, and the mixed
/// gradient coincides with the plain gradient.
///
/// # Restrictions
///
/// This struct requires `geo_ndim == space_ndim` (2D or 3D solids); shells and cables are not
/// supported.
///
/// # Coordinates
///
/// [MixedPad::new] allocates both pads but leaves their coordinates unset; the caller must fill
/// them with [Scratchpad::set_xx] or use [crate::mesh::Mesh::get_mixed_pad], which sets both
/// coordinate sets from a cell. The reduced pad holds only its `k` corner coordinates.
pub struct MixedPad {
    /// The full-order shape, which defines the Jacobian and the integration weight
    pub full: Scratchpad,

    /// The reduced shape, whose gradient uses the full shape's inverse Jacobian
    pub reduced: Scratchpad,
}

impl MixedPad {
    /// Allocates the pair; the reduced kind is `full_kind.lower_order().unwrap_or(full_kind)`
    ///
    /// The reduced pad takes the first `k = reduced.nnode()` nodes of the full one (gemlab's
    /// convention: corner nodes come first in every quadratic kind).
    ///
    /// **Note:** The coordinates are left unset. Use [Scratchpad::set_xx] on `full` and on
    /// `reduced` afterwards, or [crate::mesh::Mesh::get_mixed_pad], which does it for a cell.
    ///
    /// # Input
    ///
    /// * `space_ndim` -- space dimension (2 or 3); must equal the shape's `geo_ndim`
    /// * `full_kind` -- the kind of the full-order shape (e.g. [GeoKind::Qua8])
    ///
    /// # Returns
    ///
    /// Returns an error if `geo_ndim != space_ndim` ("MixedPad requires geo_ndim = space_ndim")
    /// or if the scratchpads cannot be allocated (see [Scratchpad::new]).
    pub fn new(space_ndim: usize, full_kind: GeoKind) -> Result<Self, StrError> {
        if full_kind.ndim() != space_ndim {
            return Err("MixedPad requires geo_ndim = space_ndim");
        }
        let reduced_kind = full_kind.lower_order().unwrap_or(full_kind);
        Ok(MixedPad {
            full: Scratchpad::new(space_ndim, full_kind)?,
            reduced: Scratchpad::new(space_ndim, reduced_kind)?,
        })
    }

    /// Fills both pads at ξ and returns the full shape's det(J)
    ///
    /// After this call:
    ///
    /// ```text
    /// full.interp, full.gradient       -- the full shape
    /// reduced.interp                   -- the reduced shape values
    /// reduced.deriv                    -- dN/dξ of the reduced shape
    /// reduced.gradient                 -- L_reduced · J_full⁻¹
    /// ```
    ///
    /// **Note:** `reduced.jacobian` and `reduced.inv_jacobian` are filled as a by-product of
    /// [Scratchpad::calc_jacobian] with the reduced shape's **own** Jacobian, which is discarded:
    /// they must not be used to differentiate the reduced field.
    ///
    /// # Input
    ///
    /// * `ksi` -- reference coordinates ξ with len ≥ `geo_ndim`, as in [Scratchpad::calc_gradient]
    ///
    /// # Returns
    ///
    /// Returns the determinant of the full shape's Jacobian matrix (the integration weight is
    /// `det(J) · w`, with `w` from [crate::integ::Gauss]), or an error if the coordinates have not
    /// been set or the Jacobian cannot be inverted.
    pub fn calc(&mut self, ksi: &[f64]) -> Result<f64, StrError> {
        // full shape: N and B = L · J⁻¹, plus det(J)
        let det_jac = self.full.calc_gradient(ksi)?;

        // reduced shape: N and L only (the reduced Jacobian is discarded)
        self.reduced.calc_interp(ksi);
        self.reduced.calc_jacobian(ksi)?;

        // mixed gradient: B = L_reduced · J_full⁻¹
        mat_mat_mul(
            &mut self.reduced.gradient,
            1.0,
            &self.reduced.deriv,
            &self.full.inv_jacobian,
            0.0,
        )
        .unwrap(); // cannot fail because the dimensions are fixed by the two scratchpads
        Ok(det_jac)
    }

    /// Returns the number of nodes of the full-order shape
    #[inline]
    pub fn nnode_full(&self) -> usize {
        self.full.nnode()
    }

    /// Returns the number of nodes of the reduced shape
    #[inline]
    pub fn nnode_reduced(&self) -> usize {
        self.reduced.nnode()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::MixedPad;
    use crate::integ::Gauss;
    use crate::mesh::{Mesh, Samples};
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::StrError;
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

    /// Computes the area of a straight-edged 2D element with the shoelace formula on its corners
    fn shoelace(pad: &Scratchpad) -> f64 {
        let n = pad.kind.lower_order().unwrap_or(pad.kind).nnode();
        let mut sum = 0.0;
        for m in 0..n {
            let next = if m + 1 == n { 0 } else { m + 1 };
            sum += pad.xxt.get(0, m) * pad.xxt.get(1, next) - pad.xxt.get(0, next) * pad.xxt.get(1, m);
        }
        0.5 * sum
    }

    /// Checks the four invariants of the mixed pair at every Gauss point of the full kind
    ///
    /// 1. `Σ N = 1` for both shapes (machine precision)
    /// 2. `Σ B_reduced = 0` per column (rigid-body translation)
    /// 3. a linear field is reproduced by `N_reduced` and its gradient by `B_reduced`
    /// 4. `Σ_i Σ_j ∫ N_i N_j dΩ = measure` and `Σ_i Σ_j ∫ N_i B_j dΩ = 0`
    fn check_invariants(mixed: &mut MixedPad, ngauss: usize, measure: f64) {
        let space_ndim = mixed.full.space_ndim();
        let n_full = mixed.nnode_full();
        let k = mixed.nnode_reduced();
        let gauss = Gauss::new_or_sized(mixed.full.kind, Some(ngauss)).unwrap();

        // nodal coordinates of the reduced shape (its nodes are the first k nodes of the full one)
        let mut xr = vec![vec![0.0; space_ndim]; k];
        for m in 0..k {
            for j in 0..space_ndim {
                xr[m][j] = mixed.full.xxt.get(j, m);
            }
        }

        // values of the linear field at the reduced nodes
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
            let mut sum_full = 0.0;
            for m in 0..n_full {
                sum_full += mixed.full.interp[m];
            }
            approx_eq(sum_full, 1.0, 1e-14);
            let mut sum_reduced = 0.0;
            for m in 0..k {
                sum_reduced += mixed.reduced.interp[m];
            }
            approx_eq(sum_reduced, 1.0, 1e-14);

            // 2. the reduced gradient sums to zero over the nodes (it differentiates a constant
            //    field to zero, which is the rigid-body translation condition)
            for j in 0..space_ndim {
                let mut sum = 0.0;
                for m in 0..k {
                    sum += mixed.reduced.gradient.get(m, j);
                }
                approx_eq(sum, 0.0, 1e-14);
            }

            // coordinates of the current point, from the full shape's own interpolation
            let mut x = vec![0.0; space_ndim];
            for m in 0..n_full {
                for j in 0..space_ndim {
                    x[j] += mixed.full.interp[m] * mixed.full.xxt.get(j, m);
                }
            }

            // 3. the linear field is reproduced at x(ξ), and its gradient is exact
            let mut value = 0.0;
            for m in 0..k {
                value += mixed.reduced.interp[m] * phi[m];
            }
            approx_eq(value, linear_field(&x), 1e-13);
            for j in 0..space_ndim {
                let mut g = 0.0;
                for m in 0..k {
                    g += mixed.reduced.gradient.get(m, j) * phi[m];
                }
                approx_eq(g, A[j], 1e-13);
            }

            // 4. accumulate the two integrals
            for i in 0..k {
                for j in 0..k {
                    norm_int += mixed.reduced.interp[i] * mixed.reduced.interp[j] * w;
                }
                for j in 0..n_full {
                    for d in 0..space_ndim {
                        mixed_int[d] += mixed.reduced.interp[i] * mixed.full.gradient.get(j, d) * w;
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

    /// Returns a mesh with a single cell and a MixedPad for it
    fn build(mesh: &Mesh) -> MixedPad {
        mesh.get_mixed_pad(0)
    }

    #[test]
    fn get_mixed_pad_qua8() {
        // unit square: Qua8 (full) + Qua4 (reduced)
        let mesh = Samples::one_qua8();
        let mut mixed = build(&mesh);
        assert_eq!(mixed.nnode_full(), 8);
        assert_eq!(mixed.nnode_reduced(), 4);
        check_invariants(&mut mixed, 9, 1.0);
    }

    #[test]
    fn get_mixed_pad_tri6() {
        // triangle with base 1 and height 0.85: Tri6 (full) + Tri3 (reduced), area = 0.425
        let mesh = Samples::one_tri6();
        let mut mixed = build(&mesh);
        assert_eq!(mixed.nnode_full(), 6);
        assert_eq!(mixed.nnode_reduced(), 3);
        check_invariants(&mut mixed, 7, 0.425);
    }

    #[test]
    fn get_mixed_pad_hex20() {
        // unit cube: Hex20 (full) + Hex8 (reduced), volume = 1
        let mesh = Samples::one_hex20(0.0, 0.0, 0.0);
        let mut mixed = build(&mesh);
        assert_eq!(mixed.nnode_full(), 20);
        assert_eq!(mixed.nnode_reduced(), 8);
        check_invariants(&mut mixed, 27, 1.0);
    }

    #[test]
    fn distorted_qua8_matches_the_reference_values() {
        // Distorted quadrilateral with straight edges: the midside nodes are the edge midpoints, so
        // the area is exactly the shoelace area of the corners (2.195), as measured in the
        // investigation recorded in the plan note.
        let corners = [(0.0, 0.0), (2.0, 0.1), (1.7, 1.3), (-0.2, 1.0)];
        let mut mixed = MixedPad::new(2, GeoKind::Qua8).unwrap();
        for m in 0..4 {
            mixed.full.set_xx(m, 0, corners[m].0);
            mixed.full.set_xx(m, 1, corners[m].1);
            mixed.reduced.set_xx(m, 0, corners[m].0);
            mixed.reduced.set_xx(m, 1, corners[m].1);
        }
        // midside nodes: 4 = middle(0,1), 5 = middle(1,2), 6 = middle(2,3), 7 = middle(3,0)
        let edges = [(0, 1), (1, 2), (2, 3), (3, 0)];
        for (i, (a, b)) in edges.iter().enumerate() {
            let m = 4 + i;
            mixed.full.set_xx(m, 0, 0.5 * (corners[*a].0 + corners[*b].0));
            mixed.full.set_xx(m, 1, 0.5 * (corners[*a].1 + corners[*b].1));
        }

        // the area comes from the corners only (straight edges)
        approx_eq(shoelace(&mixed.full), 2.195, 1e-14);

        // 9 Gauss points, as in the investigation
        check_invariants(&mut mixed, 9, 2.195);

        // measure the residuals quoted in the plan note, for the record
        let gauss = Gauss::new_or_sized(GeoKind::Qua8, Some(9)).unwrap();
        let k = mixed.nnode_reduced();
        let mut phi = vec![0.0; k];
        for m in 0..k {
            phi[m] = linear_field(&[mixed.full.xxt.get(0, m), mixed.full.xxt.get(1, m)]);
        }
        let mut norm_int = 0.0;
        let mut mixed_int = [0.0; 2];
        let (mut err_sum_n, mut err_sum_b) = (0.0, 0.0);
        let (mut err_value, mut err_grad) = (0.0, 0.0);
        for p in 0..gauss.npoint() {
            let iota = gauss.coords(p);
            let det_jac = mixed.calc(iota).unwrap();
            let w = det_jac * gauss.weight(p);

            let mut sum_n = 0.0;
            for m in 0..k {
                sum_n += mixed.reduced.interp[m];
            }
            err_sum_n = f64::max(err_sum_n, f64::abs(sum_n - 1.0));
            for j in 0..2 {
                let mut sum_b = 0.0;
                for m in 0..k {
                    sum_b += mixed.reduced.gradient.get(m, j);
                }
                err_sum_b = f64::max(err_sum_b, f64::abs(sum_b));
            }

            let mut x = [0.0; 2];
            for m in 0..mixed.nnode_full() {
                for j in 0..2 {
                    x[j] += mixed.full.interp[m] * mixed.full.xxt.get(j, m);
                }
            }
            let mut value = 0.0;
            for m in 0..k {
                value += mixed.reduced.interp[m] * phi[m];
            }
            err_value = f64::max(err_value, f64::abs(value - linear_field(&x)));
            for j in 0..2 {
                let mut g = 0.0;
                for m in 0..k {
                    g += mixed.reduced.gradient.get(m, j) * phi[m];
                }
                err_grad = f64::max(err_grad, f64::abs(g - A[j]));
            }

            for i in 0..k {
                for j in 0..k {
                    norm_int += mixed.reduced.interp[i] * mixed.reduced.interp[j] * w;
                }
                for j in 0..mixed.nnode_full() {
                    for d in 0..2 {
                        mixed_int[d] += mixed.reduced.interp[i] * mixed.full.gradient.get(j, d) * w;
                    }
                }
            }
        }
        println!("area (shoelace)      = {:.6}", shoelace(&mixed.full));
        println!("max |ΣN_reduced - 1| = {:.3e}", err_sum_n);
        println!("max |ΣB_reduced|     = {:.3e}", err_sum_b);
        println!("max field error      = {:.3e}", err_value);
        println!("max gradient error   = {:.3e}", err_grad);
        println!("∫N∘N∘ dΩ             = {:.6}", norm_int);
        println!("Σ∫N∘B dΩ             = ({:.3e}, {:.3e})", mixed_int[0], mixed_int[1]);

        // the four invariants at the same tolerances as above
        approx_eq(err_sum_n, 0.0, 1e-14); // ΣN = 1 to machine precision
        approx_eq(err_sum_b, 0.0, 1e-14);
        approx_eq(err_value, 0.0, 1e-13);
        approx_eq(err_grad, 0.0, 1e-13);
        approx_eq(norm_int, 2.195, 1e-12);
        approx_eq(mixed_int[0], 0.0, 1e-12);
        approx_eq(mixed_int[1], 0.0, 1e-12);
    }

    #[test]
    fn mixed_pad_without_lower_order_duplicates_the_kind() -> Result<(), StrError> {
        // Qua4 has no lower order: the reduced kind is the full kind and the mixed gradient is the
        // plain gradient
        let mut mixed = MixedPad::new(2, GeoKind::Qua4)?;
        assert_eq!(mixed.nnode_full(), 4);
        assert_eq!(mixed.nnode_reduced(), 4);

        // skewed quadrilateral
        let xx = [(0.0, 0.0), (2.0, 0.1), (1.7, 1.3), (-0.2, 1.0)];
        for m in 0..4 {
            mixed.full.set_xx(m, 0, xx[m].0);
            mixed.full.set_xx(m, 1, xx[m].1);
            mixed.reduced.set_xx(m, 0, xx[m].0);
            mixed.reduced.set_xx(m, 1, xx[m].1);
        }

        // compare with a plain scratchpad at the same ξ
        let at_ksi = [-0.5, 0.3];
        let det_jac = mixed.calc(&at_ksi)?;
        let mut pad = Scratchpad::new(2, GeoKind::Qua4)?;
        for m in 0..4 {
            pad.set_xx(m, 0, xx[m].0);
            pad.set_xx(m, 1, xx[m].1);
        }
        let det_jac_ref = pad.calc_gradient(&at_ksi)?;
        approx_eq(det_jac, det_jac_ref, 1e-15);
        mat_approx_eq(&mixed.full.gradient, &pad.gradient, 1e-15);
        mat_approx_eq(&mixed.reduced.gradient, &pad.gradient, 1e-15);
        for m in 0..4 {
            approx_eq(mixed.reduced.interp[m], pad.interp[m], 1e-15);
        }

        // the invariants also hold for the duplicated kind
        check_invariants(&mut mixed, 9, 2.195);
        Ok(())
    }

    #[test]
    fn mixed_pad_handles_errors() {
        // shells and cables are not supported
        assert_eq!(
            MixedPad::new(2, GeoKind::Lin2).err(),
            Some("MixedPad requires geo_ndim = space_ndim")
        );
        assert_eq!(
            MixedPad::new(3, GeoKind::Tri3).err(),
            Some("MixedPad requires geo_ndim = space_ndim")
        );
        assert_eq!(
            MixedPad::new(4, GeoKind::Hex8).err(),
            Some("MixedPad requires geo_ndim = space_ndim")
        );

        // error from Scratchpad::new is propagated (a line in 1D)
        assert_eq!(MixedPad::new(1, GeoKind::Lin2).err(), Some("space_ndim must be 2 or 3"));

        // coordinates not set
        let mut mixed = MixedPad::new(2, GeoKind::Qua8).unwrap();
        assert_eq!(
            mixed.calc(&[0.0, 0.0]).err(),
            Some("all components of the coordinates matrix must be set first")
        );
    }
}
