use super::Scratchpad;
use crate::StrError;
use russell_lab::{mat_inverse, mat_mat_mul};

/// Indicates that the determinant of the Jacobian is not available (e.g., Shells)
///
/// **Note:** This also indicates that the inverse Jacobian matrix has not been computed
pub const DET_JAC_NOT_AVAILABLE: f64 = -1.0;

impl Scratchpad {
    /// Calculates the Jacobian of the mapping from general to reference space
    ///
    /// The components of the Jacobian matrix are
    ///
    /// ```text
    ///        ∂xᵢ
    /// Jᵢⱼ := ——— = Σ X[m][i] * L[m][j]
    ///        ∂ξⱼ   m
    /// ```
    ///
    /// Thus, in matrix notation
    ///
    /// ```text
    /// jacobian := J = Xᵀ · L
    /// jacobian := Jcable = Xᵀ · L
    /// jacobian := Jshell = Xᵀ · L
    /// ```
    ///
    /// where:
    ///
    /// * `Jcable`` -- Jacobian for line in multi-dimensions (geom_ndim < space_ndim)
    /// * `Jshell`` -- Jacobian for 3D surfaces (geo_ndim = 2 and space_ndim = 3)
    ///
    /// For the `SOLID` case (`geo_ndim = space_ndim`), this function also computes the inverse Jacobian:
    ///
    /// ```text
    /// inv_jacobian := J⁻¹
    /// ```
    ///
    /// # Output
    ///
    /// * `deriv` -- derivatives of the interpolation functions (nnode); `L` matrix
    /// * `jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    /// * `inv_jacobian` -- If `geo_ndim = space_ndim` (`SOLID` case): inverse Jacobian matrix (space_ndim,space_ndim)
    /// * Returns one of the following:
    ///     * `CABLE`: (geo_ndim = 1 and space_ndim = 2 or 3), returns the norm of the Jacobian vector
    ///     * `SHELL`: (geo_ndim = 2 and space_ndim = 3), returns [DET_JAC_NOT_AVAILABLE] indicating that the
    ///        determinant of the Jacobian is not available and the inverse Jacobian has not been computed
    ///     * `SOLID`: (geo_ndim = space_ndim), returns the determinant of the Jacobian
    ///
    /// # Input
    ///
    /// * `ksi` -- reference coordinates ξ with len ≥ geo_ndim
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Scratchpad};
    /// use gemlab::StrError;
    /// use russell_lab::approx_eq;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     //  3-------------2         ξ₀   ξ₁
    ///     //  |      ξ₁     |  node    r    s
    ///     //  |      |      |     0 -1.0 -1.0
    ///     //  |      +--ξ₀  |     1  1.0 -1.0
    ///     //  |             |     2  1.0  1.0
    ///     //  |             |     3 -1.0  1.0
    ///     //  0-------------1
    ///
    ///     let a = 3.0;
    ///     let space_ndim = 2;
    ///     let mut pad = Scratchpad::new(space_ndim, GeoKind::Qua4)?;
    ///     pad.set_xx(0, 0, 0.0);
    ///     pad.set_xx(0, 1, 0.0);
    ///     pad.set_xx(1, 0, 2.0 * a);
    ///     pad.set_xx(1, 1, 0.0);
    ///     pad.set_xx(2, 0, 2.0 * a);
    ///     pad.set_xx(2, 1, a);
    ///     pad.set_xx(3, 0, 0.0);
    ///     pad.set_xx(3, 1, a);
    ///
    ///     let det_jac = pad.calc_jacobian(&[0.0, 0.0])?;
    ///     approx_eq(det_jac, a * a / 2.0, 1e-15);
    ///
    ///     // the solution is
    ///     //  ┌         ┐
    ///     //  │  a   0  │
    ///     //  │  0  a/2 │
    ///     //  └         ┘
    ///     assert_eq!(
    ///         format!("{}", pad.jacobian),
    ///         "┌         ┐\n\
    ///          │   3   0 │\n\
    ///          │   0 1.5 │\n\
    ///          └         ┘"
    ///     );
    ///     Ok(())
    /// }
    /// ```
    pub fn calc_jacobian(&mut self, ksi: &[f64]) -> Result<f64, StrError> {
        // check
        if !self.ok_xxt {
            return Err("all components of the coordinates matrix must be set first");
        }

        // matrix L: dNᵐ/dξ
        (self.fn_deriv)(&mut self.deriv, ksi);

        // matrix J: dx/dξ
        mat_mat_mul(&mut self.jacobian, 1.0, &self.xxt, &self.deriv)?;

        // inverse Jacobian and determinant/norm (or not possible)
        let (space_ndim, geo_ndim) = self.jacobian.dims();
        if geo_ndim == space_ndim {
            // SOLID case: inverse J (returns determinant)
            mat_inverse(&mut self.inv_jacobian, &self.jacobian)
        } else {
            // CABLE case: norm of Jacobian vector
            if geo_ndim == 1 {
                let mut norm_jac = 0.0;
                for i in 0..space_ndim {
                    norm_jac += self.jacobian.get(i, 0) * self.jacobian.get(i, 0);
                }
                Ok(f64::sqrt(norm_jac))
            } else {
                // SHELL case
                Ok(DET_JAC_NOT_AVAILABLE)
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::DET_JAC_NOT_AVAILABLE;
    use crate::shapes::scratchpad_testing::aux;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::{deriv_approx_eq, Matrix, Vector};

    #[test]
    fn calc_jacobian_handles_errors() {
        let mut pad = Scratchpad::new(2, GeoKind::Tri3).unwrap();
        assert_eq!(
            pad.calc_jacobian(&[0.0, 0.0]).err(),
            Some("all components of the coordinates matrix must be set first")
        );

        // bugged Jacobian matrix
        // (this would only happen if the user messes up the properties of Scratchpad directly)
        pad.set_xx(2, 1, 0.0); // setting the last component
                               // (cannot really check that all components have been set)
        pad.jacobian = Matrix::new(0, 0);
        assert_eq!(pad.calc_jacobian(&[0.0, 0.0]).err(), Some("matrices are incompatible"));
    }

    // Holds arguments for numerical differentiation of x with respect to ξ => Jacobian
    struct ArgsNumJac {
        pad: Scratchpad,  // scratchpad to send to calc_coords
        at_ksi: Vec<f64>, // at reference coord value
        ksi: Vec<f64>,    // temporary reference coord
        x: Vector,        // (space_ndim) coordinates at ξ
        i: usize,         // dimension index from 0 to space_ndim
        j: usize,         // dimension index from 0 to geo_ndim
    }

    // Computes xᵢ(ξ) with variable v := ξⱼ
    fn x_given_ksi(v: f64, args: &mut ArgsNumJac) -> f64 {
        args.ksi.copy_from_slice(&args.at_ksi);
        args.ksi[args.j] = v;
        args.pad.calc_coords(&mut args.x, &args.ksi).unwrap();
        args.x[args.i]
    }

    #[test]
    fn calc_jacobian_works() {
        // kind and tolerances
        let problem = vec![
            // Lin
            (GeoKind::Lin2, 1e-12),
            (GeoKind::Lin3, 1e-11),
            (GeoKind::Lin4, 1e-11),
            (GeoKind::Lin5, 1e-11),
            // Tri
            (GeoKind::Tri3, 1e-11),
            (GeoKind::Tri6, 1e-11),
            (GeoKind::Tri10, 1e-10),
            (GeoKind::Tri15, 1e-10),
            // Qua
            (GeoKind::Qua4, 1e-11),
            (GeoKind::Qua8, 1e-11),
            (GeoKind::Qua9, 1e-11),
            (GeoKind::Qua12, 1e-10),
            (GeoKind::Qua16, 1e-10),
            (GeoKind::Qua17, 1e-10),
            // Tet
            (GeoKind::Tet4, 1e-11),
            (GeoKind::Tet10, 1e-11),
            (GeoKind::Tet20, 1e-9),
            // Hex
            (GeoKind::Hex8, 1e-11),
            (GeoKind::Hex20, 1e-11),
            (GeoKind::Hex32, 1e-9),
        ];
        assert_eq!(problem.len(), GeoKind::VALUES.len());

        // loop over shapes
        for (kind, tol) in problem {
            println!("calc_jacobian: kind = {:?}", kind);
            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let mut pad = aux::gen_scratchpad_with_coords(space_ndim, kind);

            // set ξ within reference space
            let at_ksi = vec![0.25; geo_ndim];

            // compute Jacobian, its inverse, and determinant
            let det_jac = pad.calc_jacobian(&at_ksi).unwrap();
            assert!(det_jac > 0.0);

            // set arguments for numerical integration
            let args = &mut ArgsNumJac {
                pad: pad.clone(),
                at_ksi,
                ksi: vec![0.0; geo_ndim],
                x: Vector::new(space_ndim),
                i: 0,
                j: 0,
            };

            // check J(ξ) = dx(ξ)/dξ
            for i in 0..space_ndim {
                args.i = i;
                for j in 0..geo_ndim {
                    args.j = j;
                    // Jᵢⱼ := dxᵢ/dξⱼ
                    deriv_approx_eq(pad.jacobian.get(i, j), args.at_ksi[j], args, tol, x_given_ksi);
                }
            }
        }
    }

    #[test]
    fn calc_jacobian_special_cases_work() {
        // CABLE: line in 2d
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Lin2).unwrap();
        let l = 3.5;
        pad.set_xx(0, 0, 0.0); // node 0
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, l); // node 1
        pad.set_xx(1, 1, 0.0);
        let norm_jac_vec = pad.calc_jacobian(&[0.0]).unwrap();
        assert_eq!(norm_jac_vec, l / 2.0); // 2.0 = length of shape in the reference space

        // SHELL: triangle on a plane diagonal to the y-z plane
        let space_ndim = 3;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Tri3).unwrap();
        pad.set_xx(0, 0, 0.0); // node 0
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(0, 2, 0.0);
        pad.set_xx(1, 0, 1.0); // node 1
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(1, 2, 0.0);
        pad.set_xx(2, 0, 0.5); // node 2
        pad.set_xx(2, 1, 1.0);
        pad.set_xx(2, 2, 1.0);
        let norm_jac_vec = pad.calc_jacobian(&[0.0, 0.0]).unwrap();
        assert_eq!(norm_jac_vec, DET_JAC_NOT_AVAILABLE);
        // pad.draw_shape_simple("/tmp/gemlab/test_jacobian_tri3_in_3d.svg")
        //     .unwrap();
    }
}
