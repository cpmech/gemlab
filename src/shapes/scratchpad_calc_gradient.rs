use super::Scratchpad;
use crate::StrError;
use russell_lab::mat_mat_mul;

impl Scratchpad {
    /// Calculates the gradient of the interpolation functions
    ///
    /// **Note:** This function works with `geo_ndim == space_ndim` only.
    ///
    /// The gradient is given by:
    ///
    /// ```text
    ///             →
    /// →  →    dNᵐ(ξ)
    /// Bᵐ(ξ) = ——————
    ///            →
    ///           dx
    /// ```
    ///
    /// which can be organized in an (nnode,space_ndim) matrix `B` as follows
    ///
    /// ```text
    /// B = L · J⁻¹
    /// ```
    ///
    /// # Output
    ///
    /// * `deriv` -- interpolation functions (nnode)
    /// * `jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    /// * `inv_jacobian` -- inverse Jacobian matrix (space_ndim,space_ndim)
    /// * `gradient` -- gradient matrix (nnode,space_ndim)
    /// * Returns the determinant of the Jacobian matrix
    ///
    /// # Input
    ///
    /// * `ksi` -- reference coordinates ξ with len ≥ geo_ndim
    ///
    /// # Example
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Scratchpad};
    /// use gemlab::StrError;
    /// use russell_chk::vec_approx_eq;
    /// use russell_lab::Matrix;
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
    ///     pad.calc_gradient(&[0.0, 0.0])?;
    ///
    ///     let correct_gg = Matrix::from(&[
    ///         [-1.0 / (4.0 * a), -1.0 / (2.0 * a)],
    ///         [1.0 / (4.0 * a), -1.0 / (2.0 * a)],
    ///         [1.0 / (4.0 * a), 1.0 / (2.0 * a)],
    ///         [-1.0 / (4.0 * a), 1.0 / (2.0 * a)],
    ///     ]);
    ///     vec_approx_eq(pad.gradient.as_data(), correct_gg.as_data(), 1e-15);
    ///     Ok(())
    /// }
    /// ```
    pub fn calc_gradient(&mut self, ksi: &[f64]) -> Result<f64, StrError> {
        // check
        let (space_ndim, geo_ndim) = self.jacobian.dims();
        if geo_ndim != space_ndim {
            return Err("calc_gradient requires that geo_ndim = space_ndim");
        }

        // Jacobian matrix J: dx/dξ
        let det_jac = self.calc_jacobian(ksi)?;

        // gradient: B = L · J⁻¹
        mat_mat_mul(&mut self.gradient, 1.0, &self.deriv, &self.inv_jacobian).unwrap(); // cannot fail because the dims are checked
        Ok(det_jac)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::shapes::scratchpad_testing::aux;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_chk::deriv_approx_eq;
    use russell_lab::math::ONE_BY_3;
    use russell_lab::{vec_copy, Vector};

    #[test]
    fn calc_gradient_handles_errors() {
        let mut pad = Scratchpad::new(2, GeoKind::Lin2).unwrap();
        assert_eq!(
            pad.calc_gradient(&[0.0, 0.0]).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );

        let mut pad = Scratchpad::new(2, GeoKind::Tri3).unwrap();
        assert_eq!(
            pad.calc_gradient(&[0.0, 0.0]).err(),
            Some("all components of the coordinates matrix must be set first")
        );
    }

    // Holds arguments for numerical differentiation of N with respect to x => B (gradient) matrix
    struct ArgsNumGrad {
        pad: Scratchpad, // scratchpad to send to calc_coords
        at_x: Vector,    // at x coord value
        x: Vector,       // temporary x coord
        ksi: Vec<f64>,   // temporary reference coord
        m: usize,        // node index from 0 to nnode
        j: usize,        // dimension index from 0 to space_ndim
    }

    // Computes Nᵐ(ξ(x)) with variable v := xⱼ
    fn nn_given_x(v: f64, args: &mut ArgsNumGrad) -> f64 {
        vec_copy(&mut args.x, &args.at_x).unwrap();
        args.x[args.j] = v;
        args.pad.approximate_ksi(&mut args.ksi, &args.x, 10, 1e-14).unwrap();
        (args.pad.fn_interp)(&mut args.pad.interp, &args.ksi);
        args.pad.interp[args.m]
    }

    #[test]
    fn calc_gradient_works() {
        // kind (except Lin) and tolerances
        let problem = vec![
            // Tri
            (GeoKind::Tri3, 1e-12),
            (GeoKind::Tri6, 1e-10),
            (GeoKind::Tri10, 1e-9),
            (GeoKind::Tri15, 1e-9),
            // Qua
            (GeoKind::Qua4, 1e-11),
            (GeoKind::Qua8, 1e-10),
            (GeoKind::Qua9, 1e-10),
            (GeoKind::Qua12, 1e-9),
            (GeoKind::Qua16, 1e-9),
            (GeoKind::Qua17, 1e-9),
            // Tet
            (GeoKind::Tet4, 1e-12),
            (GeoKind::Tet10, 1e-9),
            (GeoKind::Tet20, 1e-9),
            // Hex
            (GeoKind::Hex8, 1e-11),
            (GeoKind::Hex20, 1e-10),
            (GeoKind::Hex32, 1e-9),
        ];

        // loop over shapes
        for (kind, tol) in problem {
            // println!("kind = {:?}", kind);

            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let mut pad = aux::gen_scratchpad_with_coords(space_ndim, kind);

            // set ξ within reference space
            let at_ksi = vec![ONE_BY_3; geo_ndim];

            // compute x corresponding to ξ using the isoparametric formula
            let mut at_x = Vector::new(space_ndim);
            pad.calc_coords(&mut at_x, &at_ksi).unwrap();

            // compute gradient
            let det_jac = pad.calc_gradient(&at_ksi).unwrap();
            assert!(det_jac > 0.0);

            // set arguments for numerical integration
            let args = &mut ArgsNumGrad {
                pad: pad.clone(),
                at_x,
                x: Vector::new(space_ndim),
                ksi: vec![0.0; geo_ndim],
                m: 0,
                j: 0,
            };

            // check Bᵐ(ξ(x)) = dNᵐ(ξ(x))/dx
            for m in 0..kind.nnode() {
                args.m = m;
                for j in 0..geo_ndim {
                    args.j = j;
                    // Bᵐⱼ := dNᵐ/dxⱼ
                    deriv_approx_eq(pad.gradient.get(m, j), args.at_x[j], args, tol, nn_given_x);
                }
            }
        }
    }
}
