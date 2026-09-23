use super::Scratchpad;
use crate::StrError;
use russell_tensor::Tensor1;

impl Scratchpad {
    /// Calculates the real coordinates x from reference coordinates ξ
    ///
    /// This function uses the isoparametric formula to calculate x given ξ:
    ///
    /// ```text
    /// → →         →  →
    /// x(ξ) = Σ Nᵐ(ξ) xᵐ
    ///        m
    ///
    /// x := Xᵀ ⋅ N
    /// ```
    ///
    /// # Output
    ///
    /// * `x` -- real coordinates (space_ndim); components beyond `space_ndim` are set to zero
    /// * `pad.interp` -- (nnode) interpolation functions @ ξ
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
    /// use russell_tensor::{Tensor1, t1_approx_eq};
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
    ///     let (x0, y0) = (3.0, 4.0);
    ///     let (w, h) = (10.0, 5.0);
    ///     let space_ndim = 2;
    ///     let mut pad = Scratchpad::new(space_ndim, GeoKind::Qua4)?;
    ///     pad.set_xx(0, 0, x0);
    ///     pad.set_xx(0, 1, y0);
    ///     pad.set_xx(1, 0, x0 + w);
    ///     pad.set_xx(1, 1, y0);
    ///     pad.set_xx(2, 0, x0 + w);
    ///     pad.set_xx(2, 1, y0 + h);
    ///     pad.set_xx(3, 0, x0);
    ///     pad.set_xx(3, 1, y0 + h);
    ///
    ///     let mut x = Tensor1::new();
    ///     pad.calc_coords(&mut x, &[0.0, 0.0])?;
    ///     let expected = Tensor1::from(&[x0 + w / 2.0, y0 + h / 2.0, 0.0]);
    ///     t1_approx_eq(&x, &expected, 1e-15);
    ///     Ok(())
    /// }
    /// ```
    pub fn calc_coords(&mut self, x: &mut Tensor1, ksi: &[f64]) -> Result<(), StrError> {
        if !self.ok_xxt {
            return Err("all components of the coordinates matrix must be set first");
        }
        let (space_ndim, _, nnode) = self.dims();
        self.calc_interp(ksi);
        for i in 0..space_ndim {
            x.set(i, 0.0);
            for j in 0..nnode {
                x.add(i, self.xxt.get(i, j) * self.interp.get(j));
            }
        }
        // reset unused trailing components (e.g., the z component in 2D)
        for i in space_ndim..3 {
            x.set(i, 0.0);
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::shapes::scratchpad_testing::aux;
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_lab::math::ONE_BY_3;
    use russell_tensor::{t1_approx_eq, Tensor1};

    #[test]
    fn calc_coords_handles_errors() {
        let mut x = Tensor1::new();
        let mut pad = Scratchpad::new(2, GeoKind::Tri3).unwrap();
        assert_eq!(
            pad.calc_coords(&mut x, &[0.0, 0.0]).err(),
            Some("all components of the coordinates matrix must be set first")
        );
        pad.set_xx(2, 1, 0.0); // setting the last component
                               // (cannot really check that all components have been set)
        pad.calc_coords(&mut x, &[0.0, 0.0]).unwrap();
    }

    #[test]
    fn calc_coords_works() {
        // kind, tol, tol_in for the case inside the shape
        let problem = vec![
            (GeoKind::Tri3, 1e-15, 0.35),  // linear maps are inaccurate for the circular wedge
            (GeoKind::Tri6, 1e-15, 0.013), // << quadratic mapping is inaccurate as well
            (GeoKind::Tri10, 1e-14, 1e-14),
            (GeoKind::Tri15, 1e-14, 1e-5), // << this triangle is inaccurate as well here
            (GeoKind::Qua4, 1e-15, 0.19),  // linear maps are inaccurate for the circular wedge
            (GeoKind::Qua8, 1e-15, 1e-14),
            (GeoKind::Qua17, 1e-15, 1e-15),
            (GeoKind::Tet4, 1e-15, 0.35),   // linear tetrahedron is also inaccurate here
            (GeoKind::Tet10, 1e-15, 0.013), // quadratic tetrahedron is also inaccurate here
            (GeoKind::Tet20, 1e-14, 1e-14), // cubic tetrahedron
            (GeoKind::Hex8, 1e-14, 0.19),   // bi-linear maps are inaccurate for the circular wedge
            (GeoKind::Hex20, 1e-15, 1e-15),
            (GeoKind::Hex32, 1e-15, 0.00012), // TODO: check why this tolerance is high
        ];

        // loop over shapes
        for (kind, tol, tol_in) in problem {
            println!("kind = {:?}", kind);

            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let mut pad = aux::gen_scratchpad_with_coords(space_ndim, kind);

            // loop over nodes of shape
            let nnode = kind.nnode();
            let mut x = Tensor1::new();
            let mut x_correct = Tensor1::new();
            let (ksi_min, ksi_del) = kind.ksi_min_ksi_del();
            for m in 0..nnode {
                // get ξᵐ corresponding to node m
                let ksi = kind.reference_coords(m);

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                pad.calc_coords(&mut x, ksi).unwrap();

                // compare xᵐ with generated coordinates
                aux::map_point_coords(&mut x_correct, ksi, ksi_min, ksi_del);
                t1_approx_eq(&x, &x_correct, tol);
            }

            // test again inside the reference domain
            let ksi_in = if kind.is_tri_or_tet() {
                vec![ONE_BY_3; geo_ndim]
            } else {
                vec![0.0; geo_ndim]
            };
            pad.calc_coords(&mut x, &ksi_in).unwrap();
            aux::map_point_coords(&mut x_correct, &ksi_in, ksi_min, ksi_del);
            t1_approx_eq(&x, &x_correct, tol_in);
        }
    }
}
