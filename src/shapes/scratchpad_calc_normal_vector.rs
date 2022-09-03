use super::Scratchpad;
use crate::StrError;
use russell_lab::{mat_mat_mul, Vector};

impl Scratchpad {
    /// Calculates the (unit) normal vector
    ///
    /// **Important:** This function only works with:
    ///
    /// * `CABLE` in 2D case (geo_ndim = 1 and space_ndim = 2) -- e.g., line in 2D, or
    /// * `SHELL` in 3D case (geo_ndim = 2 and space_ndim = 3) -- e.g., surface in 3D.
    ///
    /// The case `CABLE` in 3D is **not** available because we don't know the direction of the normal vector.
    ///
    /// # Output
    ///
    /// * `un` -- (space_ndim) the **unit** normal vector
    /// * `deriv` -- derivatives of the interpolation functions (nnode); `L` matrix
    /// * `jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
    /// * Returns the magnitude of the normal vector
    ///
    /// # Input
    ///
    /// * `ksi` -- reference coordinates ξ with len ≥ geo_ndim
    ///
    /// # Examples
    ///
    /// ## Line in multi-dimensions (geo_ndim = 1 and space_ndim > 1)
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Scratchpad};
    /// use gemlab::util::SQRT_2;
    /// use gemlab::StrError;
    /// use russell_lab::Vector;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     //  →  __       1   -----
    ///     //  n |.      ,'     / \
    ///     //      '.  ,'        |
    ///     //        2'          H
    ///     //      ,'            |
    ///     //    ,'   45°        |
    ///     //   0  ____________ \ /
    ///
    ///     const H: f64 = 5.0;
    ///     let space_ndim = 2;
    ///     let mut pad = Scratchpad::new(space_ndim, GeoKind::Lin3)?;
    ///     pad.set_xx(0, 0, 0.0);
    ///     pad.set_xx(0, 1, 0.0);
    ///     pad.set_xx(1, 0, H);
    ///     pad.set_xx(1, 1, H);
    ///     pad.set_xx(2, 0, H / 2.0);
    ///     pad.set_xx(2, 1, H / 2.0);
    ///
    ///     // ||n|| = L/2 (2 is the length in the natural space)
    ///     // nx = -(L/2)sin(45) = -(H√2/2) √2/2 = -H/2
    ///     // ny = +(L/2)cos(45) = +(H√2/2) √2/2 = +H/2
    ///     // unx = -(H/2)/(H√2/2) = -1/√2
    ///     // uny = +(H/2)/(H√2/2) = +1/√2
    ///     let mut un = Vector::new(2);
    ///     let mag_n = pad.calc_normal_vector(&mut un, &[0.0, 0.0])?;
    ///     assert_eq!(un.as_data(), &[-1.0 / SQRT_2, 1.0 / SQRT_2]);
    ///     assert_eq!(mag_n, H * SQRT_2 / 2.0);
    ///     Ok(())
    /// }
    /// ```
    ///
    /// ## Boundary surface (geo_ndim = 2 and space_ndim = 3)
    ///
    /// ```
    /// use gemlab::shapes::{GeoKind, Scratchpad};
    /// use gemlab::StrError;
    /// use russell_lab::Vector;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     //           .   .  .   . ,.2|
    ///     //         ' .           ,,'||
    ///     //       '   .         ,,'  ||
    ///     //     '     .       .,'    ||  →
    ///     //  .  .   . .   .  3'      ||  n
    ///     //           z     ||   ==========)
    ///     //  .        |     ||       ||
    ///     //          ,*---y || .  . ,1
    ///     //  .      x       ||    ,,'
    ///     //      ,'         ||  ,,'
    ///     //  . ,'           ||,,'
    ///     //  . . .   .   .  |0'
    ///
    ///     let space_ndim = 3;
    ///     let mut pad = Scratchpad::new(space_ndim, GeoKind::Qua4)?;
    ///     pad.set_xx(0, 0, 1.0); // node 0
    ///     pad.set_xx(0, 1, 1.0);
    ///     pad.set_xx(0, 2, 0.0);
    ///     pad.set_xx(1, 0, 0.0); // node 1
    ///     pad.set_xx(1, 1, 1.0);
    ///     pad.set_xx(1, 2, 0.0);
    ///     pad.set_xx(2, 0, 0.0); // node 2
    ///     pad.set_xx(2, 1, 1.0);
    ///     pad.set_xx(2, 2, 1.0);
    ///     pad.set_xx(3, 0, 1.0); // node 3
    ///     pad.set_xx(3, 1, 1.0);
    ///     pad.set_xx(3, 2, 1.0);
    ///
    ///     let mut un = Vector::new(3);
    ///     let mag_n = pad.calc_normal_vector(&mut un, &[0.0, 0.0, 0.0])?;
    ///     assert_eq!(un.as_data(), &[0.0, 1.0, 0.0]);
    ///     assert_eq!(mag_n, 1.0 / 4.0);
    ///     Ok(())
    /// }
    /// ```
    pub fn calc_normal_vector(&mut self, un: &mut Vector, ksi: &[f64]) -> Result<f64, StrError> {
        // check
        let (space_ndim, geo_ndim) = self.jacobian.dims();
        if space_ndim == 2 && geo_ndim != 1 {
            return Err("calc_normal_vector requires geo_ndim = 1 in 2D (CABLE in 2D)");
        }
        if space_ndim == 3 && geo_ndim != 2 {
            return Err("calc_normal_vector requires geo_ndim = 2 in 3D (SHELL in 3D)");
        }
        if un.dim() != space_ndim {
            return Err("un.dim() must be equal to space_ndim");
        }

        // matrix L: dNᵐ/dξ
        (self.fn_deriv)(&mut self.deriv, ksi);

        // matrix J: dx/dξ
        mat_mat_mul(&mut self.jacobian, 1.0, &self.xxt, &self.deriv)?;

        // CABLE: line in 2D (geo_ndim = 1 and space_ndim = 2)
        //          →
        //         dx
        // g₁(ξ) = —— = Xᵀ · L = first_column(J)
        //         dξ
        //
        // →   →    →
        // n = e₃ × g₁ = {-g₁_0, +g₁_1}
        if space_ndim == 2 {
            un[0] = -self.jacobian[1][0];
            un[1] = self.jacobian[0][0];
            let mag_n = f64::sqrt(un[0] * un[0] + un[1] * un[1]);
            if mag_n > 0.0 {
                un[0] /= mag_n;
                un[1] /= mag_n;
            }
            return Ok(mag_n);
        }

        // SHELL: surface in 3D (geo_ndim = 2 and space_ndim = 3)
        //          →
        // →  →    dx
        // g₁(ξ) = ——— = first_column(J)
        //         dξ₁
        //
        //          →
        // →  →    dx
        // g₂(ξ) = ——— = second_column(J)
        //         dξ₂
        //
        // →   →    →
        // n = g₁ × g₂
        let jj = &self.jacobian;
        un[0] = jj[1][0] * jj[2][1] - jj[2][0] * jj[1][1];
        un[1] = jj[2][0] * jj[0][1] - jj[0][0] * jj[2][1];
        un[2] = jj[0][0] * jj[1][1] - jj[1][0] * jj[0][1];
        let mag_n = f64::sqrt(un[0] * un[0] + un[1] * un[1] + un[2] * un[2]);
        if mag_n > 0.0 {
            un[0] /= mag_n;
            un[1] /= mag_n;
            un[2] /= mag_n;
        }
        Ok(mag_n)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::shapes::scratchpad_testing::aux;
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::util::{ONE_BY_3, SQRT_2, SQRT_3};
    use russell_chk::{approx_eq, vec_approx_eq};
    use russell_lab::{vector_norm, NormVec, Vector};

    #[test]
    fn calc_normal_vector_handles_errors() {
        let mut un = Vector::new(1);
        let mut pad = Scratchpad::new(2, GeoKind::Tri3).unwrap();
        assert_eq!(
            pad.calc_normal_vector(&mut un, &[0.0, 0.0]).err(),
            Some("calc_normal_vector requires geo_ndim = 1 in 2D (CABLE in 2D)")
        );
        let mut pad = Scratchpad::new(3, GeoKind::Lin2).unwrap();
        assert_eq!(
            pad.calc_normal_vector(&mut un, &[0.0, 0.0]).err(),
            Some("calc_normal_vector requires geo_ndim = 2 in 3D (SHELL in 3D)")
        );
        let mut pad = Scratchpad::new(3, GeoKind::Tri3).unwrap();
        assert_eq!(
            pad.calc_normal_vector(&mut un, &[0.0, 0.0]).err(),
            Some("un.dim() must be equal to space_ndim")
        );
    }

    #[test]
    fn calc_normal_vector_works_line() {
        // kind, tol_mag, tol_vec
        let problem = vec![
            (GeoKind::Lin2, 1e-15, 1e-15),
            (GeoKind::Lin3, 1e-15, 1e-15),
            (GeoKind::Lin4, 1e-14, 1e-14),
            (GeoKind::Lin5, 1e-15, 1e-14),
        ];

        // correct values
        const KSI_DEL: f64 = 2.0;
        let correct_magnitude = (aux::RMAX - aux::RMIN) / KSI_DEL;
        let correct_normal = vec![-f64::sin(aux::AMAX), f64::cos(aux::AMAX)];

        // lover over shapes
        let ksi = &[0.25];
        let mut un = Vector::new(2);
        for (kind, tol_mag, tol_vec) in problem {
            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let mut pad = aux::gen_scratchpad_with_coords(space_ndim, kind);

            // check
            let mag_n = pad.calc_normal_vector(&mut un, ksi).unwrap();
            approx_eq(mag_n, correct_magnitude, tol_mag);
            approx_eq(vector_norm(&un, NormVec::Euc), 1.0, tol_mag);
            vec_approx_eq(un.as_data(), &correct_normal, tol_vec);
        }
    }

    #[test]
    fn calc_normal_vector_works_surface_hex() {
        // kind, tol_mag, tol_vec
        let problem = vec![
            (GeoKind::Hex8, 1e-15, 1e-15),
            (GeoKind::Hex20, 1e-15, 1e-14),
            (GeoKind::Hex32, 1e-14, 1e-14),
        ];

        // correct values: face # 2 and 3
        const REF_AREA: f64 = 4.0;
        let area_face2_face3 = (aux::RMAX - aux::RMIN) * (aux::ZMAX - aux::ZMIN);
        let correct_magnitude_face2_face3 = area_face2_face3 / REF_AREA;
        let correct_normal_face2 = vec![f64::sin(aux::AMIN), -f64::cos(aux::AMIN), 0.0];
        let correct_normal_face3 = vec![-f64::sin(aux::AMAX), f64::cos(aux::AMAX), 0.0];

        // lover over shapes
        let ksi = &[ONE_BY_3, ONE_BY_3];
        let mut un = Vector::new(3);
        for (kind, tol_mag, tol_vec) in problem {
            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let pad = aux::gen_scratchpad_with_coords(space_ndim, kind);

            // face # 0
            let mut pad_face = aux::extract_face(0, &pad);
            pad_face.calc_normal_vector(&mut un, ksi).unwrap();
            assert!(un[0] < 0.0);
            assert!(un[1] < 0.0);
            approx_eq(un[2], 0.0, tol_vec);

            // face # 1
            let mut pad_face = aux::extract_face(1, &pad);
            pad_face.calc_normal_vector(&mut un, ksi).unwrap();
            assert!(un[0] > 0.0);
            assert!(un[1] > 0.0);
            approx_eq(un[2], 0.0, tol_vec);

            // face # 2
            let mut pad_face = aux::extract_face(2, &pad);
            let mag_n = pad_face.calc_normal_vector(&mut un, ksi).unwrap();
            approx_eq(mag_n, correct_magnitude_face2_face3, tol_mag);
            approx_eq(vector_norm(&un, NormVec::Euc), 1.0, tol_mag);
            vec_approx_eq(un.as_data(), &correct_normal_face2, tol_vec);

            // face # 3
            let mut pad_face = aux::extract_face(3, &pad);
            let mag_n = pad_face.calc_normal_vector(&mut un, ksi).unwrap();
            approx_eq(mag_n, correct_magnitude_face2_face3, tol_mag);
            approx_eq(vector_norm(&un, NormVec::Euc), 1.0, tol_mag);
            vec_approx_eq(un.as_data(), &correct_normal_face3, tol_vec);

            // face # 4
            let mut pad_face = aux::extract_face(4, &pad);
            pad_face.calc_normal_vector(&mut un, ksi).unwrap();
            approx_eq(un[0], 0.0, tol_vec);
            approx_eq(un[1], 0.0, tol_vec);
            assert!(un[2] < 0.0);

            // face # 5
            let mut pad_face = aux::extract_face(5, &pad);
            pad_face.calc_normal_vector(&mut un, ksi).unwrap();
            approx_eq(un[0], 0.0, tol_vec);
            approx_eq(un[1], 0.0, tol_vec);
            assert!(un[2] > 0.0);
        }
    }

    #[test]
    fn normals_are_outward_2d() {
        // select Tri and Qua
        let kinds = vec![
            // Tri
            GeoKind::Tri3,
            GeoKind::Tri6,
            GeoKind::Tri10,
            GeoKind::Tri15,
            // Qua
            GeoKind::Qua4,
            GeoKind::Qua8,
            GeoKind::Qua9,
            GeoKind::Qua12,
            GeoKind::Qua16,
            GeoKind::Qua17,
        ];

        // solution
        //
        //   →     Δℓ_edge   Δℓ_edge
        // ||n|| = ——————— = ———————
        //         Δξ_lin       2
        //                                                 →     2 √2
        // For the diagonal of Tri: Δℓ_edge = 2 √2, thus ||n|| = ————— = √2
        //                                                         2
        const OS2: f64 = 1.0 / SQRT_2;
        let tri_correct = vec![
            &[0.0, -1.0], // bottom
            &[OS2, OS2],  // diagonal
            &[-1.0, 0.0], // left
        ];
        let qua_correct = vec![
            &[0.0, -1.0], // bottom
            &[1.0, 0.0],  // right
            &[0.0, 1.0],  // top
            &[-1.0, 0.0], // left
        ];

        // auxiliary
        let mut un = Vector::new(2);
        let ksi_values = &[[0.0, 0.0], [ONE_BY_3, ONE_BY_3]];

        // loop over shapes
        for kind in kinds {
            // scratchpad with coordinates
            let pad = aux::gen_scratchpad_with_coords_aligned(kind);

            // loop over edges
            for e in 0..pad.kind.nedge() {
                for ksi in ksi_values {
                    let mut pad_edge = aux::extract_edge(e, &pad);
                    let mag_n = pad_edge.calc_normal_vector(&mut un, ksi).unwrap();
                    if pad.kind.is_tri_or_tet() {
                        // check triangle
                        if e == 1 {
                            approx_eq(mag_n, SQRT_2, 1e-15);
                        } else {
                            approx_eq(mag_n, 1.0, 1e-15);
                        }
                        vec_approx_eq(un.as_data(), tri_correct[e], 1e-15);
                    } else {
                        // check quadrilateral
                        approx_eq(mag_n, 1.0, 1e-15);
                        vec_approx_eq(un.as_data(), qua_correct[e], 1e-15);
                    }
                }
            }
        }
    }

    #[test]
    fn normals_are_outward_3d() {
        // select Tet and Hex
        let problem = vec![
            // Tet
            (GeoKind::Tet4, 1e-15),
            (GeoKind::Tet10, 1e-15),
            (GeoKind::Tet20, 1e-14),
            // Hex
            (GeoKind::Hex8, 1e-15),
            (GeoKind::Hex20, 1e-15),
            (GeoKind::Hex32, 1e-15),
        ];

        // solution
        // Hex with sides equal to h=2:
        //
        //   →         ΔA_face       h²
        // ||n|| = ——————————————— = —— = 1
        //         Δξ₁_qua Δξ₂_qua    4
        //
        // Tet with sides equal to h=2:
        //
        //   Faces orthogonal to the x,y,z axes:
        //
        //      →     ΔA_face   h²/2
        //    ||n|| = ——————— = ———— = 4
        //             ΔA_tri    1/2
        //
        //   Face orthogonal to the the diagonal:
        //
        //      →     ΔA_face   √3 h²/2
        //    ||n|| = ——————— = ——————— = 4 √3
        //             ΔA_tri     1/2
        const OS3: f64 = 1.0 / SQRT_3;
        let tet_correct = vec![
            &[-1.0, 0.0, 0.0], // negative-x face
            &[0.0, -1.0, 0.0], // negative-y face
            &[0.0, 0.0, -1.0], // negative-z face
            &[OS3, OS3, OS3],  // face orthogonal to the diagonal
        ];
        let hex_correct = vec![
            &[-1.0, 0.0, 0.0], // behind
            &[1.0, 0.0, 0.0],  // front
            &[0.0, -1.0, 0.0], // left
            &[0.0, 1.0, 0.0],  // right
            &[0.0, 0.0, -1.0], // bottom
            &[0.0, 0.0, 1.0],  // top
        ];

        // auxiliary
        let mut un = Vector::new(3);
        let ksi_values = &[[0.0, 0.0, 0.0], [ONE_BY_3, ONE_BY_3, ONE_BY_3]];

        // loop over shapes
        for (kind, tol) in problem {
            // scratchpad with coordinates
            let pad = aux::gen_scratchpad_with_coords_aligned(kind);

            // loop over faces
            for f in 0..pad.kind.nface() {
                for ksi in ksi_values {
                    let mut pad_face = aux::extract_face(f, &pad);
                    let mag_n = pad_face.calc_normal_vector(&mut un, ksi).unwrap();
                    if pad.kind.is_tri_or_tet() {
                        // check tetrahedron
                        if f == 3 {
                            approx_eq(mag_n, 4.0 * SQRT_3, tol);
                        } else {
                            approx_eq(mag_n, 4.0, tol);
                        }
                        vec_approx_eq(un.as_data(), tet_correct[f], tol);
                    } else {
                        // check hexahedron
                        approx_eq(mag_n, 1.0, tol);
                        vec_approx_eq(un.as_data(), hex_correct[f], tol);
                    }
                }
            }
        }
    }
}
