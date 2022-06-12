use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{mat_mat_mul, Vector};

/// Calculates the normal vector
///
/// **Important:** This function only works with:
///
/// * `CABLE` case (geo_ndim = 1 and space_ndim = 2) -- e.g., line in 2D, or
/// * `SHELL` case (geo_ndim = 2 and space_ndim = 3) -- e.g., surface in 3D.
///
/// i.e., `geo_ndim < space_ndim`.
///
/// # Output
///
/// * `normal` -- (space_ndim) the normal vector; not necessarily unitary
/// * `pad.deriv` -- derivatives of the interpolation functions (nnode); `L` matrix
/// * `pad.jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
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
/// use gemlab::shapes::{op, GeoKind, Scratchpad};
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
///     let mut normal = Vector::new(2);
///     op::calc_normal_vector(&mut normal, &mut pad, &[0.0, 0.0])?;
///     assert_eq!(normal.as_data(), &[-H / 2.0, H / 2.0]);
///     Ok(())
/// }
/// ```
///
/// ## Boundary surface (geo_ndim = 2 and space_ndim = 3)
///
/// ```
/// use gemlab::shapes::{op, GeoKind, Scratchpad};
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
///     let mut normal = Vector::new(3);
///     op::calc_normal_vector(&mut normal, &mut pad, &[0.0, 0.0, 0.0])?;
///     const A: f64 = 1.0;
///     assert_eq!(normal.as_data(), &[0.0, A / 4.0, 0.0]);
///     Ok(())
/// }
/// ```
pub fn calc_normal_vector(n: &mut Vector, pad: &mut Scratchpad, ksi: &[f64]) -> Result<(), StrError> {
    // check
    let (space_ndim, geo_ndim) = pad.jacobian.dims();
    if geo_ndim >= space_ndim {
        return Err("geo_ndim must be smaller than space_ndim");
    }
    if n.dim() != space_ndim {
        return Err("n.dim() must be equal to space_ndim");
    }

    // matrix L: dNᵐ/dξ
    (pad.fn_deriv)(&mut pad.deriv, ksi);

    // matrix J: dx/dξ
    mat_mat_mul(&mut pad.jacobian, 1.0, &pad.xxt, &pad.deriv)?;

    // CABLE: line in 2D (geo_ndim = 1 and self.space_ndim = 2)
    //          →
    //         dx
    // g₁(ξ) = —— = Xᵀ · L = first_column(J)
    //         dξ
    //
    // →   →    →
    // n = e₃ × g₁ = {-g₁_0, +g₁_1}
    if space_ndim == 2 {
        n[0] = -pad.jacobian[1][0];
        n[1] = pad.jacobian[0][0];
        return Ok(());
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
    let jj = &pad.jacobian;
    n[0] = jj[1][0] * jj[2][1] - jj[2][0] * jj[1][1];
    n[1] = jj[2][0] * jj[0][1] - jj[0][0] * jj[2][1];
    n[2] = jj[0][0] * jj[1][1] - jj[1][0] * jj[0][1];
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::calc_normal_vector;
    use crate::shapes::op::testing::aux::{self, extract_edge, extract_face};
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::util::ONE_BY_3;
    use crate::StrError;
    use russell_chk::{assert_approx_eq, assert_vec_approx_eq};
    use russell_lab::{vector_norm, NormVec, Vector};

    #[test]
    fn calc_normal_vector_handles_errors() {
        let mut n = Vector::new(1);
        let mut pad = Scratchpad::new(2, GeoKind::Tri3).unwrap();
        assert_eq!(
            calc_normal_vector(&mut n, &mut pad, &[0.0, 0.0]).err(),
            Some("geo_ndim must be smaller than space_ndim")
        );
        let mut pad = Scratchpad::new(3, GeoKind::Tri3).unwrap();
        assert_eq!(
            calc_normal_vector(&mut n, &mut pad, &[0.0, 0.0]).err(),
            Some("n.dim() must be equal to space_ndim")
        );
    }

    #[test]
    fn calc_normal_vector_works_line() -> Result<(), StrError> {
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
        let correct_normal = vec![
            -correct_magnitude * f64::sin(aux::AMAX),
            correct_magnitude * f64::cos(aux::AMAX),
        ];

        // lover over shapes
        let ksi = &[0.25];
        let mut normal = Vector::new(2);
        for (kind, tol_mag, tol_vec) in problem {
            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let mut pad = aux::gen_scratchpad_with_coords(space_ndim, kind);

            // check
            calc_normal_vector(&mut normal, &mut pad, ksi)?;
            assert_approx_eq!(vector_norm(&normal, NormVec::Euc), correct_magnitude, tol_mag);
            assert_vec_approx_eq!(normal.as_data(), &correct_normal, tol_vec);
        }
        Ok(())
    }

    #[test]
    fn calc_normal_vector_works_surface_hex() -> Result<(), StrError> {
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
        let correct_normal_face2 = vec![
            correct_magnitude_face2_face3 * f64::sin(aux::AMIN),
            -correct_magnitude_face2_face3 * f64::cos(aux::AMIN),
            0.0,
        ];
        let correct_normal_face3 = vec![
            -correct_magnitude_face2_face3 * f64::sin(aux::AMAX),
            correct_magnitude_face2_face3 * f64::cos(aux::AMAX),
            0.0,
        ];

        // lover over shapes
        let ksi = &[ONE_BY_3, ONE_BY_3];
        let mut normal = Vector::new(3);
        for (kind, tol_mag, tol_vec) in problem {
            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let pad = aux::gen_scratchpad_with_coords(space_ndim, kind);

            // face # 0
            let mut pad_face = extract_face(0, &pad);
            calc_normal_vector(&mut normal, &mut pad_face, ksi)?;
            assert!(normal[0] < 0.0);
            assert!(normal[1] < 0.0);
            assert_approx_eq!(normal[2], 0.0, tol_vec);

            // face # 1
            let mut pad_face = extract_face(1, &pad);
            calc_normal_vector(&mut normal, &mut pad_face, ksi)?;
            assert!(normal[0] > 0.0);
            assert!(normal[1] > 0.0);
            assert_approx_eq!(normal[2], 0.0, tol_vec);

            // face # 2
            let mut pad_face = extract_face(2, &pad);
            calc_normal_vector(&mut normal, &mut pad_face, ksi)?;
            assert_approx_eq!(
                vector_norm(&normal, NormVec::Euc),
                correct_magnitude_face2_face3,
                tol_mag
            );
            assert_vec_approx_eq!(normal.as_data(), &correct_normal_face2, tol_vec);

            // face # 3
            let mut pad_face = extract_face(3, &pad);
            calc_normal_vector(&mut normal, &mut pad_face, ksi)?;
            assert_approx_eq!(
                vector_norm(&normal, NormVec::Euc),
                correct_magnitude_face2_face3,
                tol_mag
            );
            assert_vec_approx_eq!(normal.as_data(), &correct_normal_face3, tol_vec);

            // face # 4
            let mut pad_face = extract_face(4, &pad);
            calc_normal_vector(&mut normal, &mut pad_face, ksi)?;
            assert_approx_eq!(normal[0], 0.0, tol_vec);
            assert_approx_eq!(normal[1], 0.0, tol_vec);
            assert!(normal[2] < 0.0);

            // face # 5
            let mut pad_face = extract_face(5, &pad);
            calc_normal_vector(&mut normal, &mut pad_face, ksi)?;
            assert_approx_eq!(normal[0], 0.0, tol_vec);
            assert_approx_eq!(normal[1], 0.0, tol_vec);
            assert!(normal[2] > 0.0);
        }
        Ok(())
    }

    #[test]
    fn normals_are_outward_2d() -> Result<(), StrError> {
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
        let tri_correct = vec![
            &[0.0, -1.0], // bottom
            &[1.0, 1.0],  // diagonal
            &[-1.0, 0.0], // left
        ];
        let qua_correct = vec![
            &[0.0, -1.0], // bottom
            &[1.0, 0.0],  // right
            &[0.0, 1.0],  // top
            &[-1.0, 0.0], // left
        ];

        // auxiliary
        let mut normal = Vector::new(2);
        let ksi_values = &[[0.0, 0.0], [ONE_BY_3, ONE_BY_3]];

        // loop over shapes
        for kind in kinds {
            // scratchpad with coordinates
            let pad = aux::gen_scratchpad_with_coords_aligned(kind);

            // loop over edges
            for e in 0..pad.kind.nedge() {
                for ksi in ksi_values {
                    let mut pad_edge = extract_edge(e, &pad);
                    calc_normal_vector(&mut normal, &mut pad_edge, ksi)?;
                    if pad.kind.is_tri_or_tet() {
                        // check triangle
                        assert_vec_approx_eq!(normal.as_data(), tri_correct[e], 1e-15);
                    } else {
                        // check quadrilateral
                        assert_vec_approx_eq!(normal.as_data(), qua_correct[e], 1e-15);
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn normals_are_outward_3d() -> Result<(), StrError> {
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
        let tet_correct = vec![
            &[-4.0, 0.0, 0.0], // negative-x face
            &[0.0, -4.0, 0.0], // negative-y face
            &[0.0, 0.0, -4.0], // negative-z face
            &[4.0, 4.0, 4.0],  // face orthogonal to the diagonal
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
        let mut normal = Vector::new(3);
        let ksi_values = &[[0.0, 0.0, 0.0], [ONE_BY_3, ONE_BY_3, ONE_BY_3]];

        // loop over shapes
        for (kind, tol) in problem {
            // scratchpad with coordinates
            let pad = aux::gen_scratchpad_with_coords_aligned(kind);

            // loop over faces
            for f in 0..pad.kind.nface() {
                for ksi in ksi_values {
                    let mut pad_face = extract_face(f, &pad);
                    calc_normal_vector(&mut normal, &mut pad_face, ksi)?;
                    if pad.kind.is_tri_or_tet() {
                        // check tetrahedron
                        assert_vec_approx_eq!(normal.as_data(), tet_correct[f], tol);
                    } else {
                        // check hexahedron
                        assert_vec_approx_eq!(normal.as_data(), hex_correct[f], tol);
                    }
                }
            }
        }
        Ok(())
    }
}
