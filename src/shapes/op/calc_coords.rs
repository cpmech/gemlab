use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{mat_vec_mul, Vector};

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
/// * `x` -- real coordinates (space_ndim)
/// * `pad.interp` -- (nnode) interpolation functions @ ξ
///
/// # Input
///
/// * `ksi` -- reference coordinates ξ with len ≥ geo_ndim
pub fn calc_coords(x: &mut Vector, pad: &mut Scratchpad, ksi: &[f64]) -> Result<(), StrError> {
    if !pad.ok_xxt {
        return Err("all components of the coordinates matrix must be set first");
    }
    let space_ndim = pad.xmin.len();
    if x.dim() != space_ndim {
        return Err("x.dim() must be equal to space_ndim");
    }
    (pad.fn_interp)(&mut pad.interp, ksi);
    mat_vec_mul(x, 1.0, &pad.xxt, &pad.interp)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::calc_coords;
    use crate::shapes::op::testing::aux;
    use crate::shapes::GeoKind;
    use crate::util::ONE_BY_3;
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;

    #[test]
    fn calc_coords_works() -> Result<(), StrError> {
        // kind, tol, tol_in for the case inside the shape
        let problem = vec![
            (GeoKind::Tri3, 1e-15, 0.35),  // linear maps are inaccurate for the circular wedge
            (GeoKind::Tri6, 1e-15, 0.013), // << quadratic mapping is inaccurate as well
            (GeoKind::Tri10, 1e-14, 1e-14),
            (GeoKind::Tri15, 1e-14, 1e-5), // << this triangle is inaccurate as well here
            (GeoKind::Qua4, 1e-15, 0.19),  // linear maps are inaccurate for the circular wedge
            (GeoKind::Qua8, 1e-15, 1e-15),
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
            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let mut pad = aux::gen_scratchpad_with_coords(space_ndim, kind);

            // loop over nodes of shape
            let nnode = kind.nnode();
            let mut x = Vector::new(space_ndim);
            let mut x_correct = Vector::new(space_ndim);
            let (ksi_min, ksi_del) = kind.ksi_min_ksi_del();
            for m in 0..nnode {
                // get ξᵐ corresponding to node m
                let ksi = kind.reference_coords(m);

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                calc_coords(&mut x, &mut pad, ksi)?;

                // compare xᵐ with generated coordinates
                aux::map_point_coords(&mut x_correct, ksi, ksi_min, ksi_del);
                assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol);
            }

            // test again inside the reference domain
            let ksi_in = if kind.is_tri_or_tet() {
                vec![ONE_BY_3; geo_ndim]
            } else {
                vec![0.0; geo_ndim]
            };
            calc_coords(&mut x, &mut pad, &ksi_in)?;
            aux::map_point_coords(&mut x_correct, &ksi_in, ksi_min, ksi_del);
            assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol_in);
        }
        Ok(())
    }
}
