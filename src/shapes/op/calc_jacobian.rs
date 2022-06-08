use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{inverse, mat_mat_mul};

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
/// jacobian := Jline = Xᵀ · L
/// jacobian := Jsurf = Xᵀ · L
/// ```
///
/// where:
///
/// * `Jline`` -- Jacobian for line in multi-dimensions (geom_ndim < space_ndim)
/// * `Jsurf`` -- Jacobian for 3D surfaces (geo_ndim = 2 and space_ndim = 3)
///
/// If `geo_ndim = space_ndim`, this function also computes the inverse Jacobian:
///
/// ```text
/// inv_jacobian := J⁻¹
/// ```
///
/// # Output
///
/// * `pad.deriv` -- derivatives of the interpolation functions (nnode); `L` matrix
/// * `pad.jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
/// * `pad.inv_jacobian` -- If `geo_ndim = space_ndim`: inverse Jacobian matrix (space_ndim,space_ndim)
/// * Returns one of the following:
///     * If `geo_ndim = space_ndim`, returns the determinant of the Jacobian
///     * If `geo_ndim = 1` and `space_ndim > 1`, returns the norm of the Jacobian vector
///     * Otherwise, returns zero
///
/// # Input
///
/// * `ksi` -- reference coordinates ξ with len ≥ geo_ndim
pub fn calc_jacobian(pad: &mut Scratchpad, ksi: &[f64]) -> Result<f64, StrError> {
    // check
    if !pad.ok_xxt {
        return Err("all components of the coordinates matrix must be set first");
    }

    // matrix L: dNᵐ/dξ
    (pad.fn_deriv)(&mut pad.deriv, ksi);

    // matrix J: dx/dξ
    mat_mat_mul(&mut pad.jacobian, 1.0, &pad.xxt, &pad.deriv)?;

    let (space_ndim, geo_ndim) = pad.jacobian.dims();
    if geo_ndim == space_ndim {
        // inverse J
        inverse(&mut pad.inv_jacobian, &pad.jacobian)
    } else {
        // norm of Jacobian vector
        if geo_ndim == 1 {
            let mut norm_jac = 0.0;
            for i in 0..space_ndim {
                norm_jac += pad.jacobian[i][0] * pad.jacobian[i][0];
            }
            return Ok(f64::sqrt(norm_jac));
        }
        Ok(0.0)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::calc_jacobian;
    use crate::shapes::op::calc_coords;
    use crate::shapes::op::testing::aux;
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::StrError;
    use russell_chk::assert_deriv_approx_eq;
    use russell_lab::Vector;

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
        calc_coords(&mut args.x, &mut args.pad, &args.ksi).unwrap();
        args.x[args.i]
    }

    #[test]
    fn calc_jacobian_works() -> Result<(), StrError> {
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
            (GeoKind::Qua9, 1e-12),
            (GeoKind::Qua12, 1e-10),
            (GeoKind::Qua16, 1e-10),
            (GeoKind::Qua17, 1e-10),
            // Tet
            (GeoKind::Tet4, 1e-12),
            (GeoKind::Tet10, 1e-11),
            (GeoKind::Tet20, 1e-10),
            // Hex
            (GeoKind::Hex8, 1e-11),
            (GeoKind::Hex20, 1e-11),
            (GeoKind::Hex32, 1e-9),
        ];
        assert_eq!(problem.len(), GeoKind::VALUES.len());

        // loop over shapes
        for (kind, tol) in problem {
            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let mut pad = aux::gen_scratchpad_with_coords(space_ndim, kind);

            // set ξ within reference space
            let at_ksi = vec![0.25; geo_ndim];

            // compute Jacobian, its inverse, and determinant
            let det_jac = calc_jacobian(&mut pad, &at_ksi)?;
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
                    assert_deriv_approx_eq!(pad.jacobian[i][j], args.at_ksi[j], x_given_ksi, args, tol);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn calc_jacobian_special_cases_work() -> Result<(), StrError> {
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Lin2)?;
        let l = 3.5;
        pad.set_xx(0, 0, 0.0); // node 0
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, l); // node 1
        pad.set_xx(1, 1, 0.0);
        let norm_jac_vec = calc_jacobian(&mut pad, &[0.0])?;
        assert_eq!(norm_jac_vec, l / 2.0); // 2.0 = length of shape in the reference space

        let space_ndim = 3;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Tri3)?;
        pad.set_xx(0, 0, 0.0); // node 0
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(0, 2, 0.0);
        pad.set_xx(1, 0, 1.0); // node 1
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(1, 2, 0.0);
        pad.set_xx(2, 0, 0.5); // node 2
        pad.set_xx(2, 1, 1.0);
        pad.set_xx(2, 2, 1.0);
        let norm_jac_vec = calc_jacobian(&mut pad, &[0.0, 0.0])?;
        assert_eq!(norm_jac_vec, 0.0);
        Ok(())
    }
}
