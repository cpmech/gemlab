use crate::shapes::{FnDeriv, Scratchpad};
use crate::StrError;
use russell_lab::{inverse, mat_mat_mul, Matrix};

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
///
/// # Panics
///
/// This function does not check for the vector/matrix dimensions. Thus a panic may occur
/// if they are incompatible, including if they are not consistent with `fn_deriv`.
pub fn calc_jacobian(pad: &mut Scratchpad, ksi: &[f64], xxt: &Matrix, fn_deriv: FnDeriv) -> Result<f64, StrError> {
    fn_deriv(&mut pad.deriv, ksi);
    mat_mat_mul(&mut pad.jacobian, 1.0, xxt, &pad.deriv)?;
    let (space_ndim, geo_ndim) = pad.jacobian.dims();
    if geo_ndim == space_ndim {
        inverse(&mut pad.inv_jacobian, &pad.jacobian)
    } else {
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
    use crate::shapes::{FnInterp, GeoKind, Scratchpad};
    use crate::StrError;
    use russell_chk::assert_deriv_approx_eq;
    use russell_lab::{Matrix, Vector};

    // Holds arguments for numerical differentiation of x with respect to ξ => Jacobian
    struct ArgsNumJac<'a> {
        fn_interp: FnInterp, // function to calc interpolation function
        xxt: &'a Matrix,     // transposed matrix of coordinates
        pad: Scratchpad,     // scratchpad to send to calc_coords
        at_ksi: Vec<f64>,    // at reference coord value
        ksi: Vec<f64>,       // temporary reference coord
        x: Vector,           // (space_ndim) coordinates at ξ
        i: usize,            // dimension index from 0 to space_ndim
        j: usize,            // dimension index from 0 to geo_ndim
    }

    // Computes xᵢ(ξ) with variable v := ξⱼ
    fn x_given_ksi(v: f64, args: &mut ArgsNumJac) -> f64 {
        args.ksi.copy_from_slice(&args.at_ksi);
        args.ksi[args.j] = v;
        calc_coords(&mut args.x, &mut args.pad, &args.ksi, &args.xxt, args.fn_interp).unwrap();
        args.x[args.i]
    }

    #[test]
    fn calc_jacobian_works() -> Result<(), StrError> {
        // kind and tolerances
        let problem = vec![
            // Lin
            (GeoKind::Lin2, 1e-12),
            (GeoKind::Lin3, 1e-12),
            (GeoKind::Lin4, 1e-12),
            (GeoKind::Lin5, 1e-11),
            // Tri
            (GeoKind::Tri3, 1e-12),
            (GeoKind::Tri6, 1e-11),
            (GeoKind::Tri10, 1e-11),
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
            (GeoKind::Tet10, 1e-12),
            (GeoKind::Tet20, 1e-10),
            // Hex
            (GeoKind::Hex8, 1e-11),
            (GeoKind::Hex20, 1e-11),
            (GeoKind::Hex32, 1e-10),
        ];
        assert_eq!(problem.len(), GeoKind::VALUES.len());

        // loop over shapes
        for (kind, tol) in problem {
            // generate coordinates matrix
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let xxt = aux::gen_coords_transp(space_ndim, kind);

            // scratchpad and derivative function
            let mut pad = Scratchpad::new(space_ndim, kind)?;
            let fn_deriv = kind.functions().1;

            // set ξ within reference space
            let at_ksi = vec![0.25; geo_ndim];

            // compute Jacobian, its inverse, and determinant
            let det_jac = calc_jacobian(&mut pad, &at_ksi, &xxt, fn_deriv)?;
            assert!(det_jac > 0.0);

            // set arguments for numerical integration
            let args = &mut ArgsNumJac {
                fn_interp: kind.functions().0,
                xxt: &xxt,
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
}
