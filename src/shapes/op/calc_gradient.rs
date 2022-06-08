use super::calc_jacobian;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::mat_mat_mul;

/// Calculates the gradient of the interpolation functions
///
/// **Note:** This function works with `geo_ndim == space_ndim` only.
///
/// The gradient is given by:
///
/// ```text
///             →
/// →  →    dNᵐ(ξ)
/// Gᵐ(ξ) = ——————
///            →
///           dx
/// ```
///
/// which can be organized in an (nnode,space_ndim) matrix `G` as follows
///
/// ```text
/// G = L · J⁻¹
/// ```
///
/// # Output
///
/// * `pad.deriv` -- interpolation functions (nnode)
/// * `pad.jacobian` -- Jacobian matrix (space_ndim,geo_ndim)
/// * `pad.inv_jacobian` -- inverse Jacobian matrix (space_ndim,space_ndim)
/// * `pad.gradient` -- gradient matrix (nnode,space_ndim)
/// * Returns the determinant of the Jacobian matrix
///
/// # Input
///
/// * `ksi` -- reference coordinates ξ with len ≥ geo_ndim
///
/// # Example
///
/// ```
/// use gemlab::shapes::{op, GeoKind, Scratchpad};
/// use gemlab::StrError;
/// use russell_chk::assert_vec_approx_eq;
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
///     op::calc_gradient(&mut pad, &[0.0, 0.0])?;
///
///     let correct_gg = Matrix::from(&[
///         [-1.0 / (4.0 * a), -1.0 / (2.0 * a)],
///         [1.0 / (4.0 * a), -1.0 / (2.0 * a)],
///         [1.0 / (4.0 * a), 1.0 / (2.0 * a)],
///         [-1.0 / (4.0 * a), 1.0 / (2.0 * a)],
///     ]);
///     assert_vec_approx_eq!(pad.gradient.as_data(), correct_gg.as_data(), 1e-15);
///     Ok(())
/// }
/// ```
pub fn calc_gradient(pad: &mut Scratchpad, ksi: &[f64]) -> Result<f64, StrError> {
    // check
    let (space_ndim, geo_ndim) = pad.jacobian.dims();
    if geo_ndim != space_ndim {
        return Err("geo_ndim must equal space_ndim");
    }

    // Jacobian matrix J: dx/dξ
    let det_jac = calc_jacobian(pad, ksi)?;

    // gradient: G = L · J⁻¹
    mat_mat_mul(&mut pad.gradient, 1.0, &pad.deriv, &pad.inv_jacobian).unwrap(); // cannot fail because the dims are checked
    Ok(det_jac)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::calc_gradient;
    use crate::shapes::op::testing::aux;
    use crate::shapes::op::{approximate_ksi, calc_coords};
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::util::ONE_BY_3;
    use crate::StrError;
    use russell_chk::assert_deriv_approx_eq;
    use russell_lab::{copy_vector, Vector};

    // Holds arguments for numerical differentiation of N with respect to x => G (gradient) matrix
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
        copy_vector(&mut args.x, &args.at_x).unwrap();
        args.x[args.j] = v;
        approximate_ksi(&mut args.ksi, &mut args.pad, &args.x, 10, 1e-14).unwrap();
        (args.pad.fn_interp)(&mut args.pad.interp, &args.ksi);
        args.pad.interp[args.m]
    }

    #[test]
    fn calc_gradient_works() -> Result<(), StrError> {
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
            (GeoKind::Qua12, 1e-10),
            (GeoKind::Qua16, 1e-10),
            (GeoKind::Qua17, 1e-10),
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
            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let mut pad = aux::gen_scratchpad_with_coords(space_ndim, kind);

            // set ξ within reference space
            let at_ksi = vec![ONE_BY_3; geo_ndim];

            // compute x corresponding to ξ using the isoparametric formula
            let mut at_x = Vector::new(space_ndim);
            calc_coords(&mut at_x, &mut pad, &at_ksi)?;

            // compute gradient
            let det_jac = calc_gradient(&mut pad, &at_ksi)?;
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

            // check Gᵐ(ξ(x)) = dNᵐ(ξ(x))/dx
            for m in 0..kind.nnode() {
                args.m = m;
                for j in 0..geo_ndim {
                    args.j = j;
                    // Gᵐⱼ := dNᵐ/dxⱼ
                    assert_deriv_approx_eq!(pad.gradient[m][j], args.at_x[j], nn_given_x, args, tol);
                }
            }
        }
        Ok(())
    }
}
