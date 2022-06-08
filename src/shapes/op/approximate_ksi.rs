use super::{calc_coords, calc_jacobian};
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{mat_vec_mul, vector_norm, NormVec, Vector};

/// Approximates the reference coordinates from given real coordinates (inverse mapping)
///
/// **Note:** This function works with `geo_ndim == space_ndim` only.
///
/// This function uses Newton iterations and the inverse of the Jacobian to compute `ξ(x)`.
///
/// # Output
///
/// * `ksi` -- ξ reference coordinates (geo_ndim=space_ndim)
/// * `pad.interp` -- interpolation functions (nnode)
/// * `pad.deriv` -- interpolation functions (nnode,geo_ndim=space_ndim)
/// * `pad.jacobian` -- Jacobian matrix (space_ndim,geo_ndim=space_ndim)
/// * `pad.inv_jacobian` -- inverse Jacobian matrix (space_ndim,space_ndim)
/// * Returns the number of iterations
///
/// # Input
///
/// * `x` -- real coordinates (space_ndim = geo_ndim)
/// * `nit_max` -- maximum number of iterations (e.g., 10)
/// * `tol` -- tolerance for the norm of the difference x - x(ξ) (e.g., 1e-14)
pub fn approximate_ksi(
    ksi: &mut [f64],
    pad: &mut Scratchpad,
    x: &Vector,
    nit_max: usize,
    tol: f64,
) -> Result<usize, StrError> {
    // check
    let (space_ndim, geo_ndim) = pad.jacobian.dims();
    if geo_ndim != space_ndim {
        return Err("geo_ndim must equal space_ndim");
    }
    if x.dim() != space_ndim {
        return Err("x.dim() must equal space_ndim");
    }
    if ksi.len() != geo_ndim {
        return Err("ksi.len() must equal geo_ndim");
    }

    // use linear interpolation to guess ksi
    let (kmin, kdel) = if pad.kind.is_tri_or_tet() {
        (0.0, 1.0) // Tri or Tet
    } else {
        (-1.0, 2.0) // Qua or Hex
    };
    for j in 0..space_ndim {
        ksi[j] = kmin + kdel * (x[j] - pad.xmin[j]) / (pad.xmax[j] - pad.xmin[j]);
    }

    // perform iterations
    let mut residual = Vector::new(space_ndim);
    let mut x_at_ksi = Vector::new(space_ndim);
    let mut delta_ksi = Vector::new(geo_ndim);
    for it in 0..nit_max {
        // check residual
        calc_coords(&mut x_at_ksi, pad, ksi)?;
        for i in 0..space_ndim {
            residual[i] = x[i] - x_at_ksi[i];
        }
        if vector_norm(&residual, NormVec::Euc) <= tol {
            return Ok(it);
        }

        // calc Jacobian
        calc_jacobian(pad, ksi)?;

        // calc ksi increment
        mat_vec_mul(&mut delta_ksi, 1.0, &pad.inv_jacobian, &residual).unwrap(); // cannot fail because all dims have been checked

        // update ksi
        for j in 0..geo_ndim {
            ksi[j] += delta_ksi[j];
        }
    }

    Err("approximate_ksi failed to converge")
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::approximate_ksi;
    use crate::shapes::op::calc_coords;
    use crate::shapes::op::testing::aux;
    use crate::shapes::GeoKind;
    use crate::shapes::Scratchpad;
    use crate::util::{ONE_BY_3, SQRT_3};
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;

    #[test]
    fn approximate_ksi_works() -> Result<(), StrError> {
        // select all kinds, except Lin
        let problem = vec![
            // Tri
            (GeoKind::Tri3, 1e-14),
            (GeoKind::Tri6, 1e-15),
            (GeoKind::Tri10, 1e-15),
            (GeoKind::Tri15, 1e-14),
            // Qua
            (GeoKind::Qua4, 1e-15),
            (GeoKind::Qua8, 1e-15),
            (GeoKind::Qua9, 1e-15),
            (GeoKind::Qua12, 1e-14),
            (GeoKind::Qua16, 1e-14),
            (GeoKind::Qua17, 1e-13),
            // Tet
            (GeoKind::Tet4, 1e-15),
            (GeoKind::Tet10, 1e-15),
            (GeoKind::Tet20, 1e-14),
            // Hex
            (GeoKind::Hex8, 1e-15),
            (GeoKind::Hex20, 1e-14),
            (GeoKind::Hex32, 1e-14),
        ];

        // loop over shapes
        for (kind, tol) in problem {
            // scratchpad with coordinates
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let mut pad = aux::gen_scratchpad_with_coords(space_ndim, kind);

            // loop over nodes of shape
            let nnode = kind.nnode();
            let mut x = Vector::new(space_ndim);
            let mut ksi = vec![0.0; geo_ndim];
            for m in 0..nnode {
                // get ξᵐ corresponding to node m
                let ksi_ref = kind.reference_coords(m);

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                calc_coords(&mut x, &mut pad, ksi_ref)?;

                // compute approximation of the inverse mapping ξᵐ(xᵐ)
                let nit = approximate_ksi(&mut ksi, &mut pad, &x, 10, 1e-14)?;

                // check (linear and bi-linear shapes converge with nit = 1)
                if kind == GeoKind::Tri3 || kind == GeoKind::Qua4 || kind == GeoKind::Tet4 || kind == GeoKind::Hex8 {
                    assert_eq!(nit, 1);
                }
                assert_vec_approx_eq!(ksi, ksi_ref, tol);
            }

            // test again inside the reference domain
            let ksi_in = if kind.is_tri_or_tet() {
                vec![ONE_BY_3; geo_ndim]
            } else {
                vec![0.0; geo_ndim]
            };
            calc_coords(&mut x, &mut pad, &ksi_in)?;
            approximate_ksi(&mut ksi, &mut pad, &x, 10, 1e-14)?;
            assert_vec_approx_eq!(ksi, ksi_in, tol);
        }
        Ok(())
    }

    #[test]
    fn approximate_ksi_works_outside() -> Result<(), StrError> {
        // Equilateral triangle
        //
        //           /
        //        2   \
        //       / \   \
        //      / ↑ \   l
        //     5  h  4   \
        //    /   ↓   \   \
        //   /         \   /
        //  0-----3-----1
        //
        //  |--s--|--s--|
        //
        //  |-----l-----|
        //
        // area = l * h / 2.0;
        let l = 5.0;
        let s = l / 2.0;
        let h = l * SQRT_3 / 2.0;
        let (x0, y0) = (3.0, 4.0);
        let (x1, y1) = (x0 + l, y0);
        let (x2, y2) = (x0 + s, y0 + h);
        let (x3, y3) = (x0 + s, y0);
        let (x4, y4) = (x0 + 1.5 * s, y0 + 0.5 * h);
        let (x5, y5) = (x0 + 0.5 * s, y0 + 0.5 * h);
        let space_ndim = 2;
        let mut pad = Scratchpad::new(space_ndim, GeoKind::Tri6)?;
        pad.set_xx(0, 0, x0);
        pad.set_xx(0, 1, y0);
        pad.set_xx(1, 0, x1);
        pad.set_xx(1, 1, y1);
        pad.set_xx(2, 0, x2);
        pad.set_xx(2, 1, y2);
        pad.set_xx(3, 0, x3);
        pad.set_xx(3, 1, y3);
        pad.set_xx(4, 0, x4);
        pad.set_xx(4, 1, y4);
        pad.set_xx(5, 0, x5);
        pad.set_xx(5, 1, y5);
        assert_eq!(
            format!("{:.2}", pad.xxt),
            "┌                               ┐\n\
             │ 3.00 8.00 5.50 5.50 6.75 4.25 │\n\
             │ 4.00 4.00 8.33 4.00 6.17 6.17 │\n\
             └                               ┘"
        );
        let mut ksi = vec![0.0; pad.kind.ndim()];
        for (nit_correct, x_data, tol) in &[
            (0, [3.0, 4.0], 1e-15),
            (0, [8.0, 4.0], 1e-15),
            (1, [5.5, 8.33], 1e-14),
            (0, [5.5, 4.0], 1e-15),
            (1, [6.75, 6.17], 1e-14),
            (1, [4.25, 6.17], 1e-14),
            (2, [10.0, 10.0], 1e-14),
            (2, [-10.0, -10.0], 1e-13),
            (6, [100.0, 100.0], 1e-12),
        ] {
            let x = Vector::from(x_data);
            let nit = approximate_ksi(&mut ksi, &mut pad, &x, 30, *tol)?;
            let mut x_out = Vector::new(2);
            calc_coords(&mut x_out, &mut pad, &ksi)?;
            assert_vec_approx_eq!(x.as_data(), x_out.as_data(), *tol);
            assert_eq!(nit, *nit_correct);
        }
        Ok(())
    }
}
