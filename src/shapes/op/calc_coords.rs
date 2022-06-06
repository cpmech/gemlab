use crate::shapes::{FnInterp, Scratchpad};
use russell_lab::{mat_vec_mul, Matrix, Vector};

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
/// * `xxt` -- is the transposed coordinates matrix `Xᵀ` shown below (space_ndim,nnode)
///
/// The transposed coordinates matrix (real space) is such that:
///
/// ```text
///      ┌                              ┐  superscript = node
///      | x⁰₀  x¹₀  x²₀  x³₀       xᴹ₀ |  subscript = dimension
/// Xᵀ = | x⁰₁  x¹₁  x²₁  x³₁  ...  xᴹ₁ |
///      | x⁰₂  x¹₂  x²₂  x³₂       xᴹ₂ |
///      └                              ┘_(space_ndim,nnode)
/// ```
///
/// where `M = nnode - 1`
///
/// # Panics
///
/// This function does not check for the vector/matrix dimensions. Thus a panic may occur
/// if they are incompatible, including if they are not consistent with `fn_interp`.
#[inline]
pub fn calc_coords(x: &mut Vector, pad: &mut Scratchpad, ksi: &[f64], xxt: &Matrix, fn_interp: FnInterp) {
    fn_interp(&mut pad.interp, ksi);
    mat_vec_mul(x, 1.0, &xxt, &pad.interp).unwrap();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::calc_coords;
    use crate::shapes::op::testing::aux;
    use crate::shapes::{GeoKind, Scratchpad};
    use crate::util::ONE_BY_3;
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;

    #[test]
    fn calc_coords_works() -> Result<(), StrError> {
        // kind, tol, tol_in for the case inside the shape
        let problem = vec![
            (GeoKind::Tri3, 1e-15, 0.45), // linear maps are inaccurate for the circular wedge
            (GeoKind::Tri6, 1e-15, 0.02), // << quadratic mapping is inaccurate as well
            (GeoKind::Tri10, 1e-14, 1e-14),
            (GeoKind::Tri15, 1e-14, 1e-4), // << this triangle is inaccurate as well here
            (GeoKind::Qua4, 1e-15, 0.14),  // linear maps are inaccurate for the circular wedge
            (GeoKind::Qua8, 1e-15, 1e-14),
            (GeoKind::Qua17, 1e-15, 1e-14),
            (GeoKind::Tet4, 1e-15, 0.45),   // linear tetrahedron is also inaccurate here
            (GeoKind::Tet10, 1e-15, 0.02),  // quadratic tetrahedron is also inaccurate here
            (GeoKind::Tet20, 1e-14, 1e-14), // cubic tetrahedron
            (GeoKind::Hex8, 1e-14, 0.14),   // bi-linear maps are inaccurate for the circular wedge
            (GeoKind::Hex20, 1e-15, 1e-14),
            (GeoKind::Hex32, 1e-15, 1e-4), // TODO: check why this tolerance is high
        ];

        // loop over shapes
        for (kind, tol, tol_in) in problem {
            // generate coordinates matrix
            let geo_ndim = kind.ndim();
            let space_ndim = usize::max(2, geo_ndim);
            let xxt = aux::gen_coords_transp(space_ndim, kind);

            // scratchpad and interpolation function
            let mut pad = Scratchpad::new(space_ndim, kind)?;
            let fn_interp = kind.functions().0;

            // loop over nodes of shape
            let nnode = kind.nnode();
            let mut x = Vector::new(space_ndim);
            let mut x_correct = Vector::new(space_ndim);
            for m in 0..nnode {
                // get ξᵐ corresponding to node m
                let ksi = kind.reference_coords(m);

                // calculate xᵐ(ξᵐ) using the isoparametric formula
                calc_coords(&mut x, &mut pad, ksi, &xxt, fn_interp);

                // compare xᵐ with generated coordinates
                aux::map_point_coords(&mut x_correct, ksi, kind);
                assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol);
            }

            // test again inside the reference domain
            let ksi_in = if kind.is_tri_or_tet() {
                vec![ONE_BY_3; geo_ndim]
            } else {
                vec![0.0; geo_ndim]
            };
            calc_coords(&mut x, &mut pad, &ksi_in, &xxt, fn_interp);
            aux::map_point_coords(&mut x_correct, &ksi_in, kind);
            assert_vec_approx_eq!(x.as_data(), x_correct.as_data(), tol_in);
        }
        Ok(())
    }
}
