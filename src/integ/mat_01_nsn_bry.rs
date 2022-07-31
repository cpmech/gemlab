use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the shape(N) times scalar(S) times shape(N) integration case 01 (boundary integral version)
///
/// Diffusion coefficients:
///
/// ```text
///       ⌠
/// Kᵐⁿ = │ Nᵐ s Nⁿ dΩ
///       ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///       nip-1    →     →      →       →
/// Kᵐⁿ ≈   Σ   Nᵐ(ιᵖ) s(ιᵖ) Nⁿ(ιᵖ) |J|(ιᵖ) wᵖ
///        p=0
/// ```
///
/// # Output
///
/// ```text
///     ┌                     ┐
///     | K⁰⁰ K⁰¹ K⁰² ··· K⁰ⁿ |  ⟸  ii0
///     | K¹⁰ K¹¹ K¹² ··· K¹ⁿ |
/// K = | K²⁰ K²¹ K²² ··· K²ⁿ |
///     |  ··  ··  ·· ···  ·· |
///     | Kᵐ⁰ Kᵐ¹ Kᵐ² ··· Kᵐⁿ |  ⟸  ii
///     └                     ┘
///        ⇑                ⇑
///       jj0               jj
/// ```
///
/// * `kk` -- A matrix containing all `Kᵐⁿ` values, one after another, and
///   sequentially placed as shown above. `m` and `n` are the indices of the nodes.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode` and `ncol(K) ≥ jj0 + nnode`
/// * `pad` -- Some members of the scratchpad will be modified.
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output matrix where to add components.
/// * `jj0` -- Stride marking the first column in the output matrix where to add components.
/// * `clear_kk` -- Fills `kk` matrix with zeros, otherwise accumulate values into `kk`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_s` -- Function `f(p)` that computes `s(x(ιᵖ))` with `0 ≤ p ≤ n_integ_point`
pub fn mat_01_nsn_bry<F>(
    kk: &mut Matrix,
    pad: &mut Scratchpad,
    ii0: usize,
    jj0: usize,
    clear_kk: bool,
    ips: IntegPointData,
    fn_s: F,
) -> Result<(), StrError>
where
    F: Fn(usize) -> Result<f64, StrError>,
{
    // check
    let (space_ndim, nnode) = pad.xxt.dims();
    if space_ndim == 2 && pad.kind.ndim() != 1 {
        return Err("in 2D, geometry ndim must be 1 (a line)");
    }
    if space_ndim == 3 && pad.kind.ndim() != 2 {
        return Err("in 3D, geometry ndim must be 2 (a surface)");
    }
    let (nrow_kk, ncol_kk) = kk.dims();
    if nrow_kk < ii0 + nnode {
        return Err("nrow(K) must be ≥ ii0 + nnode");
    }
    if ncol_kk < jj0 + nnode {
        return Err("ncol(K) must be ≥ jj0 + nnode");
    }

    // allocate auxiliary vector
    let mut un = Vector::new(space_ndim); // unit normal vector

    // clear output matrix
    if clear_kk {
        kk.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate interpolation functions and Jacobian
        (pad.fn_interp)(&mut pad.interp, iota);
        let mag_n = pad.calc_normal_vector(&mut un, iota)?;

        // calculate s
        let s = fn_s(p)?;

        // add contribution to K matrix
        let val = s * mag_n * weight;
        let nn = &pad.interp;
        for m in 0..nnode {
            for n in 0..nnode {
                kk[ii0 + m][jj0 + n] += nn[m] * val * nn[n];
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, AnalyticalTri3, IP_LIN_LEGENDRE_1, IP_TRI_INTERNAL_1};
    use crate::shapes::{GeoClass, GeoKind, Scratchpad};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Matrix;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(2, 2);
        assert_eq!(
            integ::mat_01_nsn_bry(&mut kk, &mut pad, 1, 0, false, &[], |_| Ok(0.0)).err(),
            Some("nrow(K) must be ≥ ii0 + nnode")
        );
        assert_eq!(
            integ::mat_01_nsn_bry(&mut kk, &mut pad, 0, 1, false, &[], |_| Ok(0.0)).err(),
            Some("ncol(K) must be ≥ jj0 + nnode")
        );
        // more errors
        assert_eq!(
            integ::mat_01_nsn_bry(&mut kk, &mut pad, 0, 0, false, &IP_LIN_LEGENDRE_1, |_| Ok(0.0)).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        assert_eq!(
            integ::mat_01_nsn_bry(&mut kk, &mut pad, 0, 0, false, &IP_TRI_INTERNAL_1, |_| Err("stop")).err(),
            Some("stop")
        );
    }

    /// Returns the scratchpad for a Tri3's side
    ///
    /// **Important:** `side` must be 0, 1, or 2
    fn tri3_extract_pad_side(pad: &Scratchpad, side: usize) -> Scratchpad {
        let mut pad_side = Scratchpad::new(2, GeoKind::Lin2).unwrap();
        for i in 0..2 {
            let m = pad.kind.edge_node_id(side, i);
            pad_side.set_xx(i, 0, pad.xxt[0][m]);
            pad_side.set_xx(i, 1, pad.xxt[1][m]);
        }
        pad_side
    }

    #[test]
    fn mat_01_nsn_bry_tri3_works() {
        let pad = aux::gen_pad_tri3();
        let s = 12.0;
        let ana = AnalyticalTri3::new(&pad);
        let ips = integ::points(GeoClass::Lin, 2).unwrap();

        // side # 0
        let mut pad_side = tri3_extract_pad_side(&pad, 0);
        let kk_correct = ana.mat_01_nsn_bry(0, s);
        let mut kk = Matrix::new(2, 2);
        integ::mat_01_nsn_bry(&mut kk, &mut pad_side, 0, 0, true, ips, |_| Ok(s)).unwrap();
        assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), 1e-15);

        // side # 1
        let mut pad_side = tri3_extract_pad_side(&pad, 1);
        let kk_correct = ana.mat_01_nsn_bry(1, s);
        integ::mat_01_nsn_bry(&mut kk, &mut pad_side, 0, 0, true, ips, |_| Ok(s)).unwrap();
        assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), 1e-14);

        // side # 2
        let mut pad_side = tri3_extract_pad_side(&pad, 2);
        let kk_correct = ana.mat_01_nsn_bry(2, s);
        integ::mat_01_nsn_bry(&mut kk, &mut pad_side, 0, 0, true, ips, |_| Ok(s)).unwrap();
        assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), 1e-14);
    }
}
