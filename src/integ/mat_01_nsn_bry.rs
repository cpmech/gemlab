use super::CommonArgs;
use crate::StrError;
use russell_lab::{Matrix, Vector};

/// Implements the shape(N) times scalar(S) times shape(N) integration case 01 (boundary integral version)
///
/// Callback function: `s ← f(p, un, N)`
///
/// Diffusion coefficients:
///
/// ```text
///       ⌠
/// Kᵐⁿ = │ Nᵐ s Nⁿ α dΓ
///       ⌡
///       Γₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///       nip-1    →     →      →       →
/// Kᵐⁿ ≈   Σ   Nᵐ(ιᵖ) s(ιᵖ) Nⁿ(ιᵖ) |J|(ιᵖ) wᵖ α
///        p=0
/// ```
///
/// # Results
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
/// # Arguments
///
/// * `kk` -- A matrix containing all `Kᵐⁿ` values, one after another, and
///   sequentially placed as shown above. `m` and `n` are the indices of the nodes.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode` and `ncol(K) ≥ jj0 + nnode`
/// * `args` --- Common arguments
/// * `fn_s` -- Function `f(p,un,N)→s` that computes `s(x(ιᵖ))`, given `0 ≤ p ≤ n_integ_point`,
///   the **unit** normal vector `un(x(ιᵖ))`, and shape functions N(ιᵖ).
pub fn mat_01_nsn_bry<F>(kk: &mut Matrix, args: &mut CommonArgs, mut fn_s: F) -> Result<(), StrError>
where
    F: FnMut(usize, &Vector, &Vector) -> Result<f64, StrError>,
{
    // check
    let (space_ndim, nnode) = args.pad.xxt.dims();
    let geo_ndim = args.pad.deriv.dims().1;
    if space_ndim == 2 && geo_ndim != 1 {
        return Err("in 2D, geometry ndim must be equal to 1 (a line)");
    }
    if space_ndim == 3 && geo_ndim != 2 {
        return Err("in 3D, geometry ndim must be equal to 2 (a surface)");
    }
    let (nrow_kk, ncol_kk) = kk.dims();
    let (ii0, jj0) = (args.ii0, args.jj0);
    if nrow_kk < ii0 + nnode {
        return Err("nrow(K) must be ≥ ii0 + nnode");
    }
    if ncol_kk < jj0 + nnode {
        return Err("ncol(K) must be ≥ jj0 + nnode");
    }

    // allocate auxiliary vector
    let mut un = Vector::new(space_ndim); // unit normal vector

    // clear output matrix
    if args.clear {
        kk.fill(0.0);
    }

    // loop over integration points
    for p in 0..args.ips.len() {
        // ksi coordinates and weight
        let iota = &args.ips[p];
        let weight = args.ips[p][3];

        // calculate interpolation functions and Jacobian
        (args.pad.fn_interp)(&mut args.pad.interp, iota); // N
        let mag_n = args.pad.calc_normal_vector(&mut un, iota).unwrap(); // un

        // calculate s
        let nn = &args.pad.interp;
        let s = fn_s(p, &un, nn)?;

        // calculate coefficient
        let c = if args.axisymmetric {
            let mut r = 0.0; // radius @ x(ιᵖ)
            for m in 0..nnode {
                r += nn[m] * args.pad.xxt[0][m];
            }
            s * mag_n * weight * args.alpha * r
        } else {
            s * mag_n * weight * args.alpha
        };

        // add contribution to K matrix
        for m in 0..nnode {
            for n in 0..nnode {
                kk[ii0 + m][jj0 + n] += nn[m] * c * nn[n];
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, AnalyticalTri3, CommonArgs, IP_LIN_LEGENDRE_1};
    use crate::shapes::{GeoClass, GeoKind, Scratchpad};
    use russell_chk::vec_approx_eq;
    use russell_lab::{Matrix, Vector};

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(3, 3);
        let un = Vector::new(0);
        let nn = Vector::new(0);
        let f = |_, _: &Vector, _: &Vector| Ok(0.0);
        assert_eq!(f(0, &un, &nn).unwrap(), 0.0);
        let mut args = CommonArgs::new(&mut pad, &[]);
        assert_eq!(
            integ::mat_01_nsn_bry(&mut kk, &mut args, f).err(),
            Some("in 2D, geometry ndim must be equal to 1 (a line)")
        );
        let mut pad = aux::gen_pad_tet4();
        let mut args = CommonArgs::new(&mut pad, &[]);
        let mut kk = Matrix::new(4, 4);
        assert_eq!(
            integ::mat_01_nsn_bry(&mut kk, &mut args, f).err(),
            Some("in 3D, geometry ndim must be equal to 2 (a surface)")
        );
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut args = CommonArgs::new(&mut pad, &[]);
        let mut kk = Matrix::new(2, 2);
        args.ii0 = 1;
        assert_eq!(
            integ::mat_01_nsn_bry(&mut kk, &mut args, f).err(),
            Some("nrow(K) must be ≥ ii0 + nnode")
        );
        args.ii0 = 0;
        args.jj0 = 1;
        assert_eq!(
            integ::mat_01_nsn_bry(&mut kk, &mut args, f).err(),
            Some("ncol(K) must be ≥ jj0 + nnode")
        );
        args.jj0 = 0;
        args.ips = &IP_LIN_LEGENDRE_1;
        assert_eq!(
            integ::mat_01_nsn_bry(&mut kk, &mut args, |_, _, _| Err("stop")).err(),
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
        let kk_correct = ana.mat_01_nsn_bry(0, s, false);
        let mut kk = Matrix::new(2, 2);
        let mut args = CommonArgs::new(&mut pad_side, ips);
        integ::mat_01_nsn_bry(&mut kk, &mut args, |_, _, _| Ok(s)).unwrap();
        vec_approx_eq(kk.as_data(), kk_correct.as_data(), 1e-15);

        // side # 1
        let mut pad_side = tri3_extract_pad_side(&pad, 1);
        let kk_correct = ana.mat_01_nsn_bry(1, s, false);
        let mut args = CommonArgs::new(&mut pad_side, ips);
        integ::mat_01_nsn_bry(&mut kk, &mut args, |_, _, _| Ok(s)).unwrap();
        vec_approx_eq(kk.as_data(), kk_correct.as_data(), 1e-14);

        // side # 2
        let mut pad_side = tri3_extract_pad_side(&pad, 2);
        let kk_correct = ana.mat_01_nsn_bry(2, s, false);
        let mut args = CommonArgs::new(&mut pad_side, ips);
        integ::mat_01_nsn_bry(&mut kk, &mut args, |_, _, _| Ok(s)).unwrap();
        vec_approx_eq(kk.as_data(), kk_correct.as_data(), 1e-14);
    }
}
