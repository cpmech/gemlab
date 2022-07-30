use super::IntegPointData;
use crate::shapes::Scratchpad;
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::Matrix;
use russell_tensor::Tensor4;

/// Implements the gradient(G) dot 4th-tensor(D) dot gradient(G) integration case 10 (e.g., stiffness matrix)
///
/// Stiffness tensors (note the sum over repeated lower indices):
///
/// ```text
///       ⌠                       →    →
/// Kᵐⁿ = │ Σ Σ Σ Σ Gᵐₖ Dᵢₖⱼₗ Gⁿₗ eᵢ ⊗ eⱼ dΩ
/// ▔     ⌡ i j k l
///       Ωₑ
/// ```
///
/// The numerical integration is (note the sum over repeated lower indices):
///
/// ```text
///         nip-1         →         →       →          →
/// Kᵐⁿᵢⱼ ≈   Σ   Σ Σ Gᵐₖ(ιᵖ) Dᵢₖⱼₗ(ιᵖ) Gⁿₗ(ιᵖ) |J|(ιᵖ) wᵖ
///          p=0  k l
/// ```
///
/// # Output
///
/// ```text
///     ┌                                               ┐
///     | K⁰⁰₀₀ K⁰⁰₀₁ K⁰¹₀₀ K⁰¹₀₁ K⁰²₀₀ K⁰²₀₁ ··· K⁰ⁿ₀ⱼ |  ⟸  ii0
///     | K⁰⁰₁₀ K⁰⁰₁₁ K⁰¹₁₀ K⁰¹₁₁ K⁰²₁₀ K⁰²₁₁ ··· K⁰ⁿ₁ⱼ |
///     | K¹⁰₀₀ K¹⁰₀₁ K¹¹₀₀ K¹¹₀₁ K¹²₀₀ K¹²₀₁ ··· K¹ⁿ₀ⱼ |
/// K = | K¹⁰₁₀ K¹⁰₁₁ K¹¹₁₀ K¹¹₁₁ K¹²₁₀ K¹²₁₁ ··· K¹ⁿ₁ⱼ |
///     | K²⁰₀₀ K²⁰₀₁ K²¹₀₀ K²¹₀₁ K²²₀₀ K²²₀₁ ··· K²ⁿ₀ⱼ |
///     | K²⁰₁₀ K²⁰₁₁ K²¹₁₀ K²¹₁₁ K²²₁₀ K²²₁₁ ··· K²ⁿ₁ⱼ |
///     |  ···   ···   ···   ···   ···   ···  ···  ···  |
///     | Kᵐ⁰ᵢ₀ Kᵐ⁰ᵢ₁ Kᵐ¹ᵢ₀ Kᵐ¹ᵢ₁ Kᵐ²ᵢ₀ Kᵐ²ᵢ₁ ··· Kᵐⁿᵢⱼ |  ⟸  ii := i + m ⋅ space_ndim
///     └                                               ┘
///        ⇑                                        ⇑
///       jj0                                       jj := j + n ⋅ space_ndim
///
/// m = ii / space_ndim    n = jj / space_ndim
/// i = ii % space_ndim    j = jj % space_ndim
/// ```
///
/// * `kk` -- A matrix containing all `Kᵐⁿᵢⱼ` values, one after another, and sequentially placed as shown
///   above (in 2D). `m` and `n` are the indices of the node and `i` and `j` correspond to `space_ndim`.
///   The dimensions must be `nrow(K) ≥ ii0 + nnode ⋅ space_ndim` and `ncol(K) ≥ jj0 + nnode ⋅ space_ndim`.
/// * `pad` -- Some members of the scratchpad will be modified.
///
/// # Input
///
/// * `ii0` -- Stride marking the first row in the output matrix where to add components.
/// * `jj0` -- Stride marking the first column in the output matrix where to add components.
/// * `clear_kk` -- Fills `kk` matrix with zeros, otherwise accumulate values into `kk`
/// * `ips` -- Integration points (n_integ_point)
/// * `fn_dd` -- Function f(D,p) that computes `D(x(ιᵖ))`
///
/// # Examples
///
/// ```
/// use gemlab::integ;
/// use gemlab::shapes::{GeoKind, Scratchpad};
/// use gemlab::StrError;
/// use russell_lab::{copy_matrix, Matrix};
/// use russell_tensor::LinElasticity;
///
/// fn main() -> Result<(), StrError> {
///     // shape and state
///     let space_ndim = 3;
///     let mut pad = Scratchpad::new(space_ndim, GeoKind::Tet4)?;
///     pad.set_xx(0, 0, 2.0);
///     pad.set_xx(0, 1, 3.0);
///     pad.set_xx(0, 2, 4.0);
///     pad.set_xx(1, 0, 6.0);
///     pad.set_xx(1, 1, 3.0);
///     pad.set_xx(1, 2, 2.0);
///     pad.set_xx(2, 0, 2.0);
///     pad.set_xx(2, 1, 5.0);
///     pad.set_xx(2, 2, 1.0);
///     pad.set_xx(3, 0, 4.0);
///     pad.set_xx(3, 1, 3.0);
///     pad.set_xx(3, 2, 6.0);
///
///     // constants
///     let young = 480.0;
///     let poisson = 1.0 / 3.0;
///     let two_dim = false;
///     let plane_stress = false;
///     let model = LinElasticity::new(young, poisson, two_dim, plane_stress);
///
///     // stiffness
///     let nrow = pad.kind.nnode() * space_ndim;
///     let mut kk = Matrix::new(nrow, nrow);
///     let ips = integ::default_points(pad.kind);
///     integ::mat_10_gdg(&mut kk, &mut pad, 0, 0, true, ips, |dd, _| {
///         copy_matrix(&mut dd.mat, &model.get_modulus().mat)
///     })?;
///
///     // check
///     let correct =
///     "┌                                                                         ┐\n\
///      │   745   540   120    -5    30    60  -270  -240     0  -470  -330  -180 │\n\
///      │   540  1720   270  -120   520   210  -120 -1080   -60  -300 -1160  -420 │\n\
///      │   120   270   565     0   150   175     0  -120  -270  -120  -300  -470 │\n\
///      │    -5  -120     0   145   -90   -60   -90   120     0   -50    90    60 │\n\
///      │    30   520   150   -90   220    90    60  -360   -60     0  -380  -180 │\n\
///      │    60   210   175   -60    90   145     0  -120   -90     0  -180  -230 │\n\
///      │  -270  -120     0   -90    60     0   180     0     0   180    60     0 │\n\
///      │  -240 -1080  -120   120  -360  -120     0   720     0   120   720   240 │\n\
///      │     0   -60  -270     0   -60   -90     0     0   180     0   120   180 │\n\
///      │  -470  -300  -120   -50     0     0   180   120     0   340   180   120 │\n\
///      │  -330 -1160  -300    90  -380  -180    60   720   120   180   820   360 │\n\
///      │  -180  -420  -470    60  -180  -230     0   240   180   120   360   520 │\n\
///      └                                                                         ┘";
///     assert_eq!(format!("{:.0}", kk), correct);
///     Ok(())
/// }
/// ```
pub fn mat_10_gdg<F>(
    kk: &mut Matrix,
    pad: &mut Scratchpad,
    ii0: usize,
    jj0: usize,
    clear_kk: bool,
    ips: IntegPointData,
    fn_dd: F,
) -> Result<(), StrError>
where
    F: Fn(&mut Tensor4, usize) -> Result<(), StrError>,
{
    // check
    let (space_ndim, nnode) = pad.xxt.dims();
    let (nrow_kk, ncol_kk) = kk.dims();
    if nrow_kk < ii0 + nnode * space_ndim {
        return Err("nrow(K) must be ≥ ii0 + nnode ⋅ space_ndim");
    }
    if ncol_kk < jj0 + nnode * space_ndim {
        return Err("ncol(K) must be ≥ jj0 + nnode ⋅ space_ndim");
    }

    // allocate auxiliary tensor
    let mut dd = Tensor4::new(true, space_ndim == 2);

    // clear output matrix
    if clear_kk {
        kk.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate Jacobian and Gradient
        let det_jac = pad.calc_gradient(iota)?;

        // calculate D tensor
        fn_dd(&mut dd, p)?;

        // add contribution to K matrix
        let c = det_jac * weight;
        mat_10_gdg_add_to_mat_kk(kk, ii0, jj0, &dd, c, pad);
    }
    Ok(())
}

/// Adds contribution to the K-matrix in integ_mat_10_gdg
#[inline]
#[rustfmt::skip]
fn mat_10_gdg_add_to_mat_kk(kk: &mut Matrix, ii0: usize, jj0: usize, dd: &Tensor4, c: f64, pad: &mut Scratchpad) {
    let s = SQRT_2;
    let g = &pad.gradient;
    let d = &dd.mat;
    let (space_ndim, nnode) = pad.xxt.dims();
    if space_ndim == 2 {
        for m in 0..nnode {
            for n in 0..nnode {
                kk[ii0+0+m*2][jj0+0+n*2] += c * (g[m][1]*g[n][1]*d[3][3] + s*g[m][1]*g[n][0]*d[3][0] + s*g[m][0]*g[n][1]*d[0][3] + 2.0*g[m][0]*g[n][0]*d[0][0]) / 2.0;
                kk[ii0+0+m*2][jj0+1+n*2] += c * (g[m][1]*g[n][0]*d[3][3] + s*g[m][1]*g[n][1]*d[3][1] + s*g[m][0]*g[n][0]*d[0][3] + 2.0*g[m][0]*g[n][1]*d[0][1]) / 2.0;
                kk[ii0+1+m*2][jj0+0+n*2] += c * (g[m][0]*g[n][1]*d[3][3] + s*g[m][0]*g[n][0]*d[3][0] + s*g[m][1]*g[n][1]*d[1][3] + 2.0*g[m][1]*g[n][0]*d[1][0]) / 2.0;
                kk[ii0+1+m*2][jj0+1+n*2] += c * (g[m][0]*g[n][0]*d[3][3] + s*g[m][0]*g[n][1]*d[3][1] + s*g[m][1]*g[n][0]*d[1][3] + 2.0*g[m][1]*g[n][1]*d[1][1]) / 2.0;
            }
        }
    } else {
        for m in 0..nnode {
            for n in 0..nnode {
                kk[ii0+0+m*3][jj0+0+n*3] += c * (g[m][2]*g[n][2]*d[5][5] + g[m][2]*g[n][1]*d[5][3] + s*g[m][2]*g[n][0]*d[5][0] + g[m][1]*g[n][2]*d[3][5] + g[m][1]*g[n][1]*d[3][3] + s*g[m][1]*g[n][0]*d[3][0] + s*g[m][0]*g[n][2]*d[0][5] + s*g[m][0]*g[n][1]*d[0][3] + 2.0*g[m][0]*g[n][0]*d[0][0]) / 2.0;
                kk[ii0+0+m*3][jj0+1+n*3] += c * (g[m][2]*g[n][2]*d[5][4] + g[m][2]*g[n][0]*d[5][3] + s*g[m][2]*g[n][1]*d[5][1] + g[m][1]*g[n][2]*d[3][4] + g[m][1]*g[n][0]*d[3][3] + s*g[m][1]*g[n][1]*d[3][1] + s*g[m][0]*g[n][2]*d[0][4] + s*g[m][0]*g[n][0]*d[0][3] + 2.0*g[m][0]*g[n][1]*d[0][1]) / 2.0;
                kk[ii0+0+m*3][jj0+2+n*3] += c * (g[m][2]*g[n][0]*d[5][5] + g[m][2]*g[n][1]*d[5][4] + s*g[m][2]*g[n][2]*d[5][2] + g[m][1]*g[n][0]*d[3][5] + g[m][1]*g[n][1]*d[3][4] + s*g[m][1]*g[n][2]*d[3][2] + s*g[m][0]*g[n][0]*d[0][5] + s*g[m][0]*g[n][1]*d[0][4] + 2.0*g[m][0]*g[n][2]*d[0][2]) / 2.0;
                kk[ii0+1+m*3][jj0+0+n*3] += c * (g[m][2]*g[n][2]*d[4][5] + g[m][2]*g[n][1]*d[4][3] + s*g[m][2]*g[n][0]*d[4][0] + g[m][0]*g[n][2]*d[3][5] + g[m][0]*g[n][1]*d[3][3] + s*g[m][0]*g[n][0]*d[3][0] + s*g[m][1]*g[n][2]*d[1][5] + s*g[m][1]*g[n][1]*d[1][3] + 2.0*g[m][1]*g[n][0]*d[1][0]) / 2.0;
                kk[ii0+1+m*3][jj0+1+n*3] += c * (g[m][2]*g[n][2]*d[4][4] + g[m][2]*g[n][0]*d[4][3] + s*g[m][2]*g[n][1]*d[4][1] + g[m][0]*g[n][2]*d[3][4] + g[m][0]*g[n][0]*d[3][3] + s*g[m][0]*g[n][1]*d[3][1] + s*g[m][1]*g[n][2]*d[1][4] + s*g[m][1]*g[n][0]*d[1][3] + 2.0*g[m][1]*g[n][1]*d[1][1]) / 2.0;
                kk[ii0+1+m*3][jj0+2+n*3] += c * (g[m][2]*g[n][0]*d[4][5] + g[m][2]*g[n][1]*d[4][4] + s*g[m][2]*g[n][2]*d[4][2] + g[m][0]*g[n][0]*d[3][5] + g[m][0]*g[n][1]*d[3][4] + s*g[m][0]*g[n][2]*d[3][2] + s*g[m][1]*g[n][0]*d[1][5] + s*g[m][1]*g[n][1]*d[1][4] + 2.0*g[m][1]*g[n][2]*d[1][2]) / 2.0;
                kk[ii0+2+m*3][jj0+0+n*3] += c * (g[m][0]*g[n][2]*d[5][5] + g[m][0]*g[n][1]*d[5][3] + s*g[m][0]*g[n][0]*d[5][0] + g[m][1]*g[n][2]*d[4][5] + g[m][1]*g[n][1]*d[4][3] + s*g[m][1]*g[n][0]*d[4][0] + s*g[m][2]*g[n][2]*d[2][5] + s*g[m][2]*g[n][1]*d[2][3] + 2.0*g[m][2]*g[n][0]*d[2][0]) / 2.0;
                kk[ii0+2+m*3][jj0+1+n*3] += c * (g[m][0]*g[n][2]*d[5][4] + g[m][0]*g[n][0]*d[5][3] + s*g[m][0]*g[n][1]*d[5][1] + g[m][1]*g[n][2]*d[4][4] + g[m][1]*g[n][0]*d[4][3] + s*g[m][1]*g[n][1]*d[4][1] + s*g[m][2]*g[n][2]*d[2][4] + s*g[m][2]*g[n][0]*d[2][3] + 2.0*g[m][2]*g[n][1]*d[2][1]) / 2.0;
                kk[ii0+2+m*3][jj0+2+n*3] += c * (g[m][0]*g[n][0]*d[5][5] + g[m][0]*g[n][1]*d[5][4] + s*g[m][0]*g[n][2]*d[5][2] + g[m][1]*g[n][0]*d[4][5] + g[m][1]*g[n][1]*d[4][4] + s*g[m][1]*g[n][2]*d[4][2] + s*g[m][2]*g[n][0]*d[2][5] + s*g[m][2]*g[n][1]*d[2][4] + 2.0*g[m][2]*g[n][2]*d[2][2]) / 2.0;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::integ::testing::aux;
    use crate::integ::{self, AnalyticalTet4, AnalyticalTri3, IP_LIN_LEGENDRE_1, IP_TRI_INTERNAL_1};
    use crate::shapes::{GeoKind, Scratchpad};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::{copy_matrix, Matrix};
    use russell_tensor::LinElasticity;

    #[test]
    fn capture_some_errors() {
        let mut pad = aux::gen_pad_lin2(1.0);
        let mut kk = Matrix::new(4, 4);
        assert_eq!(
            integ::mat_10_gdg(&mut kk, &mut pad, 1, 0, false, &[], |_, _| Ok(())).err(),
            Some("nrow(K) must be ≥ ii0 + nnode ⋅ space_ndim")
        );
        assert_eq!(
            integ::mat_10_gdg(&mut kk, &mut pad, 0, 1, false, &[], |_, _| Ok(())).err(),
            Some("ncol(K) must be ≥ jj0 + nnode ⋅ space_ndim")
        );
        // more errors
        assert_eq!(
            integ::mat_10_gdg(&mut kk, &mut pad, 0, 0, false, &IP_LIN_LEGENDRE_1, |_, _| Ok(())).err(),
            Some("calc_gradient requires that geo_ndim = space_ndim")
        );
        let mut pad = aux::gen_pad_tri3();
        let mut kk = Matrix::new(6, 6);
        assert_eq!(
            integ::mat_10_gdg(&mut kk, &mut pad, 0, 0, false, &IP_TRI_INTERNAL_1, |_, _| Err("stop")).err(),
            Some("stop")
        );
    }

    #[test]
    fn mat_10_gdg_works_tri3_plane_stress() {
        // Element # 0 from example 1.6 from @bhatti page 32
        // Solid bracket with thickness = 0.25
        //              1     -10                connectivity:
        // y=2.0 (-100) o'-,__                    eid : vertices
        //              |     '-,__ 3   -10         0 :  0, 2, 3
        // y=1.5 - - -  |        ,'o-,__            1 :  3, 1, 0
        //              |  1   ,'  |    '-,__ 5     2 :  2, 4, 5
        //              |    ,'    |  3   ,-'o      3 :  5, 3, 2
        //              |  ,'  0   |   ,-'   |
        //              |,'        |,-'   2  |   constraints:
        // y=0.0 (-100) o----------o---------o    -100 : fixed on x and y
        //              0          2         4
        //             x=0.0     x=2.0     x=4.0
        // @bhatti Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.

        // scratchpad
        let mut pad = Scratchpad::new(2, GeoKind::Tri3).unwrap();
        pad.set_xx(0, 0, 0.0);
        pad.set_xx(0, 1, 0.0);
        pad.set_xx(1, 0, 2.0);
        pad.set_xx(1, 1, 0.0);
        pad.set_xx(2, 0, 2.0);
        pad.set_xx(2, 1, 1.5);

        // constants
        let young = 10_000.0;
        let poisson = 0.2;
        let th = 0.25; // thickness
        let plane_stress = true;
        let model = LinElasticity::new(young, poisson, true, plane_stress);

        // stiffness
        let class = pad.kind.class();
        let (space_ndim, nnode) = pad.xxt.dims();
        let nrow = nnode * space_ndim;
        let mut kk = Matrix::new(nrow, nrow);
        let ips = integ::points(class, 1).unwrap();
        integ::mat_10_gdg(&mut kk, &mut pad, 0, 0, true, ips, |dd, _| {
            let in_array = model.get_modulus().mat.as_data();
            let out_array = dd.mat.as_mut_data();
            for i in 0..in_array.len() {
                out_array[i] = th * in_array[i];
            }
            Ok(())
        })
        .unwrap();

        // compare against results from Bhatti's book
        #[rustfmt::skip]
        let kk_bhatti = Matrix::from( &[
            [  9.765625000000001e+02,  0.000000000000000e+00, -9.765625000000001e+02,  2.604166666666667e+02,  0.000000000000000e+00, -2.604166666666667e+02],
            [  0.000000000000000e+00,  3.906250000000000e+02,  5.208333333333334e+02, -3.906250000000000e+02, -5.208333333333334e+02,  0.000000000000000e+00],
            [ -9.765625000000001e+02,  5.208333333333334e+02,  1.671006944444445e+03, -7.812500000000000e+02, -6.944444444444445e+02,  2.604166666666667e+02],
            [  2.604166666666667e+02, -3.906250000000000e+02, -7.812500000000000e+02,  2.126736111111111e+03,  5.208333333333334e+02, -1.736111111111111e+03],
            [  0.000000000000000e+00, -5.208333333333334e+02, -6.944444444444445e+02,  5.208333333333334e+02,  6.944444444444445e+02,  0.000000000000000e+00],
            [ -2.604166666666667e+02,  0.000000000000000e+00,  2.604166666666667e+02, -1.736111111111111e+03,  0.000000000000000e+00,  1.736111111111111e+03],
        ]);
        assert_vec_approx_eq!(kk.as_data(), kk_bhatti.as_data(), 1e-12);

        // analytical solution
        let ana = AnalyticalTri3::new(&pad);
        let kk_correct = ana.integ_gdg(young, poisson, plane_stress, th).unwrap();

        // compare against analytical solution
        let tolerances = [1e-12, 1e-12, 1e-12, 1e-11, 1e-12];
        let selection: Vec<_> = [1, 3, 4, 12, 16]
            .iter()
            .map(|n| integ::points(class, *n).unwrap())
            .collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_10_gdg(&mut kk, &mut pad, 0, 0, true, ips, |dd, _| {
                let in_array = model.get_modulus().mat.as_data();
                let out_array = dd.mat.as_mut_data();
                for i in 0..in_array.len() {
                    out_array[i] = th * in_array[i];
                }
                Ok(())
            })
            .unwrap();
            assert_vec_approx_eq!(kk_correct.as_data(), kk.as_data(), tol); // 1e-12
        });
    }

    #[test]
    fn mat_10_gdg_works_tet4() {
        // scratchpad
        let mut pad = aux::gen_pad_tet4();

        // constants
        let young = 480.0;
        let poisson = 1.0 / 3.0;
        let model = LinElasticity::new(young, poisson, false, false);

        // analytical solution
        let mut ana = AnalyticalTet4::new(&pad);
        let kk_correct = ana.integ_gdg(young, poisson).unwrap();

        // check
        let class = pad.kind.class();
        let (space_ndim, nnode) = pad.xxt.dims();
        let nrow = nnode * space_ndim;
        let mut kk = Matrix::new(nrow, nrow);
        let tolerances = [1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12];
        let selection: Vec<_> = [1, 4, 5, 8, 14, 15, 24]
            .iter()
            .map(|n| integ::points(class, *n).unwrap())
            .collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            integ::mat_10_gdg(&mut kk, &mut pad, 0, 0, true, ips, |dd, _| {
                copy_matrix(&mut dd.mat, &model.get_modulus().mat)
            })
            .unwrap();
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol); //1e-12
        });
    }
}