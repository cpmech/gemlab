use super::IntegPointData;
use crate::shapes::{Shape, StateOfShape};
use crate::util::SQRT_2;
use crate::StrError;
use russell_lab::Matrix;
use russell_tensor::Tensor4;

/// Implements the gradient(G) dot 4th-tensor(D) dot gradient(G) integration case (stiffness matrix)
///
/// Stiffness tensors:
///
/// ```text
///       ⌠               →    →
/// Kᵐⁿ = │ Gᵐₖ Dᵢₖⱼₗ Gⁿₗ eᵢ ⊗ eⱼ tₕ dΩ
/// ▔     ⌡
///       Ωₑ
/// ```
///
/// The numerical integration is:
///
/// ```text
///         nip-1     →         →       →          →
/// Kᵐⁿᵢⱼ ≈   Σ   Gᵐₖ(ιᵖ) Dᵢₖⱼₗ(ιᵖ) Gⁿₗ(ιᵖ) tₕ |J|(ιᵖ) wᵖ
///          p=0
/// ```
///
/// # Output
///
/// ```text
///     ┌                                               ┐
///     | K⁰⁰₀₀ K⁰⁰₀₁ K⁰¹₀₀ K⁰¹₀₁ K⁰²₀₀ K⁰²₀₁ ··· K⁰ⁿ₀ⱼ |
///     | K⁰⁰₁₀ K⁰⁰₁₁ K⁰¹₁₀ K⁰¹₁₁ K⁰²₁₀ K⁰²₁₁ ··· K⁰ⁿ₁ⱼ |
///     | K¹⁰₀₀ K¹⁰₀₁ K¹¹₀₀ K¹¹₀₁ K¹²₀₀ K¹²₀₁ ··· K¹ⁿ₀ⱼ |
/// K = | K¹⁰₁₀ K¹⁰₁₁ K¹¹₁₀ K¹¹₁₁ K¹²₁₀ K¹²₁₁ ··· K¹ⁿ₁ⱼ |
///     | K²⁰₀₀ K²⁰₀₁ K²¹₀₀ K²¹₀₁ K²²₀₀ K²²₀₁ ··· K²ⁿ₀ⱼ |
///     | K²⁰₁₀ K²⁰₁₁ K²¹₁₀ K²¹₁₁ K²²₁₀ K²²₁₁ ··· K²ⁿ₁ⱼ |
///     |  ···   ···   ···   ···   ···   ···  ···  ···  |
///     | Kᵐ⁰ᵢ₀ Kᵐ⁰ᵢ₁ Kᵐ¹ᵢ₀ Kᵐ¹ᵢ₁ Kᵐ²ᵢ₀ Kᵐ²ᵢ₁ ··· Kᵐⁿᵢⱼ |  ⟸  ii := i + m * space_ndim
///     └                                               ┘
///                                                 ⇑
///                                                 jj := j + n * space_ndim
///
/// m = ii / space_ndim    n = jj / space_ndim
/// i = ii % space_ndim    j = jj % space_ndim
/// ```
///
/// * `kk` -- A matrix containing all `Kᵐⁿᵢⱼ` values, one after another, and sequentially placed
///           as shown above (in 2D). `m` and `n` are the indices of the node and `i` and `j`
///           correspond to `space_ndim`. The dimension of `K` must be equal to
///           (`nnode * space_ndim`, `nnode * space_ndim`).
///
/// # Updated
///
/// * `state` -- Will be updated by the Shape functions
///
/// # Input
///
/// * `shape` -- Shape functions
/// * `ips` -- Integration points (n_integ_point)
/// * `th` -- tₕ the out-of-plane thickness in 2D or 1.0 otherwise (e.g., for plane-stress models)
/// * `erase_kk` -- Fills `kk` matrix with zeros, otherwise accumulate values into `kk`
/// * `fn_dd` -- Function f(D,p) corresponding to `D(x(ιᵖ))`
///
/// # Examples
///
/// ```
/// use gemlab::integ::{mat_gdg_stiffness, select_integ_points};
/// use gemlab::shapes::{GeoClass, Shape, StateOfShape};
/// use gemlab::StrError;
/// use russell_tensor::LinElasticity;
/// use russell_lab::{copy_matrix, Matrix};
///
/// fn main() -> Result<(), StrError> {
///     // shape and state
///     let space_ndim = 3;
///     let geo_ndim = 3;
///     let nnode = 4;
///     let shape = Shape::new(space_ndim, geo_ndim, nnode)?;
///     let mut state = StateOfShape::new(
///         shape.geo_ndim,
///         &[[2.0, 3.0, 4.0],
///           [6.0, 3.0, 2.0],
///           [2.0, 5.0, 1.0],
///           [4.0, 3.0, 6.0]],
///     )?;
///
///     // constants
///     let young = 480.0;
///     let poisson = 1.0 / 3.0;
///     let two_dim = false;
///     let plane_stress = false;
///     let model = LinElasticity::new(young, poisson, two_dim, plane_stress);
///
///     // stiffness
///     let nrow = shape.nnode * shape.space_ndim;
///     let mut kk = Matrix::new(nrow, nrow);
///     let ips = select_integ_points(GeoClass::Tet, 1)?;
///     mat_gdg_stiffness(&mut kk, &mut state, &shape, ips, 1.0, true, |dd, _| {
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
pub fn mat_gdg_stiffness<F>(
    kk: &mut Matrix,
    state: &mut StateOfShape,
    shape: &Shape,
    ips: IntegPointData,
    th: f64,
    erase_kk: bool,
    fn_dd: F,
) -> Result<(), StrError>
where
    F: Fn(&mut Tensor4, usize) -> Result<(), StrError>,
{
    // check
    let (nrow_kk, ncol_kk) = kk.dims();
    if nrow_kk != ncol_kk || nrow_kk != shape.nnode * shape.space_ndim {
        return Err("K.dims() must be equal to (nnode*space_ndim,nnode*space_ndim)");
    }

    // allocate auxiliary tensor
    let mut dd = Tensor4::new(true, shape.space_ndim == 2);

    // clear output matrix
    if erase_kk {
        kk.fill(0.0);
    }

    // loop over integration points
    for p in 0..ips.len() {
        // ksi coordinates and weight
        let iota = &ips[p];
        let weight = ips[p][3];

        // calculate Jacobian and Gradient
        let det_jac = shape.calc_gradient(state, iota)?;

        // calculate D tensor
        fn_dd(&mut dd, p)?;

        // add contribution to K matrix
        let c = det_jac * weight * th;
        mat_gdg_add_to_mat_kk(kk, &dd, c, shape, state);
    }
    Ok(())
}

/// Adds contribution to the K-matrix in integ_mat_10_gdg
#[inline]
#[rustfmt::skip]
fn mat_gdg_add_to_mat_kk(kk: &mut Matrix, dd: &Tensor4, c: f64, shape: &Shape, state: &mut StateOfShape) {
    let s = SQRT_2;
    let g = &state.gradient;
    let d = &dd.mat;
    if shape.space_ndim == 2 {
        for m in 0..shape.nnode {
            for n in 0..shape.nnode {
                kk[0+m*2][0+n*2] += c * (g[m][1]*g[n][1]*d[3][3] + s*g[m][1]*g[n][0]*d[3][0] + s*g[m][0]*g[n][1]*d[0][3] + 2.0*g[m][0]*g[n][0]*d[0][0]) / 2.0;
                kk[0+m*2][1+n*2] += c * (g[m][1]*g[n][0]*d[3][3] + s*g[m][1]*g[n][1]*d[3][1] + s*g[m][0]*g[n][0]*d[0][3] + 2.0*g[m][0]*g[n][1]*d[0][1]) / 2.0;
                kk[1+m*2][0+n*2] += c * (g[m][0]*g[n][1]*d[3][3] + s*g[m][0]*g[n][0]*d[3][0] + s*g[m][1]*g[n][1]*d[1][3] + 2.0*g[m][1]*g[n][0]*d[1][0]) / 2.0;
                kk[1+m*2][1+n*2] += c * (g[m][0]*g[n][0]*d[3][3] + s*g[m][0]*g[n][1]*d[3][1] + s*g[m][1]*g[n][0]*d[1][3] + 2.0*g[m][1]*g[n][1]*d[1][1]) / 2.0;
            }
        }
    } else {
        for m in 0..shape.nnode {
            for n in 0..shape.nnode {
                kk[0+m*3][0+n*3] += c * (g[m][2]*g[n][2]*d[5][5] + g[m][2]*g[n][1]*d[5][3] + s*g[m][2]*g[n][0]*d[5][0] + g[m][1]*g[n][2]*d[3][5] + g[m][1]*g[n][1]*d[3][3] + s*g[m][1]*g[n][0]*d[3][0] + s*g[m][0]*g[n][2]*d[0][5] + s*g[m][0]*g[n][1]*d[0][3] + 2.0*g[m][0]*g[n][0]*d[0][0]) / 2.0;
                kk[0+m*3][1+n*3] += c * (g[m][2]*g[n][2]*d[5][4] + g[m][2]*g[n][0]*d[5][3] + s*g[m][2]*g[n][1]*d[5][1] + g[m][1]*g[n][2]*d[3][4] + g[m][1]*g[n][0]*d[3][3] + s*g[m][1]*g[n][1]*d[3][1] + s*g[m][0]*g[n][2]*d[0][4] + s*g[m][0]*g[n][0]*d[0][3] + 2.0*g[m][0]*g[n][1]*d[0][1]) / 2.0;
                kk[0+m*3][2+n*3] += c * (g[m][2]*g[n][0]*d[5][5] + g[m][2]*g[n][1]*d[5][4] + s*g[m][2]*g[n][2]*d[5][2] + g[m][1]*g[n][0]*d[3][5] + g[m][1]*g[n][1]*d[3][4] + s*g[m][1]*g[n][2]*d[3][2] + s*g[m][0]*g[n][0]*d[0][5] + s*g[m][0]*g[n][1]*d[0][4] + 2.0*g[m][0]*g[n][2]*d[0][2]) / 2.0;
                kk[1+m*3][0+n*3] += c * (g[m][2]*g[n][2]*d[4][5] + g[m][2]*g[n][1]*d[4][3] + s*g[m][2]*g[n][0]*d[4][0] + g[m][0]*g[n][2]*d[3][5] + g[m][0]*g[n][1]*d[3][3] + s*g[m][0]*g[n][0]*d[3][0] + s*g[m][1]*g[n][2]*d[1][5] + s*g[m][1]*g[n][1]*d[1][3] + 2.0*g[m][1]*g[n][0]*d[1][0]) / 2.0;
                kk[1+m*3][1+n*3] += c * (g[m][2]*g[n][2]*d[4][4] + g[m][2]*g[n][0]*d[4][3] + s*g[m][2]*g[n][1]*d[4][1] + g[m][0]*g[n][2]*d[3][4] + g[m][0]*g[n][0]*d[3][3] + s*g[m][0]*g[n][1]*d[3][1] + s*g[m][1]*g[n][2]*d[1][4] + s*g[m][1]*g[n][0]*d[1][3] + 2.0*g[m][1]*g[n][1]*d[1][1]) / 2.0;
                kk[1+m*3][2+n*3] += c * (g[m][2]*g[n][0]*d[4][5] + g[m][2]*g[n][1]*d[4][4] + s*g[m][2]*g[n][2]*d[4][2] + g[m][0]*g[n][0]*d[3][5] + g[m][0]*g[n][1]*d[3][4] + s*g[m][0]*g[n][2]*d[3][2] + s*g[m][1]*g[n][0]*d[1][5] + s*g[m][1]*g[n][1]*d[1][4] + 2.0*g[m][1]*g[n][2]*d[1][2]) / 2.0;
                kk[2+m*3][0+n*3] += c * (g[m][0]*g[n][2]*d[5][5] + g[m][0]*g[n][1]*d[5][3] + s*g[m][0]*g[n][0]*d[5][0] + g[m][1]*g[n][2]*d[4][5] + g[m][1]*g[n][1]*d[4][3] + s*g[m][1]*g[n][0]*d[4][0] + s*g[m][2]*g[n][2]*d[2][5] + s*g[m][2]*g[n][1]*d[2][3] + 2.0*g[m][2]*g[n][0]*d[2][0]) / 2.0;
                kk[2+m*3][1+n*3] += c * (g[m][0]*g[n][2]*d[5][4] + g[m][0]*g[n][0]*d[5][3] + s*g[m][0]*g[n][1]*d[5][1] + g[m][1]*g[n][2]*d[4][4] + g[m][1]*g[n][0]*d[4][3] + s*g[m][1]*g[n][1]*d[4][1] + s*g[m][2]*g[n][2]*d[2][4] + s*g[m][2]*g[n][0]*d[2][3] + 2.0*g[m][2]*g[n][1]*d[2][1]) / 2.0;
                kk[2+m*3][2+n*3] += c * (g[m][0]*g[n][0]*d[5][5] + g[m][0]*g[n][1]*d[5][4] + s*g[m][0]*g[n][2]*d[5][2] + g[m][1]*g[n][0]*d[4][5] + g[m][1]*g[n][1]*d[4][4] + s*g[m][1]*g[n][2]*d[4][2] + s*g[m][2]*g[n][0]*d[2][5] + s*g[m][2]*g[n][1]*d[2][4] + 2.0*g[m][2]*g[n][2]*d[2][2]) / 2.0;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::mat_gdg_stiffness;
    use crate::integ::{select_integ_points, AnalyticalTet4, AnalyticalTri3};
    use crate::shapes::{GeoClass, Shape, StateOfShape};
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::{copy_matrix, Matrix};
    use russell_tensor::LinElasticity;

    #[test]
    fn capture_some_errors() {
        let shape = Shape::new(2, 1, 2).unwrap();
        let mut state = StateOfShape::new(shape.geo_ndim, &[[0.0, 0.0], [1.0, 0.0]]).unwrap();
        let mut kk = Matrix::new(2, 2);
        assert_eq!(
            mat_gdg_stiffness(&mut kk, &mut state, &shape, &[], 1.0, false, |_, _| Ok(())).err(),
            Some("K.dims() must be equal to (nnode*space_ndim,nnode*space_ndim)")
        );
    }

    #[test]
    fn mat_gdg_stiffness_works_tri3_plane_stress() -> Result<(), StrError> {
        // Element # 0 from example 1.6 from [@bhatti] page 32
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
        // [@bhatti] Bhatti, M.A. (2005) Fundamental Finite Element Analysis
        //           and Applications, Wiley, 700p.

        // shape and state
        let shape = Shape::new(2, 2, 3)?;
        let mut state = StateOfShape::new(shape.geo_ndim, &[[0.0, 0.0], [2.0, 0.0], [2.0, 1.5]])?;

        // constants
        let young = 10_000.0;
        let poisson = 0.2;
        let th = 0.25; // thickness
        let plane_stress = true;
        let model = LinElasticity::new(young, poisson, true, plane_stress);

        // stiffness
        let nrow = shape.nnode * shape.space_ndim;
        let mut kk = Matrix::new(nrow, nrow);
        let ips = select_integ_points(GeoClass::Tri, 1)?;
        mat_gdg_stiffness(&mut kk, &mut state, &shape, ips, th, true, |dd, _| {
            copy_matrix(&mut dd.mat, &model.get_modulus().mat)
        })?;

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
        let ana = AnalyticalTri3::new(&shape, &mut state);
        let kk_correct = ana.integ_stiffness(young, poisson, plane_stress, th)?;

        // compare against analytical solution
        let tolerances = [1e-12, 1e-12, 1e-12, 1e-12, 1e-11, 1e-12];
        let selection: Vec<_> = [1, 3, 1_003, 4, 12, 16]
            .iter()
            .map(|n| select_integ_points(GeoClass::Tri, *n).unwrap())
            .collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            mat_gdg_stiffness(&mut kk, &mut state, &shape, ips, th, true, |dd, _| {
                copy_matrix(&mut dd.mat, &model.get_modulus().mat)
            })
            .unwrap();
            assert_vec_approx_eq!(kk_correct.as_data(), kk.as_data(), tol); // 1e-12
        });
        Ok(())
    }

    #[test]
    fn mat_gdg_stiffness_works_tet4() -> Result<(), StrError> {
        // shape and state
        let shape = Shape::new(3, 3, 4)?;
        let mut state = StateOfShape::new(
            shape.geo_ndim,
            &[[2.0, 3.0, 4.0], [6.0, 3.0, 2.0], [2.0, 5.0, 1.0], [4.0, 3.0, 6.0]],
        )?;

        // constants
        let young = 480.0;
        let poisson = 1.0 / 3.0;
        let model = LinElasticity::new(young, poisson, false, false);

        // analytical solution
        let mut ana = AnalyticalTet4::new(&shape, &state);
        let kk_correct = ana.integ_stiffness(young, poisson)?;

        // check
        let nrow = shape.nnode * shape.space_ndim;
        let mut kk = Matrix::new(nrow, nrow);
        let tolerances = [1e-12, 1e-12, 1e-12, 1e-12, 1e-12];
        let selection: Vec<_> = [1, 4, 5, 8, 14]
            .iter()
            .map(|n| select_integ_points(GeoClass::Tet, *n).unwrap())
            .collect();
        selection.iter().zip(tolerances).for_each(|(ips, tol)| {
            // println!("nip={}, tol={:.e}", ips.len(), tol);
            mat_gdg_stiffness(&mut kk, &mut state, &shape, ips, 1.0, true, |dd, _| {
                copy_matrix(&mut dd.mat, &model.get_modulus().mat)
            })
            .unwrap();
            assert_vec_approx_eq!(kk.as_data(), kk_correct.as_data(), tol); //1e-12
        });
        Ok(())
    }
}
