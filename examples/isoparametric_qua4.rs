use gemlab::shapes::{GeoKind, Shape, StateOfShape};
use gemlab::StrError;
use russell_chk::assert_vec_approx_eq;
use russell_lab::Vector;

fn main() -> Result<(), StrError> {
    //    3-----------2         ξ₀   ξ₁
    //    |     ξ₁    |  node    r    s
    //    |     |     |     0 -1.0 -1.0
    //  h |     +--ξ₀ |     1  1.0 -1.0
    //    |           |     2  1.0  1.0
    //    |           |     3 -1.0  1.0
    //    0-----------1
    // (x0,y0)   w

    // constants
    let (x0, y0) = (3.0, 4.0);
    let (w, h) = (2.0, 1.0);

    // coordinates
    #[rustfmt::skip]
    let coords = &[
        [x0,     y0    ],
        [x0 + w, y0    ],
        [x0 + w, y0 + h],
        [x0,     y0 + h],
    ];

    // shape and state
    let shape = Shape::new(GeoKind::Qua4);
    let mut state = StateOfShape::new(shape.kind, coords)?;

    // perform interpolation
    //
    // Any coordinate within the element is calculated
    // by the following "isoparametric" formula:
    //
    // → →         →  →
    // x(ξ) = Σ Nᵐ(ξ) xᵐ
    //        m
    //
    // Let us calculate the coordinates, at the middle of the element
    // with ξ = [0.0, 0.0]

    // compute interpolation functions @ ksi_middle
    let ksi_middle = &[0.0, 0.0];
    shape.calc_interp(&mut state, ksi_middle)?;

    // perform summation
    let space_ndim = 2;
    let mut x_interpolated = Vector::new(space_ndim);
    for m in 0..shape.nnode {
        for i in 0..space_ndim {
            x_interpolated[i] += state.interp[m] * state.coords_transp[i][m];
        }
    }

    // check
    let xm = x0 + w / 2.0;
    let ym = y0 + h / 2.0;
    println!("xm = {}, ym = {}", xm, ym);
    println!("x_interpolated =\n{}", x_interpolated);
    assert_vec_approx_eq!(x_interpolated.as_data(), &[xm, ym], 1e-15);
    Ok(())
}
