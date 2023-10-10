use gemlab::shapes::{GeoKind, Scratchpad};
use gemlab::StrError;
use russell_lab::{vec_approx_eq, Vector};

fn main() -> Result<(), StrError> {
    //    3-------------2         ξ₀   ξ₁
    //    |      ξ₁     |  node    r    s
    //    |      |      |     0 -1.0 -1.0
    //  h |      +--ξ₀  |     1  1.0 -1.0
    //    |             |     2  1.0  1.0
    //    |             |     3 -1.0  1.0
    //    0-------------1
    // (x0,y0)   w

    // constants
    let (x0, y0) = (3.0, 4.0);
    let (w, h) = (2.0, 1.0);

    // scratchpad
    let space_ndim = 2;
    let mut pad = Scratchpad::new(space_ndim, GeoKind::Qua4)?;
    pad.set_xx(0, 0, x0);
    pad.set_xx(0, 1, y0);
    pad.set_xx(1, 0, x0 + w);
    pad.set_xx(1, 1, y0);
    pad.set_xx(2, 0, x0 + w);
    pad.set_xx(2, 1, y0 + h);
    pad.set_xx(3, 0, x0);
    pad.set_xx(3, 1, y0 + h);

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
    (pad.fn_interp)(&mut pad.interp, ksi_middle);

    // perform summation
    let nnode = pad.kind.nnode();
    let mut x_interpolated = Vector::new(space_ndim);
    for m in 0..nnode {
        for j in 0..space_ndim {
            x_interpolated[j] += pad.interp[m] * pad.xxt.get(j, m);
        }
    }

    // check
    let xm = x0 + w / 2.0;
    let ym = y0 + h / 2.0;
    println!("xm = {}, ym = {}", xm, ym);
    println!("x_interpolated =\n{}", x_interpolated);
    vec_approx_eq(x_interpolated.as_data(), &[xm, ym], 1e-15);
    Ok(())
}
