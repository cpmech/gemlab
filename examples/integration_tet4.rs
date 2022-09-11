use gemlab::integ;
use gemlab::shapes::{GeoKind, Scratchpad};
use gemlab::StrError;
use russell_lab::Matrix;
use russell_tensor::{copy_tensor4, LinElasticity};

fn main() -> Result<(), StrError> {
    // scratchpad
    let space_ndim = 3;
    let mut pad = Scratchpad::new(space_ndim, GeoKind::Tet4)?;
    pad.set_xx(0, 0, 2.0);
    pad.set_xx(0, 1, 3.0);
    pad.set_xx(0, 2, 4.0);
    pad.set_xx(1, 0, 6.0);
    pad.set_xx(1, 1, 3.0);
    pad.set_xx(1, 2, 2.0);
    pad.set_xx(2, 0, 2.0);
    pad.set_xx(2, 1, 5.0);
    pad.set_xx(2, 2, 1.0);
    pad.set_xx(3, 0, 4.0);
    pad.set_xx(3, 1, 3.0);
    pad.set_xx(3, 2, 6.0);

    // constants
    let young = 96.0;
    let poisson = 1.0 / 3.0;
    let two_dim = false;
    let plane_stress = false;
    let model = LinElasticity::new(young, poisson, two_dim, plane_stress);

    // stiffness
    let nnode = pad.kind.nnode();
    let nrow = nnode * space_ndim;
    let mut kk = Matrix::new(nrow, nrow);
    let ips = integ::default_points(pad.kind);
    let mut args = integ::CommonArgs::new(&mut pad, ips);
    integ::mat_10_bdb(&mut kk, &mut args, |dd, _, _, _| {
        copy_tensor4(dd, model.get_modulus()).unwrap();
        Ok(())
    })?;

    // output
    assert_eq!(
        format!("{:.0}", kk),
        "┌                                                             ┐\n\
         │  149  108   24   -1    6   12  -54  -48    0  -94  -66  -36 │\n\
         │  108  344   54  -24  104   42  -24 -216  -12  -60 -232  -84 │\n\
         │   24   54  113    0   30   35    0  -24  -54  -24  -60  -94 │\n\
         │   -1  -24    0   29  -18  -12  -18   24    0  -10   18   12 │\n\
         │    6  104   30  -18   44   18   12  -72  -12    0  -76  -36 │\n\
         │   12   42   35  -12   18   29    0  -24  -18    0  -36  -46 │\n\
         │  -54  -24    0  -18   12    0   36    0    0   36   12    0 │\n\
         │  -48 -216  -24   24  -72  -24    0  144    0   24  144   48 │\n\
         │    0  -12  -54    0  -12  -18    0    0   36    0   24   36 │\n\
         │  -94  -60  -24  -10    0    0   36   24    0   68   36   24 │\n\
         │  -66 -232  -60   18  -76  -36   12  144   24   36  164   72 │\n\
         │  -36  -84  -94   12  -36  -46    0   48   36   24   72  104 │\n\
         └                                                             ┘"
    );
    Ok(())
}
