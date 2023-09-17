use gemlab::integ::{default_points, points_coords, scalar_field};
use gemlab::prelude::*;
use gemlab::StrError;
use russell_chk::approx_eq;
use russell_lab::math::PI;

fn main() -> Result<(), StrError> {
    // generate mesh
    let (rmin, rmax) = (1.0, 3.0);
    let kind = GeoKind::Qua17;
    let mesh = Structured::quarter_ring_2d(rmin, rmax, 4, 8, kind)?;

    // allocate integration points and Scratchpad
    let ips = default_points(kind);
    let mut pad = Scratchpad::new(2, kind)?;

    // sum contribution of all cells
    let mut second_mom_inertia = 0.0; // second moment of inertia about the x-axis
    for cell in &mesh.cells {
        // set the matrix of coordinates from cell points
        mesh.set_pad(&mut pad, &cell.points);

        // calculate the coordinates of the integration points
        let x_ips = points_coords(&mut pad, ips)?;

        // perform the integration over the domain of a single cell
        second_mom_inertia += scalar_field(&mut pad, ips, |p| {
            let y = x_ips[p][1];
            Ok(y * y)
        })?;
    }

    // multiply by to obtain the 2nd moment of the whole ring
    second_mom_inertia *= 4.0;

    // compare with analytical solution
    let correct = PI * (f64::powf(rmax, 4.0) - f64::powf(rmin, 4.0)) / 4.0;
    println!(
        "second_mom_inertia = {} ({}): err = {:.2e}",
        second_mom_inertia,
        correct,
        f64::abs(second_mom_inertia - correct),
    );
    approx_eq(second_mom_inertia, correct, 1e-7);

    // draw mesh
    let mut fig = Figure::new();
    fig.cell_ids = true;
    fig.point_ids = true;
    fig.figure_size = Some((800.0, 800.0));
    mesh.draw(Some(fig), "/tmp/gemlab/example_mom_inertia_ring.svg", |_, _| {})
}
