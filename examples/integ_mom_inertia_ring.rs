use gemlab::StrError;
use gemlab::{integ, mesh, shapes, util};
use russell_chk::assert_approx_eq;

fn main() -> Result<(), StrError> {
    // generate mesh
    let (rmin, rmax) = (1.0, 3.0);
    let kind = shapes::GeoKind::Qua17;
    let mesh = mesh::Structured::quarter_ring_2d(rmin, rmax, 4, 8, kind)?;

    // allocate integration points and Scratchpad
    let ips = integ::default_points(kind);
    let mut pad = shapes::Scratchpad::new(2, kind)?;

    // sum contribution of all cells
    let mut second_mom_inertia = 0.0; // second moment of inertia about the x-axis
    for cell in &mesh.cells {
        // set the matrix of coordinates from cell points
        mesh::set_pad_coords(&mut pad, &cell.points, &mesh);

        // calculate the coordinates of the integration points
        let x_ips = integ::calc_ips_coords(&mut pad, ips)?;

        // perform the integration over the domain of a single cell
        second_mom_inertia += integ::scalar_field(&mut pad, ips, |p| {
            let y = x_ips[p][1];
            Ok(y * y)
        })?;
    }

    // multiply by to obtain the 2nd moment of the whole ring
    second_mom_inertia *= 4.0;

    // compare with analytical solution
    let correct = util::PI * (f64::powf(rmax, 4.0) - f64::powf(rmin, 4.0)) / 4.0;
    println!(
        "second_mom_inertia = {} ({}): err = {:.2e}",
        second_mom_inertia,
        correct,
        f64::abs(second_mom_inertia - correct),
    );
    assert_approx_eq!(second_mom_inertia, correct, 1e-7);

    // draw mesh
    mesh::draw_mesh(&mesh, false, "/tmp/gemlab/example_mom_inertia_ring.svg")
}
