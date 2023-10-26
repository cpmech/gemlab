use gemlab::integ::{default_points, points_coords, scalar_field};
use gemlab::prelude::*;
use gemlab::StrError;
use russell_lab::approx_eq;
use russell_lab::math::PI;

fn main() -> Result<(), StrError> {
    // generate mesh
    let r = 3.0;
    let kind = GeoKind::Qua17;
    let mesh_1 = Structured::quarter_disk_2d_a(r, 3, 3, kind, false)?;
    let mesh_2 = Structured::quarter_disk_2d_b(r / 2.0, r, 3, 3, kind, false)?;

    // allocate integration points and Scratchpad
    let ips = default_points(kind);
    let mut pad = Scratchpad::new(2, kind)?;

    // mesh 1: sum contribution of all cells
    let mut second_mom_inertia_mesh_1 = 0.0; // second moment of inertia about the x-axis
    for cell in &mesh_1.cells {
        // set the matrix of coordinates from cell points
        mesh_1.set_pad(&mut pad, &cell.points);

        // calculate the coordinates of the integration points
        let x_ips = points_coords(&mut pad, ips)?;

        // perform the integration over the domain of a single cell
        second_mom_inertia_mesh_1 += scalar_field(&mut pad, ips, |p| {
            let y = x_ips[p][1];
            Ok(y * y)
        })?;
    }
    second_mom_inertia_mesh_1 *= 4.0; // multiply by to obtain the 2nd moment of the whole disk

    // mesh 2: sum contribution of all cells
    let mut second_mom_inertia_mesh_2 = 0.0; // second moment of inertia about the x-axis
    for cell in &mesh_2.cells {
        // set the matrix of coordinates from cell points
        mesh_2.set_pad(&mut pad, &cell.points);

        // calculate the coordinates of the integration points
        let x_ips = points_coords(&mut pad, ips)?;

        // perform the integration over the domain of a single cell
        second_mom_inertia_mesh_2 += scalar_field(&mut pad, ips, |p| {
            let y = x_ips[p][1];
            Ok(y * y)
        })?;
    }
    second_mom_inertia_mesh_2 *= 4.0; // multiply by to obtain the 2nd moment of the whole disk

    // compare with analytical solution
    let correct = f64::powf(r, 4.0) * PI / 4.0; // analytical solution
    println!(
        "mesh 1: second_mom_inertia = {} ({}): err = {:.2e}",
        second_mom_inertia_mesh_1,
        correct,
        f64::abs(second_mom_inertia_mesh_1 - correct),
    );
    println!(
        "mesh 2: second_mom_inertia = {} ({}): err = {:.2e}",
        second_mom_inertia_mesh_2,
        correct,
        f64::abs(second_mom_inertia_mesh_2 - correct),
    );
    approx_eq(second_mom_inertia_mesh_1, correct, 1e-5);
    approx_eq(second_mom_inertia_mesh_2, correct, 1e-5);

    // draw meshes
    let mut fig_1 = Figure::new();
    fig_1.cell_ids = true;
    fig_1.point_ids = true;
    fig_1.figure_size = Some((800.0, 800.0));
    mesh_1.draw(Some(fig_1), "/tmp/gemlab/example_mom_inertia_disk_1.svg", |_, _| {})?;

    let mut fig_2 = Figure::new();
    fig_2.cell_ids = true;
    fig_2.point_ids = true;
    fig_2.figure_size = Some((800.0, 800.0));
    mesh_2.draw(Some(fig_2), "/tmp/gemlab/example_mom_inertia_disk_2.svg", |_, _| {})
}
