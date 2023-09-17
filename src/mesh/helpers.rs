use super::{Mesh, PointId};
use crate::shapes::Scratchpad;

/// Sets the pad's matrix of coordinates X given a list of point ids
///
/// # Panics
///
/// 1. Make sure `pad.kind.nnode() == points.len()`; otherwise a panic will occur
/// 2. This function does not check for bounds on point indices and dimensions
/// 3. Use [Mesh::check_all] to capture (some) errors
#[inline]
pub fn set_pad_coords(pad: &mut Scratchpad, points: &[PointId], mesh: &Mesh) {
    let nnode = pad.kind.nnode();
    assert_eq!(nnode, points.len());
    for m in 0..nnode {
        for j in 0..mesh.ndim {
            pad.set_xx(m, j, mesh.points[points[m]].coords[j]);
        }
    }
}

/// Returns the (min,max) point coordinates in a mesh
#[inline]
pub fn get_mesh_limits(mesh: &Mesh) -> (Vec<f64>, Vec<f64>) {
    let mut min = vec![f64::MAX; mesh.ndim];
    let mut max = vec![f64::MIN; mesh.ndim];
    for point in &mesh.points {
        for i in 0..mesh.ndim {
            min[i] = f64::min(min[i], point.coords[i]);
            max[i] = f64::max(max[i], point.coords[i]);
        }
    }
    (min, max)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{get_mesh_limits, set_pad_coords};
    use crate::mesh::Samples;
    use crate::shapes::Scratchpad;

    #[test]
    fn set_pad_coords_works_2d() {
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        let mesh = Samples::two_qua4();
        let cell = &mesh.cells[0];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        set_pad_coords(&mut pad, &cell.points, &mesh);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌         ┐\n\
             │ 0 1 1 0 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
        let cell = &mesh.cells[1];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        set_pad_coords(&mut pad, &cell.points, &mesh);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌         ┐\n\
             │ 1 2 2 1 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
    }

    #[test]
    fn set_pad_coords_works_3d() {
        //       8-------------11
        //      /.             /|
        //     / .            / |
        //    /  .           /  |
        //   /   .          /   |
        //  9-------------10    |
        //  |    .         |    |
        //  |    4---------|----7
        //  |   /.         |   /|
        //  |  / .         |  / |
        //  | /  .         | /  |
        //  |/   .         |/   |
        //  5--------------6    |
        //  |    .         |    |
        //  |    0---------|----3
        //  |   /          |   /
        //  |  /           |  /
        //  | /            | /
        //  |/             |/
        //  1--------------2
        let mesh = Samples::two_hex8();
        let cell = &mesh.cells[0];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        set_pad_coords(&mut pad, &cell.points, &mesh);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 0 0 0 0 1 1 1 1 │\n\
             └                 ┘"
        );
        let cell = &mesh.cells[1];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        set_pad_coords(&mut pad, &cell.points, &mesh);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 1 1 1 1 2 2 2 2 │\n\
             └                 ┘"
        );
    }

    #[test]
    fn get_mesh_limits_works() {
        let mesh = &Samples::two_qua4();
        let (min, max) = get_mesh_limits(&mesh);
        assert_eq!(min, &[0.0, 0.0]);
        assert_eq!(max, &[2.0, 1.0]);

        let mesh = &Samples::two_hex8();
        let (min, max) = get_mesh_limits(&mesh);
        assert_eq!(min, &[0.0, 0.0, 0.0]);
        assert_eq!(max, &[1.0, 1.0, 2.0]);
    }
}
