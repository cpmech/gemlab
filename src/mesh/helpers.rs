use super::{Mesh, PointId};
use crate::shapes::Scratchpad;

/// Sets the matrix of coordinates (transposed) Xᵀ given a list of point ids
///
/// # Panics
///
/// 1. Make sure `pad.kind.nnode() == points.len()`; otherwise a panic will occur
/// 2. This function does not check for bounds on point indices and dimensions
/// 3. Use [crate::mesh::check_all()] to capture (some) errors
#[inline]
pub fn set_xxt_from_points(pad: &mut Scratchpad, points: &Vec<PointId>, mesh: &Mesh) {
    let nnode = pad.kind.nnode();
    assert_eq!(nnode, points.len());
    for m in 0..nnode {
        for j in 0..mesh.ndim {
            pad.set_xx(m, j, mesh.points[points[m]].coords[j]);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::set_xxt_from_points;
    use crate::mesh::Samples;
    use crate::shapes::Scratchpad;
    use crate::StrError;

    #[test]
    fn set_xxt_from_points_works_2d() -> Result<(), StrError> {
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        let mesh = Samples::two_quads_horizontal();
        let cell = &mesh.cells[0];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind)?;
        set_xxt_from_points(&mut pad, &cell.points, &mesh);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌         ┐\n\
             │ 0 1 1 0 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
        let cell = &mesh.cells[1];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind)?;
        set_xxt_from_points(&mut pad, &cell.points, &mesh);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌         ┐\n\
             │ 1 2 2 1 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
        Ok(())
    }

    #[test]
    fn set_xxt_from_points_works_3d() -> Result<(), StrError> {
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
        let mesh = Samples::two_cubes_vertical();
        let cell = &mesh.cells[0];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind)?;
        set_xxt_from_points(&mut pad, &cell.points, &mesh);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 0 0 0 0 1 1 1 1 │\n\
             └                 ┘"
        );
        let cell = &mesh.cells[1];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind)?;
        set_xxt_from_points(&mut pad, &cell.points, &mesh);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 1 1 1 1 2 2 2 2 │\n\
             └                 ┘"
        );
        Ok(())
    }
}
