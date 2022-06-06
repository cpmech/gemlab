use super::{Mesh, PointId};
use crate::shapes::{GeoKind, Shape, StateOfShape};
use crate::StrError;

/// Allocates all shapes corresponding to all cells in a mesh
///
/// Returns a vector with len = number of cells.
#[inline]
pub fn allocate_cell_shapes(mesh: &Mesh) -> Vec<Shape> {
    mesh.cells.iter().map(|cell| Shape::new(cell.kind)).collect()
}

/// Allocates StateOfShape for a given set of points
#[inline]
pub fn allocate_state(mesh: &Mesh, kind: GeoKind, points: &Vec<PointId>) -> Result<StateOfShape, StrError> {
    StateOfShape::new(
        kind,
        &points
            .iter()
            .map(|id| mesh.points[*id].coords.clone())
            .collect::<Vec<_>>(),
    )
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{allocate_cell_shapes, allocate_state};
    use crate::mesh::Samples;
    use crate::shapes::GeoKind;
    use crate::StrError;

    #[test]
    fn allocate_cell_shapes_works_2d() {
        let mesh = Samples::two_quads_horizontal();
        let shapes = allocate_cell_shapes(&mesh);
        assert_eq!(shapes.len(), 2);
        assert_eq!(shapes[0].kind, GeoKind::Qua4);
        assert_eq!(shapes[1].kind, GeoKind::Qua4);
        assert_eq!(shapes[0].nnode, 4);
        assert_eq!(shapes[1].nnode, 4);
    }

    #[test]
    fn allocate_cell_shapes_works_3d() {
        let mesh = Samples::two_cubes_vertical();
        let shapes = allocate_cell_shapes(&mesh);
        assert_eq!(shapes.len(), 2);
        assert_eq!(shapes[0].kind, GeoKind::Hex8);
        assert_eq!(shapes[1].kind, GeoKind::Hex8);
        assert_eq!(shapes[0].nnode, 8);
        assert_eq!(shapes[1].nnode, 8);
    }

    #[test]
    fn allocate_state_works_2d() -> Result<(), StrError> {
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        let mesh = Samples::two_quads_horizontal();
        let state = allocate_state(&mesh, mesh.cells[0].kind, &mesh.cells[0].points)?;
        assert_eq!(
            format!("{}", state.coords_transp),
            "┌         ┐\n\
             │ 0 1 1 0 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
        let state = allocate_state(&mesh, mesh.cells[1].kind, &mesh.cells[1].points)?;
        assert_eq!(
            format!("{}", state.coords_transp),
            "┌         ┐\n\
             │ 1 2 2 1 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
        Ok(())
    }

    #[test]
    fn allocate_state_works_3d() -> Result<(), StrError> {
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
        let state = allocate_state(&mesh, mesh.cells[0].kind, &mesh.cells[0].points)?;
        assert_eq!(
            format!("{}", state.coords_transp),
            "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 0 0 0 0 1 1 1 1 │\n\
             └                 ┘"
        );

        let state = allocate_state(&mesh, mesh.cells[1].kind, &mesh.cells[1].points)?;
        assert_eq!(
            format!("{}", state.coords_transp),
            "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 1 1 1 1 2 2 2 2 │\n\
             └                 ┘"
        );
        Ok(())
    }
}
