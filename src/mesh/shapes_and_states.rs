use super::{CellId, Mesh};
use crate::shapes::{Shape, StateOfShape};
use crate::StrError;

/// Allocates all shapes corresponding to all cells in a mesh
///
/// Returns a vector with len = number of cells.
#[inline]
pub fn allocate_shapes(mesh: &Mesh) -> Result<Vec<Shape>, StrError> {
    mesh.cells
        .iter()
        .map(|cell| Shape::new(mesh.space_ndim, cell.geo_ndim, cell.points.len()))
        .collect()
}

/// Allocates StateOfShape for a given cell_id
#[inline]
pub fn allocate_state(mesh: &Mesh, cell_id: CellId) -> Result<StateOfShape, StrError> {
    let cell = &mesh.cells[cell_id];
    StateOfShape::new(
        cell.geo_ndim,
        &cell
            .points
            .iter()
            .map(|id| mesh.points[*id].coords.clone())
            .collect::<Vec<_>>(),
    )
}

/// Allocates all states corresponding to all cells/shapes in a mesh
///
/// Returns a vector with len = number of cells.
#[inline]
pub fn allocate_states(mesh: &Mesh) -> Result<Vec<StateOfShape>, StrError> {
    mesh.cells
        .iter()
        .map(|cell| {
            StateOfShape::new(
                cell.geo_ndim,
                &cell
                    .points
                    .iter()
                    .map(|id| mesh.points[*id].coords.clone())
                    .collect::<Vec<_>>(),
            )
        })
        .collect()
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{allocate_shapes, allocate_state, allocate_states};
    use crate::mesh::Samples;
    use crate::shapes::GeoKind;
    use crate::StrError;

    #[test]
    fn allocate_shapes_works_2d() -> Result<(), StrError> {
        let mesh = Samples::two_quads_horizontal();
        let shapes = allocate_shapes(&mesh)?;
        assert_eq!(shapes.len(), 2);
        assert_eq!(shapes[0].kind, GeoKind::Qua4);
        assert_eq!(shapes[1].kind, GeoKind::Qua4);
        assert_eq!(shapes[0].nnode, 4);
        assert_eq!(shapes[1].nnode, 4);
        Ok(())
    }

    #[test]
    fn allocate_shapes_works_3d() -> Result<(), StrError> {
        let mesh = Samples::two_cubes_vertical();
        let shapes = allocate_shapes(&mesh)?;
        assert_eq!(shapes.len(), 2);
        assert_eq!(shapes[0].kind, GeoKind::Hex8);
        assert_eq!(shapes[1].kind, GeoKind::Hex8);
        assert_eq!(shapes[0].nnode, 8);
        assert_eq!(shapes[1].nnode, 8);
        Ok(())
    }

    #[test]
    fn allocate_state_works_2d() -> Result<(), StrError> {
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        let mesh = Samples::two_quads_horizontal();
        let state = allocate_state(&mesh, 1)?;
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
    fn allocate_states_works_2d() -> Result<(), StrError> {
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        let mesh = Samples::two_quads_horizontal();
        let states = allocate_states(&mesh)?;
        assert_eq!(states.len(), 2);
        assert_eq!(
            format!("{}", states[0].coords_transp),
            "┌         ┐\n\
             │ 0 1 1 0 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
        assert_eq!(
            format!("{}", states[1].coords_transp),
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
        let state = allocate_state(&mesh, 1)?;
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

    #[test]
    fn allocate_states_works_3d() -> Result<(), StrError> {
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
        let states = allocate_states(&mesh)?;
        assert_eq!(states.len(), 2);
        assert_eq!(
            format!("{}", states[0].coords_transp),
            "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 0 0 0 0 1 1 1 1 │\n\
             └                 ┘"
        );
        assert_eq!(
            format!("{}", states[1].coords_transp),
            "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 1 1 1 1 2 2 2 2 │\n\
             └                 ┘"
        );
        Ok(())
    }
}
