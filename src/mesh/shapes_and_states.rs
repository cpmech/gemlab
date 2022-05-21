use super::Mesh;
use crate::shapes::{Shape, StateOfShape};
use crate::StrError;

/// Assists in allocating all shapes corresponding to all cells in a mesh
pub struct Shapes;

/// Assists in allocating all states corresponding to all shapes in a mesh
pub struct States;

impl Shapes {
    /// Allocates a new instance
    #[inline]
    pub fn new(mesh: &Mesh) -> Result<Vec<Shape>, StrError> {
        mesh.cells
            .iter()
            .map(|cell| Shape::new(mesh.space_ndim, cell.geo_ndim, cell.points.len()))
            .collect()
    }
}

impl States {
    /// Allocate a new instance
    #[inline]
    pub fn new(mesh: &Mesh, shapes: &Vec<Shape>) -> Result<Vec<StateOfShape>, StrError> {
        mesh.cells
            .iter()
            .zip(shapes)
            .map(|(cell, shape)| {
                StateOfShape::new(
                    shape.geo_ndim,
                    &cell
                        .points
                        .iter()
                        .map(|id| mesh.points[*id].coords.clone())
                        .collect::<Vec<_>>(),
                )
            })
            .collect()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Shapes, States};
    use crate::mesh::Samples;
    use crate::shapes::GeoKind;
    use crate::StrError;

    #[test]
    fn shapes_new_works() -> Result<(), StrError> {
        let mesh = Samples::two_quads_horizontal();
        let shapes = Shapes::new(&mesh)?;
        assert_eq!(shapes.len(), 2);
        assert_eq!(shapes[0].kind, GeoKind::Qua4);
        assert_eq!(shapes[1].kind, GeoKind::Qua4);
        assert_eq!(shapes[0].nnode, 4);
        assert_eq!(shapes[1].nnode, 4);
        Ok(())
    }

    #[test]
    fn states_new_works() -> Result<(), StrError> {
        let mesh = Samples::two_quads_horizontal();
        let shapes = Shapes::new(&mesh)?;
        let states = States::new(&mesh, &shapes)?;
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
}
