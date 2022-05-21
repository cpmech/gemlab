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

#[inline]
pub fn alloc_states_of_shapes(_mesh: &Mesh, shapes: &Vec<Shape>) -> Result<Vec<StateOfShape>, StrError> {
    shapes
        .iter()
        .map(|shape| {
            let state = StateOfShape::new(shape.space_ndim, shape.geo_ndim, shape.nnode);
            state
        })
        .collect()
}

// pub fn new_edge_shapes(boundaries: &Boundaries) {}
