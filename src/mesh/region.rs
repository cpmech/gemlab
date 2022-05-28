use super::{allocate_shapes, allocate_states, Boundary, Find, Mesh};
use crate::shapes::{Shape, StateOfShape};
use crate::StrError;

/// Holds all data related to a Mesh
pub struct Region {
    pub mesh: Mesh,
    pub shapes: Vec<Shape>,
    pub states: Vec<StateOfShape>,
    pub boundary: Boundary,
    pub find: Find,
}

impl Region {
    /// Allocates a new instance with a given Mesh
    pub fn with(mesh: Mesh) -> Result<Self, StrError> {
        let shapes = allocate_shapes(&mesh)?;
        let states = allocate_states(&mesh, &shapes)?;
        let boundary = Boundary::new(&mesh, &shapes)?;
        let find = Find::new(&mesh, &boundary)?;
        Ok(Region {
            mesh,
            shapes,
            states,
            boundary,
            find,
        })
    }
}
