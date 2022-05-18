#![allow(dead_code)]
#![allow(unused_imports)]

use super::{Boundary, Mesh};
use crate::util::{GridSearch, GsNdiv, GsTol};
use crate::StrError;

/// Implements functions to find points, edges, and faces on the boundary of a mesh
pub struct Find {
    grid: GridSearch,
}

impl Find {
    /// Allocates a new instance
    pub fn new(_mesh: &Mesh, boundary: &Boundary) -> Result<Self, StrError> {
        let grid = GridSearch::new(&boundary.min, &boundary.max, GsNdiv::Default, GsTol::Default)?;
        Ok(Find { grid })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::StrError;

    #[test]
    fn new_works() -> Result<(), StrError> {
        Ok(())
    }
}
