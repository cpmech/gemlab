#![allow(dead_code)]
#![allow(unused_imports)]

use super::{Boundary, Mesh};
use crate::util::{num_divisions, GridSearch};
use crate::StrError;

/// Minimum number of divisions for GridSearch
const GRID_SEARCH_NDIV_MIN: usize = 2;

/// Number of divisions for the longest direction in GridSearch
const GRID_SEARCH_NDIV_LONG: usize = 20;

/// Implements functions to find points, edges, and faces on the boundary of a mesh
pub struct Find {
    grid: GridSearch,
}

impl Find {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, boundary: &Boundary) -> Result<Self, StrError> {
        let ndiv = num_divisions(
            GRID_SEARCH_NDIV_MIN,
            GRID_SEARCH_NDIV_LONG,
            &boundary.min,
            &boundary.max,
        )?;
        let mut grid = GridSearch::new(mesh.space_ndim)?;
        grid.initialize(&ndiv, &boundary.min, &boundary.max)?;
        Ok(Find { grid })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {}
