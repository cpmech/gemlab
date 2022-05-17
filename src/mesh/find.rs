#![allow(dead_code)]
#![allow(unused_imports)]

use super::{Boundary, Mesh};
use crate::util::GridSearch;
use crate::StrError;

/// Number of divisions along the longest direction for GridSearch
///
/// ```text
/// ndiv_other = min(ndiv_max, max(ndiv_min, (ll_other/ll_long) * ndiv_long))
/// ```
const GRID_SEARCH_NDIV_LONG: usize = 20;

/// Minimum number of divisions for GridSearch
const GRID_SEARCH_NDIV_MIN: usize = 2;

/// Maximum number of divisions for GridSearch
const GRID_SEARCH_NDIV_MAX: usize = 50;

pub struct Find {
    grid: GridSearch,

    pub min: Vec<f64>,   // (space_ndim)
    pub max: Vec<f64>,   // (space_ndim)
    pub delta: Vec<f64>, // (space_ndim)
}

impl Find {
    // pub fn new(mesh: &Mesh, boundary: &Boundary) -> Result<Self, StrError> {
    //     Ok(Find { grid })
    // }
}
