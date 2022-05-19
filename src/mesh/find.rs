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
    pub fn new(mesh: &Mesh, boundary: &Boundary) -> Result<Self, StrError> {
        let mut grid = GridSearch::new(&boundary.min, &boundary.max, GsNdiv::Default, GsTol::Default)?;
        for id in &boundary.points {
            grid.insert(*id, &mesh.points[*id].coords)?;
        }
        Ok(Find { grid })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Find;
    use crate::mesh::{Boundary, Samples};
    use crate::StrError;

    #[test]
    fn new_works() -> Result<(), StrError> {
        let mesh = Samples::two_quads_horizontal();
        let boundary = Boundary::new(&mesh)?;
        let find = Find::new(&mesh, &boundary)?;
        // // do not delete this //
        // let mut plot = find.grid.plot()?;
        // plot.set_equal_axes(true).set_figure_size_points(800.0, 400.0);
        // plot.save("/tmp/gemlab/find_new_works.svg")?;
        assert_eq!(
            format!("{}", find.grid),
            "0: [0]\n\
             9: [1]\n\
             10: [1]\n\
             19: [4]\n\
             180: [3]\n\
             189: [2]\n\
             190: [2]\n\
             199: [5]\n\
             ids = [0, 1, 2, 3, 4, 5]\n\
             nitem = 6\n\
             ncontainer = 8\n\
             ndiv = [20, 10]\n"
        );
        Ok(())
    }
}
