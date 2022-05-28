use super::{allocate_shapes, allocate_states, Boundary, Find, Mesh};
use crate::shapes::{Shape, StateOfShape};
use crate::StrError;
use std::ffi::OsStr;

/// Holds all data related to a mesh including element shapes, states, boundaries, and functions to find entities
///
/// This struct is a (high-level) convenience that calls the necessary functions to generate all derived data from a mesh struct.
///
/// The Region basically calls the following functions:
///
/// ```text
/// let shapes = allocate_shapes(&mesh)?;
/// let states = allocate_states(&mesh, &shapes)?;
/// let boundary = Boundary::new(&mesh, &shapes)?;
/// let find = Find::new(&mesh, &boundary)?;
/// ```
pub struct Region {
    /// Holds the raw mesh data
    pub mesh: Mesh,

    /// Holds all shapes of all cells (len = **number of cells**)
    pub shapes: Vec<Shape>,

    /// Holds all states of shapes of all cells (len = **number of cells**)
    pub states: Vec<StateOfShape>,

    /// Holds the boundary data such as points, edges, and faces on boundary
    pub boundary: Boundary,

    /// Allows finding points, edges, and faces on the boundary by giving coordinates or keys
    pub find: Find,
}

impl Region {
    /// Allocates and prepares a new region with a given mesh
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::mesh::{At, Cell, Mesh, Point, Region};
    /// use gemlab::shapes::GeoKind;
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     //  3---------2---------5
    ///     //  |         |         |
    ///     //  |   [0]   |   [1]   |
    ///     //  |         |         |
    ///     //  0---------1---------4
    ///     let mesh = Mesh {
    ///         space_ndim: 2,
    ///         points: vec![
    ///             Point { id: 0, coords: vec![0.0, 0.0] },
    ///             Point { id: 1, coords: vec![1.0, 0.0] },
    ///             Point { id: 2, coords: vec![1.0, 1.0] },
    ///             Point { id: 3, coords: vec![0.0, 1.0] },
    ///             Point { id: 4, coords: vec![2.0, 0.0] },
    ///             Point { id: 5, coords: vec![2.0, 1.0] },
    ///         ],
    ///         cells: vec![
    ///             Cell { id: 0, attribute_id: 1, geo_ndim: 2, points: vec![0, 1, 2, 3] },
    ///             Cell { id: 1, attribute_id: 2, geo_ndim: 2, points: vec![1, 4, 5, 2] },
    ///         ],
    ///     };
    ///     let region = Region::with(mesh)?;
    ///     assert_eq!(region.mesh.space_ndim, 2);
    ///     assert_eq!(region.shapes.len(), 2);
    ///     assert_eq!(region.states.len(), 2);
    ///     assert_eq!(region.shapes[0].kind, GeoKind::Qua4);
    ///     assert_eq!(region.boundary.points.len(), 6);
    ///     assert_eq!(region.boundary.edges.len(), 6);
    ///     assert_eq!(region.boundary.faces.len(), 0);
    ///     assert_eq!(region.boundary.min, &[0.0, 0.0]);
    ///     assert_eq!(region.boundary.max, &[2.0, 1.0]);
    ///     assert_eq!(region.find.points(At::XY(0.0, 0.0))?.len(), 1);
    ///     Ok(())
    /// }
    /// ```
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

    /// Allocates and prepares a new region with a mesh defined in a text string
    #[inline]
    pub fn with_text(mesh_text: &str) -> Result<Self, StrError> {
        let mesh = Mesh::from_text(mesh_text)?;
        Region::with(mesh)
    }

    /// Allocates and prepares a new region with a mesh read from a text file
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    #[inline]
    pub fn with_text_file<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let mesh = Mesh::from_text_file(full_path)?;
        Region::with(mesh)
    }

    /// Allocates and prepares a new region with a mesh read from a binary file
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    #[inline]
    pub fn with_binary_file<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let mesh = Mesh::read(full_path)?;
        Region::with(mesh)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Mesh, Region};
    use crate::mesh::{At, Samples};
    use crate::StrError;

    #[test]
    fn with_works() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let region = Region::with(mesh)?;
        // println!("{:?}", mesh); // WRONG: mesh has been moved into region
        assert_eq!(region.mesh.space_ndim, 2);
        assert_eq!(region.shapes.len(), 2);
        assert_eq!(region.states.len(), 2);
        assert_eq!(region.boundary.points.len(), 6);
        assert_eq!(region.boundary.edges.len(), 6);
        assert_eq!(region.boundary.faces.len(), 0);
        assert_eq!(region.find.points(At::XY(0.0, 0.0))?.len(), 1);
        Ok(())
    }

    #[test]
    fn with_text_works() -> Result<(), StrError> {
        let mesh = Mesh::from_text(
            "# 1.0  3-----------2-----------5
             #      |           |           |
             #      |    [0]    |    [1]    |  [*] indicates id
             #      |    (1)    |    (2)    |  (*) indicates attribute_id
             #      |           |           |
             # 0.0  0-----------1-----------4
             #     0.0         1.0         2.0
             #
             # header
             # space_ndim npoint ncell
                     2      6     2

             # points
             # id    x   y
             0  0.0 0.0
             1  1.0 0.0
             2  1.0 1.0
             3  0.0 1.0
             4  2.0 0.0
             5  2.0 1.0

             # cells
             # id att geo_ndim nnode  point_ids...
             0   1        2     4  0 1 2 3
             1   2        2     4  1 4 5 2
             ",
        )?;
        let region = Region::with(mesh)?;
        assert_eq!(region.mesh.space_ndim, 2);
        assert_eq!(region.shapes.len(), 2);
        assert_eq!(region.states.len(), 2);
        assert_eq!(region.boundary.points.len(), 6);
        assert_eq!(region.boundary.edges.len(), 6);
        assert_eq!(region.boundary.faces.len(), 0);
        assert_eq!(region.find.points(At::XY(0.0, 0.0))?.len(), 1);
        Ok(())
    }

    #[test]
    fn with_text_file_works() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Mesh::from_text_file("./data/meshes/two_quads_horizontal.msh")?;
        let region = Region::with(mesh)?;
        assert_eq!(region.mesh.space_ndim, 2);
        assert_eq!(region.shapes.len(), 2);
        assert_eq!(region.states.len(), 2);
        assert_eq!(region.boundary.points.len(), 6);
        assert_eq!(region.boundary.edges.len(), 6);
        assert_eq!(region.boundary.faces.len(), 0);
        assert_eq!(region.find.points(At::XY(0.0, 0.0))?.len(), 1);
        Ok(())
    }

    #[test]
    fn with_binary_file_works() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let full_path = "/tmp/gemlab/test_region_two_quads_horizontal.dat";
        let mesh = Samples::two_quads_horizontal();
        mesh.write(full_path)?;
        let region = Region::with_binary_file(full_path)?;
        assert_eq!(region.mesh.space_ndim, 2);
        assert_eq!(region.shapes.len(), 2);
        assert_eq!(region.states.len(), 2);
        assert_eq!(region.boundary.points.len(), 6);
        assert_eq!(region.boundary.edges.len(), 6);
        assert_eq!(region.boundary.faces.len(), 0);
        assert_eq!(region.find.points(At::XY(0.0, 0.0))?.len(), 1);
        Ok(())
    }
}
