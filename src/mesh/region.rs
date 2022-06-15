use super::{Extract, Features, Find, MapEdge2dToCells, MapFaceToCells, Mesh};
use crate::StrError;
use std::ffi::OsStr;

/// Holds all (immutable) data related to a mesh including element shapes, boundary, interior, and functions to find features
pub struct Region {
    /// Holds the raw mesh data
    pub mesh: Mesh,

    /// Maps all edge keys to cells sharing the edge (2D only)
    pub all_2d_edges: Option<MapEdge2dToCells>,

    /// Maps all face keys to cells sharing the face (3D only)
    pub all_faces: Option<MapFaceToCells>,

    /// Holds points, edges and faces on the mesh boundary, interior, or both
    ///
    /// Depends on the Extract option.
    pub features: Features,

    /// Finds features on the boundary, interior, or both
    ///
    /// Depends on the Extract option and corresponds to `features`.
    pub find: Find,
}

impl Region {
    /// Allocates and prepares a new region with a given mesh
    ///
    /// # Input
    ///
    /// * `mesh` -- the mesh (will move to Region)
    /// * `extract` -- which features to extract?
    ///
    /// # Panics
    ///
    /// 1. This function may panic if the mesh data is inconsistent
    /// 2. You may want to call [crate::mesh::check_all()] to capture (some) errors
    pub fn with(mesh: Mesh, extract: Extract) -> Result<Self, StrError> {
        let (all_2d_edges, all_faces, features) = Features::new(&mesh, extract);
        let find = Find::new(&mesh, &features)?;
        Ok(Region {
            mesh,
            all_2d_edges,
            all_faces,
            features,
            find,
        })
    }

    /// Allocates and prepares a new region with a mesh defined in a text string
    ///
    /// # Input
    ///
    /// * `mesh_text` -- text representing the mesh
    /// * `extract` -- which features to extract?
    #[inline]
    pub fn with_text(mesh_text: &str, extract: Extract) -> Result<Self, StrError> {
        let mesh = Mesh::from_text(mesh_text)?;
        Region::with(mesh, extract)
    }

    /// Allocates and prepares a new region with a mesh read from a text file
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    /// * `extract` -- which features to extract?
    #[inline]
    pub fn with_text_file<P>(full_path: &P, extract: Extract) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let mesh = Mesh::from_text_file(full_path)?;
        Region::with(mesh, extract)
    }

    /// Allocates and prepares a new region with a mesh read from a binary file
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    /// * `extract` -- which features to extract?
    #[inline]
    pub fn with_binary_file<P>(full_path: &P, extract: Extract) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let mesh = Mesh::read(full_path)?;
        Region::with(mesh, extract)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Mesh, Region};
    use crate::mesh::{At, Extract, Samples};
    use crate::StrError;

    #[test]
    fn with_works() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let region = Region::with(mesh, Extract::Boundary)?;
        // println!("{:?}", mesh); // WRONG: mesh has been moved into region
        assert_eq!(region.mesh.ndim, 2);
        assert_eq!(region.features.points.len(), 6);
        assert_eq!(region.features.edges.len(), 6);
        assert_eq!(region.features.faces.len(), 0);
        assert_eq!(region.find.points(At::XY(0.0, 0.0))?.len(), 1);
        Ok(())
    }

    #[test]
    fn with_text_works() -> Result<(), StrError> {
        let region = Region::with_text(
            "# 1.0  3-----------2-----------5
             #      |           |           |
             #      |    [0]    |    [1]    |  [*] indicates id
             #      |    (1)    |    (2)    |  (*) indicates attribute_id
             #      |           |           |
             # 0.0  0-----------1-----------4
             #     0.0         1.0         2.0
             #
             # header
             # ndim npoint ncell
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
             # id att kind  point_ids...
                0   1 qua4  0 1 2 3
                1   2 qua4  1 4 5 2
             ",
            Extract::Boundary,
        )?;
        assert_eq!(region.mesh.ndim, 2);
        assert_eq!(region.features.points.len(), 6);
        assert_eq!(region.features.edges.len(), 6);
        assert_eq!(region.features.faces.len(), 0);
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
        let region = Region::with(mesh, Extract::Boundary)?;
        assert_eq!(region.mesh.ndim, 2);
        assert_eq!(region.features.points.len(), 6);
        assert_eq!(region.features.edges.len(), 6);
        assert_eq!(region.features.faces.len(), 0);
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
        let region = Region::with_binary_file(full_path, Extract::Boundary)?;
        assert_eq!(region.mesh.ndim, 2);
        assert_eq!(region.features.points.len(), 6);
        assert_eq!(region.features.edges.len(), 6);
        assert_eq!(region.features.faces.len(), 0);
        assert_eq!(region.find.points(At::XY(0.0, 0.0))?.len(), 1);
        Ok(())
    }
}
