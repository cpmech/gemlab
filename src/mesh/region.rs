use super::{Extract, Features, Find, MapEdge2dToCells, MapFaceToCells, Mesh};
use crate::StrError;

/// Holds all (immutable) data related to a mesh including element shapes, boundary, interior, and functions to find features
pub struct Region<'a> {
    /// Holds the raw mesh data
    pub mesh: &'a Mesh,

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

impl<'a> Region<'a> {
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
    pub fn with(mesh: &'a Mesh, extract: Extract) -> Result<Self, StrError> {
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
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Region;
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
        let region = Region::with(&mesh, Extract::Boundary)?;
        // println!("{:?}", mesh);
        assert_eq!(region.mesh.ndim, 2);
        assert_eq!(region.features.points.len(), 6);
        assert_eq!(region.features.edges.len(), 6);
        assert_eq!(region.features.faces.len(), 0);
        assert_eq!(region.find.points(At::XY(0.0, 0.0))?.len(), 1);
        Ok(())
    }
}
