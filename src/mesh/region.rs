use super::{allocate_shapes, Extract, Features, Find, MapEdge2dToCells, MapFaceToCells, Mesh};
use crate::shapes::Shape;
use crate::StrError;
use std::ffi::OsStr;

/// Holds all (immutable) data related to a mesh including element shapes, boundary, interior, and functions to find features
pub struct Region {
    /// Holds the raw mesh data
    pub mesh: Mesh,

    /// Holds all shapes of all cells (len = **number of cells**)
    pub shapes: Vec<Shape>,

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
    pub fn with(mesh: Mesh, extract: Extract) -> Result<Self, StrError> {
        let shapes = allocate_shapes(&mesh)?;
        let (all_2d_edges, all_faces, features) = Features::new(&mesh, &shapes, extract)?;
        let find = Find::new(&mesh, &features)?;
        Ok(Region {
            mesh,
            shapes,
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
mod tests {}
