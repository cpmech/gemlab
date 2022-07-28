use super::{Extract, Features, Find, MapEdge2dToCells, MapFaceToCells, Mesh};
use crate::StrError;

/// Holds all (immutable) data related to a mesh and functions to find features (points, edges, faces)
///
/// # Examples
///
/// ## Two-dimensions
///
/// ```
/// use gemlab::mesh::{At, Cell, Extract, Mesh, Point, Region};
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
///         ndim: 2,
///         #[rustfmt::skip]
///          points: vec![
///              Point { id: 0, coords: vec![0.0, 0.0] },
///              Point { id: 1, coords: vec![1.0, 0.0] },
///              Point { id: 2, coords: vec![1.0, 1.0] },
///              Point { id: 3, coords: vec![0.0, 1.0] },
///              Point { id: 4, coords: vec![2.0, 0.0] },
///              Point { id: 5, coords: vec![2.0, 1.0] },
///          ],
///         cells: vec![
///             Cell {
///                 id: 0,
///                 attribute_id: 1,
///                 kind: GeoKind::Qua4,
///                 points: vec![0, 1, 2, 3],
///             },
///             Cell {
///                 id: 1,
///                 attribute_id: 2,
///                 kind: GeoKind::Qua4,
///                 points: vec![1, 4, 5, 2],
///             },
///         ],
///     };
///
///     let region = Region::new(&mesh, Extract::Boundary)?;
///
///     // check features
///
///     let mut points: Vec<_> = region.features.points.iter().copied().collect();
///     points.sort();
///     assert_eq!(points, [0, 1, 2, 3, 4, 5]);
///
///     let mut edges: Vec<_> = region.features.edges.keys().copied().collect();
///     edges.sort();
///     assert_eq!(edges, [(0, 1), (0, 3), (1, 4), (2, 3), (2, 5), (4, 5)]);
///
///     // find features
///
///     let mut points: Vec<_> = region.find.points(At::X(2.0))?.iter().copied().collect();
///     points.sort();
///     assert_eq!(points, &[4, 5]);
///
///     let mut edges: Vec<_> = region.find.edges(At::Y(1.0))?.iter().copied().collect();
///     edges.sort();
///     assert_eq!(edges, &[(2, 3), (2, 5)]);
///
///     Ok(())
/// }
/// ```
///
/// ## Three-dimensions
///
/// ```
/// use gemlab::mesh::{At, Cell, Extract, Mesh, Point, Region};
/// use gemlab::shapes::GeoKind;
/// use gemlab::StrError;
///
/// fn main() -> Result<(), StrError> {
///     //          .4--------------7
///     //        ,' |            ,'|
///     //      ,'              ,'  |
///     //    ,'     |        ,'    |
///     //  5'==============6'      |
///     //  |               |       |
///     //  |        |      |       |
///     //  |       ,0- - - | - - - 3
///     //  |     ,'        |     ,'
///     //  |   ,'          |   ,'
///     //  | ,'            | ,'
///     //  1'--------------2'
///     #[rustfmt::skip]
///      let mesh = Mesh {
///          ndim: 3,
///          points: vec![
///              Point { id: 0, coords: vec![0.0, 0.0, 0.0] },
///              Point { id: 1, coords: vec![1.0, 0.0, 0.0] },
///              Point { id: 2, coords: vec![1.0, 1.0, 0.0] },
///              Point { id: 3, coords: vec![0.0, 1.0, 0.0] },
///              Point { id: 4, coords: vec![0.0, 0.0, 1.0] },
///              Point { id: 5, coords: vec![1.0, 0.0, 1.0] },
///              Point { id: 6, coords: vec![1.0, 1.0, 1.0] },
///              Point { id: 7, coords: vec![0.0, 1.0, 1.0] },
///          ],
///          cells: vec![
///              Cell { id: 0, attribute_id: 1, kind: GeoKind::Hex8,
///                                  points: vec![0,1,2,3, 4,5,6,7] },
///          ],
///      };
///
///     let region = Region::new(&mesh, Extract::Boundary)?;
///
///     // check features
///
///     let mut points: Vec<_> = region.features.points.iter().copied().collect();
///     points.sort();
///     assert_eq!(points, (0..8).collect::<Vec<_>>());
///
///     let mut edges: Vec<_> = region.features.edges.keys().copied().collect();
///     edges.sort();
///     #[rustfmt::skip]
///     assert_eq!(
///         edges,
///         [
///             (0, 1), (0, 3), (0, 4), (1, 2),
///             (1, 5), (2, 3), (2, 6), (3, 7),
///             (4, 5), (4, 7), (5, 6), (6, 7)
///         ]
///     );
///
///     let mut faces: Vec<_> = region.features.faces.keys().copied().collect();
///     faces.sort();
///     assert_eq!(
///         faces,
///         [
///             (0, 1, 2, 3),
///             (0, 1, 4, 5),
///             (0, 3, 4, 7),
///             (1, 2, 5, 6),
///             (2, 3, 6, 7),
///             (4, 5, 6, 7),
///         ]
///     );
///
///     // find features
///
///     let mut points: Vec<_> = region.find.points(At::XY(1.0, 1.0))?.iter().copied().collect();
///     points.sort();
///     assert_eq!(points, &[2, 6]);
///
///     let mut edges: Vec<_> = region.find.edges(At::YZ(1.0, 1.0))?.iter().copied().collect();
///     edges.sort();
///     assert_eq!(edges, &[(6, 7)]);
///
///     let mut faces: Vec<_> = region.find.faces(At::Y(1.0))?.iter().copied().collect();
///     faces.sort();
///     assert_eq!(faces, &[(2, 3, 6, 7)]);
///     Ok(())
/// }
/// ```
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
    pub fn new(mesh: &'a Mesh, extract: Extract) -> Result<Self, StrError> {
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

    #[test]
    fn with_works() {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_qua4();
        let region = Region::new(&mesh, Extract::Boundary).unwrap();
        assert_eq!(region.mesh.ndim, 2);
        assert_eq!(region.features.points.len(), 6);
        assert_eq!(region.features.edges.len(), 6);
        assert_eq!(region.features.faces.len(), 0);
        assert_eq!(region.find.points(At::XY(0.0, 0.0)).unwrap().len(), 1);
    }
}
