use super::algorithms;
use super::{CellId, Mesh, PointId};
use crate::shapes::{GeoKind, Shape};
use std::collections::{HashMap, HashSet};

/// Aliases (usize,usize) as the key of edges
///
/// **Note:** Since the local numbering scheme runs over "corners" first,
/// we can compare edges using only two points; i.e., the middle points don't matter.
pub type EdgeKey = (usize, usize);

/// Aliases (usize,usize,usize,usize) as the key of faces
///
/// **Note:** If a face has at most 3 points, the fourth entry in the key will be
/// set to the total number of points. In this way, we can compare 4-node (or more nodes)
/// faces with each other. Further, since the local numbering scheme runs over the
/// "corners" first, the middle points don't matter.
pub type FaceKey = (usize, usize, usize, usize);

/// Holds the point ids of an edge (an entity belonging to a solid cell in 2D or a face in 3D)
#[derive(Clone, Debug)]
pub struct Edge {
    /// Geometry kind
    pub kind: GeoKind,

    /// List of points defining this edge; in the right order (unsorted)
    pub points: Vec<PointId>,
}

/// Holds the point ids of a face (an entity belonging to a solid cell in 3D)
#[derive(Clone, Debug)]
pub struct Face {
    /// Geometry kind
    pub kind: GeoKind,

    /// List of points defining this face; in the right order (unsorted)
    pub points: Vec<PointId>,
}

/// Maps edges to cells sharing the edge (2D only)
///
/// Relates edge keys to `Vec<(cell_id, e)>` where:
///
/// * `cell_id` -- is he id of the cell sharing the edge
/// * `e` -- is the cell's local edge index
pub type MapEdge2dToCells = HashMap<EdgeKey, Vec<(CellId, usize)>>;

/// Maps faces to cells sharing the face (3D only)
///
/// Relates face keys to `Vec<(cell_id, f)>` where:
///
/// * `cell_id` -- is the id of the cell sharing the face
/// * `f` -- is the cell's local face index
pub type MapFaceToCells = HashMap<FaceKey, Vec<(CellId, usize)>>;

/// Holds points, edges and faces on the mesh boundary or interior
pub struct Features {
    /// Set of points on the boundary edges/faces, on the interior edges/faces, or both boundary and interior
    ///
    /// **Notes:**
    ///
    /// 1. Here, a boundary point is such that it belongs to a boundary edge or a boundary face
    /// 2. An interior point is such that it belongs to an interior edge or an interior face
    /// 3. Thus, for the interior case, we save only the points on interior edges and faces
    ///    and **not** inside cells. For example, the middle nodes of a Qua9 are not saved.
    pub points: HashSet<PointId>,

    /// Set of edges on the mesh boundary, interior, or both boundary and interior
    ///
    /// **Notes:**
    ///
    /// 1. In 2D, a boundary edge is such that it is shared by one 2D cell only (1D cells are ignored)
    /// 2. In 3D, a boundary edge belongs to a boundary face
    /// 3. In 2D, an interior edge is such that it is shared by **more** than one 2D cell (1D cells are ignored)
    /// 4. In 3D, an interior edge belongs to an interior face
    pub edges: HashMap<EdgeKey, Edge>,

    /// Set of faces on the mesh boundary, interior, or both boundary and interior
    ///
    /// **Notes:**
    ///
    /// 1. A boundary face is such that it is shared by one 3D cell only (2D cells are ignored)
    /// 2. An interior face is such that it is shared by **more** than one 3D cell (2D cells are ignored)
    pub faces: HashMap<FaceKey, Face>,

    /// The minimum coordinates of the points (space_ndim)
    pub min: Vec<f64>,

    /// The maximum coordinates of the points (space_ndim)
    pub max: Vec<f64>,
}

/// Defines what features to extract
pub enum Extract {
    /// Extracts boundary and interior features
    All,

    /// Extracts boundary features only
    Boundary,

    /// Extracts interior features only
    Interior,
}

impl Features {
    /// Extracts features
    ///
    /// **Note:** This function is to be used by [Region] and not by the end-user.
    ///
    /// **Note:** The points of rods or shells are only extracted with the All or Boundary options.
    ///
    /// # Panics
    ///
    /// 1. This function panics if `shapes.len() != mesh.cells.len()` (i.e., the shapes vector and the mesh must be compatible)
    /// 2. This function works in 2D or 3D only; otherwise, it will panic.
    pub(crate) fn new(
        mesh: &Mesh,
        shapes: &Vec<Shape>,
        extract: Extract,
    ) -> (Option<MapEdge2dToCells>, Option<MapFaceToCells>, Features) {
        assert!(mesh.ndim >= 2 && mesh.ndim <= 3);
        let two_dim = if mesh.ndim == 2 { true } else { false };
        let do_rods_and_shells = match extract {
            Extract::All => true,
            Extract::Boundary => true,
            Extract::Interior => false,
        };
        let (edges, faces, mut features) = if two_dim {
            let edges = algorithms::extract_all_2d_edges(mesh, shapes);
            let features = algorithms::extract_features_2d(mesh, shapes, &edges, extract);
            (Some(edges), None, features)
        } else {
            let faces = algorithms::extract_all_faces(mesh, shapes);
            let features = algorithms::extract_features_3d(mesh, shapes, &faces, extract);
            (None, Some(faces), features)
        };
        if do_rods_and_shells {
            algorithms::extract_rods_and_shells(mesh, &mut features);
        }
        (edges, faces, features)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Edge, Face};
    use crate::shapes::GeoKind;

    #[test]
    fn derive_works() {
        let edge = Edge {
            kind: GeoKind::Lin3,
            points: vec![10, 20, 33],
        };
        let face = Face {
            kind: GeoKind::Qua4,
            points: vec![1, 2, 3, 4],
        };
        let edge_clone = edge.clone();
        let face_clone = face.clone();
        assert_eq!(format!("{:?}", edge), "Edge { kind: Lin3, points: [10, 20, 33] }");
        assert_eq!(format!("{:?}", face), "Face { kind: Qua4, points: [1, 2, 3, 4] }");
        assert_eq!(edge_clone.points.len(), 3);
        assert_eq!(face_clone.points.len(), 4);
    }
}
