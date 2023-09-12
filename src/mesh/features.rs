use super::algorithms;
use super::{CellId, Mesh, PointId};
use crate::shapes::GeoKind;
use std::collections::{HashMap, HashSet};
use std::fmt::Write;

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

/// Holds the point ids of an edge or a face
///
/// * An edge is an entity belonging to a solid cell in 2D or a face in 3D
/// * A face is an entity belonging to a solid cell in 3D
#[derive(Clone, Debug)]
pub struct Feature {
    /// Geometry kind
    pub kind: GeoKind,

    /// List of points defining this edge; in the right order (unsorted)
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
///
/// # Examples
///
/// ## Two-dimensions
///
/// ```
/// use gemlab::mesh::{At, Cell, Extract, Features, Mesh, Point};
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
///                 attribute: 1,
///                 kind: GeoKind::Qua4,
///                 points: vec![0, 1, 2, 3],
///             },
///             Cell {
///                 id: 1,
///                 attribute: 2,
///                 kind: GeoKind::Qua4,
///                 points: vec![1, 4, 5, 2],
///             },
///         ],
///     };
///     let features = Features::new(&mesh, Extract::Boundary);
///
///     let mut points: Vec<_> = features.points.iter().copied().collect();
///     points.sort();
///     assert_eq!(points, [0, 1, 2, 3, 4, 5]);
///
///     let mut edges: Vec<_> = features.edges.keys().copied().collect();
///     edges.sort();
///     assert_eq!(edges, [(0, 1), (0, 3), (1, 4), (2, 3), (2, 5), (4, 5)]);
///     Ok(())
/// }
/// ```
///
/// ## Three-dimensions
///
/// ```
/// use gemlab::mesh::{At, Cell, Extract, Features, Mesh, Point};
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
///              Cell { id: 0, attribute: 1, kind: GeoKind::Hex8,
///                                  points: vec![0,1,2,3, 4,5,6,7] },
///          ],
///      };
///
///     let features = Features::new(&mesh, Extract::Boundary);
///
///     let mut points: Vec<_> = features.points.iter().copied().collect();
///     points.sort();
///     assert_eq!(points, (0..8).collect::<Vec<_>>());
///
///     let mut edges: Vec<_> = features.edges.keys().copied().collect();
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
///     let mut faces: Vec<_> = features.faces.keys().copied().collect();
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
///     Ok(())
/// }
/// ```
pub struct Features {
    /// Maps all edge keys to cells sharing the edge (2D only)
    pub all_2d_edges: Option<MapEdge2dToCells>,

    /// Maps all face keys to cells sharing the face (3D only)
    pub all_faces: Option<MapFaceToCells>,

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
    pub edges: HashMap<EdgeKey, Feature>,

    /// Set of faces on the mesh boundary, interior, or both boundary and interior
    ///
    /// **Notes:**
    ///
    /// 1. A boundary face is such that it is shared by one 3D cell only (2D cells are ignored)
    /// 2. An interior face is such that it is shared by **more** than one 3D cell (2D cells are ignored)
    pub faces: HashMap<FaceKey, Feature>,

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
    /// # Notes
    ///
    /// * The points of rods or shells are only extracted when either the All or Boundary option is selected
    /// * You may want to call [super::check_all()] to capture (some) errors of the mesh first
    ///
    /// # Panics
    ///
    /// * This function will panic if the mesh data is invalid. For instance, when
    ///   the cell points array doesn't contain enough points or the indices are incorrect
    pub fn new(mesh: &Mesh, extract: Extract) -> Self {
        assert!(mesh.ndim >= 2 && mesh.ndim <= 3);
        let two_dim = if mesh.ndim == 2 { true } else { false };
        let do_rods_and_shells = match extract {
            Extract::All => true,
            Extract::Boundary => true,
            Extract::Interior => false,
        };
        let mut features = if two_dim {
            let all_2d_edges = algorithms::extract_all_2d_edges(mesh);
            let (points, edges, min, max) = algorithms::extract_features_2d(mesh, &all_2d_edges, extract);
            Features {
                all_2d_edges: Some(all_2d_edges),
                all_faces: None,
                points,
                edges,
                faces: HashMap::new(),
                min,
                max,
            }
        } else {
            let all_faces = algorithms::extract_all_faces(mesh);
            let (points, edges, faces, min, max) = algorithms::extract_features_3d(mesh, &all_faces, extract);
            Features {
                all_2d_edges: None,
                all_faces: Some(all_faces),
                points,
                edges,
                faces,
                min,
                max,
            }
        };
        if do_rods_and_shells {
            algorithms::extract_rods_and_shells(mesh, &mut features);
        }
        features
    }

    /// Returns an edge or panics
    pub fn get_edge(&self, a: usize, b: usize) -> &Feature {
        self.edges.get(&(a, b)).expect("cannot find edge with given key")
    }

    /// Returns an face or panics
    pub fn get_face(&self, a: usize, b: usize, c: usize, d: usize) -> &Feature {
        self.faces.get(&(a, b, c, d)).expect("cannot find face with given key")
    }
}

/// Returns a string with the points of a collection of Feature
pub fn display_features(features: &[&Feature]) -> String {
    let mut buffer = String::new();
    for feature in features {
        let pp = &feature.points;
        if feature.kind.is_lin() {
            // edge
            write!(&mut buffer, "({},{}) ", pp[0], pp[1]).unwrap();
        } else {
            // face
            if pp.len() > 3 {
                write!(&mut buffer, "({},{},{},{}) ", pp[0], pp[1], pp[2], pp[3]).unwrap();
            } else {
                write!(&mut buffer, "({},{},{}) ", pp[0], pp[1], pp[2]).unwrap();
            }
        }
    }
    buffer
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{display_features, Extract, Feature, Features};
    use crate::mesh::Samples;
    use crate::shapes::GeoKind;

    #[test]
    fn new_and_get_methods_work() {
        //      4--------------7  1.0
        //     /.             /|
        //    / .            / |    [#] indicates id
        //   /  .           /  |    (#) indicates attribute
        //  /   .          /   |
        // 5--------------6    |          z
        // |    .         |    |          ↑
        // |    0---------|----3  0.0     o → y
        // |   /  [0]     |   /          ↙
        // |  /   (1)     |  /          x
        // | /            | /
        // |/             |/
        // 1--------------2   1.0
        let mesh = Samples::one_hex8();
        let features = Features::new(&mesh, Extract::Boundary);
        let edge = features.get_edge(4, 5);
        let face = features.get_face(0, 1, 4, 5);
        assert_eq!(edge.points, &[4, 5]);
        assert_eq!(face.points, &[0, 1, 5, 4]);
    }

    #[test]
    #[should_panic(expected = "cannot find edge with given key")]
    fn get_edge_panics_on_error() {
        let mesh = Samples::one_tri3();
        let features = Features::new(&mesh, Extract::Boundary);
        features.get_edge(4, 5);
    }

    #[test]
    #[should_panic(expected = "cannot find face with given key")]
    fn get_face_panics_on_error() {
        let mesh = Samples::one_tet4();
        let features = Features::new(&mesh, Extract::Boundary);
        features.get_face(4, 5, 6, 100);
    }

    #[test]
    fn derive_works() {
        let edge = Feature {
            kind: GeoKind::Lin3,
            points: vec![10, 20, 33],
        };
        let face = Feature {
            kind: GeoKind::Qua4,
            points: vec![1, 2, 3, 4],
        };
        let edge_clone = edge.clone();
        let face_clone = face.clone();
        assert_eq!(format!("{:?}", edge), "Feature { kind: Lin3, points: [10, 20, 33] }");
        assert_eq!(format!("{:?}", face), "Feature { kind: Qua4, points: [1, 2, 3, 4] }");
        assert_eq!(edge_clone.points.len(), 3);
        assert_eq!(face_clone.points.len(), 4);
    }

    #[test]
    fn display_features_works() {
        let edges = vec![
            Feature {
                kind: GeoKind::Lin2,
                points: vec![0, 1],
            },
            Feature {
                kind: GeoKind::Lin2,
                points: vec![4, 5],
            },
            Feature {
                kind: GeoKind::Lin3,
                points: vec![7, 11, 9],
            },
        ];
        let faces = vec![
            Feature {
                kind: GeoKind::Tri3,
                points: vec![0, 3, 2],
            },
            Feature {
                kind: GeoKind::Qua4,
                points: vec![8, 9, 10, 11],
            },
            Feature {
                kind: GeoKind::Qua4,
                points: vec![6, 7, 11, 10],
            },
        ];
        let features: Vec<&Feature> = vec![edges.iter().collect::<Vec<_>>(), faces.iter().collect::<Vec<_>>()].concat();
        let res = display_features(&features);
        assert_eq!(res, "(0,1) (4,5) (7,11) (0,3,2) (8,9,10,11) (6,7,11,10) ");
    }
}
