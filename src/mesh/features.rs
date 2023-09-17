use super::{algorithms, CellId, Extract, Mesh, PointId};
use crate::shapes::GeoKind;
use crate::util::GridSearch;
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

/// Holds the point ids of an edge or a face
///
/// * An edge is an entity belonging to a solid cell in 2D or a face in 3D
/// * A face is an entity belonging to a solid cell in 3D
#[derive(Clone, Debug)]
pub struct Feature {
    /// Geometry kind
    pub kind: GeoKind,

    /// List of points defining this edge or face; in the right (FEM) order (i.e., unsorted)
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

/// Maps a point id to edges sharing the point
///
/// Relates a point id to a unique set of EdgeKey
pub type MapPointToEdges = HashMap<PointId, HashSet<EdgeKey>>;

/// Maps a point id to faces sharing the point
///
/// Relates a point id to a unique set of FaceKey
pub type MapPointToFaces = HashMap<PointId, HashSet<FaceKey>>;

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
///     #[rustfmt::skip]
///     let mesh = Mesh {
///         ndim: 2,
///         points: vec![
///             Point { id: 0, marker: 0, coords: vec![0.0, 0.0] },
///             Point { id: 1, marker: 0, coords: vec![1.0, 0.0] },
///             Point { id: 2, marker: 0, coords: vec![1.0, 1.0] },
///             Point { id: 3, marker: 0, coords: vec![0.0, 1.0] },
///             Point { id: 4, marker: 0, coords: vec![2.0, 0.0] },
///             Point { id: 5, marker: 0, coords: vec![2.0, 1.0] },
///         ],
///         cells: vec![
///             Cell { id: 0, attribute: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
///             Cell { id: 1, attribute: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
///         ],
///     };
///
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
///              Point { id: 0, marker: 0, coords: vec![0.0, 0.0, 0.0] },
///              Point { id: 1, marker: 0, coords: vec![1.0, 0.0, 0.0] },
///              Point { id: 2, marker: 0, coords: vec![1.0, 1.0, 0.0] },
///              Point { id: 3, marker: 0, coords: vec![0.0, 1.0, 0.0] },
///              Point { id: 4, marker: 0, coords: vec![0.0, 0.0, 1.0] },
///              Point { id: 5, marker: 0, coords: vec![1.0, 0.0, 1.0] },
///              Point { id: 6, marker: 0, coords: vec![1.0, 1.0, 1.0] },
///              Point { id: 7, marker: 0, coords: vec![0.0, 1.0, 1.0] },
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
    ///
    /// Relates edge keys to `Vec<(cell_id, e)>` where:
    ///
    /// * `cell_id` -- is he id of the cell sharing the edge
    /// * `e` -- is the cell's local edge index
    pub all_2d_edges: MapEdge2dToCells,

    /// Maps all face keys to cells sharing the face (3D only)
    ///
    /// Relates face keys to `Vec<(cell_id, f)>` where:
    ///
    /// * `cell_id` -- is the id of the cell sharing the face
    /// * `f` -- is the cell's local face index
    pub all_faces: MapFaceToCells,

    /// Set of points on the boundary edges/faces, on the interior edges/faces, or both boundary and interior
    ///
    /// **Notes:**
    ///
    /// 1. Here, a boundary point is such that it belongs to a boundary edge or a boundary face
    /// 2. An interior point is such that it belongs to an interior edge or an interior face
    /// 3. Thus, for the interior case, we save only the points on interior edges and faces
    ///    and **not** inside cells. For example, the central nodes of a Qua9 are not saved.
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

    /// Space number of dimension (needed for the find functions)
    pub space_ndim: usize,

    /// Tool to quickly find points by coordinates
    pub grid: GridSearch,

    /// Maps a point id to edges sharing the point
    pub point_to_edges: MapPointToEdges,

    /// Maps a point id to faces sharing the point
    pub point_to_faces: MapPointToFaces,
}

impl Features {
    /// Extracts features
    ///
    /// # Notes
    ///
    /// * The points of rods or shells are only extracted when either the All or Boundary option is selected
    /// * You may want to call [Mesh::check_all] to capture (some) errors of the mesh first
    ///
    /// # Panics
    ///
    /// * This function will panic if the mesh data is invalid. For instance, when
    ///   the cell points array doesn't contain enough points or the indices are incorrect
    pub fn new(mesh: &Mesh, extract: Extract) -> Self {
        // options
        assert!(mesh.ndim >= 2 && mesh.ndim <= 3);
        let space_ndim = mesh.ndim;
        let do_rods_and_shells = match extract {
            Extract::All => true,
            Extract::Boundary => true,
            Extract::Interior => false,
        };

        // define variables
        let all_2d_edges: MapEdge2dToCells;
        let all_faces: MapFaceToCells;
        let mut points: HashSet<PointId>;
        let edges: HashMap<EdgeKey, Feature>;
        let faces: HashMap<FaceKey, Feature>;
        let mut min: Vec<f64>;
        let mut max: Vec<f64>;

        // extract features
        if space_ndim == 2 {
            all_2d_edges = algorithms::extract_all_2d_edges(mesh);
            all_faces = HashMap::new();
            faces = HashMap::new();
            (points, edges, min, max) = algorithms::extract_features_2d(mesh, &all_2d_edges, extract);
        } else {
            all_2d_edges = HashMap::new();
            all_faces = algorithms::extract_all_faces(mesh);
            (points, edges, faces, min, max) = algorithms::extract_features_3d(mesh, &all_faces, extract);
        };

        // handle rods and shells
        if do_rods_and_shells {
            mesh.cells.iter().for_each(|cell| {
                let geo_ndim = cell.kind.ndim();
                if geo_ndim == 1 || (geo_ndim == 2 && mesh.ndim == 3) {
                    cell.points.iter().for_each(|id| {
                        points.insert(*id);
                        for j in 0..mesh.ndim {
                            min[j] = f64::min(min[j], mesh.points[*id].coords[j]);
                            max[j] = f64::max(max[j], mesh.points[*id].coords[j]);
                        }
                    });
                }
            });
        }

        // add point ids to grid
        let mut grid = GridSearch::new(&min, &max, None, None, None).unwrap();
        for point_id in &points {
            grid.insert(*point_id, &mesh.points[*point_id].coords).unwrap();
        }

        // map point ids to edges
        let mut point_to_edges = HashMap::new();
        for (edge_key, edge) in &edges {
            for point_id in &edge.points {
                point_to_edges
                    .entry(*point_id)
                    .or_insert(HashSet::new())
                    .insert(*edge_key);
            }
        }

        // map point ids to faces
        let mut point_to_faces = HashMap::new();
        for (face_key, face) in &faces {
            for point_id in &face.points {
                point_to_faces
                    .entry(*point_id)
                    .or_insert(HashSet::new())
                    .insert(*face_key);
            }
        }

        // results
        Features {
            all_2d_edges,
            all_faces,
            points,
            edges,
            faces,
            min,
            max,
            space_ndim,
            grid,
            point_to_edges,
            point_to_faces,
        }
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

/// Returns the 2D edge key of a given cell and edge index
///
/// # Input
///
/// * `mesh` -- the mesh
/// * `cell_id` -- the cell ID
/// * `e` -- the index of the edge (see [GeoKind::edge_node_id])
///
/// # Output
///
/// Returns `(a, b)`, a **sorted** pair of point indices.
///
/// # Examples
///
/// ```
/// use gemlab::mesh::{get_edge_key_2d, Cell, Mesh, Point};
/// use gemlab::shapes::GeoKind;
/// use gemlab::StrError;
///
/// fn main() -> Result<(), StrError> {
///     //  3---------2---------5
///     //  |         |         |
///     //  |   [0]   |   [1]   |
///     //  |         |         |
///     //  0---------1---------4
///     #[rustfmt::skip]
///     let mesh = Mesh {
///         ndim: 2,
///         points: vec![
///             Point { id: 0, marker: 0, coords: vec![0.0, 0.0] },
///             Point { id: 1, marker: 0, coords: vec![1.0, 0.0] },
///             Point { id: 2, marker: 0, coords: vec![1.0, 1.0] },
///             Point { id: 3, marker: 0, coords: vec![0.0, 1.0] },
///             Point { id: 4, marker: 0, coords: vec![2.0, 0.0] },
///             Point { id: 5, marker: 0, coords: vec![2.0, 1.0] },
///         ],
///         cells: vec![
///             Cell { id: 0, attribute: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
///             Cell { id: 1, attribute: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
///         ],
///     };
///
///     assert_eq!(get_edge_key_2d(&mesh, 0, 0), (0, 1));
///     assert_eq!(get_edge_key_2d(&mesh, 0, 1), (1, 2));
///     assert_eq!(get_edge_key_2d(&mesh, 0, 2), (2, 3));
///     assert_eq!(get_edge_key_2d(&mesh, 0, 3), (0, 3));
///
///     assert_eq!(get_edge_key_2d(&mesh, 1, 0), (1, 4));
///     assert_eq!(get_edge_key_2d(&mesh, 1, 1), (4, 5));
///     assert_eq!(get_edge_key_2d(&mesh, 1, 2), (2, 5));
///     assert_eq!(get_edge_key_2d(&mesh, 1, 3), (1, 2));
///     Ok(())
/// }
/// ```
pub fn get_edge_key_2d(mesh: &Mesh, cell_id: CellId, e: usize) -> (usize, usize) {
    let cell = &mesh.cells[cell_id];
    let local_a = cell.kind.edge_node_id(e, 0);
    let local_b = cell.kind.edge_node_id(e, 1);
    let mut a = cell.points[local_a];
    let mut b = cell.points[local_b];
    if b < a {
        let temp = a;
        a = b;
        b = temp;
    }
    (a, b)
}

/// Returns all neighbors of a 2D cell
///
/// # Input
///
/// * `mesh` -- the mesh
/// * `edges` -- the edge-to-cells map
/// * `e` -- the index of the edge (see [GeoKind::edge_node_id])
///
/// # Output
///
/// Returns a vector with the pairs `(e, neigh_cell_id, neigh_e)` where:
///
/// * `e` -- the local edge of this cell through which the neighbor is in contact
/// * `neigh_cell_id` -- is the ID of the neighbor cell
/// * `neigh_e` -- is the neighbor's local edge through which this cell is in contact
///
/// # Examples
///
/// ```
/// use gemlab::mesh::{get_neighbors_2d, Cell, Extract, Features, Mesh, Point};
/// use gemlab::shapes::GeoKind;
/// use gemlab::StrError;
///
/// fn main() -> Result<(), StrError> {
///     //  3---------2---------5
///     //  |         |         |
///     //  |   [0]   |   [1]   |
///     //  |         |         |
///     //  0---------1---------4
///     #[rustfmt::skip]
///     let mesh = Mesh {
///         ndim: 2,
///         points: vec![
///             Point { id: 0, marker: 0, coords: vec![0.0, 0.0] },
///             Point { id: 1, marker: 0, coords: vec![1.0, 0.0] },
///             Point { id: 2, marker: 0, coords: vec![1.0, 1.0] },
///             Point { id: 3, marker: 0, coords: vec![0.0, 1.0] },
///             Point { id: 4, marker: 0, coords: vec![2.0, 0.0] },
///             Point { id: 5, marker: 0, coords: vec![2.0, 1.0] },
///         ],
///         cells: vec![
///             Cell { id: 0, attribute: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
///             Cell { id: 1, attribute: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
///         ],
///     };
///
///     let features = Features::new(&mesh, Extract::All);
///     let edges = features.all_2d_edges;
///
///     let neighbors = get_neighbors_2d(&mesh, &edges, 0);
///     assert_eq!(neighbors.len(), 1);
///     assert!(neighbors.contains(&(1, 1, 3)));
///
///     let neighbors = get_neighbors_2d(&mesh, &edges, 1);
///     assert_eq!(neighbors.len(), 1);
///     assert!(neighbors.contains(&(3, 0, 1)));
///
///     Ok(())
/// }
/// ```
pub fn get_neighbors_2d(mesh: &Mesh, edges: &MapEdge2dToCells, cell_id: CellId) -> Vec<(usize, CellId, usize)> {
    let cell = &mesh.cells[cell_id];
    let nedge = cell.kind.nedge();
    let mut res = Vec::new();
    for e in 0..nedge {
        let local_a = cell.kind.edge_node_id(e, 0);
        let local_b = cell.kind.edge_node_id(e, 1);
        let mut this_a = cell.points[local_a];
        let mut this_b = cell.points[local_b];
        if this_b < this_a {
            let temp = this_a;
            this_a = this_b;
            this_b = temp;
        }
        let shares = edges.get(&(this_a, this_b)).unwrap();
        for share in shares {
            if share.0 != cell_id {
                res.push((e, share.0, share.1));
            }
        }
    }
    res
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{get_edge_key_2d, get_neighbors_2d, Extract, Feature, Features};
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
    fn get_edge_key_2d_works() {
        let mesh = Samples::block_2d_four_qua12().clone();

        assert_eq!(get_edge_key_2d(&mesh, 0, 0), (0, 1));
        assert_eq!(get_edge_key_2d(&mesh, 0, 1), (1, 2));
        assert_eq!(get_edge_key_2d(&mesh, 0, 2), (2, 3));
        assert_eq!(get_edge_key_2d(&mesh, 0, 3), (0, 3));

        assert_eq!(get_edge_key_2d(&mesh, 1, 0), (1, 12));
        assert_eq!(get_edge_key_2d(&mesh, 1, 1), (12, 13));
        assert_eq!(get_edge_key_2d(&mesh, 1, 2), (2, 13));
        assert_eq!(get_edge_key_2d(&mesh, 1, 3), (1, 2));

        assert_eq!(get_edge_key_2d(&mesh, 2, 0), (2, 3));
        assert_eq!(get_edge_key_2d(&mesh, 2, 1), (2, 20));
        assert_eq!(get_edge_key_2d(&mesh, 2, 2), (20, 21));
        assert_eq!(get_edge_key_2d(&mesh, 2, 3), (3, 21));

        assert_eq!(get_edge_key_2d(&mesh, 3, 0), (2, 13));
        assert_eq!(get_edge_key_2d(&mesh, 3, 1), (13, 28));
        assert_eq!(get_edge_key_2d(&mesh, 3, 2), (20, 28));
        assert_eq!(get_edge_key_2d(&mesh, 3, 3), (2, 20));
    }

    #[test]
    fn neighbors_are_correct_2d() {
        let mesh = Samples::block_2d_four_qua12().clone();
        let features = Features::new(&mesh, Extract::All);
        let edges = features.all_2d_edges;
        let mut keys: Vec<_> = edges.keys().into_iter().collect();
        keys.sort();
        // for key in &keys {
        //     println!("{:2?} => {:?}", key, edges.get(key).unwrap());
        // }
        assert_eq!(edges.len(), 12);
        assert_eq!(
            keys,
            &[
                &(0, 1),
                &(0, 3),
                &(1, 2),
                &(1, 12),
                &(2, 3),
                &(2, 13),
                &(2, 20),
                &(3, 21),
                &(12, 13),
                &(13, 28),
                &(20, 21),
                &(20, 28)
            ]
        );
        assert_eq!(edges.get(&(0, 1)), Some(&vec![(0, 0)]));
        assert_eq!(edges.get(&(0, 3)), Some(&vec![(0, 3)]));
        assert_eq!(edges.get(&(1, 2)), Some(&vec![(0, 1), (1, 3)]));
        assert_eq!(edges.get(&(1, 12)), Some(&vec![(1, 0)]));
        assert_eq!(edges.get(&(2, 3)), Some(&vec![(0, 2), (2, 0)]));
        assert_eq!(edges.get(&(2, 13)), Some(&vec![(1, 2), (3, 0)]));
        assert_eq!(edges.get(&(2, 20)), Some(&vec![(2, 1), (3, 3)]));
        assert_eq!(edges.get(&(3, 21)), Some(&vec![(2, 3)]));
        assert_eq!(edges.get(&(12, 13)), Some(&vec![(1, 1)]));
        assert_eq!(edges.get(&(13, 28)), Some(&vec![(3, 1)]));
        assert_eq!(edges.get(&(20, 21)), Some(&vec![(2, 2)]));
        assert_eq!(edges.get(&(20, 28)), Some(&vec![(3, 2)]));
    }

    #[test]
    fn get_neighbors_2d_works() {
        let mesh = Samples::block_2d_four_qua4();
        let features = Features::new(&mesh, Extract::All);
        let edges = features.all_2d_edges;

        let neighbors = get_neighbors_2d(&mesh, &edges, 0);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(1, 1, 3)));
        assert!(neighbors.contains(&(2, 2, 0)));

        let neighbors = get_neighbors_2d(&mesh, &edges, 1);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(3, 0, 1)));
        assert!(neighbors.contains(&(2, 3, 0)));

        let neighbors = get_neighbors_2d(&mesh, &edges, 2);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(0, 0, 2)));
        assert!(neighbors.contains(&(1, 3, 3)));

        let neighbors = get_neighbors_2d(&mesh, &edges, 3);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(0, 1, 2)));
        assert!(neighbors.contains(&(3, 2, 1)));
    }
}
