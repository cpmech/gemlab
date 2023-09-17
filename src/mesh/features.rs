use super::{algorithms, At, CellId, Extract, Mesh, PointId};
use crate::shapes::GeoKind;
use crate::util::GridSearch;
use crate::StrError;
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

    /// Total number of points in the mesh (needed to generate face keys of Tri3)
    ///
    /// **readonly**
    pub num_points: usize,

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
        if mesh.ndim == 2 {
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
            space_ndim: mesh.ndim,
            num_points: mesh.points.len(),
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
    /// use gemlab::mesh::{Cell, Extract, Features, Mesh, Point};
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
    ///
    ///     let neighbors = features.get_neighbors_2d(&mesh, 0);
    ///     assert_eq!(neighbors.len(), 1);
    ///     assert!(neighbors.contains(&(1, 1, 3)));
    ///
    ///     let neighbors = features.get_neighbors_2d(&mesh, 1);
    ///     assert_eq!(neighbors.len(), 1);
    ///     assert!(neighbors.contains(&(3, 0, 1)));
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn get_neighbors_2d(&self, mesh: &Mesh, cell_id: CellId) -> Vec<(usize, CellId, usize)> {
        assert_eq!(mesh.ndim, 2);
        assert_eq!(self.space_ndim, 2);
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
            let shares = self.all_2d_edges.get(&(this_a, this_b)).unwrap();
            for share in shares {
                if share.0 != cell_id {
                    res.push((e, share.0, share.1));
                }
            }
        }
        res
    }

    /// Searches point ids
    ///
    /// # Input
    ///
    /// * `at` -- the constraint
    /// * `filter` -- function `fn(x) -> bool` that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true).
    ///   Use `|_| true` or [crate::util::any_x] to allow any point in the resulting array.
    ///   Example of filter: `|x| x[0] > 1.4 && x[0] < 1.6`
    ///
    /// # Output
    ///
    /// * If at least one point has been found, returns a **sorted** array of point ids.
    /// * Otherwise, returns an error
    pub fn search_point_ids<F>(&self, at: At, filter: F) -> Result<Vec<PointId>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        let mut point_ids: HashSet<PointId> = HashSet::new();
        match at {
            At::X(x) => {
                if self.space_ndim == 2 {
                    for id in self.grid.find_on_line(&[x, 0.0], &[x, 1.0], filter)? {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid.find_on_plane_yz(x, filter)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::Y(y) => {
                if self.space_ndim == 2 {
                    for id in self.grid.find_on_line(&[0.0, y], &[1.0, y], filter)? {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid.find_on_plane_xz(y, filter)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::Z(z) => {
                if self.space_ndim == 2 {
                    return Err("At::Z works in 3D only");
                } else {
                    for id in self.grid.find_on_plane_xy(z, filter)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::XY(x, y) => {
                if self.space_ndim == 2 {
                    if let Some(id) = self.grid.find(&[x, y])? {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid.find_on_line(&[x, y, 0.0], &[x, y, 1.0], filter)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::YZ(y, z) => {
                if self.space_ndim == 2 {
                    return Err("At::YZ works in 3D only");
                } else {
                    for id in self.grid.find_on_line(&[0.0, y, z], &[1.0, y, z], filter)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::XZ(x, z) => {
                if self.space_ndim == 2 {
                    return Err("At::XZ works in 3D only");
                } else {
                    for id in self.grid.find_on_line(&[x, 0.0, z], &[x, 1.0, z], filter)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::XYZ(x, y, z) => {
                if self.space_ndim == 2 {
                    return Err("At::XYZ works in 3D only");
                } else {
                    if let Some(id) = self.grid.find(&[x, y, z])? {
                        point_ids.insert(id);
                    }
                }
            }
            At::Circle(x, y, r) => {
                if self.space_ndim == 2 {
                    for id in self.grid.find_on_circle(&[x, y], r, filter)? {
                        point_ids.insert(id);
                    }
                } else {
                    return Err("At::Circle works in 2D only");
                }
            }
            At::Cylinder(ax, ay, az, bx, by, bz, r) => {
                if self.space_ndim == 2 {
                    return Err("At::Cylinder works in 3D only");
                } else {
                    for id in self.grid.find_on_cylinder(&[ax, ay, az], &[bx, by, bz], r, filter)? {
                        point_ids.insert(id);
                    }
                }
            }
        }
        if point_ids.len() == 0 {
            return Err("cannot find any point with given constraints/filter");
        }
        let mut ids: Vec<_> = point_ids.iter().copied().collect();
        ids.sort();
        Ok(ids)
    }

    /// Searches edge keys
    ///
    /// # Input
    ///
    /// * `at` -- the constraint
    /// * `filter` -- function `fn(x) -> bool` that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true).
    ///   Use `|_| true` or [crate::util::any_x] to allow any point in the resulting array.
    ///   Example of filter: `|x| x[0] > 1.4 && x[0] < 1.6`
    ///
    /// # Output
    ///
    /// * If at least one point has been found, returns a **sorted** array of edge keys
    /// * Otherwise, returns an error
    pub fn search_edge_keys<F>(&self, at: At, filter: F) -> Result<Vec<EdgeKey>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        let mut edge_keys: HashSet<EdgeKey> = HashSet::new();
        // find all points constrained by "at" and "filter"
        let point_ids = self.search_point_ids(at, filter)?;
        for point_id in &point_ids {
            // select all edges connected to the found points
            let edges = self.point_to_edges.get(point_id).unwrap(); // unwrap here because there should be no hanging edges
            for edge_key in edges {
                // accept edge when at least two edge points validate "At"
                if point_ids.contains(&edge_key.0) && point_ids.contains(&edge_key.1) {
                    edge_keys.insert(*edge_key);
                }
            }
        }
        if edge_keys.len() == 0 {
            return Err("cannot find any edge with given constraints/filter");
        }
        let mut keys: Vec<_> = edge_keys.iter().copied().collect();
        keys.sort();
        Ok(keys)
    }

    /// Searches face keys
    ///
    /// # Input
    ///
    /// * `at` -- the constraint
    /// * `filter` -- function `fn(x) -> bool` that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true).
    ///   Use `|_| true` or [crate::util::any_x] to allow any point in the resulting array.
    ///   Example of filter: `|x| x[0] > 1.4 && x[0] < 1.6`
    ///
    /// # Output
    ///
    /// * If at least one point has been found, returns a **sorted** array of face keys
    /// * Otherwise, returns an error
    pub fn search_face_keys<F>(&self, at: At, filter: F) -> Result<Vec<FaceKey>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        if self.space_ndim != 3 {
            return Err("cannot find face keys in 2D");
        }
        let mut face_keys: HashSet<FaceKey> = HashSet::new();
        // find all points constrained by "at" and "filter"
        let point_ids = self.search_point_ids(at, filter)?;
        for point_id in &point_ids {
            // select all faces connected to the found points
            let faces = self.point_to_faces.get(point_id).unwrap(); // unwrap here because there should be no hanging faces
            for face_key in faces {
                // accept face when at least four face points validate "At"
                let fourth_is_ok = if face_key.3 == self.num_points {
                    true
                } else {
                    point_ids.contains(&face_key.3)
                };
                if point_ids.contains(&face_key.0)
                    && point_ids.contains(&face_key.1)
                    && point_ids.contains(&face_key.2)
                    && fourth_is_ok
                {
                    face_keys.insert(*face_key);
                }
            }
        }
        if face_keys.len() == 0 {
            return Err("cannot find any face with given constraints/filter");
        }
        let mut keys: Vec<_> = face_keys.iter().copied().collect();
        keys.sort();
        Ok(keys)
    }

    /// Searches edges
    ///
    /// # Input
    ///
    /// * `at` -- the constraint
    /// * `filter` -- function `fn(x) -> bool` that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true).
    ///   Use `|_| true` or [crate::util::any_x] to allow any point in the resulting array.
    ///   Example of filter: `|x| x[0] > 1.4 && x[0] < 1.6`
    ///
    /// # Output
    ///
    /// * If at least one point has been found, returns an array such that the edge keys are **sorted**
    /// * Otherwise, returns an error
    pub fn search_edges<F>(&self, at: At, filter: F) -> Result<Vec<&Feature>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        let results: Result<Vec<_>, _> = self
            .search_edge_keys(at, filter)?
            .iter()
            .map(|key| {
                self.edges
                    .get(key)
                    .ok_or("INTERNAL ERROR: features.edges data is inconsistent")
            })
            .collect();
        results
    }

    /// Searches faces
    ///
    /// # Input
    ///
    /// * `at` -- the constraint
    /// * `filter` -- function `fn(x) -> bool` that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true).
    ///   Use `|_| true` or [crate::util::any_x] to allow any point in the resulting array.
    ///   Example of filter: `|x| x[0] > 1.4 && x[0] < 1.6`
    ///
    /// # Output
    ///
    /// * If at least one point has been found, returns an array such that the face keys are **sorted**
    /// * Otherwise, returns an error
    pub fn search_faces<F>(&self, at: At, filter: F) -> Result<Vec<&Feature>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        let results: Result<Vec<_>, _> = self
            .search_face_keys(at, filter)?
            .iter()
            .map(|key| {
                self.faces
                    .get(key)
                    .ok_or("INTERNAL ERROR: features.faces data is inconsistent")
            })
            .collect();
        results
    }

    /// Searches many edges using a list of constraints
    ///
    /// # Input
    ///
    /// * `at` -- the constraint
    /// * `filter` -- function `fn(x) -> bool` that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true).
    ///   Use `|_| true` or [crate::util::any_x] to allow any point in the resulting array.
    ///   Example of filter: `|x| x[0] > 1.4 && x[0] < 1.6`
    ///
    /// # Output
    ///
    /// * Returns edges sorted by keys
    /// * **Warning** Every `At` in the `ats` must generate at least one edge,
    ///   otherwise an error will occur.
    pub fn search_many_edges<F>(&self, ats: &[At], mut filter: F) -> Result<Vec<&Feature>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        let mut edge_keys: HashSet<EdgeKey> = HashSet::new();
        for at in ats {
            let found = self.search_edge_keys(*at, &mut filter)?;
            for edge in found {
                edge_keys.insert(edge);
            }
        }
        let mut keys: Vec<_> = edge_keys.iter().collect();
        keys.sort();
        let results: Result<Vec<_>, _> = keys
            .iter()
            .map(|key| {
                self.edges
                    .get(key)
                    .ok_or("INTERNAL ERROR: features.edges data is inconsistent")
            })
            .collect();
        results
    }

    /// Search many faces using a list of constraints
    ///
    /// # Input
    ///
    /// * `at` -- the constraint
    /// * `filter` -- function `fn(x) -> bool` that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true).
    ///   Use `|_| true` or [crate::util::any_x] to allow any point in the resulting array.
    ///   Example of filter: `|x| x[0] > 1.4 && x[0] < 1.6`
    ///
    /// # Output
    ///
    /// * Returns faces sorted by keys
    /// * **Warning** Every `At` in the `ats` must generate at least one face,
    ///   otherwise an error will occur.
    pub fn search_many_faces<F>(&self, ats: &[At], mut filter: F) -> Result<Vec<&Feature>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        let mut face_keys: HashSet<FaceKey> = HashSet::new();
        for at in ats {
            let found = self.search_face_keys(*at, &mut filter)?;
            for face in found {
                face_keys.insert(face);
            }
        }
        let mut keys: Vec<_> = face_keys.iter().collect();
        keys.sort();
        let results: Result<Vec<_>, _> = keys
            .iter()
            .map(|key| {
                self.faces
                    .get(key)
                    .ok_or("INTERNAL ERROR: features.faces data is inconsistent")
            })
            .collect();
        results
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Feature, Features};
    use crate::mesh::{At, Extract, Samples};
    use crate::shapes::GeoKind;
    use crate::util::any_x;
    use plotpy::Plot;
    use russell_lab::math::{PI, SQRT_2};

    const SAVE_FIGURE: bool = false;

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
    fn new_method_allocates_grid_2d() {
        let mesh = Samples::two_qua4();
        let feat = Features::new(&mesh, Extract::Boundary);
        assert_eq!(
            format!("{}", feat.grid),
            "0: [0]\n\
             9: [1]\n\
             19: [4]\n\
             180: [3]\n\
             189: [2]\n\
             199: [5]\n\
             ids = [0, 1, 2, 3, 4, 5]\n\
             nitem = 6\n\
             ncontainer = 6\n\
             ndiv = [20, 10]\n"
        );
        if SAVE_FIGURE {
            let mut plot = Plot::new();
            feat.grid.draw(&mut plot).unwrap();
            plot.set_equal_axes(true).set_figure_size_points(800.0, 400.0);
            plot.save("/tmp/gemlab/test_features_grid_two_qua4.svg").unwrap();
        }
    }

    #[test]
    fn new_method_allocates_grid_3d() {
        let mesh = Samples::two_hex8();
        let feat = Features::new(&mesh, Extract::Boundary);
        assert_eq!(
            format!("{}", feat.grid),
            "0: [0]\n\
             9: [1]\n\
             90: [3]\n\
             99: [2]\n\
             900: [4]\n\
             909: [5]\n\
             990: [7]\n\
             999: [6]\n\
             1900: [8]\n\
             1909: [9]\n\
             1990: [11]\n\
             1999: [10]\n\
             ids = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]\n\
             nitem = 12\n\
             ncontainer = 12\n\
             ndiv = [10, 10, 20]\n"
        );
        if SAVE_FIGURE {
            let mut plot = Plot::new();
            feat.grid.draw(&mut plot).unwrap();
            plot.set_equal_axes(true).set_figure_size_points(800.0, 800.0);
            plot.save("/tmp/gemlab/test_features_grid_two_hex8.svg").unwrap();
        }
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
    fn neighbors_are_correct_2d() {
        let mesh = Samples::block_2d_four_qua12().clone();
        let features = Features::new(&mesh, Extract::All);
        let edges = features.all_2d_edges;
        let mut keys: Vec<_> = edges.keys().into_iter().collect();
        keys.sort();
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

        let neighbors = features.get_neighbors_2d(&mesh, 0);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(1, 1, 3)));
        assert!(neighbors.contains(&(2, 2, 0)));

        let neighbors = features.get_neighbors_2d(&mesh, 1);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(3, 0, 1)));
        assert!(neighbors.contains(&(2, 3, 0)));

        let neighbors = features.get_neighbors_2d(&mesh, 2);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(0, 0, 2)));
        assert!(neighbors.contains(&(1, 3, 3)));

        let neighbors = features.get_neighbors_2d(&mesh, 3);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(0, 1, 2)));
        assert!(neighbors.contains(&(3, 2, 1)));
    }

    #[test]
    fn search_points_fails_on_wrong_input() {
        // any_x works
        assert_eq!(any_x(&vec![]), true);

        // 2d
        let mesh = Samples::two_qua4();
        let feat = Features::new(&mesh, Extract::Boundary);
        assert_eq!(
            feat.search_point_ids(At::Z(0.0), any_x).err(),
            Some("At::Z works in 3D only")
        );
        assert_eq!(
            feat.search_point_ids(At::YZ(0.0, 0.0), any_x).err(),
            Some("At::YZ works in 3D only")
        );
        assert_eq!(
            feat.search_point_ids(At::XZ(0.0, 0.0), any_x).err(),
            Some("At::XZ works in 3D only")
        );
        assert_eq!(
            feat.search_point_ids(At::XYZ(0.0, 0.0, 0.0), any_x).err(),
            Some("At::XYZ works in 3D only")
        );
        assert_eq!(
            feat.search_point_ids(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), any_x)
                .err(),
            Some("At::Cylinder works in 3D only")
        );

        // 3d
        let mesh = Samples::two_hex8();
        let boundary = Features::new(&mesh, Extract::Boundary);
        assert_eq!(
            boundary.search_point_ids(At::Circle(0.0, 0.0, 0.0), any_x).err(),
            Some("At::Circle works in 2D only")
        );
    }

    #[test]
    fn search_points_works_2d() {
        // `.       `.
        //   3--------2--------5
        //   | `.     | `.     |
        //   |   `~.  | circle |
        //   |      `.|        |
        //   0--------1--------4
        //           circle
        let mesh = Samples::two_qua4();
        let feat = Features::new(&mesh, Extract::All);
        assert_eq!(feat.search_point_ids(At::XY(0.0, 0.0), any_x).unwrap(), &[0]);
        assert_eq!(feat.search_point_ids(At::XY(2.0, 1.0), any_x).unwrap(), &[5]);
        assert_eq!(
            feat.search_point_ids(At::XY(10.0, 0.0), any_x).err(),
            Some("cannot find point because the coordinates are outside the grid")
        );
        assert_eq!(
            feat.search_point_ids(At::Circle(0.0, 0.0, 1.0), any_x).unwrap(),
            &[1, 3]
        );
        assert_eq!(
            feat.search_point_ids(At::Circle(0.0, 0.0, SQRT_2), any_x).unwrap(),
            &[2]
        );
        assert_eq!(
            feat.search_point_ids(At::Circle(0.0, 0.0, 10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
    }

    #[test]
    fn search_points_works_3d() {
        //      8-----------11  2.0
        //     /.           /|
        //    / .          / |
        //   /  .         /  |
        //  9-----------10   |
        //  |   .        |   |
        //  |   4--------|---7  1.0
        //  |  /.        |  /|
        //  | / .        | / |
        //  |/  .        |/  |
        //  5------------6   |          z
        //  |   .        |   |          ↑
        //  |   0--------|---3  0.0     o → y
        //  |  /         |  /          ↙
        //  | /          | /          x
        //  |/           |/
        //  1------------2   1.0
        // 0.0          1.0
        let mesh = Samples::two_hex8();
        let feat = Features::new(&mesh, Extract::Boundary);
        assert_eq!(feat.search_point_ids(At::X(0.0), any_x).unwrap(), &[0, 3, 4, 7, 8, 11]);
        assert_eq!(feat.search_point_ids(At::X(1.0), any_x).unwrap(), &[1, 2, 5, 6, 9, 10]);
        assert_eq!(
            feat.search_point_ids(At::X(10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(feat.search_point_ids(At::Y(0.0), any_x).unwrap(), &[0, 1, 4, 5, 8, 9]);
        assert_eq!(feat.search_point_ids(At::Y(1.0), any_x).unwrap(), &[2, 3, 6, 7, 10, 11]);
        assert_eq!(
            feat.search_point_ids(At::Y(10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(feat.search_point_ids(At::Z(0.0), any_x).unwrap(), &[0, 1, 2, 3]);
        assert_eq!(feat.search_point_ids(At::Z(1.0), any_x).unwrap(), &[4, 5, 6, 7]);
        assert_eq!(feat.search_point_ids(At::Z(2.0), any_x).unwrap(), &[8, 9, 10, 11]);
        assert_eq!(
            feat.search_point_ids(At::Z(10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(feat.search_point_ids(At::XY(0.0, 0.0), any_x).unwrap(), &[0, 4, 8]);
        assert_eq!(feat.search_point_ids(At::XY(1.0, 1.0), any_x).unwrap(), &[2, 6, 10]);
        assert_eq!(
            feat.search_point_ids(At::XY(10.0, 10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(feat.search_point_ids(At::YZ(0.0, 0.0), any_x).unwrap(), &[0, 1]);
        assert_eq!(feat.search_point_ids(At::YZ(1.0, 1.0), any_x).unwrap(), &[6, 7]);
        assert_eq!(feat.search_point_ids(At::XZ(0.0, 0.0), any_x).unwrap(), &[0, 3]);
        assert_eq!(feat.search_point_ids(At::XZ(1.0, 0.0), any_x).unwrap(), &[1, 2]);
        assert_eq!(feat.search_point_ids(At::XZ(1.0, 2.0), any_x).unwrap(), &[9, 10]);
        assert_eq!(feat.search_point_ids(At::XYZ(0.0, 0.0, 0.0), any_x).unwrap(), &[0]);
        assert_eq!(feat.search_point_ids(At::XYZ(1.0, 1.0, 2.0), any_x).unwrap(), &[10]);
        assert_eq!(
            feat.search_point_ids(At::XYZ(10.0, 0.0, 0.0), any_x).err(),
            Some("cannot find point because the coordinates are outside the grid")
        );
        assert_eq!(
            feat.search_point_ids(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0), any_x)
                .unwrap(),
            &[1, 3, 5, 7, 9, 11],
        );
        assert_eq!(
            feat.search_point_ids(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, SQRT_2), any_x)
                .unwrap(),
            &[2, 6, 10],
        );
        assert_eq!(
            feat.search_point_ids(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 10.0), any_x)
                .err(),
            Some("cannot find any point with given constraints/filter")
        );
    }

    #[test]
    fn search_edges_works_2d() {
        // 3--------2--------5
        // |        |        |
        // |        |        |
        // |        |        |
        // 0--------1--------4
        let mesh = Samples::two_qua4();
        let feat = Features::new(&mesh, Extract::Boundary);
        assert_eq!(feat.search_edge_keys(At::Y(0.0), any_x).unwrap(), &[(0, 1), (1, 4)]);
        assert_eq!(feat.search_edge_keys(At::Y(0.0), |x| x[0] <= 1.0).unwrap(), &[(0, 1)]);
        assert_eq!(feat.search_edge_keys(At::X(2.0), any_x).unwrap(), &[(4, 5)]);
        assert_eq!(feat.search_edge_keys(At::Y(1.0), any_x).unwrap(), &[(2, 3), (2, 5)]);
        assert_eq!(feat.search_edge_keys(At::Y(1.0), |x| x[0] >= 1.0).unwrap(), &[(2, 5)]);
        assert_eq!(feat.search_edge_keys(At::X(0.0), any_x).unwrap(), &[(0, 3)]);

        // internal
        assert_eq!(
            feat.search_edge_keys(At::X(1.0), any_x).err(),
            Some("cannot find any edge with given constraints/filter")
        );

        // far away
        assert_eq!(
            feat.search_edge_keys(At::X(10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );

        // high-level function
        let res = feat.search_edges(At::Y(0.0), any_x).unwrap();
        assert_eq!(res.len(), 2);
        assert_eq!(res[0].points, &[1, 0]);
        assert_eq!(res[1].points, &[4, 1]);
        assert_eq!(
            feat.search_edges(At::XYZ(0.0, 0.0, 0.0), any_x).err(),
            Some("At::XYZ works in 3D only")
        );

        // many edges
        let res = feat
            .search_many_edges(&[At::X(0.0), At::X(2.0), At::Y(0.0), At::Y(1.0)], any_x)
            .unwrap();
        assert_eq!(res.len(), 6);
        let keys: Vec<_> = res.iter().map(|r| (r.points[0], r.points[1])).collect();
        assert_eq!(keys, &[(1, 0), (0, 3), (4, 1), (3, 2), (2, 5), (5, 4)]);
    }

    #[test]
    fn search_edges_works_3d() {
        //      8-----------11  2.0
        //     /.           /|
        //    / .          / |
        //   /  .         /  |
        //  9-----------10   |
        //  |   .        |   |
        //  |   4--------|---7  1.0
        //  |  /.        |  /|
        //  | / .        | / |
        //  |/  .        |/  |
        //  5------------6   |          z
        //  |   .        |   |          ↑
        //  |   0--------|---3  0.0     o → y
        //  |  /         |  /          ↙
        //  | /          | /          x
        //  |/           |/
        //  1------------2   1.0
        // 0.0          1.0
        let mesh = Samples::two_hex8();
        let feat = Features::new(&mesh, Extract::Boundary);
        assert_eq!(
            feat.search_edge_keys(At::X(0.0), any_x).unwrap(),
            &[(0, 3), (0, 4), (3, 7), (4, 7), (4, 8), (7, 11), (8, 11)],
        );
        assert_eq!(feat.search_edge_keys(At::X(0.0), |x| x[2] == 0.0).unwrap(), &[(0, 3)],);
        assert_eq!(
            feat.search_edge_keys(At::X(1.0), any_x).unwrap(),
            &[(1, 2), (1, 5), (2, 6), (5, 6), (5, 9), (6, 10), (9, 10)],
        );
        assert_eq!(
            feat.search_edge_keys(At::X(10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            feat.search_edge_keys(At::Y(0.0), any_x).unwrap(),
            &[(0, 1), (0, 4), (1, 5), (4, 5), (4, 8), (5, 9), (8, 9)],
        );
        assert_eq!(
            feat.search_edge_keys(At::Y(0.0), |x| x[2] <= 1.0).unwrap(),
            &[(0, 1), (0, 4), (1, 5), (4, 5)],
        );
        assert_eq!(
            feat.search_edge_keys(At::Y(1.0), any_x).unwrap(),
            &[(2, 3), (2, 6), (3, 7), (6, 7), (6, 10), (7, 11), (10, 11)],
        );
        assert_eq!(
            feat.search_edge_keys(At::Y(10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            feat.search_edge_keys(At::Z(0.0), any_x).unwrap(),
            &[(0, 1), (0, 3), (1, 2), (2, 3)]
        );
        assert_eq!(
            feat.search_edge_keys(At::Z(2.0), any_x).unwrap(),
            &[(8, 9), (8, 11), (9, 10), (10, 11)],
        );
        assert_eq!(feat.search_edge_keys(At::Z(2.0), |x| x[1] == 1.0).unwrap(), &[(10, 11)],);
        assert_eq!(
            feat.search_edge_keys(At::Z(10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            feat.search_edge_keys(At::XY(0.0, 0.0), any_x).unwrap(),
            &[(0, 4), (4, 8)]
        );
        assert_eq!(
            feat.search_edge_keys(At::XY(1.0, 1.0), any_x).unwrap(),
            &[(2, 6), (6, 10)]
        );
        assert_eq!(
            feat.search_edge_keys(At::XY(1.0, 1.0), |x| x[2] <= 1.0).unwrap(),
            &[(2, 6)]
        );
        assert_eq!(
            feat.search_edge_keys(At::XY(10.0, 10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(feat.search_edge_keys(At::YZ(0.0, 0.0), any_x).unwrap(), &[(0, 1)]);
        assert_eq!(feat.search_edge_keys(At::YZ(1.0, 1.0), any_x).unwrap(), &[(6, 7)]);
        assert_eq!(
            feat.search_edge_keys(At::YZ(10.0, 10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(feat.search_edge_keys(At::XZ(0.0, 0.0), any_x).unwrap(), &[(0, 3)]);
        assert_eq!(feat.search_edge_keys(At::XZ(1.0, 0.0), any_x).unwrap(), &[(1, 2)]);
        assert_eq!(feat.search_edge_keys(At::XZ(1.0, 2.0), any_x).unwrap(), &[(9, 10)]);
        assert_eq!(feat.search_edge_keys(At::XZ(1.0, 2.0), |x| x[0] == -1.0).ok(), None);
        assert_eq!(
            feat.search_edge_keys(At::XZ(10.0, 10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            feat.search_edge_keys(At::XYZ(0.0, 0.0, 0.0), any_x).err(),
            Some("cannot find any edge with given constraints/filter")
        );
        assert_eq!(
            feat.search_edge_keys(At::XYZ(10.0, 0.0, 0.0), any_x).err(),
            Some("cannot find point because the coordinates are outside the grid")
        );
        assert_eq!(
            feat.search_edge_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0), any_x)
                .unwrap(),
            &[(1, 5), (3, 7), (5, 9), (7, 11)],
        );
        assert_eq!(
            feat.search_edge_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0), |x| x[2] >= 1.0)
                .unwrap(),
            &[(5, 9), (7, 11)],
        );
        assert_eq!(
            feat.search_edge_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, SQRT_2), any_x)
                .unwrap(),
            &[(2, 6), (6, 10)],
        );
        assert_eq!(
            feat.search_edge_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 10.0), any_x)
                .err(),
            Some("cannot find any point with given constraints/filter")
        );

        // high-level function
        let res = feat.search_edges(At::XY(0.0, 0.0), any_x).unwrap();
        assert_eq!(res.len(), 2);
        assert_eq!(res[0].points, &[0, 4]);
        assert_eq!(res[1].points, &[4, 8]);

        // many faces
        let res = feat.search_many_faces(&[At::Z(0.0), At::Z(2.0)], any_x).unwrap();
        assert_eq!(res.len(), 2);
        let keys: Vec<_> = res
            .iter()
            .map(|r| (r.points[0], r.points[1], r.points[2], r.points[3]))
            .collect();
        assert_eq!(keys, &[(0, 3, 2, 1), (8, 9, 10, 11)]);
    }

    #[test]
    fn search_faces_returns_error_in_2d() {
        // 3--------2--------5
        // |        |        |
        // |        |        |
        // |        |        |
        // 0--------1--------4
        let mesh = Samples::two_qua4();
        let feat = Features::new(&mesh, Extract::Boundary);
        assert_eq!(
            feat.search_face_keys(At::X(0.0), any_x).err(),
            Some("cannot find face keys in 2D")
        );
    }

    #[test]
    fn search_faces_works() {
        //      8-----------11  2.0
        //     /.           /|
        //    / .          / |
        //   /  .         /  |
        //  9-----------10   |
        //  |   .        |   |
        //  |   4--------|---7  1.0
        //  |  /.        |  /|
        //  | / .        | / |
        //  |/  .        |/  |
        //  5------------6   |          z
        //  |   .        |   |          ↑
        //  |   0--------|---3  0.0     o → y
        //  |  /         |  /          ↙
        //  | /          | /          x
        //  |/           |/
        //  1------------2   1.0
        // 0.0          1.0
        let mesh = Samples::two_hex8();
        let feat = Features::new(&mesh, Extract::Boundary);
        assert_eq!(
            feat.search_face_keys(At::X(0.0), any_x).unwrap(),
            &[(0, 3, 4, 7), (4, 7, 8, 11)]
        );
        assert_eq!(
            feat.search_face_keys(At::X(1.0), any_x).unwrap(),
            &[(1, 2, 5, 6), (5, 6, 9, 10)]
        );
        assert_eq!(
            feat.search_face_keys(At::X(1.0), |x| x[2] <= 1.0).unwrap(),
            &[(1, 2, 5, 6)]
        );
        assert_eq!(
            feat.search_face_keys(At::X(10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::Y(0.0), any_x).unwrap(),
            &[(0, 1, 4, 5), (4, 5, 8, 9)]
        );
        assert_eq!(
            feat.search_face_keys(At::Y(1.0), any_x).unwrap(),
            &[(2, 3, 6, 7), (6, 7, 10, 11)]
        );
        assert_eq!(
            feat.search_face_keys(At::Y(1.0), |x| x[2] >= 1.0).unwrap(),
            &[(6, 7, 10, 11)]
        );
        assert_eq!(
            feat.search_face_keys(At::Y(10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(feat.search_face_keys(At::Z(0.0), any_x).unwrap(), &[(0, 1, 2, 3)]);
        assert_eq!(feat.search_face_keys(At::Z(2.0), any_x).unwrap(), &[(8, 9, 10, 11)]);
        assert_eq!(feat.search_face_keys(At::Z(2.0), |x| x[0] <= -1.0).ok(), None);
        assert_eq!(
            feat.search_face_keys(At::Z(10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::XY(0.0, 0.0), any_x).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::XY(1.0, 1.0), any_x).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::XY(10.0, 10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::YZ(0.0, 0.0), any_x).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::YZ(1.0, 1.0), any_x).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::YZ(10.0, 10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::XZ(0.0, 0.0), any_x).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::XZ(1.0, 0.0), any_x).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::XZ(1.0, 2.0), any_x).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::XZ(10.0, 10.0), any_x).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::XYZ(0.0, 0.0, 0.0), any_x).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::XYZ(10.0, 0.0, 0.0), any_x).err(),
            Some("cannot find point because the coordinates are outside the grid")
        );
        assert_eq!(
            feat.search_face_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0), any_x)
                .err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, SQRT_2), any_x)
                .err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            feat.search_face_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 10.0), any_x)
                .err(),
            Some("cannot find any point with given constraints/filter")
        );

        // high-level function
        let res = feat.search_faces(At::Z(0.0), any_x).unwrap();
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].points, &[0, 3, 2, 1]);
        assert_eq!(
            feat.search_faces(At::Circle(0.0, 0.0, 1.0), any_x).err(),
            Some("At::Circle works in 2D only")
        );
    }

    #[test]
    fn search_works_with_ring() {
        // 2.0   14---36--,__11
        //        |          `,-..33
        // 1.75  24   [7]   22     `-,
        //        |         ,  [5]    ,8.
        // 1.5   13--35--10/        20   `*
        //        |       ,`*32    ,'      30
        // 1.25  23 [6] 21     *.7     [3]   *
        //        |     ,  [4]  , *.          5
        // 1.0   12-34-9      19    29     18' *
        //              `31. ,' [2]   *  _,     *
        //                  6.       _.4'        *
        //                   28  _.17   *   [1]  27
        //                     3'  [0]  26        *
        //                     25        *        *
        //        +             0---15---1---16---2
        //
        //                     1.0 1.25  1.5 1.75  2.0
        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        let feat = Features::new(&mesh, Extract::Boundary);
        let (r, rr) = (1.0, 2.0);
        assert_eq!(feat.search_point_ids(At::XY(1.00, 0.00), any_x).unwrap(), &[0]);
        assert_eq!(feat.search_point_ids(At::XY(1.25, 0.00), any_x).unwrap(), &[15]);
        assert_eq!(feat.search_point_ids(At::XY(1.50, 0.00), any_x).unwrap(), &[1]);
        assert_eq!(feat.search_point_ids(At::XY(1.75, 0.00), any_x).unwrap(), &[16]);
        assert_eq!(feat.search_point_ids(At::XY(2.00, 0.00), any_x).unwrap(), &[2]);
        assert_eq!(feat.search_point_ids(At::XY(0.00, 1.00), any_x).unwrap(), &[12]);
        assert_eq!(feat.search_point_ids(At::XY(0.00, 1.25), any_x).unwrap(), &[23]);
        assert_eq!(feat.search_point_ids(At::XY(0.00, 1.75), any_x).unwrap(), &[24]);
        assert_eq!(feat.search_point_ids(At::XY(0.00, 1.50), any_x).unwrap(), &[13]);
        assert_eq!(feat.search_point_ids(At::XY(0.00, 2.00), any_x).unwrap(), &[14]);
        assert_eq!(
            feat.search_point_ids(At::XY(SQRT_2 / 2.0, SQRT_2 / 2.0), any_x)
                .unwrap(),
            &[6]
        );
        assert_eq!(feat.search_point_ids(At::XY(SQRT_2, SQRT_2), any_x).unwrap(), &[8]);
        assert_eq!(
            feat.search_point_ids(At::Circle(0.0, 0.0, r), any_x).unwrap(),
            &[0, 3, 6, 9, 12, 25, 28, 31, 34],
        );
        assert_eq!(
            feat.search_point_ids(At::Circle(0.0, 0.0, r), |x| {
                let alpha = f64::atan2(x[1], x[0]) * 180.0 / PI;
                f64::abs(alpha - 45.0) < 1e-15
            })
            .unwrap(),
            &[6],
        );
        assert_eq!(
            feat.search_point_ids(At::Circle(0.0, 0.0, rr), any_x).unwrap(),
            &[2, 5, 8, 11, 14, 27, 30, 33, 36],
        );
        assert_eq!(feat.search_edge_keys(At::Y(0.0), any_x).unwrap(), &[(0, 1), (1, 2)]);
        assert_eq!(feat.search_edge_keys(At::X(0.0), any_x).unwrap(), &[(12, 13), (13, 14)]);
        assert_eq!(
            feat.search_edge_keys(At::Circle(0.0, 0.0, r), any_x).unwrap(),
            &[(0, 3), (3, 6), (6, 9), (9, 12)],
        );
        assert_eq!(
            feat.search_edge_keys(At::Circle(0.0, 0.0, rr), any_x).unwrap(),
            &[(2, 5), (5, 8), (8, 11), (11, 14)],
        );
    }
}
