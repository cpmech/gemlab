use super::algorithms::{extract_all_2d_edges, extract_all_faces, extract_features_2d, extract_features_3d};
use super::{At, CellId, Mesh, PointId};
use super::{Edge, EdgeKey, Edges, MapEdge2dToCells, MapPointToEdges};
use super::{Face, FaceKey, Faces, MapFaceToCells, MapPointToFaces};
use crate::util::GridSearch;
use crate::StrError;
use russell_lab::sort2;
use std::collections::{HashMap, HashSet};

/// Holds derived mesh features such as points, edges, and faces on the mesh boundary or the interior
///
/// # Examples
///
/// ## Two-dimensions
///
/// ```
/// use gemlab::mesh::{At, Cell, Features, Mesh, Point};
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
///         marked_edges: Vec::new(),
///         marked_faces: Vec::new(),
///     };
///
///     let features = Features::new(&mesh, false);
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
/// use gemlab::mesh::{At, Cell, Features, Mesh, Point};
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
///          marked_edges: Vec::new(),
///          marked_faces: Vec::new(),
///      };
///
///     let features = Features::new(&mesh, false);
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
pub struct Features<'a> {
    /// Holds an access to the mesh
    pub mesh: &'a Mesh,

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
    /// 5. 1D cells in 2D go to the `cables` array
    pub edges: HashMap<EdgeKey, Edge>,

    /// Set of faces on the mesh boundary, interior, or both boundary and interior
    ///
    /// **Notes:**
    ///
    /// 1. A boundary face is such that it is shared by one 3D cell only (2D cells are ignored)
    /// 2. An interior face is such that it is shared by **more** than one 3D cell (2D cells are ignored)
    /// 3. 1D cells in 3D go to the `cables` array
    /// 4. 2D cells in 3D go to the `shells` array
    pub faces: HashMap<FaceKey, Face>,

    /// Holds the ids of linear cells; GeoKind::Lin with geo_ndim = 1 in 2D or 3D
    pub cables: Vec<CellId>,

    /// Holds the ids of 2D cells (geo_ndim = 2) in 3D (space_ndim = 3)
    pub shells: Vec<CellId>,

    /// The minimum coordinates of the points (space_ndim)
    pub min: Vec<f64>,

    /// The maximum coordinates of the points (space_ndim)
    pub max: Vec<f64>,

    /// Tool to quickly search points by coordinates
    pub grid: GridSearch,

    /// Maps a point id to edges sharing the point
    pub point_to_edges: MapPointToEdges,

    /// Maps a point id to faces sharing the point
    pub point_to_faces: MapPointToFaces,
}

impl<'a> Features<'a> {
    /// Extracts features
    ///
    /// # Input
    ///
    /// * `mesh` -- the mesh
    /// * `extract_all` -- if true, will extract the boundary and the interior features.
    ///    Otherwise, if false, will extract the **boundary only**
    ///
    /// # Notes
    ///
    /// * You may want to call [Mesh::check_all()] to capture (some) errors of the mesh first
    ///
    /// # Panics
    ///
    /// * This function will panic if the mesh data is invalid. For instance, when
    ///   the cell points array doesn't contain enough points or the indices are incorrect
    pub fn new(mesh: &'a Mesh, extract_all: bool) -> Self {
        // options
        assert!(mesh.ndim >= 2 && mesh.ndim <= 3);

        // define variables
        let all_2d_edges: MapEdge2dToCells;
        let all_faces: MapFaceToCells;
        let mut points: HashSet<PointId>;
        let edges: HashMap<EdgeKey, Edge>;
        let faces: HashMap<FaceKey, Face>;
        let mut min: Vec<f64>;
        let mut max: Vec<f64>;

        // extract features
        if mesh.ndim == 2 {
            all_2d_edges = extract_all_2d_edges(mesh);
            all_faces = HashMap::new();
            faces = HashMap::new();
            (points, edges, min, max) = extract_features_2d(mesh, &all_2d_edges, extract_all);
        } else {
            all_2d_edges = HashMap::new();
            all_faces = extract_all_faces(mesh);
            (points, edges, faces, min, max) = extract_features_3d(mesh, &all_faces, extract_all);
        };

        // handle rods and shells
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

        // CABLE cells
        let cables: Vec<_> = mesh
            .cells
            .iter()
            .filter_map(|cell| if cell.kind.ndim() == 1 { Some(cell.id) } else { None })
            .collect();

        // SHELL cells
        let shells = if mesh.ndim == 3 {
            mesh.cells
                .iter()
                .filter_map(|cell| if cell.kind.ndim() == 2 { Some(cell.id) } else { None })
                .collect()
        } else {
            Vec::new()
        };

        // results
        Features {
            mesh,
            all_2d_edges,
            all_faces,
            points,
            edges,
            faces,
            cables,
            shells,
            min,
            max,
            grid,
            point_to_edges,
            point_to_faces,
        }
    }

    /// Returns an edge or panics
    pub fn get_edge(&self, a: usize, b: usize) -> &Edge {
        self.edges.get(&(a, b)).expect("cannot find edge with given key")
    }

    /// Returns an edge or panics
    pub fn get_edge_by_key(&self, key: &EdgeKey) -> &Edge {
        self.edges.get(key).expect("cannot find edge with given key")
    }

    /// Returns a face or panics
    pub fn get_face(&self, a: usize, b: usize, c: usize, d: usize) -> &Face {
        self.faces.get(&(a, b, c, d)).expect("cannot find face with given key")
    }

    /// Returns cells sharing a given (2D) edge
    ///
    /// Returns a **sorted** list of Cell IDs
    pub fn get_cells_via_2d_edge(&self, edge: &Edge) -> Vec<CellId> {
        let cells = self.all_2d_edges.get(&edge.key()).expect("cannot find 2D edge");
        let mut ids: Vec<_> = cells.iter().map(|c| c.0).collect();
        ids.sort();
        ids
    }

    /// Returns cells sharing a given face
    ///
    /// Returns a **sorted** list of Cell IDs
    pub fn get_cells_via_face(&self, face: &Face) -> Vec<CellId> {
        let cells = self.all_faces.get(&face.key()).expect("cannot find face");
        let mut ids: Vec<_> = cells.iter().map(|c| c.0).collect();
        ids.sort();
        ids
    }

    /// Returns many cells sharing a given (2D) edge
    ///
    /// Returns a **sorted** list of Cell IDs
    pub fn get_cells_via_2d_edges(&self, edges: &Edges) -> Vec<CellId> {
        let mut ids: Vec<_> = edges.all.iter().flat_map(|e| self.get_cells_via_2d_edge(e)).collect();
        ids.sort();
        ids
    }

    /// Returns many cells sharing a given face
    ///
    /// Returns a **sorted** list of Cell IDs
    pub fn get_cells_via_faces(&self, faces: &Faces) -> Vec<CellId> {
        let mut ids: Vec<_> = faces.all.iter().flat_map(|f| self.get_cells_via_face(f)).collect();
        ids.sort();
        ids
    }

    /// Returns all points (sorted) on a set of (2D) edges
    pub fn get_points_via_2d_edges(&self, edges: &Edges) -> Vec<PointId> {
        let mut point_ids: Vec<_> = edges.all.iter().flat_map(|e| e.points.clone()).collect();
        point_ids.sort();
        point_ids.dedup();
        point_ids
    }

    /// Returns all points (sorted) on a set of faces
    pub fn get_points_via_faces(&self, faces: &Faces) -> Vec<PointId> {
        let mut point_ids: Vec<_> = faces.all.iter().flat_map(|f| f.points.clone()).collect();
        point_ids.sort();
        point_ids.dedup();
        point_ids
    }

    /// Returns all neighbors of a 2D cell
    ///
    /// # Input
    ///
    /// * `mesh` -- the mesh
    /// * `edges` -- the edge-to-cells map
    /// * `e` -- the index of the edge (see [edge_node_id](crate::shapes::GeoKind::edge_node_id))
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
    /// use gemlab::mesh::{Cell, Features, Mesh, Point};
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
    ///         marked_edges: Vec::new(),
    ///         marked_faces: Vec::new(),
    ///     };
    ///
    ///     let features = Features::new(&mesh, true);
    ///
    ///     let neighbors = features.get_neighbors_2d(0);
    ///     assert_eq!(neighbors.len(), 1);
    ///     assert!(neighbors.contains(&(1, 1, 3)));
    ///
    ///     let neighbors = features.get_neighbors_2d(1);
    ///     assert_eq!(neighbors.len(), 1);
    ///     assert!(neighbors.contains(&(3, 0, 1)));
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn get_neighbors_2d(&self, cell_id: CellId) -> Vec<(usize, CellId, usize)> {
        assert_eq!(self.mesh.ndim, 2);
        let cell = &self.mesh.cells[cell_id];
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

    /// Extracts all boundary edges
    pub fn get_boundary_edges(&self) -> Vec<EdgeKey> {
        if self.mesh.ndim == 2 {
            let mut boundary_edges = Vec::new();
            for (edge_key, cells) in &self.all_2d_edges {
                let shared_by_ncell = cells.len();
                if shared_by_ncell == 1 {
                    boundary_edges.push(*edge_key);
                }
            }
            boundary_edges
        } else {
            let mut edges_keys = HashSet::new();
            for (face_key, face) in &self.faces {
                let shared_by_ncell = self.all_faces.get(face_key).unwrap().len();
                if shared_by_ncell == 1 {
                    for e in 0..face.kind.nedge() {
                        let p0 = face.points[face.kind.edge_node_id(e, 0)];
                        let p1 = face.points[face.kind.edge_node_id(e, 1)];
                        let mut edge_key = (self.mesh.points[p0].id, self.mesh.points[p1].id);
                        sort2(&mut edge_key);
                        edges_keys.insert(edge_key);
                    }
                }
            }
            edges_keys.into_iter().collect()
        }
    }

    /// Triangulates the 3D boundary faces into triangles
    ///
    /// Note: the order of the output arrays is non-deterministic because it depends
    /// on the HashMap structures used to extract the boundary features.
    ///
    /// # Output
    ///
    /// * `xx` -- x coordinates of the triangle points
    /// * `yy` -- y coordinates of the triangle points
    /// * `zz` -- z coordinates of the triangle points
    /// * `triangles` -- indices of the triangle points
    ///
    /// # Panics
    ///
    /// * Panics if the mesh is not 3D
    pub fn triangulate_3d_boundary(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<Vec<PointId>>) {
        assert_eq!(self.mesh.ndim, 3);
        let mut old_point_id_to_new_point_id = HashMap::new();
        let mut xx: Vec<f64> = Vec::new();
        let mut yy: Vec<f64> = Vec::new();
        let mut zz: Vec<f64> = Vec::new();
        let mut triangles: Vec<Vec<PointId>> = Vec::new();
        for (face_key, face) in &self.faces {
            let shared_by_ncell = self.all_faces.get(face_key).unwrap().len();
            if shared_by_ncell == 1 {
                // this is a boundary face
                let face_npoint = face.points.len();
                let mut new_face_points = Vec::with_capacity(face_npoint);
                for p in &face.points {
                    let new_point_id = old_point_id_to_new_point_id
                        .entry(*p)
                        .or_insert_with(|| {
                            let npoint_new = xx.len();
                            xx.push(self.mesh.points[*p].coords[0]);
                            yy.push(self.mesh.points[*p].coords[1]);
                            zz.push(self.mesh.points[*p].coords[2]);
                            npoint_new
                        })
                        .clone();
                    new_face_points.push(new_point_id);
                }
                triangles.push(vec![new_face_points[0], new_face_points[1], new_face_points[2]]);
                if face_npoint > 3 {
                    triangles.push(vec![new_face_points[2], new_face_points[3], new_face_points[0]]);
                }
            }
        }
        (xx, yy, zz, triangles)
    }

    /// Searches edges with a given marker
    ///
    /// Returns edges sorted by their edge keys
    pub fn search_marked_edges(&self, marker: i32) -> Edges {
        let mut all: Vec<_> = self.edges.values().filter(|edge| edge.marker == marker).collect();
        all.sort_by_key(|edge| edge.key());
        Edges { all }
    }

    /// Searches faces with a given marker
    ///
    /// Returns faces sorted by their face keys
    pub fn search_marked_faces(&self, marker: i32) -> Faces {
        let mut all: Vec<_> = self.faces.values().filter(|face| face.marker == marker).collect();
        all.sort_by_key(|face| face.key());
        Faces { all }
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
        F: FnMut(&[f64]) -> bool,
    {
        let mut point_ids: HashSet<PointId> = HashSet::new();
        match at {
            At::X(x) => {
                if self.mesh.ndim == 2 {
                    for id in self.grid.search_on_line(&[x, 0.0], &[x, 1.0], filter)? {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid.search_on_plane_yz(x, filter)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::Y(y) => {
                if self.mesh.ndim == 2 {
                    for id in self.grid.search_on_line(&[0.0, y], &[1.0, y], filter)? {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid.search_on_plane_xz(y, filter)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::Z(z) => {
                if self.mesh.ndim == 2 {
                    return Err("At::Z works in 3D only");
                } else {
                    for id in self.grid.search_on_plane_xy(z, filter)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::XY(x, y) => {
                if self.mesh.ndim == 2 {
                    if let Some(id) = self.grid.search(&[x, y])? {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid.search_on_line(&[x, y, 0.0], &[x, y, 1.0], filter)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::YZ(y, z) => {
                if self.mesh.ndim == 2 {
                    return Err("At::YZ works in 3D only");
                } else {
                    for id in self.grid.search_on_line(&[0.0, y, z], &[1.0, y, z], filter)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::XZ(x, z) => {
                if self.mesh.ndim == 2 {
                    return Err("At::XZ works in 3D only");
                } else {
                    for id in self.grid.search_on_line(&[x, 0.0, z], &[x, 1.0, z], filter)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::XYZ(x, y, z) => {
                if self.mesh.ndim == 2 {
                    return Err("At::XYZ works in 3D only");
                } else {
                    if let Some(id) = self.grid.search(&[x, y, z])? {
                        point_ids.insert(id);
                    }
                }
            }
            At::Circle(x, y, r) => {
                if self.mesh.ndim == 2 {
                    for id in self.grid.search_on_circle(&[x, y], r, filter)? {
                        point_ids.insert(id);
                    }
                } else {
                    return Err("At::Circle works in 2D only");
                }
            }
            At::Cylinder(ax, ay, az, bx, by, bz, r) => {
                if self.mesh.ndim == 2 {
                    return Err("At::Cylinder works in 3D only");
                } else {
                    for id in self.grid.search_on_cylinder(&[ax, ay, az], &[bx, by, bz], r, filter)? {
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
        F: FnMut(&[f64]) -> bool,
    {
        let mut edge_keys: HashSet<EdgeKey> = HashSet::new();
        // search all points constrained by "at" and "filter"
        let point_ids = self.search_point_ids(at, filter)?;
        for point_id in &point_ids {
            // select all edges connected to the found points
            if let Some(edges) = self.point_to_edges.get(point_id) {
                for edge_key in edges {
                    // accept edge when at least two edge points validate "At"
                    if point_ids.contains(&edge_key.0) && point_ids.contains(&edge_key.1) {
                        edge_keys.insert(*edge_key);
                    }
                }
                // the None branch means that the point is not attached to any edge; i.e., a CABLE or SHELL point
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
        F: FnMut(&[f64]) -> bool,
    {
        if self.mesh.ndim != 3 {
            return Err("cannot find face keys in 2D");
        }
        let mut face_keys: HashSet<FaceKey> = HashSet::new();
        // search all points constrained by "at" and "filter"
        let point_ids = self.search_point_ids(at, filter)?;
        for point_id in &point_ids {
            // select all faces connected to the found points
            if let Some(faces) = self.point_to_faces.get(point_id) {
                for face_key in faces {
                    // accept face when at least four face points validate "At"
                    let fourth_is_ok = if face_key.3 == usize::MAX {
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
                // the None branch means that the point is not attached to any face; i.e., a SHELL point
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
    pub fn search_edges<F>(&self, at: At, filter: F) -> Result<Edges, StrError>
    where
        F: FnMut(&[f64]) -> bool,
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
        match results {
            Ok(all) => Ok(Edges { all }),
            Err(e) => Err(e),
        }
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
    pub fn search_faces<F>(&self, at: At, filter: F) -> Result<Faces, StrError>
    where
        F: FnMut(&[f64]) -> bool,
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
        match results {
            Ok(all) => Ok(Faces { all }),
            Err(e) => Err(e),
        }
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
    pub fn search_many_edges<F>(&self, ats: &[At], mut filter: F) -> Result<Edges, StrError>
    where
        F: FnMut(&[f64]) -> bool,
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
        match results {
            Ok(all) => Ok(Edges { all }),
            Err(e) => Err(e),
        }
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
    pub fn search_many_faces<F>(&self, ats: &[At], mut filter: F) -> Result<Faces, StrError>
    where
        F: FnMut(&[f64]) -> bool,
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
        match results {
            Ok(all) => Ok(Faces { all }),
            Err(e) => Err(e),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Edge, Edges, Face, Faces, Features};
    use crate::mesh::{At, Samples};
    use crate::shapes::GeoKind;
    use crate::util::any_x;
    use plotpy::Plot;
    use russell_lab::math::{PI, SQRT_2};

    const SAVE_FIGURE: bool = false;

    #[test]
    fn new_and_basic_get_methods_work_2d() {
        //      y
        //      ^
        // 1.0  3------6------2
        //      |             |    [#] indicates id
        //      |             |    (#) indicates attribute
        //      7     [0]     5
        //      |     (1)     |
        //      |             |
        // 0.0  0------4------1 -> x
        //     0.0           1.0
        let mesh = Samples::one_qua8();
        let feat = Features::new(&mesh, false);
        let edge = feat.get_edge(2, 3);
        assert_eq!(edge.points, &[3, 2, 6]);
        assert_eq!(edge.key(), (2, 3));
    }

    #[test]
    fn new_and_basic_get_methods_work_3d_1() {
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
        let feat = Features::new(&mesh, false);
        let edge = feat.get_edge(4, 5);
        let face = feat.get_face(0, 1, 4, 5);
        assert_eq!(edge.points, &[4, 5]);
        assert_eq!(face.points, &[0, 1, 5, 4]);
        assert_eq!(edge.key(), (4, 5));
        assert_eq!(face.key(), (0, 1, 4, 5));
    }

    #[test]
    fn new_and_basic_get_methods_work_3d_2() {
        let mesh = Samples::one_tet4();
        let feat = Features::new(&mesh, false);
        let edge = feat.get_edge(1, 2);
        let face = feat.get_face(0, 1, 3, usize::MAX);
        assert_eq!(edge.points, &[1, 2]);
        assert_eq!(face.points, &[0, 1, 3]);
        assert_eq!(edge.key(), (1, 2));
        assert_eq!(face.key(), (0, 1, 3, usize::MAX));
        assert_eq!(feat.get_cells_via_face(&face), &[0]);
    }

    #[test]
    #[should_panic(expected = "cannot find edge with given key")]
    fn get_edge_panics_on_notfound_edge() {
        //                       4------------7-----------10
        //                      /.           /|            |
        //                     / .          / |            |
        //                    /  .         /  |            |
        //                   /   .        /   |            |
        //                  5------------6    |            |
        //                  |    .       |`.  |            |
        //                  |    0-------|--`.3------------9
        //                  |   /        |   /`.          /
        //                  |  /         |  /   `.       /
        //                  | /          | /      `.    /
        //                  |/           |/         `. /
        //  12-----11-------1------------2------------8
        let mesh = Samples::mixed_shapes_3d();
        let feat = Features::new(&mesh, true);
        feat.get_edge(7, 10); // hanging edge, i.e., CABLE
    }

    #[test]
    #[should_panic(expected = "cannot find face with given key")]
    fn get_face_panics_on_notfound_face() {
        //                       4------------7-----------10
        //                      /.           /|            |
        //                     / .          / |            |
        //                    /  .         /  |            |
        //                   /   .        /   |            |
        //                  5------------6    |            |
        //                  |    .       |`.  |            |
        //                  |    0-------|--`.3------------9
        //                  |   /        |   /`.          /
        //                  |  /         |  /   `.       /
        //                  | /          | /      `.    /
        //                  |/           |/         `. /
        //  12-----11-------1------------2------------8
        let mesh = Samples::mixed_shapes_3d();
        let feat = Features::new(&mesh, true);
        feat.get_face(3, 7, 9, 10); // hanging face, i.e., SHELL
    }

    #[test]
    fn new_method_allocates_grid_2d() {
        let mesh = Samples::two_qua4();
        let feat = Features::new(&mesh, false);
        println!("{}", feat.grid);
        assert_eq!(
            format!("{}", feat.grid),
            "0: [0]\n\
             49: [1]\n\
             98: [4]\n\
             4900: [3]\n\
             4949: [2]\n\
             4998: [5]\n\
             ids = [0, 1, 2, 3, 4, 5]\n\
             nitem = 6\n\
             ncontainer = 6\n\
             ndiv = [100, 50]\n"
        );
        if SAVE_FIGURE {
            let mut plot = Plot::new();
            feat.grid.draw(&mut plot, false).unwrap();
            plot.set_equal_axes(true).set_figure_size_points(800.0, 400.0);
            plot.save("/tmp/gemlab/test_features_grid_two_qua4.svg").unwrap();
        }
    }

    #[test]
    fn new_method_allocates_grid_3d() {
        let mesh = Samples::two_hex8();
        let feat = Features::new(&mesh, false);
        println!("{}", feat.grid);
        assert_eq!(
            format!("{}", feat.grid),
            "0: [0]\n\
             49: [1]\n\
             2450: [3]\n\
             2499: [2]\n\
             122500: [4]\n\
             122549: [5]\n\
             124950: [7]\n\
             124999: [6]\n\
             245000: [8]\n\
             245049: [9]\n\
             247450: [11]\n\
             247499: [10]\n\
             ids = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]\n\
             nitem = 12\n\
             ncontainer = 12\n\
             ndiv = [50, 50, 100]\n"
        );
    }

    #[test]
    #[should_panic(expected = "cannot find edge with given key")]
    fn get_edge_panics_on_error() {
        let mesh = Samples::one_tri3();
        let features = Features::new(&mesh, false);
        features.get_edge(4, 5);
    }

    #[test]
    #[should_panic(expected = "cannot find face with given key")]
    fn get_face_panics_on_error() {
        let mesh = Samples::one_tet4();
        let features = Features::new(&mesh, false);
        features.get_face(4, 5, 6, 100);
    }

    #[test]
    fn derive_works() {
        let edge = Edge {
            kind: GeoKind::Lin3,
            points: vec![10, 20, 33],
            marker: 0,
        };
        let face = Face {
            kind: GeoKind::Qua4,
            points: vec![1, 2, 3, 4],
            marker: 0,
        };
        let edge_clone = edge.clone();
        let face_clone = face.clone();
        assert_eq!(
            format!("{:?}", edge),
            "Edge { kind: Lin3, points: [10, 20, 33], marker: 0 }"
        );
        assert_eq!(
            format!("{:?}", face),
            "Face { kind: Qua4, points: [1, 2, 3, 4], marker: 0 }"
        );
        assert_eq!(edge_clone.points.len(), 3);
        assert_eq!(face_clone.points.len(), 4);
        assert_eq!(edge.key(), (10, 20));
        assert_eq!(face.key(), (1, 2, 3, 4));
        let edges = Edges { all: vec![&edge] };
        let faces = Faces { all: vec![&face] };
        let edges_clone = edges.clone();
        let faces_clone = faces.clone();
        assert_eq!(
            format!("{:?}", edges_clone),
            "Edges { all: [Edge { kind: Lin3, points: [10, 20, 33], marker: 0 }] }"
        );
        assert_eq!(
            format!("{:?}", faces_clone),
            "Faces { all: [Face { kind: Qua4, points: [1, 2, 3, 4], marker: 0 }] }"
        );
    }

    #[test]
    fn neighbors_are_correct_2d() {
        let mesh = Samples::block_2d_four_qua12().clone();
        let features = Features::new(&mesh, true);
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
        let features = Features::new(&mesh, true);

        let neighbors = features.get_neighbors_2d(0);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(1, 1, 3)));
        assert!(neighbors.contains(&(2, 2, 0)));

        let neighbors = features.get_neighbors_2d(1);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(3, 0, 1)));
        assert!(neighbors.contains(&(2, 3, 0)));

        let neighbors = features.get_neighbors_2d(2);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(0, 0, 2)));
        assert!(neighbors.contains(&(1, 3, 3)));

        let neighbors = features.get_neighbors_2d(3);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(0, 1, 2)));
        assert!(neighbors.contains(&(3, 2, 1)));
    }

    #[test]
    fn get_boundary_edges_2d_works() {
        let mesh = Samples::two_qua4();
        let features = Features::new(&mesh, true);
        let mut edge_keys = features.get_boundary_edges();
        edge_keys.sort();
        assert_eq!(edge_keys, vec![(0, 1), (0, 3), (1, 4), (2, 3), (2, 5), (4, 5)]);

        let mesh = Samples::three_tri6_arrow();
        let features = Features::new(&mesh, true);
        let mut edge_keys = features.get_boundary_edges();
        edge_keys.sort();
        assert_eq!(edge_keys, vec![(0, 1), (0, 6), (1, 6)]);
    }

    #[test]
    fn get_boundary_edges_3d_works() {
        let mesh = Samples::two_hex8();
        let features = Features::new(&mesh, true);
        let mut edge_keys = features.get_boundary_edges();
        edge_keys.sort();
        assert_eq!(
            edge_keys,
            vec![
                (0, 1),   //  0
                (0, 3),   //  1
                (0, 4),   //  2
                (1, 2),   //  3
                (1, 5),   //  4
                (2, 3),   //  5
                (2, 6),   //  6
                (3, 7),   //  7
                (4, 5),   //  8
                (4, 7),   //  9
                (4, 8),   // 10
                (5, 6),   // 11
                (5, 9),   // 12
                (6, 7),   // 13
                (6, 10),  // 14
                (7, 11),  // 15
                (8, 9),   // 16
                (8, 11),  // 17
                (9, 10),  // 18
                (10, 11), // 19
            ]
        );

        let mesh = Samples::four_hex8();
        let features = Features::new(&mesh, true);
        let mut edge_keys = features.get_boundary_edges();
        edge_keys.sort();
        assert_eq!(
            edge_keys,
            vec![
                (0, 1),   //  0
                (0, 3),   //  1
                (0, 4),   //  2
                (1, 2),   //  3
                (1, 5),   //  4
                (2, 3),   //  5
                (2, 6),   //  6
                (2, 12),  //  7
                (3, 7),   //  8
                (3, 13),  //  9
                (4, 5),   // 10
                (4, 7),   // 11
                (4, 8),   // 12
                (5, 6),   // 13
                (5, 9),   // 14
                (6, 10),  // 15
                (6, 14),  // 16
                (7, 11),  // 17
                (7, 15),  // 18
                (8, 9),   // 19
                (8, 11),  // 20
                (9, 10),  // 21
                (10, 11), // 22
                (10, 16), // 23
                (11, 17), // 24
                (12, 13), // 25
                (12, 14), // 26
                (13, 15), // 27
                (14, 15), // 28
                (14, 16), // 29
                (15, 17), // 30
                (16, 17), // 32
            ]
        )
    }

    #[test]
    fn triangulate_3d_boundary_works() {
        let key = |x, y, z| format!("{:.1}, {:.1}, {:.1}", x, y, z);
        let mesh = Samples::two_hex8();
        let features = Features::new(&mesh, true);
        let (xx, yy, zz, triangles) = features.triangulate_3d_boundary();
        let pp = vec![
            "0.0, 0.0, 0.0", //  0
            "1.0, 0.0, 0.0", //  1
            "1.0, 1.0, 0.0", //  2
            "0.0, 1.0, 0.0", //  3
            "0.0, 0.0, 1.0", //  4
            "1.0, 0.0, 1.0", //  5
            "1.0, 1.0, 1.0", //  6
            "0.0, 1.0, 1.0", //  7
            "0.0, 0.0, 2.0", //  8
            "1.0, 0.0, 2.0", //  9
            "1.0, 1.0, 2.0", // 10
            "0.0, 1.0, 2.0", // 11
        ];
        let boundary_faces = vec![
            (pp[0], pp[4], pp[7], pp[3]),   // lower, -x face
            (pp[1], pp[2], pp[6], pp[5]),   // lower, +x face
            (pp[0], pp[1], pp[5], pp[4]),   // lower, -y face
            (pp[2], pp[3], pp[7], pp[6]),   // lower, +y face
            (pp[4], pp[8], pp[11], pp[7]),  // upper, -x face
            (pp[5], pp[6], pp[10], pp[9]),  // upper, +x face
            (pp[4], pp[5], pp[9], pp[8]),   // upper, -y face
            (pp[6], pp[7], pp[11], pp[10]), // upper, +y face
            (pp[0], pp[3], pp[2], pp[1]),   // lower, -z face
            (pp[8], pp[9], pp[10], pp[11]), // upper, +z face
        ];
        assert_eq!(xx.len(), 12);
        assert_eq!(yy.len(), 12);
        assert_eq!(zz.len(), 12);
        assert_eq!(triangles.len(), boundary_faces.len() * 2); // quads => tris
        for i in 0..xx.len() {
            assert!(pp.contains(&key(xx[i], yy[i], zz[i]).as_str()));
        }
        let mut connectivity = Vec::new();
        for tri in &triangles {
            let a = key(xx[tri[0]], yy[tri[0]], zz[tri[0]]);
            let b = key(xx[tri[1]], yy[tri[1]], zz[tri[1]]);
            let c = key(xx[tri[2]], yy[tri[2]], zz[tri[2]]);
            connectivity.push((a, b, c));
        }
        // expected triangles
        for (k0, k1, k2, k3) in &boundary_faces {
            assert!(connectivity.contains(&(k0.to_string(), k1.to_string(), k2.to_string())));
            assert!(connectivity.contains(&(k2.to_string(), k3.to_string(), k0.to_string())));
        }
    }

    #[test]
    fn search_marked_points_edges_faces_work() {
        let mesh = Samples::two_hex8();
        let features = Features::new(&mesh, false);

        let res = features.search_marked_edges(-4);
        assert_eq!(res.all.iter().map(|e| e.key()).collect::<Vec<_>>(), vec![(3, 7)]);

        let res = features.search_marked_edges(-5);
        assert_eq!(
            res.all.iter().map(|e| e.key()).collect::<Vec<_>>(),
            vec![(8, 9), (10, 11)]
        );

        let res = features.search_marked_faces(-8);
        assert_eq!(res.all.iter().map(|f| f.key()).collect::<Vec<_>>(), vec![(2, 3, 6, 7)]);

        let res = features.search_marked_faces(-9);
        assert_eq!(
            res.all.iter().map(|f| f.key()).collect::<Vec<_>>(),
            vec![(8, 9, 10, 11)]
        );
    }

    #[test]
    fn search_points_fails_on_wrong_input() {
        // any_x works
        assert_eq!(any_x(&vec![]), true);

        // 2d
        let mesh = Samples::two_qua4();
        let feat = Features::new(&mesh, false);
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
        let boundary = Features::new(&mesh, false);
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
        let feat = Features::new(&mesh, true);
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
        let feat = Features::new(&mesh, false);
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
    fn search_points_works_3d_mixed() {
        //                       4------------7-----------10
        //                      /.           /|            |
        //                     / .          / |            |
        //                    /  .         /  |            |
        //                   /   .        /   |            |
        //         z        5------------6    |            |
        //         ↑        |    .       |`.  |            |
        //         o → y    |    0-------|--`.3------------9
        //        ↙         |   /        |   /`.          /
        //      x           |  /         |  /   `.       /
        //                  | /          | /      `.    /
        //                  |/           |/         `. /
        //  12-----11-------1------------2------------8
        let mesh = Samples::mixed_shapes_3d();
        let feat = Features::new(&mesh, true);
        assert_eq!(feat.search_point_ids(At::X(0.0), any_x).unwrap(), &[0, 3, 4, 7, 9, 10]);
        assert_eq!(feat.search_point_ids(At::Y(1.0), any_x).unwrap(), &[2, 3, 6, 7]);
        assert_eq!(feat.search_point_ids(At::XZ(0.0, 1.0), any_x).unwrap(), &[4, 7, 10]);
        assert_eq!(
            feat.search_point_ids(At::XZ(1.0, 0.0), any_x).unwrap(),
            &[1, 2, 8, 11, 12]
        );
        assert_eq!(
            feat.search_point_ids(At::Z(0.0), any_x).unwrap(),
            &[0, 1, 2, 3, 8, 9, 11, 12]
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
        let feat = Features::new(&mesh, false);
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
        assert_eq!(res.all.len(), 2);
        assert_eq!(res.all[0].points, &[1, 0]);
        assert_eq!(res.all[1].points, &[4, 1]);
        assert_eq!(
            feat.search_edges(At::XYZ(0.0, 0.0, 0.0), any_x).err(),
            Some("At::XYZ works in 3D only")
        );

        // many edges
        let res = feat
            .search_many_edges(&[At::X(0.0), At::X(2.0), At::Y(0.0), At::Y(1.0)], any_x)
            .unwrap();
        assert_eq!(res.all.len(), 6);
        let keys: Vec<_> = res.all.iter().map(|r| (r.points[0], r.points[1])).collect();
        assert_eq!(keys, &[(1, 0), (0, 3), (4, 1), (3, 2), (2, 5), (5, 4)]);
    }

    #[test]
    fn search_edges_works_2d_mixed() {
        // 1.0              4-----------3
        //                  |           |
        //                  |    [1]    |   [*] indicates id
        //                  |    (2)    |   (*) indicates attribute
        //                  |           |
        // 0.0  0-----------1-----------2-----------5
        //           [0]                     [2]
        //           (1)                     (1)
        let mesh = Samples::mixed_shapes_2d();
        let feat = Features::new(&mesh, false);
        assert_eq!(feat.cables, &[0, 2]);
        // note that (0,1) and (2,5) are NOT edge; they're CABLE
        assert_eq!(feat.search_edge_keys(At::Y(0.0), any_x).unwrap(), &[(1, 2)]);
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
        let feat = Features::new(&mesh, false);
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
        assert_eq!(res.all.len(), 2);
        assert_eq!(res.all[0].points, &[0, 4]);
        assert_eq!(res.all[1].points, &[4, 8]);

        // many faces
        let res = feat.search_many_faces(&[At::Z(0.0), At::Z(2.0)], any_x).unwrap();
        assert_eq!(res.all.len(), 2);
        let keys: Vec<_> = res
            .all
            .iter()
            .map(|r| (r.points[0], r.points[1], r.points[2], r.points[3]))
            .collect();
        assert_eq!(keys, &[(0, 3, 2, 1), (8, 9, 10, 11)]);
    }

    #[test]
    fn search_edges_works_3d_mixed() {
        //                       4------------7-----------10
        //                      /.           /|            |
        //                     / .          / |            |
        //                    /  .         /  |            |
        //                   /   .        /   |            |
        //         z        5------------6    |            |
        //         ↑        |    .       |`.  |            |
        //         o → y    |    0-------|--`.3------------9
        //        ↙         |   /        |   /`.          /
        //      x           |  /         |  /   `.       /
        //                  | /          | /      `.    /
        //                  |/           |/         `. /
        //  12-----11-------1------------2------------8
        let mesh = Samples::mixed_shapes_3d();
        let feat = Features::new(&mesh, true);
        assert_eq!(feat.cables, &[4]);
        assert_eq!(feat.shells, &[2, 3]);
        // note that (11,12) and (1,11) are not edges; they are CABLE
        assert_eq!(
            feat.search_edge_keys(At::XZ(1.0, 0.0), any_x).unwrap(),
            &[(1, 2), (2, 8)]
        );
        // note that (7,10) and (9,10) are not edges because they belong to SHELL (hanging faces)
        assert_eq!(
            feat.search_edge_keys(At::X(0.0), any_x).unwrap(),
            &[(0, 3), (0, 4), (3, 7), (4, 7)]
        );
    }

    #[test]
    fn search_faces_returns_error_in_2d() {
        // 3--------2--------5
        // |        |        |
        // |        |        |
        // |        |        |
        // 0--------1--------4
        let mesh = Samples::two_qua4();
        let feat = Features::new(&mesh, false);
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
        let feat = Features::new(&mesh, false);
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
        assert_eq!(res.all.len(), 1);
        assert_eq!(res.all[0].points, &[0, 3, 2, 1]);
        assert_eq!(
            feat.search_faces(At::Circle(0.0, 0.0, 1.0), any_x).err(),
            Some("At::Circle works in 2D only")
        );
    }

    #[test]
    fn search_faces_works_mixed() {
        //                       4------------7-----------10
        //                      /.           /|            |
        //                     / .          / |            |
        //                    /  .         /  |            |
        //                   /   .        /   |            |
        //         z        5------------6    |            |
        //         ↑        |    .       |`.  |            |
        //         o → y    |    0-------|--`.3------------9
        //        ↙         |   /        |   /`.          /
        //      x           |  /         |  /   `.       /
        //                  | /          | /      `.    /
        //                  |/           |/         `. /
        //  12-----11-------1------------2------------8
        let mesh = Samples::mixed_shapes_3d();
        let feat = Features::new(&mesh, true);
        assert_eq!(feat.cables, &[4]);
        assert_eq!(feat.shells, &[2, 3]);
        // note that (3,7,9,10) is not face because it is SHELL (hanging face)
        assert_eq!(feat.search_face_keys(At::X(0.0), any_x).unwrap(), &[(0, 3, 4, 7)]);
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
        let feat = Features::new(&mesh, false);
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

    #[test]
    fn display_works_1() {
        //       8-------------11  2.0
        //      /.             /|
        //     / .            / |
        //    /  .           /  |
        //   /   .          /   |
        //  9-------------10    |
        //  |    .         |    |
        //  |    4---------|----7  1.0
        //  |   /. [1]     |   /|
        //  |  / . (2)     |  / |
        //  | /  .         | /  |
        //  |/   .         |/   |
        //  5--------------6    |          z
        //  |    .         |    |          ↑
        //  |    0---------|----3  0.0     o → y
        //  |   /  [0]     |   /          ↙
        //  |  /   (1)     |  /          x
        //  | /            | /
        //  |/             |/
        //  1--------------2   1.0
        // 0.0            1.0
        let mesh = Samples::two_hex8();
        let features = Features::new(&mesh, false);
        let edges = features.search_edges(At::Z(0.0), any_x).unwrap();
        let faces = features.search_faces(At::Y(1.0), any_x).unwrap();
        assert_eq!(format!("{}", edges), "(0, 1), (0, 3), (1, 2), (2, 3)");
        assert_eq!(format!("{}", faces), "(2, 3, 6, 7), (6, 7, 10, 11)");
    }

    #[test]
    fn display_works_2() {
        let mesh = Samples::one_tet4();
        let features = Features::new(&mesh, false);
        let edges = features.search_edges(At::Z(0.0), any_x).unwrap();
        let faces = features.search_faces(At::X(0.0), any_x).unwrap();
        assert_eq!(format!("{}", edges), "(0, 1), (0, 2), (1, 2)");
        assert_eq!(format!("{}", faces), "(0, 2, 3, MAX)");
    }

    #[test]
    fn get_cells_and_points_methods_work_2d_1() {
        // 1.0              4-----------3
        //                  |           |
        //                  |    [1]    |   [*] indicates id
        //                  |    (2)    |   (*) indicates attribute
        //                  |           |
        // 0.0  0-----------1-----------2-----------5
        //           [0]                     [2]
        //           (1)                     (1)
        let mesh = Samples::mixed_shapes_2d();
        let feat = Features::new(&mesh, false);
        let edge = feat.get_edge(1, 2);
        let edges = Edges { all: vec![&edge] };
        assert_eq!(feat.get_cells_via_2d_edge(&edge), &[1]);
        assert_eq!(feat.get_cells_via_2d_edges(&edges), &[1]);
        assert_eq!(feat.get_points_via_2d_edges(&edges), &[1, 2]);
    }

    #[test]
    fn get_cells_and_points_methods_work_2d_2() {
        // 7---------------6---------------8
        // |               |               |
        // |               |               |
        // |      [2]      |      [3]      |
        // |               |               |
        // |               |               |
        // 3---------------2---------------5
        // |               |               |
        // |               |               |
        // |      [0]      |      [1]      |
        // |               |               |
        // |               |               |
        // 0---------------1---------------4
        let mesh = Samples::block_2d_four_qua4();
        let feat = Features::new(&mesh, true); // need interior edges
        let edge_a = feat.get_edge(2, 3);
        let edge_b = feat.get_edge(2, 5);
        let edges = Edges {
            all: vec![&edge_a, &edge_b],
        };
        assert_eq!(feat.get_cells_via_2d_edge(&edge_a), &[0, 2]);
        assert_eq!(feat.get_cells_via_2d_edge(&edge_b), &[1, 3]);
        assert_eq!(feat.get_cells_via_2d_edges(&edges), &[0, 1, 2, 3]);
        assert_eq!(feat.get_points_via_2d_edges(&edges), &[2, 3, 5]);
    }

    #[test]
    fn get_cells_and_points_methods_work_3d_1() {
        //                       4------------7-----------10
        //                      /.           /|            |
        //                     / .          / |            |
        //                    /  .         /  |            |
        //                   /   .        /   |            |
        //         z        5------------6    |            |
        //         ↑        |    .       |`.  |            |
        //         o → y    |    0-------|--`.3------------9
        //        ↙         |   /        |   /`.          /
        //      x           |  /         |  /   `.       /
        //                  | /          | /      `.    /
        //                  |/           |/         `. /
        //  12-----11-------1------------2------------8
        let mesh = Samples::mixed_shapes_3d();
        let feat = Features::new(&mesh, false);
        let face_a = feat.get_face(2, 3, 6, 7);
        let face_b = feat.get_face(2, 3, 6, usize::MAX);
        let faces = Faces {
            all: vec![&face_a, &face_b],
        };
        assert_eq!(feat.get_cells_via_face(&face_a), &[0]);
        assert_eq!(feat.get_cells_via_face(&face_b), &[1]);
        assert_eq!(feat.get_cells_via_faces(&faces), &[0, 1]);
        assert_eq!(feat.get_points_via_faces(&faces), &[2, 3, 6, 7]);
    }

    #[test]
    fn get_cells_and_points_methods_work_3d_2() {
        let mesh = Samples::block_3d_eight_hex8();
        let feat = Features::new(&mesh, true); // need interior faces
        let face_a = feat.get_face(2, 3, 6, 7);
        let face_b = feat.get_face(6, 7, 20, 21);
        let faces = Faces {
            all: vec![&face_a, &face_b],
        };
        assert_eq!(feat.get_cells_via_face(&face_a), &[0, 2]);
        assert_eq!(feat.get_cells_via_face(&face_b), &[4, 6]);
        assert_eq!(feat.get_cells_via_faces(&faces), &[0, 2, 4, 6]);
        assert_eq!(feat.get_points_via_faces(&faces), &[2, 3, 6, 7, 20, 21]);
    }
}
