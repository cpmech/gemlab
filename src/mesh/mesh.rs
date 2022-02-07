use super::read_mesh::{parse_mesh, read_mesh};
use super::At;
use crate::shapes::{Shape, ShapeState};
use crate::util::GridSearch;
use crate::StrError;
use russell_lab::{sort2, sort4};
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::fmt::{self, Write};

/// Aliases usize as Point ID
pub type PointId = usize;

/// Aliases usize as Cell ID
pub type CellId = usize;

/// Aliases usize as Cell's attribute ID
pub type CellAttributeId = usize;

/// Aliases (usize,usize) as the key of Edge
///
/// # Note
///
/// Since the local numbering scheme runs over "corners" first, we can compare
/// edges using only two points; i.e., the middle points don't matter.
pub type EdgeKey = (usize, usize);

/// Aliases (usize,usize,usize,usize) as the key of Face
///
/// # Note
///
/// If all faces have at most 3 points, the fourth entry in the key will be equal to the total number of points.
/// In this way, we can compare 4-node (or more nodes) faces with each other, since that the local numbering
/// scheme runs over the "corners" first; i.e., the middle points don't matter.
pub type FaceKey = (usize, usize, usize, usize);

/// Holds point data
pub struct Point {
    /// Identification number which equals the index of the point in the mesh
    ///
    /// **raw data**
    pub id: PointId,

    /// Point coordinates (2D or 3D)
    ///
    /// **raw data**
    pub coords: Vec<f64>,

    /// Set of boundary edges sharing this point
    ///
    /// (derived property)
    pub shared_by_boundary_edges: HashSet<EdgeKey>,

    /// Set of boundary faces sharing this point
    ///
    /// (derived property)
    pub shared_by_boundary_faces: HashSet<FaceKey>,
}

/// Holds cell (aka geometric shape, polygon, polyhedra) data
pub struct Cell {
    /// Identification number which equals the index of the cell in the mesh
    ///
    /// **raw data**
    pub id: CellId,

    /// Attribute identification number
    ///
    /// **raw data**
    pub attribute_id: CellAttributeId,

    /// Space dimension of this cell
    ///
    /// The cell's ndim may be different than the space dimension of the mesh.
    /// For example, a 1D line in the 2D or 3D space or a 2D triangle in the 3D space.
    ///
    /// **raw data**
    pub geo_ndim: usize,

    /// List of points defining this cell (nodes); in the right order (unsorted)
    ///
    /// Note: The list of nodes must follow a **counter-clockwise order**.
    ///
    /// **raw data**
    pub points: Vec<PointId>,

    /// Shape object for numerical integration and other calculations
    ///
    /// (derived property)
    pub shape: Shape,
}

/// Holds edge data (derived data structure)
pub struct Edge {
    /// List of points defining this edge; in the right order (unsorted)
    pub points: Vec<PointId>,

    /// Set of 2D cells sharing this edge (to find the boundary)
    ///
    /// **2D mesh only**
    pub shared_by_2d_cells: HashSet<CellId>,

    /// Set of boundary faces sharing this edge
    pub shared_by_boundary_faces: HashSet<FaceKey>,

    /// Shape object for numerical integration and other calculations
    pub shape: Shape,
}

/// Holds face data (derived data structure)
pub struct Face {
    /// List of points defining this face; in the right order (unsorted)
    pub points: Vec<PointId>,

    /// Set of cells sharing this face
    pub shared_by_cells: HashSet<CellId>,

    /// Shape object for numerical integration and other calculations
    pub shape: Shape,
}

/// Holds mesh data
pub struct Mesh {
    /// Space dimension of the mesh
    ///
    /// The mesh's ndim may be different that an cell's ndim.
    /// For example, a 3D mesh may contain 1D lines or 2D triangles.
    ///
    /// **raw data**
    pub space_ndim: usize,

    /// All points in the mesh
    ///
    /// **raw data**
    pub points: Vec<Point>,

    /// All cells (aka geometric shape, polygon, polyhedra) in the mesh
    ///
    /// **raw data**
    pub cells: Vec<Cell>,

    /// Set of points on the boundaries
    ///
    /// Note: a boundary point belongs to a boundary edge or a boundary face
    ///
    /// (derived property)
    pub boundary_points: HashSet<PointId>,

    /// Set of edges on the boundaries
    ///
    /// Note: In 2D, a boundary edge is such that it's shared by one 2D cell only (ignore 1D cells)
    ///
    /// Note: In 3D, a boundary edge belongs to a boundary face
    ///
    /// (derived property)
    pub boundary_edges: HashMap<EdgeKey, Edge>,

    /// Set of faces on the boundaries
    ///
    /// Note: A boundary face is such that it's shared by one 3D cell only
    ///
    /// (derived property)
    pub boundary_faces: HashMap<FaceKey, Face>,

    /// Min coordinates
    ///
    /// (derived property)
    pub min: Vec<f64>,

    /// Max coordinates
    ///
    /// (derived property)
    pub max: Vec<f64>,

    /// Allows searching boundary points using their coordinates
    grid_boundary_points: GridSearch,

    /// Indicates whether the derived variables and maps have been computed or not
    derived_props_computed: bool,
}

impl Mesh {
    /// Returns a new empty mesh
    pub(super) fn new(space_ndim: usize) -> Result<Self, StrError> {
        if space_ndim < 2 || space_ndim > 3 {
            return Err("space_ndim must be 2 or 3");
        }
        Ok(Mesh {
            space_ndim,
            points: Vec::new(),
            cells: Vec::new(),
            boundary_points: HashSet::new(),
            boundary_edges: HashMap::new(),
            boundary_faces: HashMap::new(),
            min: Vec::new(),
            max: Vec::new(),
            grid_boundary_points: GridSearch::new(space_ndim)?,
            derived_props_computed: false,
        })
    }

    /// Parses raw mesh data from a text string and computes derived properties
    ///
    /// # Note
    ///
    /// This function calls `compute_derived_props` already.
    pub fn from_text(raw_mesh_data: &str) -> Result<Self, StrError> {
        let mut mesh = parse_mesh(raw_mesh_data)?;
        mesh.compute_derived_props()?;
        Ok(mesh)
    }

    /// Reads raw mesh data from text file and computes derived properties
    ///
    /// # Note
    ///
    /// This function calls `compute_derived_props` already.
    pub fn from_text_file<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let mut mesh = read_mesh(full_path)?;
        mesh.compute_derived_props()?;
        Ok(mesh)
    }

    /// Computes derived properties such as boundaries and limits
    ///
    /// # Note
    ///
    /// This function is automatically called by `from_text` and `from_text_file`.
    pub fn compute_derived_props(&mut self) -> Result<(), StrError> {
        self.derived_props_computed = false;
        self.boundary_points.clear();
        self.boundary_edges.clear();
        self.boundary_faces.clear();
        if self.space_ndim == 2 {
            self.compute_derived_props_2d()?;
        } else {
            self.compute_derived_props_3d()?;
        }
        self.compute_limits()?;
        self.grid_boundary_points
            .initialize(&vec![10; self.space_ndim], &self.min, &self.max)?;
        for point_id in &self.boundary_points {
            self.grid_boundary_points
                .insert(*point_id, &self.points[*point_id].coords)?;
        }
        self.derived_props_computed = true;
        Ok(())
    }

    /// Finds boundary points in the mesh
    ///
    /// # Input
    ///
    /// * `at` -- the location constraint
    ///
    /// # Output
    ///
    /// * Returns a **sorted** list of point ids (boundary points only)
    ///
    /// # Note
    ///
    /// `compute_derived_props` must be called first, otherwise this function returns an error
    ///
    /// # Warning
    ///
    /// This function cannot be used concurrently.
    pub fn find_boundary_points(&mut self, at: At) -> Result<Vec<PointId>, StrError> {
        /*
        NOTE: we use "&mut self" here because grid_boundary_points requires it
        */
        if !self.derived_props_computed {
            return Err("compute_derived_props must be called first");
        }
        let mut point_ids: HashSet<PointId> = HashSet::new();
        match at {
            At::X(x) => {
                if self.space_ndim == 2 {
                    for id in self.grid_boundary_points.find_on_line(&[x, 0.0], &[x, 1.0])? {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid_boundary_points.find_on_plane_yz(x)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::Y(y) => {
                if self.space_ndim == 2 {
                    for id in self.grid_boundary_points.find_on_line(&[0.0, y], &[1.0, y])? {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid_boundary_points.find_on_plane_xz(y)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::Z(z) => {
                if self.space_ndim == 2 {
                    return Err("At::Z works in 3D only");
                } else {
                    for id in self.grid_boundary_points.find_on_plane_xy(z)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::XY(x, y) => {
                if self.space_ndim == 2 {
                    if let Some(id) = self.grid_boundary_points.find(&[x, y])? {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid_boundary_points.find_on_line(&[x, y, 0.0], &[x, y, 1.0])? {
                        point_ids.insert(id);
                    }
                }
            }
            At::YZ(y, z) => {
                if self.space_ndim == 2 {
                    return Err("At::YZ works in 3D only");
                } else {
                    for id in self.grid_boundary_points.find_on_line(&[0.0, y, z], &[1.0, y, z])? {
                        point_ids.insert(id);
                    }
                }
            }
            At::XZ(x, z) => {
                if self.space_ndim == 2 {
                    return Err("At::XZ works in 3D only");
                } else {
                    for id in self.grid_boundary_points.find_on_line(&[x, 0.0, z], &[x, 1.0, z])? {
                        point_ids.insert(id);
                    }
                }
            }
            At::XYZ(x, y, z) => {
                if self.space_ndim == 2 {
                    return Err("At::XYZ works in 3D only");
                } else {
                    if let Some(id) = self.grid_boundary_points.find(&[x, y, z])? {
                        point_ids.insert(id);
                    }
                }
            }
            At::Circle(x, y, r) => {
                if self.space_ndim == 2 {
                    for id in self.grid_boundary_points.find_on_circle(&[x, y], r)? {
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
                    for id in self
                        .grid_boundary_points
                        .find_on_cylinder(&[ax, ay, az], &[bx, by, bz], r)?
                    {
                        point_ids.insert(id);
                    }
                }
            }
        }
        let mut ids: Vec<_> = point_ids.into_iter().collect();
        ids.sort();
        Ok(ids)
    }

    /// Finds boundary edges in the mesh
    ///
    /// # Input
    ///
    /// * `at` -- the location constraint
    ///
    /// # Output
    ///
    /// * Returns a **sorted** list of edge key ids (boundary edges only)
    ///
    /// # Note
    ///
    /// `compute_derived_props` must be called first, otherwise this function returns an error
    ///
    /// # Warning
    ///
    /// This function cannot be used concurrently.
    pub fn find_boundary_edges(&mut self, at: At) -> Result<Vec<EdgeKey>, StrError> {
        /*
        NOTE: we use "&mut self" here because grid_boundary_points requires it
        */
        if !self.derived_props_computed {
            return Err("compute_derived_props must be called first");
        }
        let mut edge_keys: HashSet<EdgeKey> = HashSet::new();
        // find all points near the geometric feature
        let point_ids = &self.find_boundary_points(at)?;
        for point_id in point_ids {
            // loop over all boundary edges touching this point
            let point = &self.points[*point_id];
            for edge_key in &point.shared_by_boundary_edges {
                // check if two edge points pass through the geometric feature
                if point_ids.contains(&edge_key.0) && point_ids.contains(&edge_key.1) {
                    if self.boundary_edges.contains_key(&edge_key) {
                        // update set
                        edge_keys.insert(*edge_key);
                    }
                }
            }
        }
        let mut keys: Vec<_> = edge_keys.into_iter().collect();
        keys.sort();
        Ok(keys)
    }

    /// Computes derived properties of 2D mesh
    fn compute_derived_props_2d(&mut self) -> Result<(), StrError> {
        // maps all edge keys to (cell_id, e) where e is the cell's local edge index
        let mut all_edges: HashMap<EdgeKey, Vec<(CellId, usize)>> = HashMap::new();

        // loop over all cells
        for cell in &mut self.cells {
            let nnode = cell.points.len();

            // handle 1D shapes
            if cell.geo_ndim != 2 {
                for m in 0..nnode {
                    self.boundary_points.insert(cell.points[m]);
                }
                continue; // skip 1D line in 2D because it's not a boundary edge
            }

            // set the cell node coordinates in its shape object
            for m in 0..nnode {
                for j in 0..self.space_ndim {
                    cell.shape.set_node(m, j, self.points[cell.points[m]].coords[j])?;
                }
            }

            // check if the determinant of Jacobian is positive => counterclockwise nodes
            let mut state = ShapeState::new(cell.shape.space_ndim, cell.shape.geo_ndim, cell.shape.nnode)?;
            let det_jac = cell.shape.calc_jacobian(&mut state, &[0.0, 0.0, 0.0])?;
            if det_jac < 0.0 {
                return Err("a cell has incorrect ordering of nodes");
            }

            // set information about all edges
            for e in 0..cell.shape.nedge {
                // define edge key (sorted point ids)
                let mut edge_key: EdgeKey = (
                    cell.points[cell.shape.get_edge_node_id(e, 0)],
                    cell.points[cell.shape.get_edge_node_id(e, 1)],
                );
                sort2(&mut edge_key);

                // configure edge_key => (cell_id, e)
                let edge_data = all_edges.entry(edge_key).or_insert(Vec::new());
                edge_data.push((cell.id, e));
            }
        }

        // loop over all edges
        for (edge_key, edge_data) in &all_edges {
            // skip inner edges
            if edge_data.len() != 1 {
                continue;
            }

            // edge data
            let (cell_id, e) = edge_data[0];
            let cell = &self.cells[cell_id];
            let edge_nnode = cell.shape.edge_nnode;
            let mut edge_points: Vec<PointId> = vec![0; edge_nnode];
            let mut edge_shape = Shape::new(self.space_ndim, 1, edge_nnode)?;

            // loop over all edge nodes
            for i in 0..edge_nnode {
                // configure edge nodes and shape object
                edge_points[i] = cell.points[cell.shape.get_edge_node_id(e, i)];
                for j in 0..self.space_ndim {
                    edge_shape.set_node(i, j, self.points[edge_points[i]].coords[j])?;
                }

                // set boundary points
                self.points[edge_points[i]].shared_by_boundary_edges.insert(*edge_key);
                self.boundary_points.insert(edge_points[i]);
            }

            // append new boundary edge
            self.boundary_edges.insert(
                *edge_key,
                Edge {
                    points: edge_points,
                    shared_by_2d_cells: HashSet::from([cell_id]),
                    shared_by_boundary_faces: HashSet::new(),
                    shape: edge_shape,
                },
            );
        }
        Ok(())
    }

    /// Computes derived properties of 3D mesh
    fn compute_derived_props_3d(&mut self) -> Result<(), StrError> {
        // maps all face keys to (cell_id, f) where f is the cell's local face index
        let mut all_faces: HashMap<FaceKey, Vec<(CellId, usize)>> = HashMap::new();

        // loop over all cells
        for cell in &mut self.cells {
            let nnode = cell.points.len();

            // handle 1D and 2D shapes
            if cell.geo_ndim != 3 {
                for m in 0..nnode {
                    self.boundary_points.insert(cell.points[m]);
                }
                continue; // skip 1D line or 2D shape in 3D because they don't have faces
            }

            // set the cell node coordinates in its shape object
            for m in 0..nnode {
                for j in 0..self.space_ndim {
                    cell.shape.set_node(m, j, self.points[cell.points[m]].coords[j])?;
                }
            }

            // check if the determinant of Jacobian is positive => counterclockwise nodes
            let mut state = ShapeState::new(cell.shape.space_ndim, cell.shape.geo_ndim, cell.shape.nnode)?;
            let det_jac = cell.shape.calc_jacobian(&mut state, &[0.0, 0.0, 0.0])?;
            if det_jac < 0.0 {
                return Err("a cell has incorrect ordering of nodes");
            }

            // set information about all faces
            for f in 0..cell.shape.nface {
                // define face key (sorted ids)
                let mut face_key: FaceKey = if cell.shape.face_nnode > 3 {
                    (
                        cell.points[cell.shape.get_face_node_id(f, 0)],
                        cell.points[cell.shape.get_face_node_id(f, 1)],
                        cell.points[cell.shape.get_face_node_id(f, 2)],
                        cell.points[cell.shape.get_face_node_id(f, 3)],
                    )
                } else {
                    (
                        cell.points[cell.shape.get_face_node_id(f, 0)],
                        cell.points[cell.shape.get_face_node_id(f, 1)],
                        cell.points[cell.shape.get_face_node_id(f, 2)],
                        self.points.len(),
                    )
                };
                sort4(&mut face_key);

                // configure face_key => (cell_id, f)
                let face_data = all_faces.entry(face_key).or_insert(Vec::new());
                face_data.push((cell.id, f));
            }
        }

        // sort face keys just so the next loop is deterministic
        let mut face_keys: Vec<_> = all_faces.keys().collect();
        face_keys.sort();

        // loop over all faces
        for face_key in face_keys {
            let face_data = all_faces.get(face_key).unwrap();
            // skip inner faces
            if face_data.len() != 1 {
                continue;
            }

            // face data
            let (cell_id, f) = face_data[0];
            let cell = &self.cells[cell_id];
            let face_nnode = cell.shape.face_nnode;
            let mut face_points: Vec<PointId> = vec![0; face_nnode];
            let mut face_shape = Shape::new(self.space_ndim, 2, face_nnode)?;

            // loop over all face nodes
            for i in 0..face_nnode {
                // configure face nodes and shape object
                face_points[i] = cell.points[cell.shape.get_face_node_id(f, i)];
                for j in 0..self.space_ndim {
                    face_shape.set_node(i, j, self.points[face_points[i]].coords[j])?;
                }

                // set boundary points
                self.points[face_points[i]].shared_by_boundary_faces.insert(*face_key);
                self.boundary_points.insert(face_points[i]);
            }

            // loop over all face edges
            for e in 0..face_shape.nedge {
                // define edge key (sorted point ids)
                let mut edge_key: EdgeKey = (
                    face_points[face_shape.get_edge_node_id(e, 0)],
                    face_points[face_shape.get_edge_node_id(e, 1)],
                );
                sort2(&mut edge_key);

                // handle boundary edge
                match self.boundary_edges.get_mut(&edge_key) {
                    Some(edge) => {
                        // set boundary edge information
                        edge.shared_by_boundary_faces.insert(*face_key);
                    }
                    None => {
                        // edge data
                        let edge_nnode = face_shape.edge_nnode;
                        let mut edge_points: Vec<PointId> = vec![0; edge_nnode];
                        let mut edge_shape = Shape::new(self.space_ndim, 1, edge_nnode)?;

                        // loop over all edge nodes
                        for i in 0..edge_nnode {
                            // configure edge nodes and shape object
                            edge_points[i] = face_points[face_shape.get_edge_node_id(e, i)];
                            for j in 0..self.space_ndim {
                                edge_shape.set_node(i, j, self.points[edge_points[i]].coords[j])?;
                            }

                            // set boundary points
                            self.points[edge_points[i]].shared_by_boundary_edges.insert(edge_key);
                        }

                        // append new boundary edge
                        self.boundary_edges.insert(
                            edge_key,
                            Edge {
                                points: edge_points,
                                shared_by_2d_cells: HashSet::new(),
                                shared_by_boundary_faces: HashSet::from([*face_key]),
                                shape: edge_shape,
                            },
                        );
                    }
                }
            }

            // append new boundary face
            self.boundary_faces.insert(
                *face_key,
                Face {
                    points: face_points,
                    shared_by_cells: HashSet::from([cell_id]),
                    shape: face_shape,
                },
            );
        }
        Ok(())
    }

    /// Computes the range of coordinates
    fn compute_limits(&mut self) -> Result<(), StrError> {
        self.min = vec![f64::MAX; self.space_ndim];
        self.max = vec![f64::MIN; self.space_ndim];
        for point in &self.points {
            for i in 0..self.space_ndim {
                if point.coords[i] < self.min[i] {
                    self.min[i] = point.coords[i];
                }
                if point.coords[i] > self.max[i] {
                    self.max[i] = point.coords[i];
                }
            }
        }
        for i in 0..self.space_ndim {
            if self.min[i] >= self.max[i] {
                return Err("mesh limits are invalid");
            }
        }
        Ok(())
    }

    /// Returns a string with information about all points
    pub(super) fn string_points(&self) -> String {
        let mut buf = String::new();
        for point in &self.points {
            let mut shared_by_boundary_edges: Vec<_> = point.shared_by_boundary_edges.iter().collect();
            let mut shared_by_boundary_faces: Vec<_> = point.shared_by_boundary_faces.iter().collect();
            shared_by_boundary_edges.sort();
            shared_by_boundary_faces.sort();
            write!(
                &mut buf,
                "i:{} x:{:?} e:{:?} f:{:?}\n",
                point.id, point.coords, shared_by_boundary_edges, shared_by_boundary_faces,
            )
            .unwrap();
        }
        buf
    }

    /// Returns a string with information about all cells
    pub(super) fn string_cells(&self) -> String {
        let mut buf = String::new();
        for cell in &self.cells {
            write!(
                &mut buf,
                "i:{} a:{} g:{} p:{:?}\n",
                cell.id, cell.attribute_id, cell.geo_ndim, cell.points,
            )
            .unwrap();
        }
        buf
    }

    /// Returns a string with information about all boundary points
    pub(super) fn string_boundary_points(&self) -> String {
        let mut buf = String::new();
        let mut point_indices: Vec<_> = self.boundary_points.iter().collect();
        point_indices.sort();
        write!(&mut buf, "{:?}\n", point_indices).unwrap();
        buf
    }

    /// Returns a string with information about all boundary edges
    pub(super) fn string_boundary_edges(&self) -> String {
        let mut buf = String::new();
        let mut keys_and_edges: Vec<_> = self.boundary_edges.keys().zip(self.boundary_edges.values()).collect();
        keys_and_edges.sort_by(|left, right| left.0.cmp(&right.0));
        for (key, edge) in keys_and_edges {
            let mut shared_by_2d_cells: Vec<_> = edge.shared_by_2d_cells.iter().collect();
            let mut shared_by_boundary_faces: Vec<_> = edge.shared_by_boundary_faces.iter().collect();
            shared_by_2d_cells.sort();
            shared_by_boundary_faces.sort();
            write!(
                &mut buf,
                "k:({},{}) p:{:?} c:{:?} f:{:?}\n",
                key.0, key.1, edge.points, shared_by_2d_cells, shared_by_boundary_faces,
            )
            .unwrap();
        }
        buf
    }

    /// Returns a string with information about all boundary faces
    pub(super) fn string_boundary_faces(&self) -> String {
        let mut buf = String::new();
        let mut keys_and_faces: Vec<_> = self.boundary_faces.keys().zip(self.boundary_faces.values()).collect();
        keys_and_faces.sort_by(|left, right| left.0.cmp(&right.0));
        for (key, face) in keys_and_faces {
            let mut shared_by_cells: Vec<_> = face.shared_by_cells.iter().collect();
            shared_by_cells.sort();
            write!(
                &mut buf,
                "k:({},{},{},{}) p:{:?} c:{:?}\n",
                key.0, key.1, key.2, key.3, face.points, shared_by_cells
            )
            .unwrap();
        }
        buf
    }
}

impl fmt::Display for Mesh {
    /// Prints mesh data (may be large)
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SUMMARY\n")?;
        write!(f, "=======\n")?;
        write!(f, "space_ndim = {}\n", self.space_ndim)?;
        write!(f, "npoint = {}\n", self.points.len())?;
        write!(f, "ncell = {}\n", self.cells.len())?;
        write!(f, "n_boundary_point = {}\n", self.boundary_points.len())?;
        write!(f, "n_boundary_edge = {}\n", self.boundary_edges.len())?;
        write!(f, "n_boundary_face = {}\n", self.boundary_faces.len())?;

        write!(f, "\nPOINTS\n")?;
        write!(f, "======\n")?;
        write!(f, "{}", self.string_points())?;

        write!(f, "\nCELLS\n")?;
        write!(f, "=====\n")?;
        write!(f, "{}", self.string_cells())?;

        write!(f, "\nBOUNDARY POINTS\n")?;
        write!(f, "===============\n")?;
        write!(f, "{}", self.string_boundary_points())?;

        write!(f, "\nBOUNDARY EDGES\n")?;
        write!(f, "==============\n")?;
        write!(f, "{}", self.string_boundary_edges())?;

        write!(f, "\nBOUNDARY FACES\n")?;
        write!(f, "==============\n")?;
        write!(f, "{}", self.string_boundary_faces())?;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Point;
    use crate::mesh::{At, Mesh};
    use crate::shapes::ShapeState;
    use crate::util::SQRT_2;
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;
    use std::collections::HashSet;

    #[test]
    fn new_fails_on_wrong_input() {
        assert_eq!(Mesh::new(1).err(), Some("space_ndim must be 2 or 3"));
        assert_eq!(Mesh::new(4).err(), Some("space_ndim must be 2 or 3"));
    }

    #[test]
    fn new_works() -> Result<(), StrError> {
        let mesh = Mesh::new(2)?;
        assert_eq!(mesh.space_ndim, 2);
        assert_eq!(mesh.points.len(), 0);
        assert_eq!(mesh.cells.len(), 0);
        assert_eq!(mesh.boundary_points.len(), 0);
        assert_eq!(mesh.boundary_edges.len(), 0);
        assert_eq!(mesh.boundary_faces.len(), 0);
        assert_eq!(mesh.min.len(), 0);
        assert_eq!(mesh.max.len(), 0);
        assert_eq!(
            format!("{}", mesh.grid_boundary_points),
            "ids = []\n\
             nitem = 0\n\
             ncontainer = 0\n"
        );
        assert_eq!(mesh.derived_props_computed, false);
        Ok(())
    }

    #[test]
    fn from_text_fails_on_invalid_data() -> Result<(), StrError> {
        assert_eq!(
            Mesh::from_text("").err(),
            Some("text string is empty or header is missing")
        );
        Ok(())
    }

    #[test]
    fn from_text_file_fails_on_invalid_data() -> Result<(), StrError> {
        assert_eq!(Mesh::from_text_file("").err(), Some("cannot open file"));
        Ok(())
    }

    #[test]
    fn compute_limits_fails_on_wrong_data() -> Result<(), StrError> {
        let mut mesh = Mesh::new(2)?;
        mesh.points.push(Point {
            id: 0,
            coords: vec![0.0, 0.0],
            shared_by_boundary_edges: HashSet::new(),
            shared_by_boundary_faces: HashSet::new(),
        });
        assert_eq!(mesh.compute_limits().err(), Some("mesh limits are invalid"));
        Ok(())
    }

    #[test]
    fn from_text_fails_on_wrong_jacobian_2d() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let res = Mesh::from_text(
            r"# header
            # space_ndim npoint ncell
                       2      6     2
            
            # points
            # id   x   y
               0 0.0 0.0
               1 1.0 0.0
               2 1.0 1.0
               3 0.0 1.0
               4 2.0 0.0
               5 2.0 1.0
            
            # cells
            # id att geo_ndim nnode  (wrong) point_ids...
               0   1        2     4  0 3 2 1
               1   0        2     4  1 2 5 4",
        );
        assert_eq!(res.err(), Some("a cell has incorrect ordering of nodes"));
        Ok(())
    }

    #[test]
    fn from_text_and_strings_work() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let mesh = Mesh::from_text(
            r"# header
            # space_ndim npoint ncell
                       2      6     2
            
            # points
            # id   x   y
               0 0.0 0.0
               1 1.0 0.0
               2 1.0 1.0
               3 0.0 1.0
               4 2.0 0.0
               5 2.0 1.0
            
            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        2     4  0 1 2 3
               1   0        2     4  1 4 5 2",
        )?;
        assert_eq!(mesh.space_ndim, 2);
        assert_eq!(mesh.points.len(), 6);
        assert_eq!(mesh.cells.len(), 2);
        assert_eq!(mesh.boundary_points.len(), 6);
        assert_eq!(mesh.boundary_edges.len(), 6);
        assert_eq!(mesh.boundary_faces.len(), 0);
        assert_eq!(mesh.min, &[0.0, 0.0]);
        assert_eq!(mesh.max, &[2.0, 1.0]);
        assert_eq!(
            format!("{}", mesh.grid_boundary_points),
            "0: [0]\n\
             4: [1]\n\
             5: [1]\n\
             9: [4]\n\
             90: [3]\n\
             94: [2]\n\
             95: [2]\n\
             99: [5]\n\
             ids = [0, 1, 2, 3, 4, 5]\n\
             nitem = 6\n\
             ncontainer = 8\n"
        );
        assert_eq!(mesh.derived_props_computed, true);
        assert_eq!(
            format!("{}", mesh.string_points()),
            "i:0 x:[0.0, 0.0] e:[(0, 1), (0, 3)] f:[]\n\
             i:1 x:[1.0, 0.0] e:[(0, 1), (1, 4)] f:[]\n\
             i:2 x:[1.0, 1.0] e:[(2, 3), (2, 5)] f:[]\n\
             i:3 x:[0.0, 1.0] e:[(0, 3), (2, 3)] f:[]\n\
             i:4 x:[2.0, 0.0] e:[(1, 4), (4, 5)] f:[]\n\
             i:5 x:[2.0, 1.0] e:[(2, 5), (4, 5)] f:[]\n"
        );
        assert_eq!(
            format!("{}", mesh.string_cells()),
            "i:0 a:1 g:2 p:[0, 1, 2, 3]\n\
             i:1 a:0 g:2 p:[1, 4, 5, 2]\n"
        );
        assert_eq!(format!("{}", mesh.string_boundary_points()), "[0, 1, 2, 3, 4, 5]\n");
        assert_eq!(
            format!("{}", mesh.string_boundary_edges()),
            "k:(0,1) p:[1, 0] c:[0] f:[]\n\
             k:(0,3) p:[0, 3] c:[0] f:[]\n\
             k:(1,4) p:[4, 1] c:[1] f:[]\n\
             k:(2,3) p:[3, 2] c:[0] f:[]\n\
             k:(2,5) p:[2, 5] c:[1] f:[]\n\
             k:(4,5) p:[5, 4] c:[1] f:[]\n"
        );
        assert_eq!(format!("{}", mesh.string_boundary_faces()), "");
        Ok(())
    }

    #[test]
    fn from_text_file_fails_on_wrong_jacobian_3d() -> Result<(), StrError> {
        let res = Mesh::from_text_file("./data/meshes/bad_wrong_jacobian.msh");
        assert_eq!(res.err(), Some("a cell has incorrect ordering of nodes"));
        Ok(())
    }

    #[test]
    fn from_text_file_fails_on_wrong_nodes_3d() -> Result<(), StrError> {
        let res = Mesh::from_text_file("./data/meshes/bad_wrong_nodes.msh");
        assert_eq!(res.err(), Some("cannot compute inverse due to zero determinant"));
        Ok(())
    }

    #[test]
    fn from_text_file_and_display_work() -> Result<(), StrError> {
        //
        //       8-------------11
        //      /.             /|
        //     / .            / |
        //    /  .           /  |
        //   /   .          /   |
        //  9-------------10    |
        //  |    .         |    |
        //  |    4---------|----7
        //  |   /.         |   /|
        //  |  / .         |  / |
        //  | /  .         | /  |
        //  |/   .         |/   |
        //  5--------------6    |
        //  |    .         |    |
        //  |    0---------|----3
        //  |   /          |   /
        //  |  /           |  /
        //  | /            | /
        //  |/             |/
        //  1--------------2
        //
        let mesh = Mesh::from_text_file("./data/meshes/ok2.msh")?;
        assert_eq!(mesh.space_ndim, 3);
        assert_eq!(mesh.points.len(), 12);
        assert_eq!(mesh.cells.len(), 2);
        assert_eq!(mesh.boundary_points.len(), 12);
        assert_eq!(mesh.boundary_edges.len(), 20);
        assert_eq!(mesh.boundary_faces.len(), 10);
        assert_eq!(mesh.min, &[0.0, 0.0, 0.0]);
        assert_eq!(mesh.max, &[1.0, 1.0, 2.0]);
        assert_eq!(
            format!("{}", mesh.grid_boundary_points),
            "0: [0]\n\
             9: [1]\n\
             90: [3]\n\
             99: [2]\n\
             400: [4]\n\
             409: [5]\n\
             490: [7]\n\
             499: [6]\n\
             500: [4]\n\
             509: [5]\n\
             590: [7]\n\
             599: [6]\n\
             900: [8]\n\
             909: [9]\n\
             990: [11]\n\
             999: [10]\n\
             ids = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]\n\
             nitem = 12\n\
             ncontainer = 16\n"
        );
        assert_eq!(mesh.derived_props_computed, true);
        assert_eq!(
            format!("{}", mesh.string_points()),
            "i:0 x:[0.0, 0.0, 0.0] e:[(0, 1), (0, 3), (0, 4)] f:[(0, 1, 2, 3), (0, 1, 4, 5), (0, 3, 4, 7)]\n\
             i:1 x:[1.0, 0.0, 0.0] e:[(0, 1), (1, 2), (1, 5)] f:[(0, 1, 2, 3), (0, 1, 4, 5), (1, 2, 5, 6)]\n\
             i:2 x:[1.0, 1.0, 0.0] e:[(1, 2), (2, 3), (2, 6)] f:[(0, 1, 2, 3), (1, 2, 5, 6), (2, 3, 6, 7)]\n\
             i:3 x:[0.0, 1.0, 0.0] e:[(0, 3), (2, 3), (3, 7)] f:[(0, 1, 2, 3), (0, 3, 4, 7), (2, 3, 6, 7)]\n\
             i:4 x:[0.0, 0.0, 1.0] e:[(0, 4), (4, 5), (4, 7), (4, 8)] f:[(0, 1, 4, 5), (0, 3, 4, 7), (4, 5, 8, 9), (4, 7, 8, 11)]\n\
             i:5 x:[1.0, 0.0, 1.0] e:[(1, 5), (4, 5), (5, 6), (5, 9)] f:[(0, 1, 4, 5), (1, 2, 5, 6), (4, 5, 8, 9), (5, 6, 9, 10)]\n\
             i:6 x:[1.0, 1.0, 1.0] e:[(2, 6), (5, 6), (6, 7), (6, 10)] f:[(1, 2, 5, 6), (2, 3, 6, 7), (5, 6, 9, 10), (6, 7, 10, 11)]\n\
             i:7 x:[0.0, 1.0, 1.0] e:[(3, 7), (4, 7), (6, 7), (7, 11)] f:[(0, 3, 4, 7), (2, 3, 6, 7), (4, 7, 8, 11), (6, 7, 10, 11)]\n\
             i:8 x:[0.0, 0.0, 2.0] e:[(4, 8), (8, 9), (8, 11)] f:[(4, 5, 8, 9), (4, 7, 8, 11), (8, 9, 10, 11)]\n\
             i:9 x:[1.0, 0.0, 2.0] e:[(5, 9), (8, 9), (9, 10)] f:[(4, 5, 8, 9), (5, 6, 9, 10), (8, 9, 10, 11)]\n\
             i:10 x:[1.0, 1.0, 2.0] e:[(6, 10), (9, 10), (10, 11)] f:[(5, 6, 9, 10), (6, 7, 10, 11), (8, 9, 10, 11)]\n\
             i:11 x:[0.0, 1.0, 2.0] e:[(7, 11), (8, 11), (10, 11)] f:[(4, 7, 8, 11), (6, 7, 10, 11), (8, 9, 10, 11)]\n"
        );
        assert_eq!(
            format!("{}", mesh.string_cells()),
            "i:0 a:1 g:3 p:[0, 1, 2, 3, 4, 5, 6, 7]\n\
             i:1 a:0 g:3 p:[4, 5, 6, 7, 8, 9, 10, 11]\n"
        );
        assert_eq!(
            format!("{}", mesh.string_boundary_points()),
            "[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]\n"
        );
        assert_eq!(
            format!("{}", mesh.string_boundary_edges()),
            "k:(0,1) p:[0, 1] c:[] f:[(0, 1, 2, 3), (0, 1, 4, 5)]\n\
             k:(0,3) p:[3, 0] c:[] f:[(0, 1, 2, 3), (0, 3, 4, 7)]\n\
             k:(0,4) p:[0, 4] c:[] f:[(0, 1, 4, 5), (0, 3, 4, 7)]\n\
             k:(1,2) p:[1, 2] c:[] f:[(0, 1, 2, 3), (1, 2, 5, 6)]\n\
             k:(1,5) p:[5, 1] c:[] f:[(0, 1, 4, 5), (1, 2, 5, 6)]\n\
             k:(2,3) p:[2, 3] c:[] f:[(0, 1, 2, 3), (2, 3, 6, 7)]\n\
             k:(2,6) p:[6, 2] c:[] f:[(1, 2, 5, 6), (2, 3, 6, 7)]\n\
             k:(3,7) p:[3, 7] c:[] f:[(0, 3, 4, 7), (2, 3, 6, 7)]\n\
             k:(4,5) p:[4, 5] c:[] f:[(0, 1, 4, 5), (4, 5, 8, 9)]\n\
             k:(4,7) p:[7, 4] c:[] f:[(0, 3, 4, 7), (4, 7, 8, 11)]\n\
             k:(4,8) p:[4, 8] c:[] f:[(4, 5, 8, 9), (4, 7, 8, 11)]\n\
             k:(5,6) p:[5, 6] c:[] f:[(1, 2, 5, 6), (5, 6, 9, 10)]\n\
             k:(5,9) p:[9, 5] c:[] f:[(4, 5, 8, 9), (5, 6, 9, 10)]\n\
             k:(6,7) p:[6, 7] c:[] f:[(2, 3, 6, 7), (6, 7, 10, 11)]\n\
             k:(6,10) p:[10, 6] c:[] f:[(5, 6, 9, 10), (6, 7, 10, 11)]\n\
             k:(7,11) p:[7, 11] c:[] f:[(4, 7, 8, 11), (6, 7, 10, 11)]\n\
             k:(8,9) p:[8, 9] c:[] f:[(4, 5, 8, 9), (8, 9, 10, 11)]\n\
             k:(8,11) p:[11, 8] c:[] f:[(4, 7, 8, 11), (8, 9, 10, 11)]\n\
             k:(9,10) p:[9, 10] c:[] f:[(5, 6, 9, 10), (8, 9, 10, 11)]\n\
             k:(10,11) p:[10, 11] c:[] f:[(6, 7, 10, 11), (8, 9, 10, 11)]\n"
        );
        assert_eq!(
            format!("{}", mesh.string_boundary_faces()),
            "k:(0,1,2,3) p:[0, 3, 2, 1] c:[0]\n\
             k:(0,1,4,5) p:[0, 1, 5, 4] c:[0]\n\
             k:(0,3,4,7) p:[0, 4, 7, 3] c:[0]\n\
             k:(1,2,5,6) p:[1, 2, 6, 5] c:[0]\n\
             k:(2,3,6,7) p:[2, 3, 7, 6] c:[0]\n\
             k:(4,5,8,9) p:[4, 5, 9, 8] c:[1]\n\
             k:(4,7,8,11) p:[4, 8, 11, 7] c:[1]\n\
             k:(5,6,9,10) p:[5, 6, 10, 9] c:[1]\n\
             k:(6,7,10,11) p:[6, 7, 11, 10] c:[1]\n\
             k:(8,9,10,11) p:[8, 9, 10, 11] c:[1]\n"
        );
        Ok(())
    }

    #[test]
    fn normals_are_correct_2d() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let mut mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;

        // the norm of the normal vector should be equal to 0.5 = edge_length / 2.0
        // where 2.0 corresponds to the edge_length in the reference system
        let l = 0.5; // norm of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 4)], [0.0, -l]),
            // right
            (vec![(4, 5)], [l, 0.0]),
            // top
            (vec![(2, 3), (2, 5)], [0.0, l]),
            // left
            (vec![(0, 3)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let mut normal = Vector::new(mesh.space_ndim);
        let ksi = &[0.0, 0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let edge = mesh.boundary_edges.get_mut(edge_key).unwrap();
                let mut state = ShapeState::new(edge.shape.space_ndim, edge.shape.geo_ndim, edge.shape.nnode)?;
                edge.shape.calc_boundary_normal(&mut state, &mut normal, ksi)?;
                assert_vec_approx_eq!(normal.as_data(), solution, 1e-15);
            }
        }
        Ok(())
    }

    #[test]
    fn display_works() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;
        assert_eq!(
            format!("{}", mesh),
            "SUMMARY\n\
             =======\n\
             space_ndim = 2\n\
             npoint = 6\n\
             ncell = 2\n\
             n_boundary_point = 6\n\
             n_boundary_edge = 6\n\
             n_boundary_face = 0\n\
             \n\
             POINTS\n\
             ======\n\
             i:0 x:[0.0, 0.0] e:[(0, 1), (0, 3)] f:[]\n\
             i:1 x:[1.0, 0.0] e:[(0, 1), (1, 4)] f:[]\n\
             i:2 x:[1.0, 1.0] e:[(2, 3), (2, 5)] f:[]\n\
             i:3 x:[0.0, 1.0] e:[(0, 3), (2, 3)] f:[]\n\
             i:4 x:[2.0, 0.0] e:[(1, 4), (4, 5)] f:[]\n\
             i:5 x:[2.0, 1.0] e:[(2, 5), (4, 5)] f:[]\n\
             \n\
             CELLS\n\
             =====\n\
             i:0 a:1 g:2 p:[0, 1, 2, 3]\n\
             i:1 a:0 g:2 p:[1, 4, 5, 2]\n\
             \n\
             BOUNDARY POINTS\n\
             ===============\n\
             [0, 1, 2, 3, 4, 5]\n\
             \n\
             BOUNDARY EDGES\n\
             ==============\n\
             k:(0,1) p:[1, 0] c:[0] f:[]\n\
             k:(0,3) p:[0, 3] c:[0] f:[]\n\
             k:(1,4) p:[4, 1] c:[1] f:[]\n\
             k:(2,3) p:[3, 2] c:[0] f:[]\n\
             k:(2,5) p:[2, 5] c:[1] f:[]\n\
             k:(4,5) p:[5, 4] c:[1] f:[]\n\
             \n\
             BOUNDARY FACES\n\
             ==============\n"
        );
        Ok(())
    }

    #[test]
    fn find_boundary_fails_on_wrong_input() -> Result<(), StrError> {
        let mut mesh = Mesh::new(2)?;
        assert_eq!(
            mesh.find_boundary_points(At::XY(0.0, 0.0)).err(),
            Some("compute_derived_props must be called first")
        );
        assert_eq!(
            mesh.find_boundary_edges(At::XY(0.0, 0.0)).err(),
            Some("compute_derived_props must be called first")
        );
        let mut mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;
        assert_eq!(
            mesh.find_boundary_points(At::Z(0.0)).err(),
            Some("At::Z works in 3D only")
        );
        assert_eq!(
            mesh.find_boundary_points(At::YZ(0.0, 0.0)).err(),
            Some("At::YZ works in 3D only")
        );
        assert_eq!(
            mesh.find_boundary_points(At::XZ(0.0, 0.0)).err(),
            Some("At::XZ works in 3D only")
        );
        assert_eq!(
            mesh.find_boundary_points(At::XYZ(0.0, 0.0, 0.0)).err(),
            Some("At::XYZ works in 3D only")
        );
        assert_eq!(
            mesh.find_boundary_points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
                .err(),
            Some("At::Cylinder works in 3D only")
        );
        let mut mesh = Mesh::from_text_file("./data/meshes/ok2.msh")?;
        assert_eq!(
            mesh.find_boundary_points(At::Circle(0.0, 0.0, 0.0)).err(),
            Some("At::Circle works in 2D only")
        );
        Ok(())
    }

    #[test]
    fn find_boundary_works_2d() -> Result<(), StrError> {
        // `.       `.
        //   3--------2--------5
        //   | `.     | `.     |
        //   |   `~.  |        |
        //   |      `.|        |
        //   0--------1--------4
        //
        let mut mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;
        assert_eq!(mesh.find_boundary_points(At::XY(0.0, 0.0))?, &[0]);
        assert_eq!(mesh.find_boundary_points(At::XY(2.0, 1.0))?, &[5]);
        assert_eq!(
            mesh.find_boundary_points(At::XY(10.0, 0.0)).err(),
            Some("point is outside the grid")
        );
        assert_eq!(mesh.find_boundary_edges(At::Y(0.0))?, &[(0, 1), (1, 4)]);
        assert_eq!(mesh.find_boundary_edges(At::X(2.0))?, &[(4, 5)]);
        assert_eq!(mesh.find_boundary_edges(At::Y(1.0))?, &[(2, 3), (2, 5)]);
        assert_eq!(mesh.find_boundary_edges(At::X(0.0))?, &[(0, 3)]);
        assert_eq!(mesh.find_boundary_edges(At::X(10.0))?, &[]);
        assert_eq!(mesh.find_boundary_points(At::Circle(0.0, 0.0, 1.0))?, &[1, 3]);
        assert_eq!(mesh.find_boundary_points(At::Circle(0.0, 0.0, SQRT_2))?, &[2]);
        assert_eq!(mesh.find_boundary_points(At::Circle(0.0, 0.0, 10.0))?, &[]);
        Ok(())
    }

    #[test]
    fn find_boundary_works_3d() -> Result<(), StrError> {
        //
        //       8-------------11
        //      /.             /|
        //     / .            / |
        //    /  .           /  |
        //   /   .          /   |
        //  9-------------10    |
        //  |    .         |    |
        //  |    4---------|----7
        //  |   /.         |   /|
        //  |  / .         |  / |
        //  | /  .         | /  |
        //  |/   .         |/   |
        //  5--------------6    |
        //  |    .         |    |
        //  |    0---------|----3
        //  |   /          |   /
        //  |  /           |  /
        //  | /            | /
        //  |/             |/
        //  1--------------2
        //
        let mut mesh = Mesh::from_text_file("./data/meshes/ok2.msh")?;
        assert_eq!(mesh.find_boundary_points(At::X(0.0))?, &[0, 3, 4, 7, 8, 11]);
        assert_eq!(mesh.find_boundary_points(At::X(1.0))?, &[1, 2, 5, 6, 9, 10]);
        assert_eq!(mesh.find_boundary_points(At::X(10.0))?, &[]);
        assert_eq!(mesh.find_boundary_points(At::Y(0.0))?, &[0, 1, 4, 5, 8, 9]);
        assert_eq!(mesh.find_boundary_points(At::Y(1.0))?, &[2, 3, 6, 7, 10, 11]);
        assert_eq!(mesh.find_boundary_points(At::Y(10.0))?, &[]);
        assert_eq!(mesh.find_boundary_points(At::Z(0.0))?, &[0, 1, 2, 3]);
        assert_eq!(mesh.find_boundary_points(At::Z(1.0))?, &[4, 5, 6, 7]);
        assert_eq!(mesh.find_boundary_points(At::Z(2.0))?, &[8, 9, 10, 11]);
        assert_eq!(mesh.find_boundary_points(At::Z(10.0))?, &[]);
        assert_eq!(mesh.find_boundary_points(At::XY(0.0, 0.0))?, &[0, 4, 8]);
        assert_eq!(mesh.find_boundary_points(At::XY(1.0, 1.0))?, &[2, 6, 10]);
        assert_eq!(mesh.find_boundary_points(At::XY(10.0, 10.0))?, &[]);
        assert_eq!(mesh.find_boundary_points(At::YZ(0.0, 0.0))?, &[0, 1]);
        assert_eq!(mesh.find_boundary_points(At::YZ(1.0, 1.0))?, &[6, 7]);
        assert_eq!(mesh.find_boundary_points(At::XZ(0.0, 0.0))?, &[0, 3]);
        assert_eq!(mesh.find_boundary_points(At::XZ(1.0, 0.0))?, &[1, 2]);
        assert_eq!(mesh.find_boundary_points(At::XZ(1.0, 2.0))?, &[9, 10]);
        assert_eq!(mesh.find_boundary_points(At::XYZ(0.0, 0.0, 0.0))?, &[0]);
        assert_eq!(mesh.find_boundary_points(At::XYZ(1.0, 1.0, 2.0))?, &[10]);
        assert_eq!(
            mesh.find_boundary_points(At::XYZ(10.0, 0.0, 0.0)).err(),
            Some("point is outside the grid")
        );
        assert_eq!(
            mesh.find_boundary_points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0))?,
            &[1, 3, 5, 7, 9, 11]
        );
        assert_eq!(
            mesh.find_boundary_points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, SQRT_2))?,
            &[2, 6, 10]
        );
        assert_eq!(
            mesh.find_boundary_points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 10.0))?,
            &[]
        );
        Ok(())
    }

    #[test]
    fn from_text_file_works_with_mixed() -> Result<(), StrError> {
        //
        //          4--------3
        //          |        |
        //          |        |
        //          |        |
        //  0-------1--------2
        //
        let mesh = Mesh::from_text_file("./data/meshes/ok_mixed_shapes2.msh")?;
        assert_eq!(
            format!("{}", mesh),
            "SUMMARY\n\
             =======\n\
             space_ndim = 2\n\
             npoint = 5\n\
             ncell = 2\n\
             n_boundary_point = 5\n\
             n_boundary_edge = 4\n\
             n_boundary_face = 0\n\
             \n\
             POINTS\n\
             ======\n\
             i:0 x:[0.0, 0.0] e:[] f:[]\n\
             i:1 x:[1.0, 0.0] e:[(1, 2), (1, 4)] f:[]\n\
             i:2 x:[2.0, 0.0] e:[(1, 2), (2, 3)] f:[]\n\
             i:3 x:[2.0, 1.0] e:[(2, 3), (3, 4)] f:[]\n\
             i:4 x:[1.0, 1.0] e:[(1, 4), (3, 4)] f:[]\n\
             \n\
             CELLS\n\
             =====\n\
             i:0 a:1 g:1 p:[0, 1]\n\
             i:1 a:2 g:2 p:[1, 2, 3, 4]\n\
             \n\
             BOUNDARY POINTS\n\
             ===============\n\
             [0, 1, 2, 3, 4]\n\
             \n\
             BOUNDARY EDGES\n\
             ==============\n\
             k:(1,2) p:[2, 1] c:[1] f:[]\n\
             k:(1,4) p:[1, 4] c:[1] f:[]\n\
             k:(2,3) p:[3, 2] c:[1] f:[]\n\
             k:(3,4) p:[4, 3] c:[1] f:[]\n\
             \n\
             BOUNDARY FACES\n\
             ==============\n"
        );

        //
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
        //
        let mesh = Mesh::from_text_file("./data/meshes/ok_mixed_shapes.msh")?;
        assert_eq!(
            format!("{}", mesh),
            "SUMMARY\n\
             =======\n\
             space_ndim = 3\n\
             npoint = 13\n\
             ncell = 5\n\
             n_boundary_point = 13\n\
             n_boundary_edge = 16\n\
             n_boundary_face = 10\n\
             \n\
             POINTS\n\
             ======\n\
             i:0 x:[0.0, 0.0, 0.0] e:[(0, 1), (0, 3), (0, 4)] f:[(0, 1, 2, 3), (0, 1, 4, 5), (0, 3, 4, 7)]\n\
             i:1 x:[1.0, 0.0, 0.0] e:[(0, 1), (1, 2), (1, 5)] f:[(0, 1, 2, 3), (0, 1, 4, 5), (1, 2, 5, 6)]\n\
             i:2 x:[1.0, 1.0, 0.0] e:[(1, 2), (2, 3), (2, 6), (2, 8)] f:[(0, 1, 2, 3), (1, 2, 5, 6), (2, 3, 6, 7), (2, 3, 6, 13), (2, 3, 8, 13), (2, 6, 8, 13)]\n\
             i:3 x:[0.0, 1.0, 0.0] e:[(0, 3), (2, 3), (3, 6), (3, 7), (3, 8)] f:[(0, 1, 2, 3), (0, 3, 4, 7), (2, 3, 6, 7), (2, 3, 6, 13), (2, 3, 8, 13), (3, 6, 8, 13)]\n\
             i:4 x:[0.0, 0.0, 1.0] e:[(0, 4), (4, 5), (4, 7)] f:[(0, 1, 4, 5), (0, 3, 4, 7), (4, 5, 6, 7)]\n\
             i:5 x:[1.0, 0.0, 1.0] e:[(1, 5), (4, 5), (5, 6)] f:[(0, 1, 4, 5), (1, 2, 5, 6), (4, 5, 6, 7)]\n\
             i:6 x:[1.0, 1.0, 1.0] e:[(2, 6), (3, 6), (5, 6), (6, 7), (6, 8)] f:[(1, 2, 5, 6), (2, 3, 6, 7), (2, 3, 6, 13), (2, 6, 8, 13), (3, 6, 8, 13), (4, 5, 6, 7)]\n\
             i:7 x:[0.0, 1.0, 1.0] e:[(3, 7), (4, 7), (6, 7)] f:[(0, 3, 4, 7), (2, 3, 6, 7), (4, 5, 6, 7)]\n\
             i:8 x:[1.0, 2.0, 0.0] e:[(2, 8), (3, 8), (6, 8)] f:[(2, 3, 8, 13), (2, 6, 8, 13), (3, 6, 8, 13)]\n\
             i:9 x:[0.0, 2.0, 0.0] e:[] f:[]\n\
             i:10 x:[0.0, 2.0, 1.0] e:[] f:[]\n\
             i:11 x:[1.0, -0.5, 0.0] e:[] f:[]\n\
             i:12 x:[1.0, -1.0, 0.0] e:[] f:[]\n\
             \n\
             CELLS\n\
             =====\n\
             i:0 a:1 g:3 p:[0, 1, 2, 3, 4, 5, 6, 7]\n\
             i:1 a:2 g:3 p:[2, 8, 3, 6]\n\
             i:2 a:2 g:2 p:[3, 9, 10, 7]\n\
             i:3 a:2 g:2 p:[8, 9, 3]\n\
             i:4 a:3 g:1 p:[1, 11, 12]\n\
             \n\
             BOUNDARY POINTS\n\
             ===============\n\
             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]\n\
             \n\
             BOUNDARY EDGES\n\
             ==============\n\
             k:(0,1) p:[0, 1] c:[] f:[(0, 1, 2, 3), (0, 1, 4, 5)]\n\
             k:(0,3) p:[3, 0] c:[] f:[(0, 1, 2, 3), (0, 3, 4, 7)]\n\
             k:(0,4) p:[0, 4] c:[] f:[(0, 1, 4, 5), (0, 3, 4, 7)]\n\
             k:(1,2) p:[1, 2] c:[] f:[(0, 1, 2, 3), (1, 2, 5, 6)]\n\
             k:(1,5) p:[5, 1] c:[] f:[(0, 1, 4, 5), (1, 2, 5, 6)]\n\
             k:(2,3) p:[2, 3] c:[] f:[(0, 1, 2, 3), (2, 3, 6, 7), (2, 3, 6, 13), (2, 3, 8, 13)]\n\
             k:(2,6) p:[6, 2] c:[] f:[(1, 2, 5, 6), (2, 3, 6, 7), (2, 3, 6, 13), (2, 6, 8, 13)]\n\
             k:(2,8) p:[2, 8] c:[] f:[(2, 3, 8, 13), (2, 6, 8, 13)]\n\
             k:(3,6) p:[3, 6] c:[] f:[(2, 3, 6, 13), (3, 6, 8, 13)]\n\
             k:(3,7) p:[3, 7] c:[] f:[(0, 3, 4, 7), (2, 3, 6, 7)]\n\
             k:(3,8) p:[8, 3] c:[] f:[(2, 3, 8, 13), (3, 6, 8, 13)]\n\
             k:(4,5) p:[4, 5] c:[] f:[(0, 1, 4, 5), (4, 5, 6, 7)]\n\
             k:(4,7) p:[7, 4] c:[] f:[(0, 3, 4, 7), (4, 5, 6, 7)]\n\
             k:(5,6) p:[5, 6] c:[] f:[(1, 2, 5, 6), (4, 5, 6, 7)]\n\
             k:(6,7) p:[6, 7] c:[] f:[(2, 3, 6, 7), (4, 5, 6, 7)]\n\
             k:(6,8) p:[6, 8] c:[] f:[(2, 6, 8, 13), (3, 6, 8, 13)]\n\
             \n\
             BOUNDARY FACES\n\
             ==============\n\
             k:(0,1,2,3) p:[0, 3, 2, 1] c:[0]\n\
             k:(0,1,4,5) p:[0, 1, 5, 4] c:[0]\n\
             k:(0,3,4,7) p:[0, 4, 7, 3] c:[0]\n\
             k:(1,2,5,6) p:[1, 2, 6, 5] c:[0]\n\
             k:(2,3,6,7) p:[2, 3, 7, 6] c:[0]\n\
             k:(2,3,6,13) p:[2, 6, 3] c:[1]\n\
             k:(2,3,8,13) p:[2, 3, 8] c:[1]\n\
             k:(2,6,8,13) p:[2, 8, 6] c:[1]\n\
             k:(3,6,8,13) p:[8, 3, 6] c:[1]\n\
             k:(4,5,6,7) p:[4, 5, 6, 7] c:[0]\n"
         );
        Ok(())
    }
}
