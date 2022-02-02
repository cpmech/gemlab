use super::read_mesh::{parse_mesh, read_mesh};
use super::At;
use crate::shapes::Shape;
use crate::util::GridSearch;
use crate::StrError;
use russell_lab::{sort2, sort4};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::fmt::{self, Write as FmtWrite};
use std::fs;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

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
#[derive(Clone, Debug, Deserialize, Serialize)]
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
#[derive(Clone, Debug, Deserialize, Serialize)]
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
    /// **raw data**
    pub points: Vec<PointId>,
}

/// Holds edge data (derived data structure)
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Edge {
    /// List of points defining this edge; in the right order (unsorted)
    pub points: Vec<PointId>,

    /// Set of 2D cells sharing this edge (to find the boundary)
    ///
    /// **2D mesh only**
    pub shared_by_2d_cells: HashSet<CellId>,

    /// Set of boundary faces sharing this edge
    pub shared_by_boundary_faces: HashSet<FaceKey>,
}

/// Holds face data (derived data structure)
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Face {
    /// List of points defining this face; in the right order (unsorted)
    pub points: Vec<PointId>,

    /// Set of cells sharing this face
    pub shared_by_cells: HashSet<CellId>,
}

/// Holds mesh data
#[derive(Deserialize, Serialize)]
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

    /// Returns a new mesh with pre-allocated (empty) Point and Cell vectors
    pub(super) fn new_sized(space_ndim: usize, npoint: usize, ncell: usize) -> Result<Self, StrError> {
        if space_ndim < 2 || space_ndim > 3 {
            return Err("space_ndim must be 2 or 3");
        }
        if npoint < 2 {
            return Err("npoint must be greater than or equal to 2");
        }
        if ncell < 1 {
            return Err("ncell must be greater than or equal to 1");
        }
        let zero_point = Point {
            id: 0,
            coords: vec![0.0; space_ndim],
            shared_by_boundary_edges: HashSet::new(),
            shared_by_boundary_faces: HashSet::new(),
        };
        let zero_cell = Cell {
            id: 0,
            attribute_id: 0,
            geo_ndim: 0,
            points: Vec::new(),
        };
        Ok(Mesh {
            space_ndim,
            points: vec![zero_point; npoint],
            cells: vec![zero_cell; ncell],
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
    pub fn from_text(raw_mesh_data: &str) -> Result<Self, StrError> {
        let mut mesh = parse_mesh(raw_mesh_data)?;
        mesh.compute_derived_props()?;
        Ok(mesh)
    }

    /// Reads raw mesh data from text file and computes derived properties
    pub fn from_text_file<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let mut mesh = read_mesh(full_path)?;
        mesh.compute_derived_props()?;
        Ok(mesh)
    }

    /// Reads a binary file containing the mesh and computed properties
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn read<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let mut file = File::open(&path).map_err(|_| "no file found")?;
        let metadata = fs::metadata(&path).map_err(|_| "unable to read metadata")?;
        let mut serialized = vec![0; metadata.len() as usize];
        file.read(&mut serialized).expect("buffer overflow");
        Mesh::decode(&serialized)
    }

    /// Writes a binary file with the mesh and computed properties
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn save<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }
        let serialized = self.encode()?;
        let mut file = File::create(&path).map_err(|_| "cannot create file")?;
        file.write_all(&serialized).map_err(|_| "cannot write file")?;
        Ok(())
    }

    /// Computes derived properties such as boundaries and limits
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
                    for id in self.grid_boundary_points.find_on_line(&[x, 0.0], &[x, 1.0]).unwrap() {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid_boundary_points.find_on_plane_yz(x).unwrap() {
                        point_ids.insert(id);
                    }
                }
            }
            At::Y(y) => {
                if self.space_ndim == 2 {
                    for id in self.grid_boundary_points.find_on_line(&[0.0, y], &[1.0, y]).unwrap() {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid_boundary_points.find_on_plane_xz(y).unwrap() {
                        point_ids.insert(id);
                    }
                }
            }
            At::Z(z) => {
                if self.space_ndim == 2 {
                    return Err("At::Z works in 3D only");
                } else {
                    for id in self.grid_boundary_points.find_on_plane_xy(z).unwrap() {
                        point_ids.insert(id);
                    }
                }
            }
            At::XY(x, y) => {
                if self.space_ndim == 2 {
                    if let Some(id) = self.grid_boundary_points.find(&[x, y]).unwrap() {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self
                        .grid_boundary_points
                        .find_on_line(&[x, y, 0.0], &[x, y, 1.0])
                        .unwrap()
                    {
                        point_ids.insert(id);
                    }
                }
            }
            At::YZ(y, z) => {
                if self.space_ndim == 2 {
                    return Err("At::YZ works in 3D only");
                } else {
                    for id in self
                        .grid_boundary_points
                        .find_on_line(&[0.0, y, z], &[1.0, y, z])
                        .unwrap()
                    {
                        point_ids.insert(id);
                    }
                }
            }
            At::XZ(x, z) => {
                if self.space_ndim == 2 {
                    return Err("At::XZ works in 3D only");
                } else {
                    for id in self
                        .grid_boundary_points
                        .find_on_line(&[x, 0.0, z], &[x, 1.0, z])
                        .unwrap()
                    {
                        point_ids.insert(id);
                    }
                }
            }
            At::XYZ(x, y, z) => {
                if self.space_ndim == 2 {
                    return Err("At::XYZ works in 3D only");
                } else {
                    if let Some(id) = self.grid_boundary_points.find(&[x, y, z]).unwrap() {
                        point_ids.insert(id);
                    }
                }
            }
            At::Circle(x, y, r) => {
                if self.space_ndim == 2 {
                    for id in self.grid_boundary_points.find_on_circle(&[x, y], r).unwrap() {
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
                        .find_on_cylinder(&[ax, ay, az], &[bx, by, bz], r)
                        .unwrap()
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
                        // // update groups map
                        // if self.groups_boundary_edges.contains_key(group) {
                        //     let keys = self.groups_boundary_edges.get_mut(group).unwrap();
                        //     keys.insert(key.clone());
                        // // new groups map
                        // } else {
                        //     let mut keys = HashSet::new();
                        //     keys.insert(key.clone());
                        //     self.groups_boundary_edges.insert(group.to_string(), keys);
                        // }
                    }
                }
            }
        }
        let mut keys: Vec<_> = edge_keys.into_iter().collect();
        keys.sort();
        Ok(keys)
    }

    /// TODO
    pub fn extract_coords(&self, shape: &mut Shape, cell_id: usize) -> Result<(), StrError> {
        let npoint = self.cells[cell_id].points.len();
        for m in 0..npoint {
            let point_id = self.cells[cell_id].points[m];
            for i in 0..self.space_ndim {
                shape.set_node(m, i, self.points[point_id].coords[i])?;
            }
        }
        Ok(())
    }

    // ======================================================================================================
    // --- private ------------------------------------------------------------------------------------------
    // ======================================================================================================

    /// Computes derived properties of 2D mesh
    fn compute_derived_props_2d(&mut self) -> Result<(), StrError> {
        // auxiliary
        let mut all_shapes: HashMap<(usize, usize), Shape> = HashMap::new();
        let mut all_edges: HashMap<EdgeKey, Edge> = HashMap::new();
        let space_ndim = self.space_ndim;

        // loop over 2D cells
        for cell in &self.cells {
            // check ndim
            let geo_ndim = cell.geo_ndim;
            if geo_ndim != 2 {
                continue; // e.g., 1D line in 2D space
            }

            // shape
            let npoint = cell.points.len();
            let shape_key = (geo_ndim, npoint);
            if !all_shapes.contains_key(&shape_key) {
                all_shapes.insert(shape_key, Shape::new(space_ndim, geo_ndim, npoint)?);
            }
            let shape = all_shapes.get(&shape_key).unwrap();

            // edges (new derived data)
            for e in 0..shape.nedge {
                // collect edge point ids
                let mut edge_points = vec![0; shape.edge_nnode];
                for i in 0..shape.edge_nnode {
                    let local_point_id = shape.get_edge_node_id(e, i);
                    edge_points[i] = cell.points[local_point_id];
                }

                // define key (sorted ids)
                let mut edge_key: EdgeKey = (edge_points[0], edge_points[1]);
                sort2(&mut edge_key);

                // existing edge
                if all_edges.contains_key(&edge_key) {
                    let edge = all_edges.get_mut(&edge_key).unwrap();
                    edge.shared_by_2d_cells.insert(cell.id);

                // new edge
                } else {
                    let mut shared_by_2d_cells = HashSet::new();
                    shared_by_2d_cells.insert(cell.id);
                    all_edges.insert(
                        edge_key,
                        Edge {
                            points: edge_points,
                            shared_by_2d_cells,
                            shared_by_boundary_faces: HashSet::new(),
                        },
                    );
                }
            }
        }

        // find boundary edges and set boundary points
        for (edge_key, edge) in &all_edges {
            if edge.shared_by_2d_cells.len() == 1 {
                self.boundary_edges.insert(*edge_key, edge.clone());
                for point_id in &edge.points {
                    self.boundary_points.insert(*point_id);
                    self.points[*point_id].shared_by_boundary_edges.insert(*edge_key);
                }
            }
        }
        Ok(())
    }

    /// Computes derived properties of 3D mesh
    fn compute_derived_props_3d(&mut self) -> Result<(), StrError> {
        // auxiliary
        let mut all_shapes: HashMap<(usize, usize), Shape> = HashMap::new();
        let mut all_faces: HashMap<FaceKey, Face> = HashMap::new();
        let space_ndim = self.space_ndim;

        // loop over 3D cells
        for cell in &self.cells {
            // check ndim
            let geo_ndim = cell.geo_ndim;
            if geo_ndim != 3 {
                continue; // e.g., 1D line in 3D space or 2D quad in 3D space
            }

            // shape
            let npoint = cell.points.len();
            let shape_key = (geo_ndim, npoint);
            if !all_shapes.contains_key(&shape_key) {
                all_shapes.insert(shape_key, Shape::new(space_ndim, geo_ndim, npoint)?);
            }
            let shape = all_shapes.get(&shape_key).unwrap();

            // faces (new derived data)
            for f in 0..shape.nface {
                // collect face point ids
                let mut face_points = vec![0; shape.face_nnode];
                for i in 0..shape.face_nnode {
                    let local_point_id = shape.get_face_node_id(f, i);
                    face_points[i] = cell.points[local_point_id];
                }

                // define key (sorted ids)
                let mut face_key: FaceKey = if face_points.len() > 3 {
                    (face_points[0], face_points[1], face_points[2], face_points[3])
                } else {
                    (face_points[0], face_points[1], face_points[2], self.points.len())
                };
                sort4(&mut face_key);

                // existing face
                if all_faces.contains_key(&face_key) {
                    let face = all_faces.get_mut(&face_key).unwrap();
                    face.shared_by_cells.insert(cell.id);

                // new face
                } else {
                    let mut shared_by_cells = HashSet::new();
                    shared_by_cells.insert(cell.id);
                    all_faces.insert(
                        face_key,
                        Face {
                            points: face_points,
                            shared_by_cells,
                        },
                    );
                }
            }
        }

        // sort face keys just so the next loop is deterministic
        let mut face_keys: Vec<_> = all_faces.keys().collect();
        face_keys.sort();

        // find boundary faces and set boundary points and boundary edges
        for face_key in face_keys {
            // skip interior faces
            let face = all_faces.get(&face_key).unwrap();
            if face.shared_by_cells.len() != 1 {
                continue;
            }

            // set boundary face
            self.boundary_faces.insert(*face_key, face.clone());

            // set boundary points
            for point_id in &face.points {
                self.boundary_points.insert(*point_id);
                self.points[*point_id].shared_by_boundary_faces.insert(*face_key);
            }

            // face shape
            let face_ndim = 2;
            let face_npoint = face.points.len();
            let face_shape_key = (face_ndim, face_npoint);
            if !all_shapes.contains_key(&face_shape_key) {
                all_shapes.insert(face_shape_key, Shape::new(space_ndim, face_ndim, face_npoint)?);
            }
            let face_shape = all_shapes.get(&face_shape_key).unwrap();

            // boundary edges
            // let face_nedge = face_shape.get_nedge();
            // let face_edge_npoint = face_shape.get_edge_npoint();
            for e in 0..face_shape.nedge {
                // collect edge point ids
                let mut edge_points = vec![0; face_shape.edge_nnode];
                for i in 0..face_shape.edge_nnode {
                    let local_point_id = face_shape.get_edge_node_id(e, i);
                    edge_points[i] = face.points[local_point_id];
                }

                // define key (sorted ids)
                let mut edge_key: EdgeKey = (edge_points[0], edge_points[1]);
                sort2(&mut edge_key);

                // set point information
                for point_id in &edge_points {
                    self.points[*point_id].shared_by_boundary_edges.insert(edge_key);
                }

                // existing edge
                if self.boundary_edges.contains_key(&edge_key) {
                    let edge = self.boundary_edges.get_mut(&edge_key).unwrap();
                    edge.shared_by_boundary_faces.insert(*face_key);

                // new edge
                } else {
                    let mut shared_by_boundary_faces = HashSet::new();
                    shared_by_boundary_faces.insert(*face_key);
                    self.boundary_edges.insert(
                        edge_key,
                        Edge {
                            points: edge_points,
                            shared_by_2d_cells: HashSet::new(),
                            shared_by_boundary_faces,
                        },
                    );
                }
            }
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

    /// Encodes this mesh to a binary object
    fn encode(&self) -> Result<Vec<u8>, StrError> {
        let mut serialized = Vec::new();
        let mut serializer = rmp_serde::Serializer::new(&mut serialized);
        self.serialize(&mut serializer).map_err(|_| "serialize failed")?;
        Ok(serialized)
    }

    /// Decodes a binary object to Mesh
    fn decode(serialized: &Vec<u8>) -> Result<Self, StrError> {
        let mut deserializer = rmp_serde::Deserializer::new(&serialized[..]);
        let res = Deserialize::deserialize(&mut deserializer).map_err(|_| "cannot deserialize data")?;
        Ok(res)
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
        write!(f, "ndim = {}\n", self.space_ndim)?;
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
    use super::*;

    #[test]
    fn new_fails_on_wrong_input() {
        assert_eq!(Mesh::new(1).err(), Some("space_ndim must be 2 or 3"));
        assert_eq!(Mesh::new(4).err(), Some("space_ndim must be 2 or 3"));
    }

    #[test]
    fn new_sized_fails_on_wrong_input() {
        assert_eq!(Mesh::new_sized(1, 2, 1).err(), Some("space_ndim must be 2 or 3"));
        assert_eq!(Mesh::new_sized(4, 2, 1).err(), Some("space_ndim must be 2 or 3"));
        assert_eq!(
            Mesh::new_sized(2, 1, 1).err(),
            Some("npoint must be greater than or equal to 2")
        );
        assert_eq!(
            Mesh::new_sized(2, 2, 0).err(),
            Some("ncell must be greater than or equal to 1")
        );
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
    fn new_sized_works() -> Result<(), StrError> {
        let mesh = Mesh::new_sized(2, 3, 1)?;
        assert_eq!(mesh.space_ndim, 2);
        assert_eq!(mesh.points.len(), 3);
        assert_eq!(mesh.cells.len(), 1);
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
            # ndim npoint ncell
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
            # id att geo_ndim npoint point_ids...
               0   1        2      4 0 1 2 3
               1   0        2      4 1 4 5 2",
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
            "k:(0,1) p:[0, 1] c:[0] f:[]\n\
             k:(0,3) p:[3, 0] c:[0] f:[]\n\
             k:(1,4) p:[1, 4] c:[1] f:[]\n\
             k:(2,3) p:[2, 3] c:[0] f:[]\n\
             k:(2,5) p:[5, 2] c:[1] f:[]\n\
             k:(4,5) p:[4, 5] c:[1] f:[]\n"
        );
        assert_eq!(format!("{}", mesh.string_boundary_faces()), "");
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
            "k:(0,1) p:[1, 0] c:[] f:[(0, 1, 2, 3), (0, 1, 4, 5)]\n\
             k:(0,3) p:[0, 3] c:[] f:[(0, 1, 2, 3), (0, 3, 4, 7)]\n\
             k:(0,4) p:[4, 0] c:[] f:[(0, 1, 4, 5), (0, 3, 4, 7)]\n\
             k:(1,2) p:[2, 1] c:[] f:[(0, 1, 2, 3), (1, 2, 5, 6)]\n\
             k:(1,5) p:[1, 5] c:[] f:[(0, 1, 4, 5), (1, 2, 5, 6)]\n\
             k:(2,3) p:[3, 2] c:[] f:[(0, 1, 2, 3), (2, 3, 6, 7)]\n\
             k:(2,6) p:[2, 6] c:[] f:[(1, 2, 5, 6), (2, 3, 6, 7)]\n\
             k:(3,7) p:[7, 3] c:[] f:[(0, 3, 4, 7), (2, 3, 6, 7)]\n\
             k:(4,5) p:[5, 4] c:[] f:[(0, 1, 4, 5), (4, 5, 8, 9)]\n\
             k:(4,7) p:[4, 7] c:[] f:[(0, 3, 4, 7), (4, 7, 8, 11)]\n\
             k:(4,8) p:[8, 4] c:[] f:[(4, 5, 8, 9), (4, 7, 8, 11)]\n\
             k:(5,6) p:[6, 5] c:[] f:[(1, 2, 5, 6), (5, 6, 9, 10)]\n\
             k:(5,9) p:[5, 9] c:[] f:[(4, 5, 8, 9), (5, 6, 9, 10)]\n\
             k:(6,7) p:[7, 6] c:[] f:[(2, 3, 6, 7), (6, 7, 10, 11)]\n\
             k:(6,10) p:[6, 10] c:[] f:[(5, 6, 9, 10), (6, 7, 10, 11)]\n\
             k:(7,11) p:[11, 7] c:[] f:[(4, 7, 8, 11), (6, 7, 10, 11)]\n\
             k:(8,9) p:[9, 8] c:[] f:[(4, 5, 8, 9), (8, 9, 10, 11)]\n\
             k:(8,11) p:[8, 11] c:[] f:[(4, 7, 8, 11), (8, 9, 10, 11)]\n\
             k:(9,10) p:[10, 9] c:[] f:[(5, 6, 9, 10), (8, 9, 10, 11)]\n\
             k:(10,11) p:[11, 10] c:[] f:[(6, 7, 10, 11), (8, 9, 10, 11)]\n"
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
    fn serialize_works() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let m1 = Mesh::from_text(
            r"# header
            # ndim npoint ncell
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
            # id att geo_ndim npoint point_ids...
               0   1        2      4 0 1 2 3
               1   0        2      4 1 4 5 2",
        )?;
        m1.save("/tmp/gemlab/test.msh")?;
        let m2 = Mesh::read("/tmp/gemlab/test.msh")?;
        assert_eq!(
            format!("{}", m2.grid_boundary_points),
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
        assert_eq!(m2.derived_props_computed, true);
        assert_eq!(
            format!("{}", m2.string_points()),
            "i:0 x:[0.0, 0.0] e:[(0, 1), (0, 3)] f:[]\n\
             i:1 x:[1.0, 0.0] e:[(0, 1), (1, 4)] f:[]\n\
             i:2 x:[1.0, 1.0] e:[(2, 3), (2, 5)] f:[]\n\
             i:3 x:[0.0, 1.0] e:[(0, 3), (2, 3)] f:[]\n\
             i:4 x:[2.0, 0.0] e:[(1, 4), (4, 5)] f:[]\n\
             i:5 x:[2.0, 1.0] e:[(2, 5), (4, 5)] f:[]\n"
        );
        assert_eq!(
            format!("{}", m2.string_cells()),
            "i:0 a:1 g:2 p:[0, 1, 2, 3]\n\
             i:1 a:0 g:2 p:[1, 4, 5, 2]\n"
        );
        assert_eq!(format!("{}", m2.string_boundary_points()), "[0, 1, 2, 3, 4, 5]\n");
        assert_eq!(
            format!("{}", m2.string_boundary_edges()),
            "k:(0,1) p:[0, 1] c:[0] f:[]\n\
             k:(0,3) p:[3, 0] c:[0] f:[]\n\
             k:(1,4) p:[1, 4] c:[1] f:[]\n\
             k:(2,3) p:[2, 3] c:[0] f:[]\n\
             k:(2,5) p:[5, 2] c:[1] f:[]\n\
             k:(4,5) p:[4, 5] c:[1] f:[]\n"
        );
        assert_eq!(format!("{}", m2.string_boundary_faces()), "");
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
        let mesh = Mesh::from_text(
            r"# header
            # ndim npoint ncell
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
            # id att geo_ndim npoint point_ids...
               0   1        2      4 0 1 2 3
               1   0        2      4 1 4 5 2",
        )?;
        assert_eq!(
            format!("{}", mesh),
            "SUMMARY\n\
             =======\n\
             ndim = 2\n\
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
             k:(0,1) p:[0, 1] c:[0] f:[]\n\
             k:(0,3) p:[3, 0] c:[0] f:[]\n\
             k:(1,4) p:[1, 4] c:[1] f:[]\n\
             k:(2,3) p:[2, 3] c:[0] f:[]\n\
             k:(2,5) p:[5, 2] c:[1] f:[]\n\
             k:(4,5) p:[4, 5] c:[1] f:[]\n\
             \n\
             BOUNDARY FACES\n\
             ==============\n"
        );
        Ok(())
    }

    #[test]
    fn find_points_and_edges_work() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let mut m1 = Mesh::from_text(
            r"# header
            # ndim npoint ncell
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
            # id att geo_ndim npoint point_ids...
               0   1        2      4 0 1 2 3
               1   0        2      4 1 4 5 2",
        )?;

        let origin = m1.find_boundary_points(At::XY(0.0, 0.0))?;
        let top_right = m1.find_boundary_points(At::XY(2.0, 1.0))?;
        let bottom = m1.find_boundary_edges(At::Y(0.0))?;
        let right = m1.find_boundary_edges(At::X(2.0))?;
        let top = m1.find_boundary_edges(At::Y(1.0))?;
        let left = m1.find_boundary_edges(At::X(0.0))?;

        assert_eq!(origin, &[0]);
        assert_eq!(top_right, &[5]);
        assert_eq!(bottom, &[(0, 1), (1, 4)]);
        assert_eq!(right, &[(4, 5)]);
        assert_eq!(top, &[(2, 3), (2, 5)]);
        assert_eq!(left, &[(0, 3)]);
        Ok(())
    }
}
