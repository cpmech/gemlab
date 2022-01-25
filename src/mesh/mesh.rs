use super::read_mesh::{parse_mesh, read_mesh};
use super::At;
use crate::shapes::Shape;
use crate::util::GridSearch;
use crate::StrError;
use russell_lab::{sort2, sort3};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::fmt;
use std::fs;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

/// Aliases usize as Point ID
pub type PointId = usize;

/// Aliases usize as Cell ID
pub type CellId = usize;

/// Aliases (usize,usize) as the key of Edge
pub type EdgeKey = (usize, usize);

/// Aliases (usize,usize,usize) as the key of Face
pub type FaceKey = (usize, usize, usize);

/// Aliases usize as the group of boundary Point, Edge or Face
pub type Group = String;

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
    pub attribute_id: usize,

    /// Space dimension of this cell
    ///
    /// The cell's ndim may be different than the space dimension of the mesh.
    /// For example, a 1D line in the 2D or 3D space or a 2D triangle in the 3D space.
    ///
    /// **raw data**
    pub geo_ndim: usize,

    /// List of points defining this cell; in the right order (unsorted)
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

    /// Collects all boundary point groups
    ///
    /// (derived property)
    pub groups_boundary_points: HashMap<Group, HashSet<PointId>>,

    /// Collects all boundary edge groups
    ///
    /// (derived property)
    pub groups_boundary_edges: HashMap<Group, HashSet<EdgeKey>>,

    /// Collects all boundary face groups
    ///
    /// (derived property)
    pub groups_boundary_faces: HashMap<Group, HashSet<FaceKey>>,

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
    pub(super) fn new(ndim: usize) -> Result<Self, StrError> {
        if ndim < 2 || ndim > 3 {
            return Err("ndim must be 2 or 3");
        }
        Ok(Mesh {
            space_ndim: ndim,
            points: Vec::new(),
            cells: Vec::new(),
            boundary_points: HashSet::new(),
            boundary_edges: HashMap::new(),
            boundary_faces: HashMap::new(),
            groups_boundary_points: HashMap::new(),
            groups_boundary_edges: HashMap::new(),
            groups_boundary_faces: HashMap::new(),
            min: Vec::new(),
            max: Vec::new(),
            grid_boundary_points: GridSearch::new(ndim)?,
            derived_props_computed: false,
        })
    }

    /// Returns a new mesh with pre-allocated (empty) Point and Cell vectors
    pub(super) fn new_sized(ndim: usize, npoint: usize, ncell: usize) -> Result<Self, StrError> {
        if ndim < 2 || ndim > 3 {
            return Err("ndim must be 2 or 3");
        }
        if npoint < 2 {
            return Err("npoint must be greater than or equal to 2");
        }
        if ncell < 1 {
            return Err("ncell must be greater than or equal to 1");
        }
        let zero_point = Point {
            id: 0,
            coords: vec![0.0; ndim],
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
            space_ndim: ndim,
            points: vec![zero_point; npoint],
            cells: vec![zero_cell; ncell],
            boundary_points: HashSet::new(),
            boundary_edges: HashMap::new(),
            boundary_faces: HashMap::new(),
            groups_boundary_points: HashMap::new(),
            groups_boundary_edges: HashMap::new(),
            groups_boundary_faces: HashMap::new(),
            min: Vec::new(),
            max: Vec::new(),
            grid_boundary_points: GridSearch::new(ndim)?,
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
        self.derived_props_computed = true;
        Ok(())
    }

    /// Clears all boundary groups
    pub fn clear_boundary_groups(&mut self) {
        self.groups_boundary_points.clear();
        self.groups_boundary_edges.clear();
        self.groups_boundary_faces.clear();
    }

    /// Sets group of points on the boundary
    pub fn set_points(&mut self, group: &str, at: At) -> Result<&mut Self, StrError> {
        if !self.derived_props_computed {
            return Err("compute_derived_props must be called first");
        }
        // find all points near the geometric feature
        for id in self.find_points(&at)? {
            // update groups map
            if self.groups_boundary_points.contains_key(group) {
                let ids = self.groups_boundary_points.get_mut(group).unwrap();
                ids.insert(id);
            // new groups map
            } else {
                let mut ids = HashSet::new();
                ids.insert(id);
                self.groups_boundary_points.insert(group.to_string(), ids);
            }
        }
        Ok(self)
    }

    /// Sets group of edges on the boundary
    pub fn set_edges(&mut self, group: &str, at: At) -> Result<&mut Self, StrError> {
        if !self.derived_props_computed {
            return Err("compute_derived_props must be called first");
        }
        // find all points near the geometric feature
        let indices = &self.find_points(&at)?;
        for id in indices {
            // loop over all boundary edges touching this point
            let point = &self.points[*id];
            for key in &point.shared_by_boundary_edges {
                // check if two edge points pass through the geometric feature
                if indices.contains(&key.0) && indices.contains(&key.1) {
                    if self.boundary_edges.contains_key(&key) {
                        // update groups map
                        if self.groups_boundary_edges.contains_key(group) {
                            let keys = self.groups_boundary_edges.get_mut(group).unwrap();
                            keys.insert(key.clone());
                        // new groups map
                        } else {
                            let mut keys = HashSet::new();
                            keys.insert(key.clone());
                            self.groups_boundary_edges.insert(group.to_string(), keys);
                        }
                    }
                }
            }
        }
        Ok(self)
    }

    pub fn get_boundary_point_ids_sorted(&self, group: &str) -> Vec<PointId> {
        if let Some(points) = self.groups_boundary_points.get(group) {
            let mut ids: Vec<_> = points.iter().map(|id| *id).collect();
            ids.sort();
            return ids;
        }
        Vec::new()
    }

    pub fn get_boundary_edge_keys_sorted(&self, group: &str) -> Vec<EdgeKey> {
        if let Some(edges) = self.groups_boundary_edges.get(group) {
            let mut keys: Vec<_> = edges.iter().map(|key| *key).collect();
            keys.sort();
            return keys;
        }
        Vec::new()
    }

    pub fn extract_coords(&self, shape: &mut Shape, cell_id: usize) -> Result<(), StrError> {
        let npoint = self.cells[cell_id].points.len();
        for m in 0..npoint {
            let point_id = self.cells[cell_id].points[m];
            for i in 0..self.space_ndim {
                shape.set_point(m, i, self.points[point_id].coords[i])?;
            }
        }
        Ok(())
    }

    // ======================================================================================================
    // --- private ------------------------------------------------------------------------------------------
    // ======================================================================================================

    /// Finds points in the mesh
    fn find_points(&mut self, at: &At) -> Result<HashSet<PointId>, StrError> {
        if !self.derived_props_computed {
            return Err("compute_derived_props must be called first");
        }

        let mut points: HashSet<PointId> = HashSet::new();
        match at {
            At::X(x) => {
                if self.space_ndim == 2 {
                    for id in self.grid_boundary_points.find_on_line(&[*x, 0.0], &[*x, 1.0]).unwrap() {
                        points.insert(id);
                    }
                } else {
                    for id in self.grid_boundary_points.find_on_plane_yz(*x).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::Y(y) => {
                if self.space_ndim == 2 {
                    for id in self.grid_boundary_points.find_on_line(&[0.0, *y], &[1.0, *y]).unwrap() {
                        points.insert(id);
                    }
                } else {
                    for id in self.grid_boundary_points.find_on_plane_xz(*y).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::Z(z) => {
                if self.space_ndim == 2 {
                    return Err("At::Z works in 3D only");
                } else {
                    for id in self.grid_boundary_points.find_on_plane_xy(*z).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::XY(x, y) => {
                if self.space_ndim == 2 {
                    if let Some(id) = self.grid_boundary_points.find(&[*x, *y]).unwrap() {
                        points.insert(id);
                    }
                } else {
                    for id in self
                        .grid_boundary_points
                        .find_on_line(&[*x, *y, 0.0], &[*x, *y, 1.0])
                        .unwrap()
                    {
                        points.insert(id);
                    }
                }
            }
            At::YZ(y, z) => {
                if self.space_ndim == 2 {
                    return Err("At::YZ works in 3D only");
                } else {
                    for id in self
                        .grid_boundary_points
                        .find_on_line(&[0.0, *y, *z], &[1.0, *y, *z])
                        .unwrap()
                    {
                        points.insert(id);
                    }
                }
            }
            At::XZ(x, z) => {
                if self.space_ndim == 2 {
                    return Err("At::XZ works in 3D only");
                } else {
                    for id in self
                        .grid_boundary_points
                        .find_on_line(&[*x, 0.0, *z], &[*x, 1.0, *z])
                        .unwrap()
                    {
                        points.insert(id);
                    }
                }
            }
            At::XYZ(x, y, z) => {
                if self.space_ndim == 2 {
                    return Err("At::XYZ works in 3D only");
                } else {
                    if let Some(id) = self.grid_boundary_points.find(&[*x, *y, *z]).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::Circle(x, y, r) => {
                if self.space_ndim == 2 {
                    for id in self.grid_boundary_points.find_on_circle(&[*x, *y], *r).unwrap() {
                        points.insert(id);
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
                        .find_on_cylinder(&[*ax, *ay, *az], &[*bx, *by, *bz], *r)
                        .unwrap()
                    {
                        points.insert(id);
                    }
                }
            }
        }
        Ok(points)
    }

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
                let mut edge_points = vec![0; shape.edge_npoint];
                for i in 0..shape.edge_npoint {
                    let local_point_id = shape.get_edge_point_id(e, i);
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
        for (key, edge) in &all_edges {
            if edge.shared_by_2d_cells.len() == 1 {
                self.boundary_edges.insert(*key, edge.clone());
                for id in &edge.points {
                    self.boundary_points.insert(*id);
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
                let mut face_points = vec![0; shape.face_npoint];
                for i in 0..shape.face_npoint {
                    let local_point_id = shape.get_face_point_id(f, i);
                    face_points[i] = cell.points[local_point_id];
                }

                // define key (sorted ids)
                let mut face_key: FaceKey = (face_points[0], face_points[1], face_points[2]);
                sort3(&mut face_key);

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

        // find boundary faces and set boundary points and boundary edges
        for (face_key, face) in &all_faces {
            // skip interior faces
            if face.shared_by_cells.len() != 1 {
                continue;
            }

            // set boundary face
            self.boundary_faces.insert(*face_key, face.clone());

            // set boundary points
            for id in &face.points {
                self.boundary_points.insert(*id);
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
                let mut edge_points = vec![0; face_shape.edge_npoint];
                for i in 0..face_shape.edge_npoint {
                    let local_point_id = face_shape.get_edge_point_id(e, i);
                    edge_points[i] = face.points[local_point_id];
                }

                // define key (sorted ids)
                let mut edge_key: EdgeKey = (edge_points[0], edge_points[1]);
                sort2(&mut edge_key);

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
}

impl fmt::Display for Mesh {
    /// Prints mesh data (may be large)
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // summary
        write!(f, "ndim = {}\n", self.space_ndim).unwrap();
        write!(f, "npoint = {}\n", self.points.len()).unwrap();
        write!(f, "ncell = {}\n", self.cells.len()).unwrap();
        write!(f, "n_boundary_point = {}\n", self.boundary_points.len()).unwrap();
        write!(f, "n_boundary_edge = {}\n", self.boundary_edges.len()).unwrap();
        write!(f, "n_boundary_face = {}\n", self.boundary_faces.len()).unwrap();

        // points: i=index, g=group, x=coordinates, e=shared_by_boundary_edges, f=shared_by_boundary_faces
        write!(f, "\npoints\n").unwrap();
        for point in &self.points {
            let mut shared_by_boundary_edges: Vec<_> = point.shared_by_boundary_edges.iter().collect();
            let mut shared_by_boundary_faces: Vec<_> = point.shared_by_boundary_faces.iter().collect();
            shared_by_boundary_edges.sort();
            shared_by_boundary_faces.sort();
            write!(
                f,
                "i:{} x:{:?} e:{:?} f:{:?}\n",
                point.id, point.coords, shared_by_boundary_edges, shared_by_boundary_faces,
            )
            .unwrap();
        }

        // cells: i=index, a=attribute_id, p=points
        write!(f, "\ncells\n").unwrap();
        for cell in &self.cells {
            write!(
                f,
                "i:{} a:{} n:{} p:{:?}\n",
                cell.id, cell.attribute_id, cell.geo_ndim, cell.points,
            )
            .unwrap();
        }

        // boundary points: i=index
        write!(f, "\nboundary_points\n").unwrap();
        let mut point_indices: Vec<_> = self.boundary_points.iter().collect();
        point_indices.sort();
        for index in point_indices {
            write!(f, "{} ", index).unwrap();
        }
        if self.boundary_points.len() > 0 {
            write!(f, "\n").unwrap();
        }

        // boundary edges: k=key(point,point), p=point_ids c=shared_by_2d_cells f=shared_by_boundary_faces
        write!(f, "\nboundary_edges\n").unwrap();
        let mut keys_and_edges: Vec<_> = self.boundary_edges.keys().zip(self.boundary_edges.values()).collect();
        keys_and_edges.sort_by(|left, right| left.0.cmp(&right.0));
        for (key, edge) in keys_and_edges {
            let mut shared_by_2d_cells: Vec<_> = edge.shared_by_2d_cells.iter().collect();
            let mut shared_by_boundary_faces: Vec<_> = edge.shared_by_boundary_faces.iter().collect();
            shared_by_2d_cells.sort();
            shared_by_boundary_faces.sort();
            write!(
                f,
                "k:({},{}) p:{:?} c:{:?} f:{:?}\n",
                key.0, key.1, edge.points, shared_by_2d_cells, shared_by_boundary_faces,
            )
            .unwrap();
        }

        // boundary_faces: k=key(point,point,point), p=point_ids, c=shared_by_cells
        write!(f, "\nboundary_faces\n").unwrap();
        let mut keys_and_faces: Vec<_> = self.boundary_faces.keys().zip(self.boundary_faces.values()).collect();
        keys_and_faces.sort_by(|left, right| left.0.cmp(&right.0));
        for (key, face) in keys_and_faces {
            let mut shared_by_cells: Vec<_> = face.shared_by_cells.iter().collect();
            shared_by_cells.sort();
            write!(
                f,
                "k:({},{},{}) p:{:?} c:{:?}\n",
                key.0, key.1, key.2, face.points, shared_by_cells
            )
            .unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_works() {
        // todo
    }

    #[test]
    fn serialize_works() {
        // todo
    }

    #[test]
    fn compute_derived_props_works() -> Result<(), StrError> {
        let mesh = Mesh::from_text(
            r"
            2 4 1
            0  0.0 0.0
            1  1.0 0.0
            2  1.0 1.0
            3  0.0 1.0
            0 0  2 4  0 1 2 3
        ",
        )?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 2\n\
             npoint = 4\n\
             ncell = 1\n\
             n_boundary_point = 4\n\
             n_boundary_edge = 4\n\
             n_boundary_face = 0\n\
             \n\
             points\n\
             i:0 x:[0.0, 0.0] e:[] f:[]\n\
             i:1 x:[1.0, 0.0] e:[] f:[]\n\
             i:2 x:[1.0, 1.0] e:[] f:[]\n\
             i:3 x:[0.0, 1.0] e:[] f:[]\n\
             \n\
             cells\n\
             i:0 a:0 n:2 p:[0, 1, 2, 3]\n\
             \n\
             boundary_points\n\
             0 1 2 3 \n\
             \n\
             boundary_edges\n\
             k:(0,1) p:[0, 1] c:[0] f:[]\n\
             k:(0,3) p:[3, 0] c:[0] f:[]\n\
             k:(1,2) p:[1, 2] c:[0] f:[]\n\
             k:(2,3) p:[2, 3] c:[0] f:[]\n\
             \n\
             boundary_faces\n"
        );
        Ok(())
    }

    #[test]
    fn display_works() {
        // todo
    }

    #[test]
    fn set_group_works() -> Result<(), StrError> {
        let mut mesh = parse_mesh(
            r"
            2 4 1
            0 1 0.0 0.0
            1 1 1.0 0.0
            2 1 1.0 1.0
            3 1 0.0 1.0
            0 1  2 4  0 1 2 3
        ",
        )?;

        mesh.set_points("origin", At::XY(0.0, 0.0))?
            .set_edges("left", At::X(0.0))?
            .set_edges("right", At::X(1.0))?
            .set_edges("bottom", At::Y(0.0))?
            .set_edges("top", At::Y(1.0))?;

        println!("{}", mesh);

        Ok(())
    }
}
