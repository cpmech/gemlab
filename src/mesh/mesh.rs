use super::{At, Geo};
use crate::shapes::{kind_from_ndim_npoint, new_shape, Kind, Shape};
use crate::util::GridSearch;
use crate::StrError;
use russell_lab::{sort2, sort3};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::fmt;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Aliases usize as the group of Point, Edge, Face, or Cell
pub type Group = usize;

/// Aliases usize as the index (==id) of Point and Cell
pub type Index = usize;

/// Aliases (usize,usize) as the key of Edge
pub type EdgeKey = (usize, usize);

/// Aliases (usize,usize,usize) as the key of Face
pub type FaceKey = (usize, usize, usize);

/// Holds point data
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Point {
    /// Identification number which equals the index of the point in the mesh
    ///
    /// **raw data**
    pub id: Index,

    /// Group in which this point belongs to
    ///
    /// **raw data**
    pub group: Group,

    /// Point coordinates (2D or 3D)
    ///
    /// **raw data**
    pub coords: Vec<f64>,

    /// Set of cells sharing this point
    pub shared_by_cells: HashSet<Index>,

    /// Set of edges sharing this point
    pub shared_by_edges: HashSet<EdgeKey>,

    /// Set of faces sharing this point
    pub shared_by_faces: HashSet<FaceKey>,
}

/// Holds edge data (derived data structure)
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Edge {
    /// Group in which this edge belongs to
    pub group: Group,

    /// List of points defining this edge; in the right order (unsorted)
    pub points: Vec<Index>,

    /// Set of cells sharing this edge
    pub shared_by_cells: HashSet<Index>,
}

/// Holds face data (derived data structure)
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Face {
    /// Group in which this face belongs to
    pub group: Group,

    /// List of points defining this face; in the right order (unsorted)
    pub points: Vec<Index>,

    /// Set of cells sharing this edge
    pub shared_by_cells: HashSet<Index>,
}

/// Holds cell (aka geometric shape, polygon, polyhedra) data
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Cell {
    /// Identification number which equals the index of the cell in the mesh
    ///
    /// **raw data**
    pub id: Index,

    /// Group in which this face belongs to
    ///
    /// **raw data**
    pub group: Group,

    /// Space dimension of this cell
    ///
    /// # Note
    ///
    /// The cell's ndim may be different than the space dimension of the mesh.
    /// For example, a 1D line in the 2D or 3D space or a 2D triangle in the 3D space.
    ///
    /// **raw data**
    pub ndim: usize,

    /// List of points defining this cell; in the right order (unsorted)
    ///
    /// **raw data**
    pub points: Vec<Index>,

    /// Set of edges defining the boundary of this cell; in the right order
    pub boundary_edges: HashSet<EdgeKey>,

    /// Set of faces defining the boundary of this cell; in the right order
    pub boundary_faces: HashSet<FaceKey>,
}

/// Holds mesh data
#[derive(Deserialize, Serialize)]
pub struct Mesh {
    /// Space dimension of the mesh
    ///
    /// # Note
    ///
    /// The mesh's ndim may be different that an cell's ndim.
    /// For example, a 3D mesh may contain 1D lines or 2D triangles.
    ///
    /// **raw data**
    pub ndim: usize,

    /// All points in the mesh
    ///
    /// **raw data**
    pub points: Vec<Point>,

    /// All cells (aka geometric shape, polygon, polyhedra) in the mesh
    ///
    /// **raw data**
    pub cells: Vec<Cell>,

    /// Set of points on the boundaries
    pub boundary_points: HashSet<Index>,

    /// Set of edges on the boundaries
    pub boundary_edges: HashMap<EdgeKey, Edge>,

    /// Set of faces on the boundaries
    pub boundary_faces: HashMap<FaceKey, Face>,

    /// Maps group index to a set of points
    pub point_groups: HashMap<Group, HashSet<Index>>,

    /// Maps group index to a set of edges
    pub edge_groups: HashMap<Group, HashSet<EdgeKey>>,

    /// Maps group index to a set of faces
    pub face_groups: HashMap<Group, HashSet<FaceKey>>,

    /// Maps group index to a set of cells
    pub cell_groups: HashMap<Group, HashSet<Index>>,

    /// Min coordinates
    pub min: Vec<f64>,

    /// Max coordinates
    pub max: Vec<f64>,

    /// Allows searching using point coordinates
    pub grid: GridSearch,

    /// Indicates whether the derived variables and maps have been computed or not
    derived_props_computed: bool,
}

impl Mesh {
    pub(crate) fn new(ndim: usize) -> Result<Self, StrError> {
        if ndim < 2 || ndim > 3 {
            return Err("ndim must be 2 or 3");
        }
        Ok(Mesh {
            ndim,
            points: Vec::new(),
            cells: Vec::new(),
            boundary_points: HashSet::new(),
            boundary_edges: HashMap::new(),
            boundary_faces: HashMap::new(),
            point_groups: HashMap::new(),
            cell_groups: HashMap::new(),
            edge_groups: HashMap::new(),
            face_groups: HashMap::new(),
            min: Vec::new(),
            max: Vec::new(),
            grid: GridSearch::new(ndim)?,
            derived_props_computed: false,
        })
    }

    pub(crate) fn new_sized(ndim: usize, npoint: usize, ncell: usize) -> Result<Self, StrError> {
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
            group: 0,
            coords: vec![0.0; ndim],
            shared_by_cells: HashSet::new(),
            shared_by_edges: HashSet::new(),
            shared_by_faces: HashSet::new(),
        };
        let zero_cell = Cell {
            id: 0,
            group: 0,
            ndim: 0,
            points: Vec::new(),
            boundary_edges: HashSet::new(),
            boundary_faces: HashSet::new(),
        };
        Ok(Mesh {
            ndim,
            points: vec![zero_point; npoint],
            cells: vec![zero_cell; ncell],
            boundary_points: HashSet::new(),
            boundary_edges: HashMap::new(),
            boundary_faces: HashMap::new(),
            point_groups: HashMap::new(),
            cell_groups: HashMap::new(),
            edge_groups: HashMap::new(),
            face_groups: HashMap::new(),
            min: Vec::new(),
            max: Vec::new(),
            grid: GridSearch::new(ndim)?,
            derived_props_computed: false,
        })
    }

    pub fn compute_derived_props(&mut self) -> Result<(), StrError> {
        self.derived_props_computed = false;

        // auxiliary maps
        let mut all_shapes: HashMap<Kind, Box<dyn Shape>> = HashMap::new();
        let mut all_edges: HashMap<EdgeKey, Edge> = HashMap::new();
        let mut all_faces: HashMap<FaceKey, Face> = HashMap::new();

        // loop over cells
        for cell in &self.cells {
            // kind and shape
            let ndim = cell.ndim;
            let npoint = cell.points.len();
            let kind = match kind_from_ndim_npoint(ndim, npoint) {
                Some(v) => v,
                None => return Err("cannot find the Kind of cell"),
            };
            if !all_shapes.contains_key(&kind) {
                all_shapes.insert(kind, new_shape(kind));
            }
            let shape = all_shapes.get(&kind).unwrap();

            // points
            for m in 0..npoint {
                let point_id = cell.points[m];
                self.points[point_id].shared_by_cells.insert(cell.id);
            }

            // edges
            let nedge = shape.get_nedge();
            let edge_npoint = shape.get_edge_npoint();
            for e in 0..nedge {
                // collect point ids
                let mut edge_points = vec![0; edge_npoint];
                for i in 0..edge_npoint {
                    let local_point_id = shape.get_edge_local_point_id(e, i);
                    edge_points[i] = cell.points[local_point_id];
                }

                // define key (sorted ids)
                let mut edge_key: EdgeKey = (edge_points[0], edge_points[1]);
                sort2(&mut edge_key);

                // existing item
                if all_edges.contains_key(&edge_key) {
                    let edge = all_edges.get_mut(&edge_key).unwrap();
                    edge.shared_by_cells.insert(cell.id);

                // new item
                } else {
                    let mut shared_by_cells = HashSet::new();
                    shared_by_cells.insert(cell.id);
                    all_edges.insert(
                        edge_key,
                        Edge {
                            group: cell.group,
                            points: edge_points,
                            shared_by_cells,
                        },
                    );
                }
            }

            // faces
            let nface = shape.get_nface();
            let face_npoint = shape.get_face_npoint();
            for f in 0..nface {
                // collect point ids
                let mut face_points = vec![0; face_npoint];
                for i in 0..face_npoint {
                    let local_point_id = shape.get_face_local_point_id(f, i);
                    face_points[i] = cell.points[local_point_id];
                }

                // define key (sorted ids)
                let mut face_key: FaceKey = (face_points[0], face_points[1], face_points[2]);
                sort3(&mut face_key);

                // existing item
                if all_faces.contains_key(&face_key) {
                    let face = all_faces.get_mut(&face_key).unwrap();
                    face.shared_by_cells.insert(cell.id);

                // new item
                } else {
                    let mut shared_by_cells = HashSet::new();
                    shared_by_cells.insert(cell.id);
                    all_faces.insert(
                        face_key,
                        Face {
                            group: cell.group,
                            points: face_points,
                            shared_by_cells,
                        },
                    );
                }
            }
        }

        println!("ALL EDGES");
        println!("{:?}", all_edges);
        println!("\nALL FACES");
        println!("{:?}", all_faces);

        // limits
        self.min = vec![f64::MAX; self.ndim];
        self.max = vec![f64::MIN; self.ndim];
        for point in &self.points {
            for i in 0..self.ndim {
                if point.coords[i] < self.min[i] {
                    self.min[i] = point.coords[i];
                }
                if point.coords[i] > self.max[i] {
                    self.max[i] = point.coords[i];
                }
            }
        }
        for i in 0..self.ndim {
            if self.min[i] >= self.max[i] {
                return Err("mesh limits are invalid");
            }
        }

        self.derived_props_computed = true;
        Ok(())
    }

    pub fn from_serialized(serialized: &Vec<u8>) -> Result<Self, StrError> {
        let mut deserializer = rmp_serde::Deserializer::new(&serialized[..]);
        let res = Deserialize::deserialize(&mut deserializer).map_err(|_| "cannot deserialize data")?;
        Ok(res)
    }

    pub fn get_serialized(&self) -> Result<Vec<u8>, StrError> {
        let mut serialized = Vec::new();
        let mut serializer = rmp_serde::Serializer::new(&mut serialized);
        self.serialize(&mut serializer).map_err(|_| "serialize failed")?;
        Ok(serialized)
    }

    pub fn write_serialized<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }
        let serialized = self.get_serialized()?;
        let mut file = File::create(path).map_err(|_| "cannot create file")?;
        file.write_all(&serialized).map_err(|_| "cannot write file")?;
        Ok(())
    }

    pub fn set_group(&mut self, group: Group, geo: Geo) -> Result<(), StrError> {
        if !self.derived_props_computed {
            return Err("compute_derived_props must be called first");
        }

        match geo {
            Geo::Point(at) => {
                // find all points near the geometric feature
                for index in self.find_points(&at)? {
                    // update map
                    if self.point_groups.contains_key(&group) {
                        let indices = self.point_groups.get_mut(&group).unwrap();
                        indices.insert(index);
                    // new map
                    } else {
                        let mut indices = HashSet::new();
                        indices.insert(index);
                        self.point_groups.insert(group, indices);
                    }
                }
            }
            Geo::Edge(at) => {
                // find all points near the geometric feature
                let indices = &self.find_points(&at)?;
                for index in indices {
                    // loop over all edges touching this point
                    let point = &self.points[*index];
                    for key in &point.shared_by_edges {
                        // check if two edge points pass through the geometric feature
                        if indices.contains(&key.0) && indices.contains(&key.1) {
                            // update map
                            if self.edge_groups.contains_key(&group) {
                                let keys = self.edge_groups.get_mut(&group).unwrap();
                                keys.insert(key.clone());
                            // new map
                            } else {
                                let mut keys = HashSet::new();
                                keys.insert(key.clone());
                                self.edge_groups.insert(group, keys);
                            }
                        }
                    }
                }
            }
            Geo::Face(_) => {
                // todo
            }
        };
        Ok(())
    }

    fn find_points(&mut self, at: &At) -> Result<HashSet<Index>, StrError> {
        if !self.derived_props_computed {
            return Err("compute_derived_props must be called first");
        }

        let mut points: HashSet<Index> = HashSet::new();
        match at {
            At::X(x) => {
                if self.ndim == 2 {
                    for id in self.grid.find_on_line(&[*x, 0.0], &[*x, 1.0]).unwrap() {
                        points.insert(id);
                    }
                } else {
                    for id in self.grid.find_on_plane_yz(*x).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::Y(y) => {
                if self.ndim == 2 {
                    for id in self.grid.find_on_line(&[0.0, *y], &[1.0, *y]).unwrap() {
                        points.insert(id);
                    }
                } else {
                    for id in self.grid.find_on_plane_xz(*y).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::Z(z) => {
                if self.ndim == 2 {
                    return Err("At::Z works in 3D only");
                } else {
                    for id in self.grid.find_on_plane_xy(*z).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::XY(x, y) => {
                if self.ndim == 2 {
                    if let Some(id) = self.grid.find(&[*x, *y]).unwrap() {
                        points.insert(id);
                    }
                } else {
                    for id in self.grid.find_on_line(&[*x, *y, 0.0], &[*x, *y, 1.0]).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::YZ(y, z) => {
                if self.ndim == 2 {
                    return Err("At::YZ works in 3D only");
                } else {
                    for id in self.grid.find_on_line(&[0.0, *y, *z], &[1.0, *y, *z]).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::XZ(x, z) => {
                if self.ndim == 2 {
                    return Err("At::XZ works in 3D only");
                } else {
                    for id in self.grid.find_on_line(&[*x, 0.0, *z], &[*x, 1.0, *z]).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::XYZ(x, y, z) => {
                if self.ndim == 2 {
                    return Err("At::XYZ works in 3D only");
                } else {
                    if let Some(id) = self.grid.find(&[*x, *y, *z]).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::Circle(x, y, r) => {
                if self.ndim == 2 {
                    for id in self.grid.find_on_circle(&[*x, *y], *r).unwrap() {
                        points.insert(id);
                    }
                } else {
                    return Err("At::Circle works in 2D only");
                }
            }
            At::Cylinder(ax, ay, az, bx, by, bz, r) => {
                if self.ndim == 2 {
                    return Err("At::Cylinder works in 3D only");
                } else {
                    for id in self
                        .grid
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
}

impl fmt::Display for Mesh {
    /// Prints mesh data (may be large)
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // summary
        write!(f, "ndim = {}\n", self.ndim).unwrap();
        write!(f, "npoint = {}\n", self.points.len()).unwrap();
        write!(f, "ncell = {}\n", self.cells.len()).unwrap();
        write!(f, "n_boundary_point = {}\n", self.boundary_points.len()).unwrap();
        write!(f, "n_boundary_edge = {}\n", self.boundary_edges.len()).unwrap();
        write!(f, "n_boundary_face = {}\n", self.boundary_faces.len()).unwrap();

        // points: i=index, g=group, x=coordinates, c=shared_by_cell_ids
        write!(f, "\npoints\n").unwrap();
        for point in &self.points {
            let mut shared_by_cells: Vec<_> = point.shared_by_cells.iter().collect();
            shared_by_cells.sort();
            write!(
                f,
                "i:{} g:{} x:{:?} c:{:?}\n",
                point.id, point.group, point.coords, shared_by_cells
            )
            .unwrap();
        }

        // cells: i=index, g=group, p=points, e=boundary_edge_ids, f=boundary_face_ids
        write!(f, "\ncells\n").unwrap();
        for cell in &self.cells {
            let mut boundary_edges: Vec<_> = cell.boundary_edges.iter().collect();
            let mut boundary_faces: Vec<_> = cell.boundary_faces.iter().collect();
            boundary_edges.sort();
            boundary_faces.sort();
            write!(
                f,
                "i:{} g:{} n:{} p:{:?} e:{:?} f:{:?}\n",
                cell.id, cell.group, cell.ndim, cell.points, boundary_edges, boundary_faces
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

        // boundary edges: i=index, g=group, k=key(point,point), p=point_ids c=shared_by_cell_ids
        write!(f, "\nboundary_edges\n").unwrap();
        let mut keys_and_edges: Vec<_> = self.boundary_edges.keys().zip(self.boundary_edges.values()).collect();
        keys_and_edges.sort_by(|left, right| left.0.cmp(&right.0));
        for (key, edge) in keys_and_edges {
            let mut shared_by_cells: Vec<_> = edge.shared_by_cells.iter().collect();
            shared_by_cells.sort();
            write!(
                f,
                "g:{} k:({},{}) p:{:?} c:{:?}\n",
                edge.group, key.0, key.1, edge.points, shared_by_cells
            )
            .unwrap();
        }

        // boundary_faces: i=index, g=group, k=key(point,point,point), p=point_ids, c=shared_by_cell_ids
        write!(f, "\nboundary_faces\n").unwrap();
        let mut keys_and_faces: Vec<_> = self.boundary_faces.keys().zip(self.boundary_faces.values()).collect();
        keys_and_faces.sort_by(|left, right| left.0.cmp(&right.0));
        for (key, face) in keys_and_faces {
            let mut shared_by_cells: Vec<_> = face.shared_by_cells.iter().collect();
            shared_by_cells.sort();
            write!(
                f,
                "g:{} k:({},{},{}) p:{:?} c:{:?}\n",
                face.group, key.0, key.1, key.2, face.points, shared_by_cells
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
    use crate::mesh::parse_mesh;

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
        mesh.compute_derived_props()?;
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

        mesh.set_group(1, Geo::Point(At::XY(0.0, 0.0)))?;
        mesh.set_group(2, Geo::Point(At::X(0.0)))?;
        mesh.set_group(3, Geo::Point(At::X(1.0)))?;
        mesh.set_group(4, Geo::Point(At::Y(0.0)))?;
        mesh.set_group(5, Geo::Point(At::Y(1.0)))?;

        println!("{}", mesh);

        Ok(())
    }
}
