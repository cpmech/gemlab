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

pub type Group = usize;
pub type Index = usize;
pub type CellId = usize;
pub type EdgeKey = (usize, usize);
pub type FaceKey = (usize, usize, usize);

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Point {
    pub id: Index,
    pub group: Group,
    pub coords: Vec<f64>,
    pub shared_by_cells: HashSet<Index>,
    pub shared_by_edges: HashSet<EdgeKey>,
    pub shared_by_faces: HashSet<FaceKey>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Edge {
    pub group: Group,
    pub points: Vec<Index>, // in the right order (unsorted)
    pub shared_by_cells: HashSet<Index>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Face {
    pub group: Group,
    pub points: Vec<Index>, // in the right order (unsorted)
    pub shared_by_cells: HashSet<Index>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Cell {
    pub id: Index,
    pub group: Group,
    pub ndim: usize,
    pub points: Vec<Index>,
    pub boundary_edges: HashSet<EdgeKey>,
    pub boundary_faces: HashSet<FaceKey>,
}

#[derive(Deserialize, Serialize)]
pub struct Mesh {
    pub ndim: usize,
    pub points: Vec<Point>,
    pub cells: Vec<Cell>,
    pub boundary_points: HashSet<Index>,
    pub boundary_edges: HashMap<EdgeKey, Edge>,
    pub boundary_faces: HashMap<FaceKey, Face>,
    pub point_groups: HashMap<String, HashSet<Index>>,
    pub cell_groups: HashMap<String, HashSet<Index>>,
    pub edge_groups: HashMap<String, HashSet<EdgeKey>>,
    pub face_groups: HashMap<String, HashSet<FaceKey>>,
    pub min: Vec<f64>,
    pub max: Vec<f64>,
    pub grid: GridSearch,
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
        })
    }

    pub(crate) fn new_zeroed(ndim: usize, npoint: usize, ncell: usize) -> Result<Self, StrError> {
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
        })
    }

    pub(crate) fn compute_derived_props(&mut self) -> Result<(), StrError> {
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

    pub fn set_group(&mut self, name: &str, geo: Geo) -> Result<(), StrError> {
        match geo {
            Geo::Point(at) => {
                // find all points near the geometric feature
                for index in self.find_points(&at) {
                    // new map
                    if self.point_groups.contains_key(name) {
                        let group = self.point_groups.get_mut(name).unwrap();
                        group.insert(index);
                    // update map
                    } else {
                        let mut indices = HashSet::new();
                        indices.insert(index);
                        self.point_groups.insert(name.to_string(), indices);
                    }
                }
            }
            Geo::Edge(at) => {
                // find all points near the geometric feature
                let points = self.find_points(&at);
                for index in &points {
                    // loop over all edges touching this point
                    let point = &self.points[*index];
                    for key in &point.shared_by_edges {
                        // check if two edge points pass through the geometric feature
                        if points.contains(&key.0) && points.contains(&key.1) {
                            // new map
                            if self.edge_groups.contains_key(name) {
                                let group = self.edge_groups.get_mut(name).unwrap();
                                group.insert(key.clone());
                            // update map
                            } else {
                                let mut keys = HashSet::new();
                                keys.insert(key.clone());
                                self.edge_groups.insert(name.to_string(), keys);
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

    fn find_points(&mut self, at: &At) -> HashSet<Index> {
        let mut points: HashSet<Index> = HashSet::new();
        match at {
            At::XY(x, y) => {
                if self.ndim == 2 {
                    if let Some(id) = self.grid.find(&[*x, *y]).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::XYZ(x, y, z) => {
                if self.ndim == 3 {
                    if let Some(id) = self.grid.find(&[*x, *y, *z]).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::Horizontal(x, y) => {
                if self.ndim == 2 {
                    let ids = self.grid.find_along_segment(&[*x, *y], &[*x + 1.0, *y]).unwrap();
                    for id in ids {
                        points.insert(id);
                    }
                }
            }
            At::Vertical(x, y) => {
                if self.ndim == 2 {
                    for id in self.grid.find_along_segment(&[*x, *y], &[*x, *y + 1.0]).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::Circle(x, y, r) => {
                if self.ndim == 2 {
                    for id in self.grid.find_along_circumference(&[*x, *y], *r).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::Circle3D(x, y, z, r) => {
                if self.ndim == 3 {
                    for id in self.grid.find_along_circumference(&[*x, *y, *z], *r).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::CylinderX(x, y, z, r) => {
                if self.ndim == 3 {
                    for id in self.grid.find_along_cylinder_x(&[*x, *y, *z], *r).unwrap() {
                        points.insert(id);
                    }
                }
            }
            At::CylinderY(x, y, z, r) => {
                println!("{} {} {} {}", x, y, z, r);
            }
            At::CylinderZ(x, y, z, r) => {
                println!("{} {} {} {}", x, y, z, r);
            }
            At::PlaneNormalX(x, y, z) => {
                println!("{} {} {}", x, y, z);
            }
            At::PlaneNormalY(x, y, z) => {
                println!("{} {} {}", x, y, z);
            }
            At::PlaneNormalZ(x, y, z) => {
                println!("{} {} {}", x, y, z);
            }
        }
        points
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
    use std::collections::{HashMap, HashSet};

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
        let mut mesh = Mesh {
            ndim: 2,
            points: vec![
                Point {
                    id: 0,
                    group: 1,
                    coords: vec![0.0, 0.0],
                    shared_by_cells: HashSet::new(),
                    shared_by_edges: HashSet::new(),
                    shared_by_faces: HashSet::new(),
                },
                Point {
                    id: 1,
                    group: 1,
                    coords: vec![1.0, 0.0],
                    shared_by_cells: HashSet::new(),
                    shared_by_edges: HashSet::new(),
                    shared_by_faces: HashSet::new(),
                },
                Point {
                    id: 2,
                    group: 1,
                    coords: vec![1.0, 1.0],
                    shared_by_cells: HashSet::new(),
                    shared_by_edges: HashSet::new(),
                    shared_by_faces: HashSet::new(),
                },
                Point {
                    id: 3,
                    group: 1,
                    coords: vec![0.0, 1.0],
                    shared_by_cells: HashSet::new(),
                    shared_by_edges: HashSet::new(),
                    shared_by_faces: HashSet::new(),
                },
            ],
            cells: vec![Cell {
                id: 0,
                group: 1,
                ndim: 2,
                points: vec![0, 1, 2, 3],
                boundary_edges: HashSet::new(),
                boundary_faces: HashSet::new(),
            }],
            boundary_points: HashSet::new(),
            boundary_edges: HashMap::new(),
            boundary_faces: HashMap::new(),
            point_groups: HashMap::new(),
            cell_groups: HashMap::new(),
            edge_groups: HashMap::new(),
            face_groups: HashMap::new(),
            min: Vec::new(),
            max: Vec::new(),
            grid: GridSearch::new(2)?,
        };
        mesh.compute_derived_props()?;
        Ok(())
    }

    #[test]
    fn display_works() {
        // todo
    }

    #[test]
    fn set_group_works() -> Result<(), StrError> {
        return Ok(());

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

        mesh.set_group("origin", Geo::Point(At::XY(0.0, 0.0)))?;
        mesh.set_group("origin", Geo::Point(At::Horizontal(0.0, 0.0)))?;
        mesh.set_group("origin", Geo::Point(At::Vertical(0.0, 0.0)))?;
        mesh.set_group("origin", Geo::Point(At::XYZ(0.0, 0.0, 0.0)))?;
        mesh.set_group("bedrock", Geo::Edge(At::Horizontal(0.0, 0.0)))?;
        mesh.set_group("left", Geo::Edge(At::Vertical(0.0, 0.0)))?;
        mesh.set_group("right", Geo::Edge(At::Vertical(0.0, 0.0)))?;
        mesh.set_group("tunnel", Geo::Face(At::Circle(0.0, 0.0, 0.5)))?;
        mesh.set_group("tunnel", Geo::Face(At::Circle3D(0.0, 0.0, 0.0, 0.5)))?;
        mesh.set_group("tunnel", Geo::Face(At::CylinderX(0.0, 0.0, 0.0, 0.5)))?;
        mesh.set_group("bedrock", Geo::Face(At::PlaneNormalZ(0.0, 0.0, 0.0)))?;

        println!("{}", mesh);

        Ok(())
    }
}
