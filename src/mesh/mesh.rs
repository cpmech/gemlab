use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fmt;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Point {
    pub id: usize,
    pub group: usize,
    pub coords: Vec<f64>,
    pub shared_by_cell_ids: Vec<usize>,
    // pub shared_by_boundary_edge_ids: Vec<usize>,
    // pub shared_by_boundary_face_ids: Vec<usize>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Edge {
    pub id: usize,
    pub group: usize,
    pub point_ids: Vec<usize>, // in the right order (unsorted)
    pub shared_by_cell_ids: Vec<usize>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Face {
    pub id: usize,
    pub group: usize,
    pub point_ids: Vec<usize>, // in the right order (unsorted)
    pub shared_by_cell_ids: Vec<usize>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Cell {
    pub id: usize,
    pub group: usize,
    pub point_ids: Vec<usize>,
    pub boundary_edge_ids: Vec<usize>,
    pub boundary_face_ids: Vec<usize>,
}

#[derive(Deserialize, Serialize)]
pub struct Mesh {
    pub ndim: usize,
    pub points: Vec<Point>,
    pub cells: Vec<Cell>,
    pub boundary_points: HashMap<usize, bool>,
    pub boundary_edges: HashMap<(usize, usize), Edge>, // the key is sorted
    pub boundary_faces: HashMap<(usize, usize, usize), Face>, // the key is sorted
}

impl Mesh {
    pub fn new(ndim: usize) -> Result<Self, &'static str> {
        if ndim < 2 || ndim > 3 {
            return Err("ndim must be 2 or 3");
        }
        Ok(Mesh {
            ndim,
            points: Vec::new(),
            cells: Vec::new(),
            boundary_points: HashMap::new(),
            boundary_edges: HashMap::new(),
            boundary_faces: HashMap::new(),
        })
    }

    pub fn new_zeroed(ndim: usize, npoint: usize, ncell: usize) -> Result<Self, &'static str> {
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
            shared_by_cell_ids: Vec::new(),
        };
        let zero_cell = Cell {
            id: 0,
            group: 0,
            point_ids: Vec::new(),
            boundary_edge_ids: Vec::new(),
            boundary_face_ids: Vec::new(),
        };
        Ok(Mesh {
            ndim,
            points: vec![zero_point; npoint],
            cells: vec![zero_cell; ncell],
            boundary_points: HashMap::new(),
            boundary_edges: HashMap::new(),
            boundary_faces: HashMap::new(),
        })
    }

    pub fn from_serialized(serialized: &Vec<u8>) -> Result<Self, &'static str> {
        let mut deserializer = rmp_serde::Deserializer::new(&serialized[..]);
        let res = Deserialize::deserialize(&mut deserializer).map_err(|_| "cannot deserialize data")?;
        Ok(res)
    }

    pub fn get_serialized(&self) -> Result<Vec<u8>, &'static str> {
        let mut serialized = Vec::new();
        let mut serializer = rmp_serde::Serializer::new(&mut serialized);
        self.serialize(&mut serializer).map_err(|_| "serialize failed")?;
        Ok(serialized)
    }

    pub fn write_serialized<P>(&self, full_path: &P) -> Result<(), &'static str>
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

    pub fn get_limits(&self) -> (Vec<f64>, Vec<f64>) {
        let mut min = vec![f64::MAX; self.ndim];
        let mut max = vec![f64::MIN; self.ndim];
        for point in &self.points {
            for i in 0..self.ndim {
                if point.coords[i] < min[i] {
                    min[i] = point.coords[i];
                }
                if point.coords[i] > max[i] {
                    max[i] = point.coords[i];
                }
            }
        }
        return (min, max);
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
            write!(
                f,
                "i:{} g:{} x:{:?} c:{:?}\n",
                point.id, point.group, point.coords, point.shared_by_cell_ids
            )
            .unwrap();
        }

        // cells: i=index, g=group, p=points, e=boundary_edge_ids, f=boundary_face_ids
        write!(f, "\ncells\n").unwrap();
        for cell in &self.cells {
            write!(
                f,
                "i:{} g:{} p:{:?} e:{:?} f:{:?}\n",
                cell.id, cell.group, cell.point_ids, cell.boundary_edge_ids, cell.boundary_face_ids
            )
            .unwrap();
        }

        // boundary points: i=index
        write!(f, "\nboundary_points\n").unwrap();
        let mut point_indices: Vec<_> = self.boundary_points.keys().collect();
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
        keys_and_edges.sort_by(|left, right| left.1.id.cmp(&right.1.id));
        for (key, edge) in keys_and_edges {
            write!(
                f,
                "i:{} g:{} k:({},{}) p:{:?} c:{:?}\n",
                edge.id, edge.group, key.0, key.1, edge.point_ids, edge.shared_by_cell_ids
            )
            .unwrap();
        }

        // boundary_faces: i=index, g=group, k=key(point,point,point), p=point_ids, c=shared_by_cell_ids
        write!(f, "\nboundary_faces\n").unwrap();
        let mut keys_and_faces: Vec<_> = self.boundary_faces.keys().zip(self.boundary_faces.values()).collect();
        keys_and_faces.sort_by(|left, right| left.1.id.cmp(&right.1.id));
        for (key, face) in keys_and_faces {
            write!(
                f,
                "i:{} g:{} k:({},{},{}) p:{:?} c:{:?}\n",
                face.id, face.group, key.0, key.1, key.2, face.point_ids, face.shared_by_cell_ids
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
    fn from_string_works() -> Result<(), &'static str> {
        let data = "ndim = 2\n\
             npoint = 9\n\
             ncell = 4\n\
             n_boundary_point = 8\n\
             n_boundary_edge = 8\n\
             n_boundary_face = 0\n\
             \n\
             points\n\
             i:0 g:1 x:[0.0, 0.0] c:[0]\n\
             i:1 g:1 x:[1.0, 0.0] c:[0, 1]\n\
             i:2 g:1 x:[1.0, 1.0] c:[0, 1, 2, 3]\n\
             i:3 g:1 x:[0.0, 1.0] c:[0, 2]\n\
             i:4 g:1 x:[2.0, 0.0] c:[1]\n\
             i:5 g:1 x:[2.0, 1.0] c:[1, 3]\n\
             i:6 g:1 x:[1.0, 2.0] c:[2, 3]\n\
             i:7 g:1 x:[0.0, 2.0] c:[2]\n\
             i:8 g:1 x:[2.0, 2.0] c:[3]\n\
             \n\
             cells\n\
             i:0 g:1 p:[0, 1, 2, 3] e:[0, 1] f:[]\n\
             i:1 g:1 p:[1, 4, 5, 2] e:[2, 3] f:[]\n\
             i:2 g:1 p:[3, 2, 6, 7] e:[4, 5] f:[]\n\
             i:3 g:1 p:[2, 5, 8, 6] e:[6, 7] f:[]\n\
             \n\
             boundary_points\n\
             0 1 3 4 5 6 7 8 \n\
             \n\
             boundary_edges\n\
             i:0 g:1 k:(0,1) p:[0, 1] c:[0]\n\
             i:1 g:1 k:(0,3) p:[3, 0] c:[0]\n\
             i:2 g:1 k:(1,4) p:[1, 4] c:[1]\n\
             i:3 g:1 k:(4,5) p:[4, 5] c:[1]\n\
             i:4 g:1 k:(6,7) p:[6, 7] c:[2]\n\
             i:5 g:1 k:(3,7) p:[7, 3] c:[2]\n\
             i:6 g:1 k:(5,8) p:[5, 8] c:[3]\n\
             i:7 g:1 k:(6,8) p:[8, 6] c:[3]\n\
             \n\
             boundary_faces\n"
            .to_string();
        // Mesh::from_string(&data)?;
        Ok(())
    }

    #[test]
    fn display_works() {
        // todo
    }
}
