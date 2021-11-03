use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Point {
    pub id: usize,
    pub group: usize,
    pub coords: Vec<f64>,
    pub shared_by_cell_ids: Vec<usize>,
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
    pub fn new(ndim: usize) -> Self {
        Mesh {
            ndim,
            points: Vec::new(),
            cells: Vec::new(),
            boundary_points: HashMap::new(),
            boundary_edges: HashMap::new(),
            boundary_faces: HashMap::new(),
        }
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

    pub fn is_point_on_boundary(&self, id: usize) -> bool {
        match self.boundary_points.get(&id) {
            Some(_) => true,
            None => false,
        }
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
    // use super::*;

    #[test]
    fn new_works() {
        // todo
    }

    #[test]
    fn serialize_works() {
        // todo
    }

    #[test]
    fn display_works() {
        // todo
    }
}
