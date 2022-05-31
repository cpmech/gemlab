use super::{At, EdgeKey, FaceKey, Features, Mesh, PointId};
use crate::util::{GridSearch, GsNdiv, GsTol};
use crate::StrError;
use std::collections::{HashMap, HashSet};

/// Implements functions to find mesh features (points, edges, faces)
///
/// # Warning
///
/// All public members are **read-only**.
pub struct Find {
    /// Space number of dimension (needed for the find functions)
    space_ndim: usize,

    /// Total number of points in the mesh (needed to generate face keys of Tri3)
    num_points: usize,

    /// Tool to quickly find points by coordinates
    grid: GridSearch,

    /// Maps a point id to edges sharing the point
    pub point_to_edges: HashMap<PointId, HashSet<EdgeKey>>,

    /// Maps a point id to faces sharing the point
    pub point_to_faces: HashMap<PointId, HashSet<FaceKey>>,
}

impl Find {
    /// Allocates a new instance
    ///
    /// **Note:** This function is to be used by [Region] and not by the end-user.
    ///
    /// # Panics
    ///
    /// If the `features` does not correspond to `mesh`, a panic will occur.
    pub(crate) fn new(mesh: &Mesh, features: &Features) -> Result<Self, StrError> {
        // expand limits a little bit to accommodate imprecisions on the (min,max) values
        const PCT: f64 = 1.0 / 100.0; // 1% is enough
        let (mut min, mut max) = (vec![0.0; mesh.space_ndim], vec![0.0; mesh.space_ndim]);
        for i in 0..mesh.space_ndim {
            let del = features.max[i] - features.min[i];
            min[i] = features.min[i] - PCT * del;
            max[i] = features.max[i] + PCT * del;
        }

        // add point ids to grid
        let mut grid = GridSearch::new(&min, &max, GsNdiv::Default, GsTol::Default)?;
        for point_id in &features.points {
            grid.insert(*point_id, &mesh.points[*point_id].coords)?;
        }

        // map point ids to edges
        let mut point_to_edges: HashMap<PointId, HashSet<EdgeKey>> = HashMap::new();
        for (edge_key, edge) in &features.edges {
            for point_id in &edge.points {
                point_to_edges
                    .entry(*point_id)
                    .or_insert(HashSet::new())
                    .insert(*edge_key);
            }
        }

        // map point ids to faces
        let mut point_to_faces: HashMap<PointId, HashSet<FaceKey>> = HashMap::new();
        for (face_key, face) in &features.faces {
            for point_id in &face.points {
                point_to_faces
                    .entry(*point_id)
                    .or_insert(HashSet::new())
                    .insert(*face_key);
            }
        }

        // done
        Ok(Find {
            space_ndim: mesh.space_ndim,
            num_points: mesh.points.len(),
            grid,
            point_to_edges,
            point_to_faces,
        })
    }

    /// Finds points
    ///
    /// # Input
    ///
    /// * `at` -- the location constraint
    ///
    /// # Output
    ///
    /// * Returns a set of point ids.
    ///   You may sort point ids using the following code snippet:
    ///
    /// ``` text
    /// let mut ids: Vec<_> = point_ids.iter().copied().collect();
    /// ids.sort();
    /// ```
    pub fn points(&self, at: At) -> Result<HashSet<PointId>, StrError> {
        let mut point_ids: HashSet<PointId> = HashSet::new();
        match at {
            At::X(x) => {
                if self.space_ndim == 2 {
                    for id in self.grid.find_on_line(&[x, 0.0], &[x, 1.0])? {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid.find_on_plane_yz(x)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::Y(y) => {
                if self.space_ndim == 2 {
                    for id in self.grid.find_on_line(&[0.0, y], &[1.0, y])? {
                        point_ids.insert(id);
                    }
                } else {
                    for id in self.grid.find_on_plane_xz(y)? {
                        point_ids.insert(id);
                    }
                }
            }
            At::Z(z) => {
                if self.space_ndim == 2 {
                    return Err("At::Z works in 3D only");
                } else {
                    for id in self.grid.find_on_plane_xy(z)? {
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
                    for id in self.grid.find_on_line(&[x, y, 0.0], &[x, y, 1.0])? {
                        point_ids.insert(id);
                    }
                }
            }
            At::YZ(y, z) => {
                if self.space_ndim == 2 {
                    return Err("At::YZ works in 3D only");
                } else {
                    for id in self.grid.find_on_line(&[0.0, y, z], &[1.0, y, z])? {
                        point_ids.insert(id);
                    }
                }
            }
            At::XZ(x, z) => {
                if self.space_ndim == 2 {
                    return Err("At::XZ works in 3D only");
                } else {
                    for id in self.grid.find_on_line(&[x, 0.0, z], &[x, 1.0, z])? {
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
                    for id in self.grid.find_on_circle(&[x, y], r)? {
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
                    for id in self.grid.find_on_cylinder(&[ax, ay, az], &[bx, by, bz], r)? {
                        point_ids.insert(id);
                    }
                }
            }
        }
        Ok(point_ids)
    }

    /// Finds edges
    ///
    /// # Input
    ///
    /// * `at` -- the location constraint
    ///
    /// # Output
    ///
    /// * Returns a set of edge keys.
    ///   You may sort the edge keys using the following code snippet:
    ///
    /// ``` text
    /// let mut keys: Vec<_> = edge_keys.iter().copied().collect();
    /// keys.sort();
    /// ```
    pub fn edges(&self, at: At) -> Result<HashSet<EdgeKey>, StrError> {
        let mut edge_keys: HashSet<EdgeKey> = HashSet::new();
        // find all points constrained by "at"
        let point_ids = self.points(at)?;
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
        Ok(edge_keys)
    }

    /// Finds faces
    ///
    /// # Input
    ///
    /// * `at` -- the location constraint
    ///
    /// # Output
    ///
    /// * Returns a set of face keys.
    ///   You may sort the face keys using the following code snippet:
    ///
    /// ``` text
    /// let mut keys: Vec<_> = face_keys.iter().copied().collect();
    /// keys.sort();
    /// ```
    pub fn faces(&self, at: At) -> Result<HashSet<FaceKey>, StrError> {
        let mut face_keys: HashSet<FaceKey> = HashSet::new();
        if self.space_ndim != 3 {
            return Ok(face_keys);
        }
        // find all points constrained by "at"
        let point_ids = self.points(at)?;
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
        Ok(face_keys)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {}
