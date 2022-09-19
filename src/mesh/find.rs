use super::{EdgeKey, Extract, FaceKey, Feature, Features, Mesh, PointId};
use crate::util::GridSearch;
use crate::StrError;
use std::collections::{HashMap, HashSet};

/// Defines the location of points
#[derive(Clone, Copy, Debug)]
pub enum At {
    /// Fixed x, any (y,z). 2D (vertical-line) or 3D (yz-plane)
    X(f64),

    /// Fixed y, any (x,z). 2D (horizontal-line) or 3D (xz-plane)
    Y(f64),

    /// Fixed z, any (y,z). 3D only (xy-plane)
    Z(f64),

    /// Fixed (x,y), any z. 2D (point) or 3D (z-line)
    XY(f64, f64),

    /// Fixed (y,z), any x. 3D only (x-line)
    YZ(f64, f64),

    /// Fixed (x,z), any y. 3D only (y-line)
    XZ(f64, f64),

    /// At (x,y,z). 3D only (point)
    XYZ(f64, f64, f64),

    /// Circumference of circle at (x,y,radius). 2D only
    Circle(f64, f64, f64),

    /// Surface of cylinder given two axes. 3D only
    ///
    /// Holds: (axis_1_x, axis_1_y, axis_1_z, axis_2_x, axis_2_y, axis_2_z, radius).
    Cylinder(f64, f64, f64, f64, f64, f64, f64),
}

/// Implements functions to find mesh features (points, edges, faces)
pub struct Find {
    /// Holds features
    ///
    /// **readonly**
    pub features: Features,

    /// Space number of dimension (needed for the find functions)
    ///
    /// **readonly**
    pub space_ndim: usize,

    /// Total number of points in the mesh (needed to generate face keys of Tri3)
    ///
    /// **readonly**
    pub num_points: usize,

    /// Tool to quickly find points by coordinates
    ///
    /// **readonly**
    pub grid: GridSearch,

    /// Maps a point id to edges sharing the point
    ///
    /// **readonly**
    pub point_to_edges: HashMap<PointId, HashSet<EdgeKey>>,

    /// Maps a point id to faces sharing the point
    ///
    /// **readonly**
    pub point_to_faces: HashMap<PointId, HashSet<FaceKey>>,
}

impl Find {
    /// Allocates a new instance
    ///
    /// The default extracted features option is `Extract::Boundary`
    pub fn new(mesh: &Mesh, extract: Option<Extract>) -> Self {
        // extract features
        let features = match extract {
            Some(e) => Features::new(mesh, e),
            None => Features::new(mesh, Extract::Boundary),
        };

        // add point ids to grid
        let mut grid = GridSearch::new(&features.min, &features.max, None, None, None).unwrap();
        for point_id in &features.points {
            grid.insert(*point_id, &mesh.points[*point_id].coords).unwrap();
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
        Find {
            features,
            space_ndim: mesh.ndim,
            num_points: mesh.points.len(),
            grid,
            point_to_edges,
            point_to_faces,
        }
    }

    /// Finds points ids
    ///
    /// # Input
    ///
    /// * `at` -- the "main" constraint
    /// * `filter` -- fn(x) -> bool that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true)
    ///
    /// # Output
    ///
    /// * If at least one point has been found, returns a **sorted** array of point ids.
    /// * Otherwise, returns an error
    pub fn point_ids<F>(&self, at: At, filter: F) -> Result<Vec<PointId>, StrError>
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

    /// Finds edge keys
    ///
    ///
    /// # Input
    ///
    /// * `at` -- the "main" constraint
    /// * `filter` -- fn(x) -> bool that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true)
    ///
    /// # Output
    ///
    /// * If at least one point has been found, returns a **sorted** array of edge keys
    /// * Otherwise, returns an error
    pub fn edge_keys<F>(&self, at: At, filter: F) -> Result<Vec<EdgeKey>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        let mut edge_keys: HashSet<EdgeKey> = HashSet::new();
        // find all points constrained by "at" and "filter"
        let point_ids = self.point_ids(at, filter)?;
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

    /// Finds face keys
    ///
    /// # Input
    ///
    /// * `at` -- the "main" constraint
    /// * `filter` -- fn(x) -> bool that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true)
    ///
    /// # Output
    ///
    /// * If at least one point has been found, returns a **sorted** array of face keys
    /// * Otherwise, returns an error
    pub fn face_keys<F>(&self, at: At, filter: F) -> Result<Vec<FaceKey>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        if self.space_ndim != 3 {
            return Err("cannot find face keys in 2D");
        }
        let mut face_keys: HashSet<FaceKey> = HashSet::new();
        // find all points constrained by "at" and "filter"
        let point_ids = self.point_ids(at, filter)?;
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

    /// Finds edges
    ///
    /// # Input
    ///
    /// * `at` -- the "main" constraint
    /// * `filter` -- fn(x) -> bool that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true)
    ///
    /// # Output
    ///
    /// * If at least one point has been found, returns an array such that the edge keys are **sorted**
    /// * Otherwise, returns an error
    pub fn edges<F>(&self, at: At, filter: F) -> Result<Vec<&Feature>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        let results: Result<Vec<_>, _> = self
            .edge_keys(at, filter)?
            .iter()
            .map(|key| {
                self.features
                    .edges
                    .get(key)
                    .ok_or("INTERNAL ERROR: features.edges data is inconsistent")
            })
            .collect();
        results
    }

    /// Finds faces
    ///
    /// # Input
    ///
    /// * `at` -- the "main" constraint
    /// * `filter` -- fn(x) -> bool that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true)
    ///
    /// # Output
    ///
    /// * If at least one point has been found, returns an array such that the face keys are **sorted**
    /// * Otherwise, returns an error
    pub fn faces<F>(&self, at: At, filter: F) -> Result<Vec<&Feature>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        let results: Result<Vec<_>, _> = self
            .face_keys(at, filter)?
            .iter()
            .map(|key| {
                self.features
                    .faces
                    .get(key)
                    .ok_or("INTERNAL ERROR: features.faces data is inconsistent")
            })
            .collect();
        results
    }

    /// Finds (many) edges using a list of constraints
    ///
    /// # Input
    ///
    /// * `at` -- the "main" constraint
    /// * `filter` -- fn(x) -> bool that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true)
    ///
    /// # Output
    ///
    /// * Returns edges sorted by keys
    /// * **Warning** Every `At` in the `ats` must generate at least one edge,
    ///   otherwise an error will occur.
    pub fn many_edges<F>(&self, ats: &[At], mut filter: F) -> Result<Vec<&Feature>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        let mut edge_keys: HashSet<EdgeKey> = HashSet::new();
        for at in ats {
            let found = self.edge_keys(*at, &mut filter)?;
            for edge in found {
                edge_keys.insert(edge);
            }
        }
        let mut keys: Vec<_> = edge_keys.iter().collect();
        keys.sort();
        let results: Result<Vec<_>, _> = keys
            .iter()
            .map(|key| {
                self.features
                    .edges
                    .get(key)
                    .ok_or("INTERNAL ERROR: features.edges data is inconsistent")
            })
            .collect();
        results
    }

    /// Finds (many) faces using a list of constraints
    ///
    /// # Input
    ///
    /// * `at` -- the "main" constraint
    /// * `filter` -- fn(x) -> bool that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true)
    ///
    /// # Output
    ///
    /// * Returns faces sorted by keys
    /// * **Warning** Every `At` in the `ats` must generate at least one face,
    ///   otherwise an error will occur.
    pub fn many_faces<F>(&self, ats: &[At], mut filter: F) -> Result<Vec<&Feature>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        let mut face_keys: HashSet<FaceKey> = HashSet::new();
        for at in ats {
            let found = self.face_keys(*at, &mut filter)?;
            for face in found {
                face_keys.insert(face);
            }
        }
        let mut keys: Vec<_> = face_keys.iter().collect();
        keys.sort();
        let results: Result<Vec<_>, _> = keys
            .iter()
            .map(|key| {
                self.features
                    .faces
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
    use super::Find;
    use crate::mesh::{At, Extract, Samples};
    use crate::util::any;
    use russell_lab::math::{PI, SQRT_2};

    #[allow(unused_imports)]
    use plotpy::Plot;

    // DO NOT DELETE the following lines
    // fn plot_grid_two_qua4(find: &Find) {
    //     let mut plot = Plot::new();
    //     find.grid.draw(&mut plot).unwrap();
    //     plot.set_equal_axes(true).set_figure_size_points(800.0, 400.0);
    //     plot.save("/tmp/gemlab/test_find_with_two_qua4.svg").unwrap();
    // }

    // DO NOT DELETE the following lines
    // fn plot_grid_two_hex8(find: &Find) {
    //     let mut plot = Plot::new();
    //     find.grid.draw(&mut plot).unwrap();
    //     plot.set_equal_axes(true).set_figure_size_points(2048.0, 2048.0);
    //     plot.save("/tmp/gemlab/test_find_with_two_hex8.svg").unwrap();
    // }

    #[test]
    fn derive_methods_work() {
        let at = At::X(0.0);
        let at_clone = at.clone();
        assert_eq!(format!("{:?}", at_clone), "X(0.0)");
    }

    #[test]
    fn new_works() {
        let mesh = Samples::two_qua4();
        let find = Find::new(&mesh, None);
        assert_eq!(
            format!("{}", find.grid),
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
        // plot_grid_two_qua4(&find);
    }

    #[test]
    fn find_points_fails_on_wrong_input() {
        assert_eq!(any(&vec![]), true);
        // 2d
        let mesh = Samples::two_qua4();
        let find = Find::new(&mesh, None);
        assert_eq!(find.point_ids(At::Z(0.0), any).err(), Some("At::Z works in 3D only"));
        assert_eq!(
            find.point_ids(At::YZ(0.0, 0.0), any).err(),
            Some("At::YZ works in 3D only")
        );
        assert_eq!(
            find.point_ids(At::XZ(0.0, 0.0), any).err(),
            Some("At::XZ works in 3D only")
        );
        assert_eq!(
            find.point_ids(At::XYZ(0.0, 0.0, 0.0), any).err(),
            Some("At::XYZ works in 3D only")
        );
        assert_eq!(
            find.point_ids(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), any)
                .err(),
            Some("At::Cylinder works in 3D only")
        );

        // 3d
        let mesh = Samples::two_hex8();
        let find = Find::new(&mesh, None);
        assert_eq!(
            find.point_ids(At::Circle(0.0, 0.0, 0.0), any).err(),
            Some("At::Circle works in 2D only")
        );
    }

    #[test]
    fn find_points_works_2d() {
        // `.       `.
        //   3--------2--------5
        //   | `.     | `.     |
        //   |   `~.  | circle |
        //   |      `.|        |
        //   0--------1--------4
        //           circle
        let mesh = Samples::two_qua4();
        let find = Find::new(&mesh, Some(Extract::All));
        assert_eq!(find.point_ids(At::XY(0.0, 0.0), any).unwrap(), &[0]);
        assert_eq!(find.point_ids(At::XY(2.0, 1.0), any).unwrap(), &[5]);
        assert_eq!(
            find.point_ids(At::XY(10.0, 0.0), any).err(),
            Some("cannot find point because the coordinates are outside the grid")
        );
        assert_eq!(find.point_ids(At::Circle(0.0, 0.0, 1.0), any).unwrap(), &[1, 3]);
        assert_eq!(find.point_ids(At::Circle(0.0, 0.0, SQRT_2), any).unwrap(), &[2]);
        assert_eq!(
            find.point_ids(At::Circle(0.0, 0.0, 10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
    }

    #[test]
    fn find_points_works_3d() {
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
        let find = Find::new(&mesh, None);
        // plot_grid_two_hex8(&find);
        assert_eq!(find.point_ids(At::X(0.0), any).unwrap(), &[0, 3, 4, 7, 8, 11]);
        assert_eq!(find.point_ids(At::X(1.0), any).unwrap(), &[1, 2, 5, 6, 9, 10]);
        assert_eq!(
            find.point_ids(At::X(10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(find.point_ids(At::Y(0.0), any).unwrap(), &[0, 1, 4, 5, 8, 9]);
        assert_eq!(find.point_ids(At::Y(1.0), any).unwrap(), &[2, 3, 6, 7, 10, 11]);
        assert_eq!(
            find.point_ids(At::Y(10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(find.point_ids(At::Z(0.0), any).unwrap(), &[0, 1, 2, 3]);
        assert_eq!(find.point_ids(At::Z(1.0), any).unwrap(), &[4, 5, 6, 7]);
        assert_eq!(find.point_ids(At::Z(2.0), any).unwrap(), &[8, 9, 10, 11]);
        assert_eq!(
            find.point_ids(At::Z(10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(find.point_ids(At::XY(0.0, 0.0), any).unwrap(), &[0, 4, 8]);
        assert_eq!(find.point_ids(At::XY(1.0, 1.0), any).unwrap(), &[2, 6, 10]);
        assert_eq!(
            find.point_ids(At::XY(10.0, 10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(find.point_ids(At::YZ(0.0, 0.0), any).unwrap(), &[0, 1]);
        assert_eq!(find.point_ids(At::YZ(1.0, 1.0), any).unwrap(), &[6, 7]);
        assert_eq!(find.point_ids(At::XZ(0.0, 0.0), any).unwrap(), &[0, 3]);
        assert_eq!(find.point_ids(At::XZ(1.0, 0.0), any).unwrap(), &[1, 2]);
        assert_eq!(find.point_ids(At::XZ(1.0, 2.0), any).unwrap(), &[9, 10]);
        assert_eq!(find.point_ids(At::XYZ(0.0, 0.0, 0.0), any).unwrap(), &[0]);
        assert_eq!(find.point_ids(At::XYZ(1.0, 1.0, 2.0), any).unwrap(), &[10]);
        assert_eq!(
            find.point_ids(At::XYZ(10.0, 0.0, 0.0), any).err(),
            Some("cannot find point because the coordinates are outside the grid")
        );
        assert_eq!(
            find.point_ids(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0), any)
                .unwrap(),
            &[1, 3, 5, 7, 9, 11],
        );
        assert_eq!(
            find.point_ids(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, SQRT_2), any)
                .unwrap(),
            &[2, 6, 10],
        );
        assert_eq!(
            find.point_ids(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 10.0), any)
                .err(),
            Some("cannot find any point with given constraints/filter")
        );
    }

    #[test]
    fn find_edges_works_2d() {
        // 3--------2--------5
        // |        |        |
        // |        |        |
        // |        |        |
        // 0--------1--------4
        let mesh = Samples::two_qua4();
        let find = Find::new(&mesh, None);
        assert_eq!(find.edge_keys(At::Y(0.0), any).unwrap(), &[(0, 1), (1, 4)]);
        assert_eq!(find.edge_keys(At::Y(0.0), |x| x[0] <= 1.0).unwrap(), &[(0, 1)]);
        assert_eq!(find.edge_keys(At::X(2.0), any).unwrap(), &[(4, 5)]);
        assert_eq!(find.edge_keys(At::Y(1.0), any).unwrap(), &[(2, 3), (2, 5)]);
        assert_eq!(find.edge_keys(At::Y(1.0), |x| x[0] >= 1.0).unwrap(), &[(2, 5)]);
        assert_eq!(find.edge_keys(At::X(0.0), any).unwrap(), &[(0, 3)]);

        // internal
        assert_eq!(
            find.edge_keys(At::X(1.0), any).err(),
            Some("cannot find any edge with given constraints/filter")
        );

        // far away
        assert_eq!(
            find.edge_keys(At::X(10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );

        // high-level function
        let res = find.edges(At::Y(0.0), any).unwrap();
        assert_eq!(res.len(), 2);
        assert_eq!(res[0].points, &[1, 0]);
        assert_eq!(res[1].points, &[4, 1]);
        assert_eq!(
            find.edges(At::XYZ(0.0, 0.0, 0.0), any).err(),
            Some("At::XYZ works in 3D only")
        );

        // many edges
        let res = find
            .many_edges(&[At::X(0.0), At::X(2.0), At::Y(0.0), At::Y(1.0)], any)
            .unwrap();
        assert_eq!(res.len(), 6);
        let keys: Vec<_> = res.iter().map(|r| (r.points[0], r.points[1])).collect();
        assert_eq!(keys, &[(1, 0), (0, 3), (4, 1), (3, 2), (2, 5), (5, 4)]);
    }

    #[test]
    fn find_edges_works_3d() {
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
        let find = Find::new(&mesh, None);
        assert_eq!(
            find.edge_keys(At::X(0.0), any).unwrap(),
            &[(0, 3), (0, 4), (3, 7), (4, 7), (4, 8), (7, 11), (8, 11)],
        );
        assert_eq!(find.edge_keys(At::X(0.0), |x| x[2] == 0.0).unwrap(), &[(0, 3)],);
        assert_eq!(
            find.edge_keys(At::X(1.0), any).unwrap(),
            &[(1, 2), (1, 5), (2, 6), (5, 6), (5, 9), (6, 10), (9, 10)],
        );
        assert_eq!(
            find.edge_keys(At::X(10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            find.edge_keys(At::Y(0.0), any).unwrap(),
            &[(0, 1), (0, 4), (1, 5), (4, 5), (4, 8), (5, 9), (8, 9)],
        );
        assert_eq!(
            find.edge_keys(At::Y(0.0), |x| x[2] <= 1.0).unwrap(),
            &[(0, 1), (0, 4), (1, 5), (4, 5)],
        );
        assert_eq!(
            find.edge_keys(At::Y(1.0), any).unwrap(),
            &[(2, 3), (2, 6), (3, 7), (6, 7), (6, 10), (7, 11), (10, 11)],
        );
        assert_eq!(
            find.edge_keys(At::Y(10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            find.edge_keys(At::Z(0.0), any).unwrap(),
            &[(0, 1), (0, 3), (1, 2), (2, 3)]
        );
        assert_eq!(
            find.edge_keys(At::Z(2.0), any).unwrap(),
            &[(8, 9), (8, 11), (9, 10), (10, 11)],
        );
        assert_eq!(find.edge_keys(At::Z(2.0), |x| x[1] == 1.0).unwrap(), &[(10, 11)],);
        assert_eq!(
            find.edge_keys(At::Z(10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(find.edge_keys(At::XY(0.0, 0.0), any).unwrap(), &[(0, 4), (4, 8)]);
        assert_eq!(find.edge_keys(At::XY(1.0, 1.0), any).unwrap(), &[(2, 6), (6, 10)]);
        assert_eq!(find.edge_keys(At::XY(1.0, 1.0), |x| x[2] <= 1.0).unwrap(), &[(2, 6)]);
        assert_eq!(
            find.edge_keys(At::XY(10.0, 10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(find.edge_keys(At::YZ(0.0, 0.0), any).unwrap(), &[(0, 1)]);
        assert_eq!(find.edge_keys(At::YZ(1.0, 1.0), any).unwrap(), &[(6, 7)]);
        assert_eq!(
            find.edge_keys(At::YZ(10.0, 10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(find.edge_keys(At::XZ(0.0, 0.0), any).unwrap(), &[(0, 3)]);
        assert_eq!(find.edge_keys(At::XZ(1.0, 0.0), any).unwrap(), &[(1, 2)]);
        assert_eq!(find.edge_keys(At::XZ(1.0, 2.0), any).unwrap(), &[(9, 10)]);
        assert_eq!(find.edge_keys(At::XZ(1.0, 2.0), |x| x[0] == -1.0).ok(), None);
        assert_eq!(
            find.edge_keys(At::XZ(10.0, 10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            find.edge_keys(At::XYZ(0.0, 0.0, 0.0), any).err(),
            Some("cannot find any edge with given constraints/filter")
        );
        assert_eq!(
            find.edge_keys(At::XYZ(10.0, 0.0, 0.0), any).err(),
            Some("cannot find point because the coordinates are outside the grid")
        );
        assert_eq!(
            find.edge_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0), any)
                .unwrap(),
            &[(1, 5), (3, 7), (5, 9), (7, 11)],
        );
        assert_eq!(
            find.edge_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0), |x| x[2] >= 1.0)
                .unwrap(),
            &[(5, 9), (7, 11)],
        );
        assert_eq!(
            find.edge_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, SQRT_2), any)
                .unwrap(),
            &[(2, 6), (6, 10)],
        );
        assert_eq!(
            find.edge_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 10.0), any)
                .err(),
            Some("cannot find any point with given constraints/filter")
        );

        // high-level function
        let res = find.edges(At::XY(0.0, 0.0), any).unwrap();
        assert_eq!(res.len(), 2);
        assert_eq!(res[0].points, &[0, 4]);
        assert_eq!(res[1].points, &[4, 8]);

        // many faces
        let res = find.many_faces(&[At::Z(0.0), At::Z(2.0)], any).unwrap();
        assert_eq!(res.len(), 2);
        let keys: Vec<_> = res
            .iter()
            .map(|r| (r.points[0], r.points[1], r.points[2], r.points[3]))
            .collect();
        assert_eq!(keys, &[(0, 3, 2, 1), (8, 9, 10, 11)]);
    }

    #[test]
    fn find_faces_returns_error_in_2d() {
        // 3--------2--------5
        // |        |        |
        // |        |        |
        // |        |        |
        // 0--------1--------4
        let mesh = Samples::two_qua4();
        let find = Find::new(&mesh, None);
        assert_eq!(
            find.face_keys(At::X(0.0), any).err(),
            Some("cannot find face keys in 2D")
        );
    }

    #[test]
    fn find_faces_works() {
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
        let find = Find::new(&mesh, None);
        assert_eq!(find.face_keys(At::X(0.0), any).unwrap(), &[(0, 3, 4, 7), (4, 7, 8, 11)]);
        assert_eq!(find.face_keys(At::X(1.0), any).unwrap(), &[(1, 2, 5, 6), (5, 6, 9, 10)]);
        assert_eq!(find.face_keys(At::X(1.0), |x| x[2] <= 1.0).unwrap(), &[(1, 2, 5, 6)]);
        assert_eq!(
            find.face_keys(At::X(10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(find.face_keys(At::Y(0.0), any).unwrap(), &[(0, 1, 4, 5), (4, 5, 8, 9)]);
        assert_eq!(
            find.face_keys(At::Y(1.0), any).unwrap(),
            &[(2, 3, 6, 7), (6, 7, 10, 11)]
        );
        assert_eq!(find.face_keys(At::Y(1.0), |x| x[2] >= 1.0).unwrap(), &[(6, 7, 10, 11)]);
        assert_eq!(
            find.face_keys(At::Y(10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(find.face_keys(At::Z(0.0), any).unwrap(), &[(0, 1, 2, 3)]);
        assert_eq!(find.face_keys(At::Z(2.0), any).unwrap(), &[(8, 9, 10, 11)]);
        assert_eq!(find.face_keys(At::Z(2.0), |x| x[0] <= -1.0).ok(), None);
        assert_eq!(
            find.face_keys(At::Z(10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::XY(0.0, 0.0), any).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::XY(1.0, 1.0), any).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::XY(10.0, 10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::YZ(0.0, 0.0), any).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::YZ(1.0, 1.0), any).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::YZ(10.0, 10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::XZ(0.0, 0.0), any).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::XZ(1.0, 0.0), any).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::XZ(1.0, 2.0), any).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::XZ(10.0, 10.0), any).err(),
            Some("cannot find any point with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::XYZ(0.0, 0.0, 0.0), any).err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::XYZ(10.0, 0.0, 0.0), any).err(),
            Some("cannot find point because the coordinates are outside the grid")
        );
        assert_eq!(
            find.face_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0), any)
                .err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, SQRT_2), any)
                .err(),
            Some("cannot find any face with given constraints/filter")
        );
        assert_eq!(
            find.face_keys(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 10.0), any)
                .err(),
            Some("cannot find any point with given constraints/filter")
        );

        // high-level function
        let res = find.faces(At::Z(0.0), any).unwrap();
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].points, &[0, 3, 2, 1]);
        assert_eq!(
            find.faces(At::Circle(0.0, 0.0, 1.0), any).err(),
            Some("At::Circle works in 2D only")
        );
    }

    #[test]
    fn find_works_with_ring() {
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
        let find = Find::new(&mesh, None);
        let (r, rr) = (1.0, 2.0);
        assert_eq!(find.point_ids(At::XY(1.00, 0.00), any).unwrap(), &[0]);
        assert_eq!(find.point_ids(At::XY(1.25, 0.00), any).unwrap(), &[15]);
        assert_eq!(find.point_ids(At::XY(1.50, 0.00), any).unwrap(), &[1]);
        assert_eq!(find.point_ids(At::XY(1.75, 0.00), any).unwrap(), &[16]);
        assert_eq!(find.point_ids(At::XY(2.00, 0.00), any).unwrap(), &[2]);
        assert_eq!(find.point_ids(At::XY(0.00, 1.00), any).unwrap(), &[12]);
        assert_eq!(find.point_ids(At::XY(0.00, 1.25), any).unwrap(), &[23]);
        assert_eq!(find.point_ids(At::XY(0.00, 1.75), any).unwrap(), &[24]);
        assert_eq!(find.point_ids(At::XY(0.00, 1.50), any).unwrap(), &[13]);
        assert_eq!(find.point_ids(At::XY(0.00, 2.00), any).unwrap(), &[14]);
        assert_eq!(find.point_ids(At::XY(SQRT_2 / 2.0, SQRT_2 / 2.0), any).unwrap(), &[6]);
        assert_eq!(find.point_ids(At::XY(SQRT_2, SQRT_2), any).unwrap(), &[8]);
        assert_eq!(
            find.point_ids(At::Circle(0.0, 0.0, r), any).unwrap(),
            &[0, 3, 6, 9, 12, 25, 28, 31, 34],
        );
        assert_eq!(
            find.point_ids(At::Circle(0.0, 0.0, r), |x| {
                let alpha = f64::atan2(x[1], x[0]) * 180.0 / PI;
                f64::abs(alpha - 45.0) < 1e-15
            })
            .unwrap(),
            &[6],
        );
        assert_eq!(
            find.point_ids(At::Circle(0.0, 0.0, rr), any).unwrap(),
            &[2, 5, 8, 11, 14, 27, 30, 33, 36],
        );
        assert_eq!(find.edge_keys(At::Y(0.0), any).unwrap(), &[(0, 1), (1, 2)]);
        assert_eq!(find.edge_keys(At::X(0.0), any).unwrap(), &[(12, 13), (13, 14)]);
        assert_eq!(
            find.edge_keys(At::Circle(0.0, 0.0, r), any).unwrap(),
            &[(0, 3), (3, 6), (6, 9), (9, 12)],
        );
        assert_eq!(
            find.edge_keys(At::Circle(0.0, 0.0, rr), any).unwrap(),
            &[(2, 5), (5, 8), (8, 11), (11, 14)],
        );
    }
}
