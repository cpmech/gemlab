use super::{EdgeKey, FaceKey, Features, Mesh, PointId};
use crate::util::{GridSearch, GsNdiv, GsTol};
use crate::StrError;
use std::collections::{HashMap, HashSet};

/// Defines the location of points
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
    /// Space number of dimension (needed for the find functions)
    space_ndim: usize,

    /// Total number of points in the mesh (needed to generate face keys of Tri3)
    num_points: usize,

    /// Tool to quickly find points by coordinates
    grid: GridSearch,

    /// Maps a point id to edges sharing the point
    point_to_edges: HashMap<PointId, HashSet<EdgeKey>>,

    /// Maps a point id to faces sharing the point
    point_to_faces: HashMap<PointId, HashSet<FaceKey>>,
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
        let (mut min, mut max) = (vec![0.0; mesh.ndim], vec![0.0; mesh.ndim]);
        for i in 0..mesh.ndim {
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
            space_ndim: mesh.ndim,
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
mod tests {
    use super::Find;
    use crate::mesh::algorithms::{extract_all_2d_edges, extract_all_faces, extract_features_2d, extract_features_3d};
    use crate::mesh::{allocate_cell_shapes, At, Extract, Samples};
    use crate::util::SQRT_2;
    use crate::StrError;
    use std::collections::HashSet;

    #[allow(dead_code)]
    fn plot_grid_two_quads_horizontal(find: &Find) -> Result<(), StrError> {
        let mut plot = find.grid.plot()?;
        plot.set_equal_axes(true).set_figure_size_points(800.0, 400.0);
        plot.save("/tmp/gemlab/find_with_two_quads_horizontal.svg")
    }

    #[allow(dead_code)]
    fn plot_grid_two_cubes_vertical(find: &Find) -> Result<(), StrError> {
        let mut plot = find.grid.plot()?;
        plot.set_equal_axes(true).set_figure_size_points(2048.0, 2048.0);
        plot.save("/tmp/gemlab/find_with_two_cubes_vertical.svg")
    }

    #[test]
    fn new_works() -> Result<(), StrError> {
        let mesh = Samples::two_quads_horizontal();
        let shapes = allocate_cell_shapes(&mesh);
        let edges = extract_all_2d_edges(&mesh, &shapes);
        let boundary = extract_features_2d(&mesh, &shapes, &edges, Extract::Boundary);
        let find = Find::new(&mesh, &boundary)?;
        // plot_grid_two_quads_horizontal(&find)?;
        assert_eq!(
            format!("{}", find.grid),
            "0: [0]\n\
             9: [1]\n\
             10: [1]\n\
             19: [4]\n\
             180: [3]\n\
             189: [2]\n\
             190: [2]\n\
             199: [5]\n\
             ids = [0, 1, 2, 3, 4, 5]\n\
             nitem = 6\n\
             ncontainer = 8\n\
             ndiv = [20, 10]\n"
        );
        Ok(())
    }

    #[test]
    fn find_points_fails_on_wrong_input() -> Result<(), StrError> {
        // 2d
        let mesh = Samples::two_quads_horizontal();
        let shapes = allocate_cell_shapes(&mesh);
        let edges = extract_all_2d_edges(&mesh, &shapes);
        let boundary = extract_features_2d(&mesh, &shapes, &edges, Extract::Boundary);
        let find = Find::new(&mesh, &boundary)?;
        assert_eq!(find.points(At::Z(0.0)).err(), Some("At::Z works in 3D only"));
        assert_eq!(find.points(At::YZ(0.0, 0.0)).err(), Some("At::YZ works in 3D only"));
        assert_eq!(find.points(At::XZ(0.0, 0.0)).err(), Some("At::XZ works in 3D only"));
        assert_eq!(
            find.points(At::XYZ(0.0, 0.0, 0.0)).err(),
            Some("At::XYZ works in 3D only")
        );
        assert_eq!(
            find.points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)).err(),
            Some("At::Cylinder works in 3D only")
        );

        // 3d
        let mesh = Samples::two_cubes_vertical();
        let shapes = allocate_cell_shapes(&mesh);
        let faces = extract_all_faces(&mesh, &shapes);
        let boundary = extract_features_3d(&mesh, &shapes, &faces, Extract::Boundary);
        let find = Find::new(&mesh, &boundary)?;
        assert_eq!(
            find.points(At::Circle(0.0, 0.0, 0.0)).err(),
            Some("At::Circle works in 2D only")
        );
        Ok(())
    }

    fn check<T>(found: &HashSet<T>, correct: &[T])
    where
        T: Copy + Ord + std::fmt::Debug,
    {
        let mut ids: Vec<T> = found.iter().copied().collect();
        ids.sort();
        assert_eq!(ids, correct);
    }

    #[test]
    fn find_points_works_2d() -> Result<(), StrError> {
        // `.       `.
        //   3--------2--------5
        //   | `.     | `.     |
        //   |   `~.  | circle |
        //   |      `.|        |
        //   0--------1--------4
        //           circle
        let mesh = Samples::two_quads_horizontal();
        let shapes = allocate_cell_shapes(&mesh);
        let edges = extract_all_2d_edges(&mesh, &shapes);
        let boundary = extract_features_2d(&mesh, &shapes, &edges, Extract::Boundary);
        let find = Find::new(&mesh, &boundary)?;
        check(&find.points(At::XY(0.0, 0.0))?, &[0]);
        check(&find.points(At::XY(2.0, 1.0))?, &[5]);
        assert_eq!(
            find.points(At::XY(10.0, 0.0)).err(),
            Some("cannot find point with coordinates outside the grid")
        );
        check(&find.points(At::Circle(0.0, 0.0, 1.0))?, &[1, 3]);
        check(&find.points(At::Circle(0.0, 0.0, SQRT_2))?, &[2]);
        check(&find.points(At::Circle(0.0, 0.0, 10.0))?, &[]);
        Ok(())
    }

    #[test]
    fn find_points_works_3d() -> Result<(), StrError> {
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
        let mesh = Samples::two_cubes_vertical();
        let shapes = allocate_cell_shapes(&mesh);
        let faces = extract_all_faces(&mesh, &shapes);
        let boundary = extract_features_3d(&mesh, &shapes, &faces, Extract::Boundary);
        let find = Find::new(&mesh, &boundary)?;
        // plot_grid_two_cubes_vertical(&find)?;
        check(&find.points(At::X(0.0))?, &[0, 3, 4, 7, 8, 11]);
        check(&find.points(At::X(1.0))?, &[1, 2, 5, 6, 9, 10]);
        check(&find.points(At::X(10.0))?, &[]);
        check(&find.points(At::Y(0.0))?, &[0, 1, 4, 5, 8, 9]);
        check(&find.points(At::Y(1.0))?, &[2, 3, 6, 7, 10, 11]);
        check(&find.points(At::Y(10.0))?, &[]);
        check(&find.points(At::Z(0.0))?, &[0, 1, 2, 3]);
        check(&find.points(At::Z(1.0))?, &[4, 5, 6, 7]);
        check(&find.points(At::Z(2.0))?, &[8, 9, 10, 11]);
        check(&find.points(At::Z(10.0))?, &[]);
        check(&find.points(At::XY(0.0, 0.0))?, &[0, 4, 8]);
        check(&find.points(At::XY(1.0, 1.0))?, &[2, 6, 10]);
        check(&find.points(At::XY(10.0, 10.0))?, &[]);
        check(&find.points(At::YZ(0.0, 0.0))?, &[0, 1]);
        check(&find.points(At::YZ(1.0, 1.0))?, &[6, 7]);
        check(&find.points(At::XZ(0.0, 0.0))?, &[0, 3]);
        check(&find.points(At::XZ(1.0, 0.0))?, &[1, 2]);
        check(&find.points(At::XZ(1.0, 2.0))?, &[9, 10]);
        check(&find.points(At::XYZ(0.0, 0.0, 0.0))?, &[0]);
        check(&find.points(At::XYZ(1.0, 1.0, 2.0))?, &[10]);
        assert_eq!(
            find.points(At::XYZ(10.0, 0.0, 0.0)).err(),
            Some("cannot find point with coordinates outside the grid")
        );
        check(
            &find.points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0))?,
            &[1, 3, 5, 7, 9, 11],
        );
        check(
            &find.points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, SQRT_2))?,
            &[2, 6, 10],
        );
        check(&find.points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 10.0))?, &[]);
        Ok(())
    }

    #[test]
    fn find_edges_works_2d() -> Result<(), StrError> {
        // 3--------2--------5
        // |        |        |
        // |        |        |
        // |        |        |
        // 0--------1--------4
        let mesh = Samples::two_quads_horizontal();
        let shapes = allocate_cell_shapes(&mesh);
        let edges = extract_all_2d_edges(&mesh, &shapes);
        let boundary = extract_features_2d(&mesh, &shapes, &edges, Extract::Boundary);
        let find = Find::new(&mesh, &boundary)?;
        check(&find.edges(At::Y(0.0))?, &[(0, 1), (1, 4)]);
        check(&find.edges(At::X(2.0))?, &[(4, 5)]);
        check(&find.edges(At::Y(1.0))?, &[(2, 3), (2, 5)]);
        check(&find.edges(At::X(0.0))?, &[(0, 3)]);
        check(&find.edges(At::X(1.0))?, &[]); // internal
        check(&find.edges(At::X(10.0))?, &[]); // far away
        Ok(())
    }

    #[test]
    fn find_edges_works_3d() -> Result<(), StrError> {
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
        let mesh = Samples::two_cubes_vertical();
        let shapes = allocate_cell_shapes(&mesh);
        let faces = extract_all_faces(&mesh, &shapes);
        let boundary = extract_features_3d(&mesh, &shapes, &faces, Extract::Boundary);
        let find = Find::new(&mesh, &boundary)?;
        check(
            &find.edges(At::X(0.0))?,
            &[(0, 3), (0, 4), (3, 7), (4, 7), (4, 8), (7, 11), (8, 11)],
        );
        check(
            &find.edges(At::X(1.0))?,
            &[(1, 2), (1, 5), (2, 6), (5, 6), (5, 9), (6, 10), (9, 10)],
        );
        check(&find.edges(At::X(10.0))?, &[]);
        check(
            &find.edges(At::Y(0.0))?,
            &[(0, 1), (0, 4), (1, 5), (4, 5), (4, 8), (5, 9), (8, 9)],
        );
        check(
            &find.edges(At::Y(1.0))?,
            &[(2, 3), (2, 6), (3, 7), (6, 7), (6, 10), (7, 11), (10, 11)],
        );
        check(&find.edges(At::Y(10.0))?, &[]);
        check(&find.edges(At::Z(0.0))?, &[(0, 1), (0, 3), (1, 2), (2, 3)]);
        check(&find.edges(At::Z(2.0))?, &[(8, 9), (8, 11), (9, 10), (10, 11)]);
        check(&find.edges(At::Z(10.0))?, &[]);
        check(&find.edges(At::XY(0.0, 0.0))?, &[(0, 4), (4, 8)]);
        check(&find.edges(At::XY(1.0, 1.0))?, &[(2, 6), (6, 10)]);
        check(&find.edges(At::XY(10.0, 10.0))?, &[]);
        check(&find.edges(At::YZ(0.0, 0.0))?, &[(0, 1)]);
        check(&find.edges(At::YZ(1.0, 1.0))?, &[(6, 7)]);
        check(&find.edges(At::YZ(10.0, 10.0))?, &[]);
        check(&find.edges(At::XZ(0.0, 0.0))?, &[(0, 3)]);
        check(&find.edges(At::XZ(1.0, 0.0))?, &[(1, 2)]);
        check(&find.edges(At::XZ(1.0, 2.0))?, &[(9, 10)]);
        check(&find.edges(At::XZ(10.0, 10.0))?, &[]);
        check(&find.edges(At::XYZ(0.0, 0.0, 0.0))?, &[]);
        assert_eq!(
            find.edges(At::XYZ(10.0, 0.0, 0.0)).err(),
            Some("cannot find point with coordinates outside the grid")
        );
        check(
            &find.edges(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0))?,
            &[(1, 5), (3, 7), (5, 9), (7, 11)],
        );
        check(
            &find.edges(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, SQRT_2))?,
            &[(2, 6), (6, 10)],
        );
        check(&find.edges(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 10.0))?, &[]);
        Ok(())
    }

    #[test]
    fn find_faces_returns_empty_in_2d() -> Result<(), StrError> {
        // 3--------2--------5
        // |        |        |
        // |        |        |
        // |        |        |
        // 0--------1--------4
        let mesh = Samples::two_quads_horizontal();
        let shapes = allocate_cell_shapes(&mesh);
        let edges = extract_all_2d_edges(&mesh, &shapes);
        let boundary = extract_features_2d(&mesh, &shapes, &edges, Extract::Boundary);
        let find = Find::new(&mesh, &boundary)?;
        assert_eq!(find.faces(At::X(0.0))?.len(), 0);
        Ok(())
    }

    #[test]
    fn find_faces_works() -> Result<(), StrError> {
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
        let mesh = Samples::two_cubes_vertical();
        let shapes = allocate_cell_shapes(&mesh);
        let faces = extract_all_faces(&mesh, &shapes);
        let boundary = extract_features_3d(&mesh, &shapes, &faces, Extract::Boundary);
        let find = Find::new(&mesh, &boundary)?;
        check(&find.faces(At::X(0.0))?, &[(0, 3, 4, 7), (4, 7, 8, 11)]);
        check(&find.faces(At::X(1.0))?, &[(1, 2, 5, 6), (5, 6, 9, 10)]);
        check(&find.faces(At::X(10.0))?, &[]);
        check(&find.faces(At::Y(0.0))?, &[(0, 1, 4, 5), (4, 5, 8, 9)]);
        check(&find.faces(At::Y(1.0))?, &[(2, 3, 6, 7), (6, 7, 10, 11)]);
        check(&find.faces(At::Y(10.0))?, &[]);
        check(&find.faces(At::Z(0.0))?, &[(0, 1, 2, 3)]);
        check(&find.faces(At::Z(2.0))?, &[(8, 9, 10, 11)]);
        check(&find.faces(At::Z(10.0))?, &[]);
        check(&find.faces(At::XY(0.0, 0.0))?, &[]);
        check(&find.faces(At::XY(1.0, 1.0))?, &[]);
        check(&find.faces(At::XY(10.0, 10.0))?, &[]);
        check(&find.faces(At::YZ(0.0, 0.0))?, &[]);
        check(&find.faces(At::YZ(1.0, 1.0))?, &[]);
        check(&find.faces(At::YZ(10.0, 10.0))?, &[]);
        check(&find.faces(At::XZ(0.0, 0.0))?, &[]);
        check(&find.faces(At::XZ(1.0, 0.0))?, &[]);
        check(&find.faces(At::XZ(1.0, 2.0))?, &[]);
        check(&find.faces(At::XZ(10.0, 10.0))?, &[]);
        check(&find.faces(At::XYZ(0.0, 0.0, 0.0))?, &[]);
        assert_eq!(
            find.faces(At::XYZ(10.0, 0.0, 0.0)).err(),
            Some("cannot find point with coordinates outside the grid")
        );
        check(&find.faces(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0))?, &[]);
        check(&find.faces(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, SQRT_2))?, &[]);
        check(&find.faces(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 10.0))?, &[]);
        Ok(())
    }

    #[test]
    fn find_works_with_ring() -> Result<(), StrError> {
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
        let shapes = allocate_cell_shapes(&mesh);
        let edges = extract_all_2d_edges(&mesh, &shapes);
        let boundary = extract_features_2d(&mesh, &shapes, &edges, Extract::Boundary);
        let find = Find::new(&mesh, &boundary)?;
        let (r, rr) = (1.0, 2.0);
        check(&find.points(At::XY(1.00, 0.00))?, &[0]);
        check(&find.points(At::XY(1.25, 0.00))?, &[15]);
        check(&find.points(At::XY(1.50, 0.00))?, &[1]);
        check(&find.points(At::XY(1.75, 0.00))?, &[16]);
        check(&find.points(At::XY(2.00, 0.00))?, &[2]);
        check(&find.points(At::XY(0.00, 1.00))?, &[12]);
        check(&find.points(At::XY(0.00, 1.25))?, &[23]);
        check(&find.points(At::XY(0.00, 1.75))?, &[24]);
        check(&find.points(At::XY(0.00, 1.50))?, &[13]);
        check(&find.points(At::XY(0.00, 2.00))?, &[14]);
        check(&find.points(At::XY(SQRT_2 / 2.0, SQRT_2 / 2.0))?, &[6]);
        check(&find.points(At::XY(SQRT_2, SQRT_2))?, &[8]);
        check(
            &find.points(At::Circle(0.0, 0.0, r))?,
            &[0, 3, 6, 9, 12, 25, 28, 31, 34],
        );
        check(
            &find.points(At::Circle(0.0, 0.0, rr))?,
            &[2, 5, 8, 11, 14, 27, 30, 33, 36],
        );
        check(&find.edges(At::Y(0.0))?, &[(0, 1), (1, 2)]);
        check(&find.edges(At::X(0.0))?, &[(12, 13), (13, 14)]);
        check(
            &find.edges(At::Circle(0.0, 0.0, r))?,
            &[(0, 3), (3, 6), (6, 9), (9, 12)],
        );
        check(
            &find.edges(At::Circle(0.0, 0.0, rr))?,
            &[(2, 5), (5, 8), (8, 11), (11, 14)],
        );
        Ok(())
    }
}
