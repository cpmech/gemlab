use super::{At, EdgeKey, ExtractedFeatures, FaceKey, Mesh, PointId};
use crate::util::{GridSearch, GsNdiv, GsTol};
use crate::StrError;
use std::collections::{HashMap, HashSet};

/// Implements functions to find points, edges, and faces on the boundary of a mesh
///
/// # Examples
///
/// ## Two-dimensional
///
/// ```
/// use gemlab::mesh::{all_edges_2d, allocate_shapes, At, Cell, Extract, ExtractedFeatures, Find, Mesh, Point};
/// use gemlab::StrError;
///
/// fn main() -> Result<(), StrError> {
///     //  3---------2---------5
///     //  |         |         |
///     //  |   [0]   |   [1]   |
///     //  |         |         |
///     //  0---------1---------4
///     let mesh = Mesh {
///         space_ndim: 2,
///         points: vec![
///             Point { id: 0, coords: vec![0.0, 0.0] },
///             Point { id: 1, coords: vec![1.0, 0.0] },
///             Point { id: 2, coords: vec![1.0, 1.0] },
///             Point { id: 3, coords: vec![0.0, 1.0] },
///             Point { id: 4, coords: vec![2.0, 0.0] },
///             Point { id: 5, coords: vec![2.0, 1.0] },
///         ],
///         cells: vec![
///             Cell { id: 0, attribute_id: 1, geo_ndim: 2, points: vec![0, 1, 2, 3] },
///             Cell { id: 1, attribute_id: 2, geo_ndim: 2, points: vec![1, 4, 5, 2] },
///         ],
///     };
///
///     let shapes = allocate_shapes(&mesh)?;
///     let edges = all_edges_2d(&mesh, &shapes)?;
///     let boundary = ExtractedFeatures::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;
///     let find = Find::new(&mesh, &boundary)?;
///
///     let mut points: Vec<_> = find.points(At::X(2.0))?.iter().copied().collect();
///     points.sort();
///     assert_eq!(points, &[4, 5]);
///
///     let mut edges: Vec<_> = find.edges(At::Y(1.0))?.iter().copied().collect();
///     edges.sort();
///     assert_eq!(edges, &[(2, 3), (2, 5)]);
///     Ok(())
/// }
/// ```
///
/// ## Three-dimensional
///
/// ```
/// use gemlab::mesh::{all_faces_3d, allocate_shapes, At, Cell, Extract, ExtractedFeatures, Find, Mesh, Point};
/// use gemlab::StrError;
///
/// fn main() -> Result<(), StrError> {
///     //          .4--------------7
///     //        ,' |            ,'|
///     //      ,'              ,'  |
///     //    ,'     |        ,'    |
///     //  5'==============6'      |
///     //  |               |       |
///     //  |        |      |       |
///     //  |       ,0- - - | - - - 3
///     //  |     ,'        |     ,'
///     //  |   ,'          |   ,'
///     //  | ,'            | ,'
///     //  1'--------------2'
///     let mesh = Mesh {
///         space_ndim: 3,
///         points: vec![
///             Point { id: 0, coords: vec![0.0, 0.0, 0.0] },
///             Point { id: 1, coords: vec![1.0, 0.0, 0.0] },
///             Point { id: 2, coords: vec![1.0, 1.0, 0.0] },
///             Point { id: 3, coords: vec![0.0, 1.0, 0.0] },
///             Point { id: 4, coords: vec![0.0, 0.0, 1.0] },
///             Point { id: 5, coords: vec![1.0, 0.0, 1.0] },
///             Point { id: 6, coords: vec![1.0, 1.0, 1.0] },
///             Point { id: 7, coords: vec![0.0, 1.0, 1.0] },
///         ],
///         cells: vec![
///             Cell { id: 0, attribute_id: 1, geo_ndim: 3, points: vec![0,1,2,3, 4,5,6,7] },
///         ],
///     };
///
///     let shapes = allocate_shapes(&mesh)?;
///     let faces = all_faces_3d(&mesh, &shapes)?;
///     let boundary = ExtractedFeatures::extract(&mesh, &shapes, None, Some(&faces), Extract::Boundary)?;
///     let find = Find::new(&mesh, &boundary)?;
///
///     let mut points: Vec<_> = find.points(At::XY(1.0, 1.0))?.iter().copied().collect();
///     points.sort();
///     assert_eq!(points, &[2, 6]);
///
///     let mut edges: Vec<_> = find.edges(At::YZ(1.0, 1.0))?.iter().copied().collect();
///     edges.sort();
///     assert_eq!(edges, &[(6, 7)]);
///
///     let mut faces: Vec<_> = find.faces(At::Y(1.0))?.iter().copied().collect();
///     faces.sort();
///     assert_eq!(faces, &[(2, 3, 6, 7)]);
///     Ok(())
/// }
/// ```
pub struct Find {
    /// Space number of dimension
    space_ndim: usize,

    /// Total number of points in the mesh
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
    /// # Panics
    ///
    /// If the `boundary` does not correspond to `mesh`, a panic will occur.
    pub fn new(mesh: &Mesh, boundary: &ExtractedFeatures) -> Result<Self, StrError> {
        // expand limits a little bit to accommodate imprecisions on the (min,max) values
        const PCT: f64 = 1.0 / 100.0; // 1% is enough
        let (mut min, mut max) = (vec![0.0; mesh.space_ndim], vec![0.0; mesh.space_ndim]);
        for i in 0..mesh.space_ndim {
            let del = boundary.max[i] - boundary.min[i];
            min[i] = boundary.min[i] - PCT * del;
            max[i] = boundary.max[i] + PCT * del;
        }

        // add point ids to grid
        let mut grid = GridSearch::new(&min, &max, GsNdiv::Default, GsTol::Default)?;
        for point_id in &boundary.points {
            grid.insert(*point_id, &mesh.points[*point_id].coords)?;
        }
        // map point ids to edges
        let mut point_to_edges: HashMap<PointId, HashSet<EdgeKey>> = HashMap::new();
        for (edge_key, edge) in &boundary.edges {
            for point_id in &edge.points {
                point_to_edges
                    .entry(*point_id)
                    .or_insert(HashSet::new())
                    .insert(*edge_key);
            }
        }
        // map point ids to faces
        let mut point_to_faces: HashMap<PointId, HashSet<FaceKey>> = HashMap::new();
        for (face_key, face) in &boundary.faces {
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

    /// Finds boundary points
    ///
    /// # Input
    ///
    /// * `at` -- the location constraint
    ///
    /// # Output
    ///
    /// * Returns a set of point ids (boundary points only). You may sort
    ///   the point ids using the following code snipped:
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

    /// Finds boundary edges
    ///
    /// # Input
    ///
    /// * `at` -- the location constraint
    ///
    /// # Output
    ///
    /// * Returns a set of edge keys (boundary edges only). You may sort
    ///   the edge keys using the following code snipped:
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

    /// Finds boundary faces
    ///
    /// # Input
    ///
    /// * `at` -- the location constraint
    ///
    /// # Output
    ///
    /// * Returns a set of face keys (boundary faces only). You may sort
    ///   the face keys using the following code snipped:
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
    use crate::mesh::{all_edges_2d, all_faces_3d, allocate_shapes, At, Extract, ExtractedFeatures, Samples};
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
        let shapes = allocate_shapes(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let boundary = ExtractedFeatures::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;
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
        let shapes = allocate_shapes(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let boundary = ExtractedFeatures::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;
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
        let shapes = allocate_shapes(&mesh)?;
        let faces = all_faces_3d(&mesh, &shapes)?;
        let boundary = ExtractedFeatures::extract(&mesh, &shapes, None, Some(&faces), Extract::Boundary)?;
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
        let shapes = allocate_shapes(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let boundary = ExtractedFeatures::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;
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
        let shapes = allocate_shapes(&mesh)?;
        let faces = all_faces_3d(&mesh, &shapes)?;
        let boundary = ExtractedFeatures::extract(&mesh, &shapes, None, Some(&faces), Extract::Boundary)?;
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
        let shapes = allocate_shapes(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let boundary = ExtractedFeatures::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;
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
        let shapes = allocate_shapes(&mesh)?;
        let faces = all_faces_3d(&mesh, &shapes)?;
        let boundary = ExtractedFeatures::extract(&mesh, &shapes, None, Some(&faces), Extract::Boundary)?;
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
        let shapes = allocate_shapes(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let boundary = ExtractedFeatures::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;
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
        let shapes = allocate_shapes(&mesh)?;
        let faces = all_faces_3d(&mesh, &shapes)?;
        let boundary = ExtractedFeatures::extract(&mesh, &shapes, None, Some(&faces), Extract::Boundary)?;
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
        let shapes = allocate_shapes(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let boundary = ExtractedFeatures::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;
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
