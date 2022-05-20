use super::{At, Boundary, EdgeKey, FaceKey, Mesh, PointId};
use crate::util::{GridSearch, GsNdiv, GsTol};
use crate::StrError;
use std::collections::{HashMap, HashSet};

/// Implements functions to find points, edges, and faces on the boundary of a mesh
pub struct Find {
    space_ndim: usize,
    grid: GridSearch,
    point_to_edges: HashMap<PointId, HashSet<EdgeKey>>,
    point_to_faces: HashMap<PointId, HashSet<FaceKey>>,
}

impl Find {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, boundary: &Boundary) -> Result<Self, StrError> {
        // add point ids to grid
        let mut grid = GridSearch::new(&boundary.min, &boundary.max, GsNdiv::Default, GsTol::Default)?;
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
            grid,
            point_to_edges,
            point_to_faces,
        })
    }

    /// Finds boundary points in the mesh
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
    /// let mut ids: Vec<_> = point_ids.iter().collect();
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

    /// Finds boundary edges in the mesh
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
    /// let mut keys: Vec<_> = edge_keys.iter().collect();
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
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Find;
    use crate::mesh::{At, Boundary, Samples};
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
        let boundary = Boundary::new(&mesh)?;
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
        let find = Find::new(&mesh, &Boundary::new(&mesh)?)?;
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
        let find = Find::new(&mesh, &Boundary::new(&mesh)?)?;
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
        let find = Find::new(&mesh, &Boundary::new(&mesh)?)?;
        check(&find.points(At::XY(0.0, 0.0))?, &[0]);
        check(&find.points(At::XY(2.0, 1.0))?, &[5]);
        assert_eq!(find.points(At::XY(10.0, 0.0)).err(), Some("point is outside the grid"));
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
        let find = Find::new(&mesh, &Boundary::new(&mesh)?)?;
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
            Some("point is outside the grid")
        );
        check(
            &find.points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0))?,
            &[1, 3, 5, 7, 9, 11],
        );
        check(
            &find.points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, SQRT_2))?,
            &[2, 6, 10],
        );
        check(
            &find.points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 10.0))?,
            &[] as &[usize],
        );
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
        let find = Find::new(&mesh, &Boundary::new(&mesh)?)?;
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
        let find = Find::new(&mesh, &Boundary::new(&mesh)?)?;
        check(
            &find.edges(At::X(0.0))?,
            &[(0, 3), (0, 4), (3, 7), (4, 7), (4, 8), (7, 11), (8, 11)],
        );
        Ok(())
    }
}
