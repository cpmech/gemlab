/// Number of divisions along the longest direction for GridSearch
///
/// ```text
/// ndiv_other = min(ndiv_max, max(ndiv_min, (ll_other/ll_long) * ndiv_long))
/// ```
const GRID_SEARCH_NDIV_LONG: usize = 20;

/// Minimum number of divisions for GridSearch
const GRID_SEARCH_NDIV_MIN: usize = 2;

/// Maximum number of divisions for GridSearch
const GRID_SEARCH_NDIV_MAX: usize = 50;

/// Finds boundary points in the mesh
///
/// # Input
///
/// * `at` -- the location constraint
///
/// # Output
///
/// * Returns a **sorted** list of point ids (boundary points only)
///
/// # Note
///
/// `compute_derived_props` must be called first, otherwise this function returns an error
pub fn find_boundary_points(&self, at: At) -> Result<Vec<PointId>, StrError> {
    if !self.derived_props_computed {
        return Err("compute_derived_props must be called first");
    }
    let mut point_ids: HashSet<PointId> = HashSet::new();
    match at {
        At::X(x) => {
            if self.space_ndim == 2 {
                for id in self.grid_boundary_points.find_on_line(&[x, 0.0], &[x, 1.0])? {
                    point_ids.insert(id);
                }
            } else {
                for id in self.grid_boundary_points.find_on_plane_yz(x)? {
                    point_ids.insert(id);
                }
            }
        }
        At::Y(y) => {
            if self.space_ndim == 2 {
                for id in self.grid_boundary_points.find_on_line(&[0.0, y], &[1.0, y])? {
                    point_ids.insert(id);
                }
            } else {
                for id in self.grid_boundary_points.find_on_plane_xz(y)? {
                    point_ids.insert(id);
                }
            }
        }
        At::Z(z) => {
            if self.space_ndim == 2 {
                return Err("At::Z works in 3D only");
            } else {
                for id in self.grid_boundary_points.find_on_plane_xy(z)? {
                    point_ids.insert(id);
                }
            }
        }
        At::XY(x, y) => {
            if self.space_ndim == 2 {
                if let Some(id) = self.grid_boundary_points.find(&[x, y])? {
                    point_ids.insert(id);
                }
            } else {
                for id in self.grid_boundary_points.find_on_line(&[x, y, 0.0], &[x, y, 1.0])? {
                    point_ids.insert(id);
                }
            }
        }
        At::YZ(y, z) => {
            if self.space_ndim == 2 {
                return Err("At::YZ works in 3D only");
            } else {
                for id in self.grid_boundary_points.find_on_line(&[0.0, y, z], &[1.0, y, z])? {
                    point_ids.insert(id);
                }
            }
        }
        At::XZ(x, z) => {
            if self.space_ndim == 2 {
                return Err("At::XZ works in 3D only");
            } else {
                for id in self.grid_boundary_points.find_on_line(&[x, 0.0, z], &[x, 1.0, z])? {
                    point_ids.insert(id);
                }
            }
        }
        At::XYZ(x, y, z) => {
            if self.space_ndim == 2 {
                return Err("At::XYZ works in 3D only");
            } else {
                if let Some(id) = self.grid_boundary_points.find(&[x, y, z])? {
                    point_ids.insert(id);
                }
            }
        }
        At::Circle(x, y, r) => {
            if self.space_ndim == 2 {
                for id in self.grid_boundary_points.find_on_circle(&[x, y], r)? {
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
                for id in self
                    .grid_boundary_points
                    .find_on_cylinder(&[ax, ay, az], &[bx, by, bz], r)?
                {
                    point_ids.insert(id);
                }
            }
        }
    }
    let mut ids: Vec<_> = point_ids.into_iter().collect();
    ids.sort();
    Ok(ids)
}

/// Finds boundary edges in the mesh
///
/// # Input
///
/// * `at` -- the location constraint
///
/// # Output
///
/// * Returns a **sorted** list of edge key ids (boundary edges only)
///
/// # Note
///
/// `compute_derived_props` must be called first, otherwise this function returns an error
pub fn find_boundary_edges(&self, at: At) -> Result<Vec<EdgeKey>, StrError> {
    if !self.derived_props_computed {
        return Err("compute_derived_props must be called first");
    }
    let mut edge_keys: HashSet<EdgeKey> = HashSet::new();
    // find all points near the geometric feature
    let point_ids = &self.find_boundary_points(at)?;
    for point_id in point_ids {
        // loop over all boundary edges touching this point
        let point = &self.points[*point_id];
        for edge_key in &point.shared_by_boundary_edges {
            // check if two edge points pass through the geometric feature
            if point_ids.contains(&edge_key.0) && point_ids.contains(&edge_key.1) {
                if self.boundary_edges.contains_key(&edge_key) {
                    edge_keys.insert(*edge_key);
                }
            }
        }
    }
    let mut keys: Vec<_> = edge_keys.into_iter().collect();
    keys.sort();
    Ok(keys)
}

/// Finds boundary faces in the mesh
///
/// # Input
///
/// * `at` -- the location constraint
///
/// # Output
///
/// * Returns a **sorted** list of face key ids (boundary faces only)
///
/// # Note
///
/// `compute_derived_props` must be called first, otherwise this function returns an error
pub fn find_boundary_faces(&self, at: At) -> Result<Vec<FaceKey>, StrError> {
    if !self.derived_props_computed {
        return Err("compute_derived_props must be called first");
    }
    let mut face_keys: HashSet<FaceKey> = HashSet::new();
    // find all points near the geometric feature
    let point_ids = &self.find_boundary_points(at)?;
    for point_id in point_ids {
        // loop over all boundary faces touching this point
        let point = &self.points[*point_id];
        for face_key in &point.shared_by_boundary_faces {
            // check if the fourth point_id in the face_key is ok
            let fourth_is_ok = if face_key.3 == self.points.len() {
                true
            } else {
                point_ids.contains(&face_key.3)
            };
            // check if the face points pass through the geometric feature
            if fourth_is_ok
                && point_ids.contains(&face_key.0)
                && point_ids.contains(&face_key.1)
                && point_ids.contains(&face_key.2)
            {
                if self.boundary_faces.contains_key(&face_key) {
                    face_keys.insert(*face_key);
                }
            }
        }
    }
    let mut keys: Vec<_> = face_keys.into_iter().collect();
    keys.sort();
    Ok(keys)
}

#[test]
fn grid_search_initialization_works_2d() -> Result<(), StrError> {
    //
    // 3.1   6---------12
    //       |    [6]   |       SOLID
    // 3.0   5---------11       <-- layer separation
    //       | [5]  ,-' |       POROUS
    //       |   ,-'    |
    // 2.5   4.-'       |
    //       | '.   [4] |       L
    //       | [3].     |       A
    // 2.0   3.    '.   |       Y
    //       | '--__ '. |       E
    //       |      '--10  1.8  R
    //       |          |       2
    //       |   [2]    |
    //       |          |
    // 1.0   2----------9       <-- layer separation
    //       |          |       L
    //       |    [1]   |       A
    // 0.5   1.__       |       Y
    //       |   '--..  |       E
    //       |  [0]   '-8  0.2  R
    // 0.0   0----------7       1
    //
    //      0.0        1.0
    //
    let mesh = Mesh::from_text_file("./data/meshes/column_distorted_tris_quads.msh")?;

    assert_eq!(mesh.space_ndim, 2);
    assert_eq!(mesh.points.len(), 13);
    assert_eq!(mesh.cells.len(), 7);
    assert_eq!(mesh.boundary_points.len(), 13);
    assert_eq!(mesh.boundary_edges.len(), 13);
    assert_eq!(mesh.boundary_faces.len(), 0);
    assert_eq!(mesh.coords_min, &[0.0, 0.0]);
    assert_eq!(mesh.coords_max, &[1.0, 3.1]);

    // DO NOT REMOVE THE CODE BELOW
    // println!("{}", mesh.grid_boundary_points);
    // let mut plot = mesh.grid_boundary_points.plot()?;
    // plot.set_figure_size_points(600.0, 1200.0)
    //     .set_equal_axes(true)
    //     .save("/tmp/gemlab/grid_search_initialization_works_2d.svg")?;

    assert_eq!(
        format!("{}", mesh.grid_boundary_points),
        "0: [0]\n\
             5: [7]\n\
             11: [8]\n\
             18: [1]\n\
             36: [2]\n\
             41: [9]\n\
             71: [10]\n\
             72: [3]\n\
             96: [4]\n\
             114: [5, 6]\n\
             119: [11, 12]\n\
             ids = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]\n\
             nitem = 13\n\
             ncontainer = 11\n\
             ndiv = [6, 20]\n"
    );
    Ok(())
}

#[test]
fn find_boundary_fails_on_wrong_input() -> Result<(), StrError> {
    let mesh = Mesh::new(2)?;
    assert_eq!(
        mesh.find_boundary_points(At::XY(0.0, 0.0)).err(),
        Some("compute_derived_props must be called first")
    );
    assert_eq!(
        mesh.find_boundary_edges(At::XY(0.0, 0.0)).err(),
        Some("compute_derived_props must be called first")
    );
    let mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;
    assert_eq!(
        mesh.find_boundary_points(At::Z(0.0)).err(),
        Some("At::Z works in 3D only")
    );
    assert_eq!(
        mesh.find_boundary_points(At::YZ(0.0, 0.0)).err(),
        Some("At::YZ works in 3D only")
    );
    assert_eq!(
        mesh.find_boundary_points(At::XZ(0.0, 0.0)).err(),
        Some("At::XZ works in 3D only")
    );
    assert_eq!(
        mesh.find_boundary_points(At::XYZ(0.0, 0.0, 0.0)).err(),
        Some("At::XYZ works in 3D only")
    );
    assert_eq!(
        mesh.find_boundary_points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
            .err(),
        Some("At::Cylinder works in 3D only")
    );
    let mesh = Mesh::from_text_file("./data/meshes/ok2.msh")?;
    assert_eq!(
        mesh.find_boundary_points(At::Circle(0.0, 0.0, 0.0)).err(),
        Some("At::Circle works in 2D only")
    );
    let mesh = Mesh::new(3)?;
    assert_eq!(
        mesh.find_boundary_faces(At::Z(0.0)).err(),
        Some("compute_derived_props must be called first")
    );
    Ok(())
}

#[test]
fn find_boundary_works_2d() -> Result<(), StrError> {
    // `.       `.
    //   3--------2--------5
    //   | `.     | `.     |
    //   |   `~.  |        |
    //   |      `.|        |
    //   0--------1--------4
    //
    let mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;
    assert_eq!(mesh.find_boundary_points(At::XY(0.0, 0.0))?, &[0]);
    assert_eq!(mesh.find_boundary_points(At::XY(2.0, 1.0))?, &[5]);
    assert_eq!(
        mesh.find_boundary_points(At::XY(10.0, 0.0)).err(),
        Some("point is outside the grid")
    );
    assert_eq!(mesh.find_boundary_edges(At::Y(0.0))?, &[(0, 1), (1, 4)]);
    assert_eq!(mesh.find_boundary_edges(At::X(2.0))?, &[(4, 5)]);
    assert_eq!(mesh.find_boundary_edges(At::Y(1.0))?, &[(2, 3), (2, 5)]);
    assert_eq!(mesh.find_boundary_edges(At::X(0.0))?, &[(0, 3)]);
    assert_eq!(mesh.find_boundary_edges(At::X(10.0))?, &[]);
    assert_eq!(mesh.find_boundary_points(At::Circle(0.0, 0.0, 1.0))?, &[1, 3]);
    assert_eq!(mesh.find_boundary_points(At::Circle(0.0, 0.0, SQRT_2))?, &[2]);
    assert_eq!(mesh.find_boundary_points(At::Circle(0.0, 0.0, 10.0))?, &[] as &[usize]);
    Ok(())
}

#[test]
fn find_boundary_works_3d() -> Result<(), StrError> {
    //
    //       8-------------11
    //      /.             /|
    //     / .            / |
    //    /  .           /  |
    //   /   .          /   |
    //  9-------------10    |
    //  |    .         |    |
    //  |    4---------|----7
    //  |   /.         |   /|
    //  |  / .         |  / |
    //  | /  .         | /  |
    //  |/   .         |/   |
    //  5--------------6    |
    //  |    .         |    |
    //  |    0---------|----3
    //  |   /          |   /
    //  |  /           |  /
    //  | /            | /
    //  |/             |/
    //  1--------------2
    //
    let mesh = Mesh::from_text_file("./data/meshes/ok2.msh")?;
    assert_eq!(mesh.find_boundary_points(At::X(0.0))?, &[0, 3, 4, 7, 8, 11]);
    assert_eq!(mesh.find_boundary_points(At::X(1.0))?, &[1, 2, 5, 6, 9, 10]);
    assert_eq!(mesh.find_boundary_points(At::X(10.0))?, &[] as &[usize]);
    assert_eq!(mesh.find_boundary_points(At::Y(0.0))?, &[0, 1, 4, 5, 8, 9]);
    assert_eq!(mesh.find_boundary_points(At::Y(1.0))?, &[2, 3, 6, 7, 10, 11]);
    assert_eq!(mesh.find_boundary_points(At::Y(10.0))?, &[] as &[usize]);
    assert_eq!(mesh.find_boundary_points(At::Z(0.0))?, &[0, 1, 2, 3]);
    assert_eq!(mesh.find_boundary_points(At::Z(1.0))?, &[4, 5, 6, 7]);
    assert_eq!(mesh.find_boundary_points(At::Z(2.0))?, &[8, 9, 10, 11]);
    assert_eq!(mesh.find_boundary_points(At::Z(10.0))?, &[] as &[usize]);
    assert_eq!(mesh.find_boundary_points(At::XY(0.0, 0.0))?, &[0, 4, 8]);
    assert_eq!(mesh.find_boundary_points(At::XY(1.0, 1.0))?, &[2, 6, 10]);
    assert_eq!(mesh.find_boundary_points(At::XY(10.0, 10.0))?, &[] as &[usize]);
    assert_eq!(mesh.find_boundary_points(At::YZ(0.0, 0.0))?, &[0, 1]);
    assert_eq!(mesh.find_boundary_points(At::YZ(1.0, 1.0))?, &[6, 7]);
    assert_eq!(mesh.find_boundary_points(At::XZ(0.0, 0.0))?, &[0, 3]);
    assert_eq!(mesh.find_boundary_points(At::XZ(1.0, 0.0))?, &[1, 2]);
    assert_eq!(mesh.find_boundary_points(At::XZ(1.0, 2.0))?, &[9, 10]);
    assert_eq!(mesh.find_boundary_points(At::XYZ(0.0, 0.0, 0.0))?, &[0]);
    assert_eq!(mesh.find_boundary_points(At::XYZ(1.0, 1.0, 2.0))?, &[10]);
    assert_eq!(
        mesh.find_boundary_points(At::XYZ(10.0, 0.0, 0.0)).err(),
        Some("point is outside the grid")
    );
    assert_eq!(
        mesh.find_boundary_points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0))?,
        &[1, 3, 5, 7, 9, 11]
    );
    assert_eq!(
        mesh.find_boundary_points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, SQRT_2))?,
        &[2, 6, 10]
    );
    assert_eq!(
        mesh.find_boundary_points(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 10.0))?,
        &[] as &[usize]
    );
    assert_eq!(mesh.find_boundary_faces(At::X(0.0))?, &[(0, 3, 4, 7), (4, 7, 8, 11)]);
    assert_eq!(mesh.find_boundary_faces(At::X(1.0))?, &[(1, 2, 5, 6), (5, 6, 9, 10)]);
    assert_eq!(mesh.find_boundary_faces(At::Y(0.0))?, &[(0, 1, 4, 5), (4, 5, 8, 9)]);
    assert_eq!(mesh.find_boundary_faces(At::Y(1.0))?, &[(2, 3, 6, 7), (6, 7, 10, 11)]);
    assert_eq!(mesh.find_boundary_faces(At::Z(0.0))?, &[(0, 1, 2, 3)]);
    assert_eq!(mesh.find_boundary_faces(At::Z(2.0))?, &[(8, 9, 10, 11)]);
    Ok(())
}
