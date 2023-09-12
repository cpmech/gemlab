use super::{get_mesh_limits, set_pad_coords, EdgeKey, FaceKey, Feature, Mesh};
use crate::shapes::DET_JAC_NOT_AVAILABLE;
use crate::shapes::{geo_case, GeoCase, Scratchpad};
use crate::util::GridSearch;
use crate::StrError;
use russell_chk::approx_eq;
use russell_lab::math::ONE_BY_3;
use russell_lab::Vector;
use std::collections::HashMap;

/// Checks if the IDs of points and cells are sequential and within bounds
///
/// This function checks that:
///
/// 1. The position of a point in the points vector matches the id of the point
/// 2. The position of a cell in the cells vector matches the id of the cell
/// 3. The number of points of a cell matches the [crate::shapes::GeoKind]
/// 4. The points of a cell are within the range of available points; i.e., 0 and npoint
pub fn check_ids_and_kind(mesh: &Mesh) -> Result<(), StrError> {
    let npoint = mesh.points.len();
    for i in 0..npoint {
        if mesh.points[i].id != i {
            return Err("incorrect point id found; ids must be sequential");
        }
    }
    for i in 0..mesh.cells.len() {
        if mesh.cells[i].id != i {
            return Err("incorrect cell id found; ids must be sequential");
        }
        if mesh.cells[i].points.len() != mesh.cells[i].kind.nnode() {
            return Err("number of cell points does not correspond to kind.nnode()");
        }
        for p in &mesh.cells[i].points {
            if *p >= npoint {
                return Err("the id of a point specified in the points list of a cell is out of bounds");
            }
        }
    }
    Ok(())
}

/// Checks if the determinants of the Jacobian of all cells are non-negative
pub fn check_jacobian(mesh: &Mesh) -> Result<(), StrError> {
    let ksi = [ONE_BY_3, ONE_BY_3, ONE_BY_3];
    for cell in &mesh.cells {
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind)?;
        set_pad_coords(&mut pad, &cell.points, mesh);
        let det_jac = pad.calc_jacobian(&ksi)?;
        if geo_case(cell.kind.ndim(), mesh.ndim) == GeoCase::Shell {
            assert_eq!(det_jac, DET_JAC_NOT_AVAILABLE);
        } else {
            if det_jac < 0.0 {
                return Err("negative determinant of Jacobian found");
            }
        }
    }
    Ok(())
}

/// This is a convenience function that calls check ids and jacobian functions
///
/// Calls [check_ids_and_kind()] and [check_jacobian()]
pub fn check_all(mesh: &Mesh) -> Result<(), StrError> {
    check_ids_and_kind(mesh)?;
    check_jacobian(mesh)
}

/// Checks if unit normal vectors of some 2D edges are correct
///
/// Note: the solutions map holds the magnitude of the normal, followed by the unit normal.
pub fn check_2d_edge_normals(
    mesh: &Mesh,
    edges: &HashMap<EdgeKey, Feature>,
    solutions: &HashMap<EdgeKey, (f64, [f64; 2])>,
    tolerance: f64,
) -> Result<(), StrError> {
    let ksi = &[0.0, 0.0];
    let mut un = Vector::new(mesh.ndim);
    for (edge_key, (correct_mag_n, correct_un)) in solutions {
        let edge = edges.get(edge_key).ok_or("cannot find edge_key in edges map")?;
        let mut pad = Scratchpad::new(mesh.ndim, edge.kind)?;
        set_pad_coords(&mut pad, &edge.points, mesh);
        let mag_n = pad.calc_normal_vector(&mut un, ksi)?;
        approx_eq(mag_n, *correct_mag_n, tolerance);
        for i in 0..mesh.ndim {
            if f64::abs(un[i] - correct_un[i]) > tolerance {
                return Err("wrong 2d edge unit normal vector found");
            }
        }
    }
    Ok(())
}

/// Checks if unit normal vectors of some faces are correct
///
/// Note: the solutions map holds the magnitude of the normal, followed by the unit normal.
pub fn check_face_normals(
    mesh: &Mesh,
    faces: &HashMap<FaceKey, Feature>,
    solutions: &HashMap<FaceKey, (f64, [f64; 3])>,
    tolerance: f64,
) -> Result<(), StrError> {
    let ksi = &[0.0, 0.0, 0.0];
    let mut un = Vector::new(mesh.ndim);
    for (face_key, (correct_mag_n, correct_un)) in solutions {
        let face = faces.get(face_key).ok_or("cannot find face_key in faces map")?;
        let mut pad = Scratchpad::new(mesh.ndim, face.kind)?;
        set_pad_coords(&mut pad, &face.points, mesh);
        let mag_n = pad.calc_normal_vector(&mut un, ksi)?;
        approx_eq(mag_n, *correct_mag_n, tolerance);
        for i in 0..mesh.ndim {
            if f64::abs(un[i] - correct_un[i]) > tolerance {
                return Err("wrong face unit normal vector found");
            }
        }
    }
    Ok(())
}

/// Checks if there are overlapping points
pub fn check_overlapping_points(mesh: &Mesh, tol: f64) -> Result<(), StrError> {
    let (xmin, xmax) = get_mesh_limits(mesh);
    let mut grid = GridSearch::new(&xmin, &xmax, Some(5), Some(tol), None)?;
    for point in &mesh.points {
        grid.insert(point.id, &point.coords)?;
    }
    for point in &mesh.points {
        match grid.find(&point.coords)? {
            Some(id) => {
                if id != point.id {
                    println!("found overlapping points: {} => {}", id, point.id);
                    return Err("found overlapping points");
                }
            }
            None => (),
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{check_2d_edge_normals, check_ids_and_kind, check_jacobian, check_overlapping_points};
    use crate::mesh::{check_face_normals, Cell, Feature, Mesh, Point};
    use crate::shapes::GeoKind;
    use std::collections::HashMap;

    #[test]
    fn check_ids_and_kind_works() {
        //  3-----------2-----------5
        //  |           |           |
        //  |    [0]    |    [1]    |
        //  |    (1)    |    (2)    |
        //  |           |           |
        //  0-----------1-----------4
        #[rustfmt::skip]
        let mut mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0] },
                Point { id: 3, coords: vec![0.0, 1.0] },
                Point { id: 4, coords: vec![2.0, 0.0] },
                Point { id: 5, coords: vec![2.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
            ],
        };
        check_ids_and_kind(&mesh).expect("should not fail");

        mesh.points[0].id = 6;
        assert_eq!(
            check_ids_and_kind(&mesh).err(),
            Some("incorrect point id found; ids must be sequential")
        );
        mesh.points[0].id = 0;

        assert_eq!(check_ids_and_kind(&mesh).expect("should not fail"), ());

        mesh.cells[0].id = 2;
        assert_eq!(
            check_ids_and_kind(&mesh).err(),
            Some("incorrect cell id found; ids must be sequential")
        );
        mesh.cells[0].id = 0;

        mesh.cells[0].points[0] = 8;
        assert_eq!(
            check_ids_and_kind(&mesh).err(),
            Some("the id of a point specified in the points list of a cell is out of bounds")
        );

        mesh.cells[0].points = vec![0, 1, 2];
        assert_eq!(
            check_ids_and_kind(&mesh).err(),
            Some("number of cell points does not correspond to kind.nnode()")
        );
    }

    #[test]
    fn check_jacobian_works() {
        //  3-----------2-----------5
        //  |           |           |
        //  |    [0]    |    [1]    |
        //  |    (1)    |    (2)    |
        //  |           |           |
        //  0-----------1-----------4
        #[rustfmt::skip]
        let mut mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0] },
                Point { id: 3, coords: vec![0.0, 1.0] },
                Point { id: 4, coords: vec![2.0, 0.0] },
                Point { id: 5, coords: vec![2.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
            ],
        };
        check_jacobian(&mesh).expect("should not fail");

        mesh.cells[0].points[1] = 3;
        mesh.cells[0].points[2] = 3;
        assert_eq!(
            check_jacobian(&mesh).err(),
            Some("cannot compute inverse due to zero determinant")
        );

        mesh.cells[0].points = vec![0, 3, 2, 1];
        assert_eq!(
            check_jacobian(&mesh).err(),
            Some("negative determinant of Jacobian found")
        );
    }

    #[test]
    fn check_2d_edge_normals_works() {
        //
        //             ^         ^
        //             |         |
        //             |         |
        //        3-------->2-------->5
        //        ^         |         |
        // <===== |   [0]   |   [1]   | =====>
        //        |         |         V
        //        0<--------1<--------4
        //             |         |
        //             |         |
        //             V         V
        //
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0] },
                Point { id: 3, coords: vec![0.0, 1.0] },
                Point { id: 4, coords: vec![2.0, 0.0] },
                Point { id: 5, coords: vec![2.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
            ],
        };

        // the magnitude (l) of the normal vector should be equal to
        // 0.5 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
        let l = 0.5; // magnitude of normal vector

        #[rustfmt::skip]
        let mut edges = HashMap::from([
            ((0, 3), Feature { kind: GeoKind::Lin2, points: vec![0, 3] }),
            ((2, 3), Feature { kind: GeoKind::Lin2, points: vec![3, 2] }),
            ((2, 5), Feature { kind: GeoKind::Lin2, points: vec![2, 5] }),
            ((4, 5), Feature { kind: GeoKind::Lin2, points: vec![5, 4] }),
            ((1, 4), Feature { kind: GeoKind::Lin2, points: vec![4, 1] }),
            ((0, 1), Feature { kind: GeoKind::Lin2, points: vec![1, 0] }),
        ]);
        let solutions = HashMap::from([
            ((0, 3), (l, [-1.0, 0.0])),
            ((2, 3), (l, [0.0, 1.0])),
            ((2, 5), (l, [0.0, 1.0])),
            ((4, 5), (l, [1.0, 0.0])),
            ((1, 4), (l, [0.0, -1.0])),
            ((0, 1), (l, [0.0, -1.0])),
        ]);
        check_2d_edge_normals(&mesh, &edges, &solutions, 1e-15).expect("should not fail");

        let points = &mut edges.get_mut(&(0, 3)).unwrap().points;
        points[0] = 3;
        points[1] = 0;
        assert_eq!(
            check_2d_edge_normals(&mesh, &edges, &solutions, 1e-15).err(),
            Some("wrong 2d edge unit normal vector found")
        );

        let solutions = HashMap::from([((10, 20), (l, [0.0, 1.0]))]);
        assert_eq!(
            check_2d_edge_normals(&mesh, &edges, &solutions, 1e-15).err(),
            Some("cannot find edge_key in edges map")
        );
    }

    #[test]
    fn check_face_normals_works() {
        //       4--------------7  1.0
        //      /.             /|
        //     / .            / |
        //    /  .           /  |
        //   /   .          /   |
        //  5--------------6    |          z
        //  |    .         |    |          ↑
        //  |    0---------|----3  0.0     o → y
        //  |   /          |   /          ↙
        //  |  /           |  /          x
        //  | /            | /
        //  |/             |/
        //  1--------------2   1.0
        // 0.0            1.0
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 3,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0, 0.0] },
                Point { id: 3, coords: vec![0.0, 1.0, 0.0] },
                Point { id: 4, coords: vec![0.0, 0.0, 1.0] },
                Point { id: 5, coords: vec![1.0, 0.0, 1.0] },
                Point { id: 6, coords: vec![1.0, 1.0, 1.0] },
                Point { id: 7, coords: vec![0.0, 1.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Hex8, points: vec![0,1,2,3, 4,5,6,7] },
            ],
        };

        // the magnitude (l) of the normal vector should be equal to
        // 0.25 = face_area / 4.0 where 4.0 corresponds to the face_area in the reference system
        let l = 0.25; // magnitude of normal vector

        #[rustfmt::skip]
        let mut faces = HashMap::from([
            ((0, 3, 4, 7), Feature { kind: GeoKind::Qua4, points: vec![0, 4, 7, 3] }),
            ((1, 2, 5, 6), Feature { kind: GeoKind::Qua4, points: vec![1, 2, 6, 5] }),
            ((0, 1, 4, 5), Feature { kind: GeoKind::Qua4, points: vec![0, 1, 5, 4] }),
            ((2, 3, 6, 7), Feature { kind: GeoKind::Qua4, points: vec![2, 3, 7, 6] }),
            ((0, 1, 2, 3), Feature { kind: GeoKind::Qua4, points: vec![0, 3, 2, 1] }),
            ((4, 5, 6, 7), Feature { kind: GeoKind::Qua4, points: vec![4, 5, 6, 7] }),
        ]);
        let solutions = HashMap::from([
            ((0, 3, 4, 7), (l, [-1.0, 0.0, 0.0])),
            ((1, 2, 5, 6), (l, [1.0, 0.0, 0.0])),
            ((0, 1, 4, 5), (l, [0.0, -1.0, 0.0])),
            ((2, 3, 6, 7), (l, [0.0, 1.0, 0.0])),
            ((0, 1, 2, 3), (l, [0.0, 0.0, -1.0])),
            ((4, 5, 6, 7), (l, [0.0, 0.0, 1.0])),
        ]);
        // check_face_normals(&mesh, &faces, &solutions, 1e-15).expect("should not fail");

        let points = &mut faces.get_mut(&(0, 3, 4, 7)).unwrap().points;
        points[1] = 3;
        points[3] = 4;
        assert_eq!(
            check_face_normals(&mesh, &faces, &solutions, 1e-15).err(),
            Some("wrong face unit normal vector found")
        );

        let solutions = HashMap::from([((10, 20, 30, 40), (l, [0.0, 0.0, 1.0]))]);
        assert_eq!(
            check_face_normals(&mesh, &faces, &solutions, 1e-15).err(),
            Some("cannot find face_key in faces map")
        );
    }

    #[test]
    fn check_overlapping_points_works() {
        #[rustfmt::skip]
        let mut mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0] },
                Point { id: 3, coords: vec![0.0, 1.0] },
                Point { id: 4, coords: vec![2.0, 0.0] },
                Point { id: 5, coords: vec![2.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
            ],
        };
        check_overlapping_points(&mesh, 1e-2).expect("should not fail");

        mesh.points[1].coords[0] = 1e-3;
        assert_eq!(
            check_overlapping_points(&mesh, 1e-2).err(),
            Some("found overlapping points")
        );
    }
}
