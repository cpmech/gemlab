use super::{allocate_state, Edge, EdgeKey, Face, FaceKey, Mesh};
use crate::shapes::Shape;
use crate::StrError;
use russell_lab::Vector;
use std::collections::HashMap;

/// Checks if the IDs of points and cells are sequential
///
/// This function checks that:
///
/// * the index of a point in the points vector matches the id of the point
/// * the index of a cell in the cells vector matches the id of the cell
pub fn check_ids_and_kind(mesh: &Mesh) -> Result<(), StrError> {
    for i in 0..mesh.points.len() {
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
    }
    Ok(())
}

/// Checks if the determinant of the Jacobian of all cells are non-negative
pub fn check_jacobian(mesh: &Mesh) -> Result<(), StrError> {
    let ksi = [0.0, 0.0, 0.0];
    for cell in &mesh.cells {
        let shape = Shape::new(cell.kind);
        let mut state = allocate_state(mesh, cell.kind, &cell.points)?;
        let det_jac = shape.calc_jacobian(&mut state, &ksi)?;
        if det_jac < 0.0 {
            return Err("negative determinant of Jacobian found");
        }
    }
    Ok(())
}

/// Checks if normal vectors of some 2D edges are correct
pub fn check_2d_edge_normals(
    mesh: &Mesh,
    edges: &HashMap<EdgeKey, Edge>,
    solutions: &HashMap<EdgeKey, [f64; 2]>,
    tolerance: f64,
) -> Result<(), StrError> {
    let ksi = &[0.0, 0.0];
    let mut normal = Vector::new(mesh.ndim);
    for (edge_key, solution) in solutions {
        let edge = edges.get(edge_key).ok_or("cannot find edge_key in edges map")?;
        let shape = Shape::new(edge.kind);
        let mut state = allocate_state(mesh, shape.kind, &edge.points)?;
        shape.calc_boundary_normal(&mut normal, &mut state, ksi)?;
        for i in 0..mesh.ndim {
            if f64::abs(normal[i] - solution[i]) > tolerance {
                return Err("wrong normal vector found");
            }
        }
    }
    Ok(())
}

/// Checks if normal vectors of some faces are correct
pub fn check_face_normals(
    mesh: &Mesh,
    faces: &HashMap<FaceKey, Face>,
    solutions: &HashMap<FaceKey, [f64; 3]>,
    tolerance: f64,
) -> Result<(), StrError> {
    let ksi = &[0.0, 0.0, 0.0];
    let mut normal = Vector::new(mesh.ndim);
    for (face_key, solution) in solutions {
        let face = faces.get(face_key).ok_or("cannot find face_key in faces map")?;
        let shape = Shape::new(face.kind);
        let mut state = allocate_state(mesh, face.kind, &face.points)?;
        shape.calc_boundary_normal(&mut normal, &mut state, ksi)?;
        for i in 0..mesh.ndim {
            if f64::abs(normal[i] - solution[i]) > tolerance {
                return Err("wrong normal vector found");
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{check_2d_edge_normals, check_ids_and_kind, check_jacobian};
    use crate::mesh::{check_face_normals, Cell, Edge, EdgeKey, Face, FaceKey, Mesh, Point};
    use crate::shapes::GeoKind;
    use std::collections::HashMap;

    #[test]
    fn check_ids_works() {
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
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute_id: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
            ],
        };
        assert_eq!(check_ids_and_kind(&mesh).expect("should not fail"), ());

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
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute_id: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
            ],
        };
        assert_eq!(check_jacobian(&mesh).expect("should not fail"), ());

        mesh.cells[0].points[2] = 3;
        mesh.cells[0].points[3] = 2;
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
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute_id: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
            ],
        };

        // the magnitude (l) of the normal vector should be equal to
        // 0.5 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
        let l = 0.5; // magnitude of normal vector

        #[rustfmt::skip]
        let mut edges: HashMap<EdgeKey, Edge> = HashMap::from([
            ((0, 3), Edge { kind: GeoKind::Lin2, points: vec![0, 3] }),
            ((2, 3), Edge { kind: GeoKind::Lin2, points: vec![3, 2] }),
            ((2, 5), Edge { kind: GeoKind::Lin2, points: vec![2, 5] }),
            ((4, 5), Edge { kind: GeoKind::Lin2, points: vec![5, 4] }),
            ((1, 4), Edge { kind: GeoKind::Lin2, points: vec![4, 1] }),
            ((0, 1), Edge { kind: GeoKind::Lin2, points: vec![1, 0] }),
        ]);
        let solutions: HashMap<EdgeKey, [f64; 2]> = HashMap::from([
            ((0, 3), [-l, 0.0]),
            ((2, 3), [0.0, l]),
            ((2, 5), [0.0, l]),
            ((4, 5), [l, 0.0]),
            ((1, 4), [0.0, -l]),
            ((0, 1), [0.0, -l]),
        ]);
        assert_eq!(check_2d_edge_normals(&mesh, &edges, &solutions, 1e-15).expect("ok"), ());

        let points = &mut edges.get_mut(&(0, 3)).unwrap().points;
        points[0] = 3;
        points[1] = 0;
        assert_eq!(
            check_2d_edge_normals(&mesh, &edges, &solutions, 1e-15).err(),
            Some("wrong normal vector found")
        );

        let solutions: HashMap<EdgeKey, [f64; 2]> = HashMap::from([((10, 20), [0.0, l])]);
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
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Hex8, points: vec![0,1,2,3, 4,5,6,7] },
            ],
        };

        // the magnitude (l) of the normal vector should be equal to
        // 0.25 = face_area / 4.0 where 4.0 corresponds to the face_area in the reference system
        let l = 0.25; // magnitude of normal vector

        #[rustfmt::skip]
        let mut faces: HashMap<FaceKey, Face> = HashMap::from([
            ((0, 3, 4, 7), Face { kind: GeoKind::Qua4, points: vec![0, 4, 7, 3] }),
            ((1, 2, 5, 6), Face { kind: GeoKind::Qua4, points: vec![1, 2, 6, 5] }),
            ((0, 1, 4, 5), Face { kind: GeoKind::Qua4, points: vec![0, 1, 5, 4] }),
            ((2, 3, 6, 7), Face { kind: GeoKind::Qua4, points: vec![2, 3, 7, 6] }),
            ((0, 1, 2, 3), Face { kind: GeoKind::Qua4, points: vec![0, 3, 2, 1] }),
            ((4, 5, 6, 7), Face { kind: GeoKind::Qua4, points: vec![4, 5, 6, 7] }),
        ]);
        let solutions: HashMap<FaceKey, [f64; 3]> = HashMap::from([
            ((0, 3, 4, 7), [-l, 0.0, 0.0]),
            ((1, 2, 5, 6), [l, 0.0, 0.0]),
            ((0, 1, 4, 5), [0.0, -l, 0.0]),
            ((2, 3, 6, 7), [0.0, l, 0.0]),
            ((0, 1, 2, 3), [0.0, 0.0, -l]),
            ((4, 5, 6, 7), [0.0, 0.0, l]),
        ]);
        assert_eq!(check_face_normals(&mesh, &faces, &solutions, 1e-15).expect("ok"), ());

        let points = &mut faces.get_mut(&(0, 3, 4, 7)).unwrap().points;
        points[0] = 4;
        points[1] = 0;
        assert_eq!(
            check_face_normals(&mesh, &faces, &solutions, 1e-15).err(),
            Some("wrong normal vector found")
        );

        let solutions: HashMap<FaceKey, [f64; 3]> = HashMap::from([((10, 20, 30, 40), [0.0, 0.0, l])]);
        assert_eq!(
            check_face_normals(&mesh, &faces, &solutions, 1e-15).err(),
            Some("cannot find face_key in faces map")
        );
    }
}
