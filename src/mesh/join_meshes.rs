use super::{Cell, Features, Mesh, Point};
use crate::util::GridSearch;
use crate::StrError;
use russell_lab::{sort2, sort4};
use std::collections::HashMap;

/// Joins two meshes by comparing the coordinates on the boundary of the first mesh
///
/// # Notes
///
/// 1. The meshes must have the same space_ndim.
/// 2. The boundary of mesh A is used to find overlapping points with mesh B.
/// 3. The boundary markers of mesh A have precedence over those of mesh B.
///
/// **Important:** This function does not guarantee the "mesh compatibility" requirements
/// for finite element analyses.
fn join_two_meshes(a: &Mesh, b: &Mesh) -> Result<Mesh, StrError> {
    // check
    if a.ndim != b.ndim {
        return Err("meshes must have the same ndim");
    }

    // find the boundary of mesh A
    let boundary_a = Features::new(a, false);

    // allocate and prepare a GridSearch for mesh A
    let mut grid_a = GridSearch::new(&boundary_a.min, &boundary_a.max, None, None, None)?;
    for m in 0..a.points.len() {
        grid_a.insert(a.points[m].id, &a.points[m].coords)?;
    }

    // create a new mesh with all mesh A data, except the marked edges/faces
    let mut mesh = a.clone();

    // renumber the points of mesh B and add them to the new mesh (if not present yet)
    let mut new_point_id = a.points.len();
    let mut map_old_to_new_point_id_b = vec![0; b.points.len()];
    for m in 0..b.points.len() {
        let x = &b.points[m].coords;
        let maybe_point_id_a = if grid_a.is_outside(x) { None } else { grid_a.search(x)? };
        let id = match maybe_point_id_a {
            Some(point_id_a) => point_id_a,
            None => {
                mesh.points.push(Point {
                    id: new_point_id,
                    marker: b.points[m].marker,
                    coords: x.clone(),
                });
                new_point_id += 1;
                new_point_id - 1
            }
        };
        map_old_to_new_point_id_b[m] = id;
    }

    // insert cells of mesh B into new mesh
    let mut new_cell_id = a.cells.len();
    for cell in &b.cells {
        mesh.cells.push(Cell {
            id: new_cell_id,
            attribute: cell.attribute,
            kind: cell.kind,
            points: cell.points.iter().map(|id| map_old_to_new_point_id_b[*id]).collect(),
        });
        new_cell_id += 1;
    }

    // create a map of boundary markers from mesh A
    let mut marked_edges_map = HashMap::new();
    let mut marked_faces_map = HashMap::new();
    a.marked_edges.iter().for_each(|(marker, p1, p2)| {
        let mut edge_key = (*p1, *p2);
        sort2(&mut edge_key);
        marked_edges_map.insert(edge_key, *marker);
    });
    a.marked_faces.iter().for_each(|(marker, p1, p2, p3, p4)| {
        let mut face_key = (*p1, *p2, *p3, *p4);
        sort4(&mut face_key);
        marked_faces_map.insert(face_key, *marker);
    });

    // insert marked edges/faces from mesh B into the boundary maps
    b.marked_edges.iter().for_each(|(marker, p1, p2)| {
        let p1new = map_old_to_new_point_id_b[*p1];
        let p2new = map_old_to_new_point_id_b[*p2];
        let mut edge_key = (p1new, p2new);
        sort2(&mut edge_key);
        if marked_edges_map.get(&edge_key).is_none() {
            marked_edges_map.insert(edge_key, *marker);
        }
    });
    b.marked_faces.iter().for_each(|(marker, p1, p2, p3, p4)| {
        let p1new = map_old_to_new_point_id_b[*p1];
        let p2new = map_old_to_new_point_id_b[*p2];
        let p3new = map_old_to_new_point_id_b[*p3];
        let p4new = map_old_to_new_point_id_b[*p4];
        let mut face_key = (p1new, p2new, p3new, p4new);
        sort4(&mut face_key);
        if marked_faces_map.get(&face_key).is_none() {
            marked_faces_map.insert(face_key, *marker);
        }
    });

    // rebuild marked edges/faces in the new mesh
    mesh.marked_edges.clear();
    mesh.marked_faces.clear();
    if marked_edges_map.len() > 0 {
        let mut edge_keys: Vec<_> = marked_edges_map.keys().cloned().collect();
        edge_keys.sort(); // sort just for deterministic behavior (important for tests)
        for edge_key in &edge_keys {
            let marker = marked_edges_map[edge_key];
            mesh.marked_edges.push((marker, edge_key.0, edge_key.1));
        }
    }
    if marked_faces_map.len() > 0 {
        let mut face_keys: Vec<_> = marked_faces_map.keys().cloned().collect();
        face_keys.sort(); // sort just for deterministic behavior (important for tests)
        for face_key in &face_keys {
            let marker = marked_faces_map[face_key];
            mesh.marked_faces
                .push((marker, face_key.0, face_key.1, face_key.2, face_key.3));
        }
    }

    // done
    Ok(mesh)
}

/// Joins meshes by comparing the coordinates on the boundary
///
/// # Notes
///
/// 1. The meshes must have the same space_ndim.
/// 2. The boundary of mesh[i] is used to find overlapping points with mesh[i+1].
/// 3. The boundary markers of mesh[i] have precedence over those of mesh[i+1].
///
/// **Important:** This function does not guarantee the "mesh compatibility" requirements
/// for finite element analyses.
pub fn join_meshes(meshes: &[&Mesh]) -> Result<Mesh, StrError> {
    if meshes.len() < 2 {
        return Err("meshes.len() must be at least 2");
    }
    let mut new_mesh = join_two_meshes(meshes[0], meshes[1])?;
    for i in 2..meshes.len() {
        let temp_mesh = join_two_meshes(&new_mesh, meshes[i])?;
        new_mesh = temp_mesh;
    }
    Ok(new_mesh)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{join_meshes, join_two_meshes};
    use crate::mesh::Samples;

    #[test]
    fn join_two_meshes_handles_errors() {
        let a = Samples::two_qua4();
        let b = Samples::two_hex8();
        assert_eq!(join_two_meshes(&a, &b).err(), Some("meshes must have the same ndim"));
    }

    #[test]
    fn join_two_meshes_works_2d() {
        //       -300        -400
        //   3-----------2-----------5
        //   |           |           |
        //   |    [0]    |    [1]    |
        //   |    (1)    |    (2)    |
        //   |           |           |
        //   0-----------1-----------4
        //
        //             UNION
        //
        //       -100        -200
        //   3-----------2-----------5
        //   |           |           |
        //   |    [0]    |    [1]    |
        //   |    (1)    |    (2)    |
        //   |           |           |
        //   0-----------1-----------4
        let a = Samples::two_qua4();
        let mut b = Samples::two_qua4();

        // shift B-mesh up
        for m in 0..b.points.len() {
            b.points[m].coords[1] += 1.0;
        }
        // change boundary markers of B-mesh
        b.marked_edges[0].0 = -300;
        b.marked_edges[1].0 = -400;

        //      y
        //      ↑    -300       -400
        // 2.0  7-----------6-----------8
        //      |           |           |
        //      |    [2]    |    [3]    |
        //      |    (1)    |    (2)    |
        //      |    -100   |   -200    |
        // 1.0  3-----------2-----------5
        //      |           |           |
        //      |    [0]    |    [1]    |
        //      |    (1)    |    (2)    |
        //      |           |           |
        // 0.0  0-----------1-----------4  → x
        //     0.0         1.0         2.0
        let mesh = join_two_meshes(&a, &b).unwrap();
        // println!("{}", mesh);
        mesh.check_ids_and_kind().unwrap();
        mesh.check_jacobian().unwrap();
        mesh.check_overlapping_points(0.01).unwrap();
        assert_eq!(mesh.points[0].coords, &[0.0, 0.0]);
        assert_eq!(mesh.points[1].coords, &[1.0, 0.0]);
        assert_eq!(mesh.points[2].coords, &[1.0, 1.0]);
        assert_eq!(mesh.points[3].coords, &[0.0, 1.0]);
        assert_eq!(mesh.points[4].coords, &[2.0, 0.0]);
        assert_eq!(mesh.points[5].coords, &[2.0, 1.0]);
        assert_eq!(mesh.points[6].coords, &[1.0, 2.0]);
        assert_eq!(mesh.points[7].coords, &[0.0, 2.0]);
        assert_eq!(mesh.points[8].coords, &[2.0, 2.0]);
        assert_eq!(mesh.cells[0].points, &[0, 1, 2, 3]);
        assert_eq!(mesh.cells[1].points, &[1, 4, 5, 2]);
        assert_eq!(mesh.cells[2].points, &[3, 2, 6, 7]);
        assert_eq!(mesh.cells[3].points, &[2, 5, 8, 6]);
        assert_eq!(mesh.cells[0].attribute, 1);
        assert_eq!(mesh.cells[1].attribute, 2);
        assert_eq!(mesh.cells[2].attribute, 1);
        assert_eq!(mesh.cells[3].attribute, 2);
        assert_eq!(
            mesh.marked_edges,
            vec![(-100, 2, 3), (-200, 2, 5), (-300, 6, 7), (-400, 6, 8)]
        );
    }

    #[test]
    fn join_two_meshes_works_3d() {
        //       8-------------11                   8-------------11
        //      /.             /|                  /.             /|
        // {-5}/ .        {-5}/ |            {-5} / .        {-5}/ |
        //    /  .   {-9}    /  |                /  .   {-9}    /  |
        //   /   .          /   |{123}          /   .          /   |{123}
        //  9-------------10    |              9-------------10    |
        //  |    .         |    |              |    .         |    |
        //  |    4---------|----7              |    4---------|----7
        //  |   /. [1]     |   /|              |   /. [1]     |   /|
        //  |  / . (2)     |  / |    UNION     |  / . (2)     |  / |
        //  | /  .         | /  |              | /  .         | /  |
        //  |/   .         |/   |{-4}          |/   .         |/   |{-4}
        //  5--------------6    |              5--------------6    |
        //  |    .         |{-8}|              |    .         |{-8}|
        //  |    0---------|----3              |    0---------|----3
        //  |   /  [0]     |   /               |   /  [0]     |   /
        //  |  /   (1)     |  /                |  /   (1)     |  /
        //  | /            | /                 | /            | /
        //  |/             |/                  |/             |/
        //  1--------------2                   1--------------2   1.0
        let a = Samples::two_hex8();
        let mut b = Samples::two_hex8();

        // shift mesh B along y (to the right)
        for m in 0..b.points.len() {
            b.points[m].coords[1] += 1.0;
        }

        //       8-------------11-------------17
        //      /.             /.             /|
        // {-5}/ .        {-5}/ .        {-5}/ |
        //    /  .   {-9}    /  .    {-9}   /  |
        //   /   .          /   .{123}     /   |{123}
        //  9=============10=============16    |
        //  |    .         |    .         |    |
        //  |    4---------|----7---------|---15
        //  |   /. [1]     |   /. [3]     |   /|
        //  |  / . (2)     |  / . (2)     |  / |
        //  | /  .         | /  .         | /  |
        //  |/   .         |/   .{-4}     |/   |{-4}
        //  5--------------6-------------14    |
        //  |    .         |{-8}.         |{-8}|
        //  |    0---------|----3---------|---13
        //  |   /  [0]     |   /  [2]     |   /
        //  |  /   (1)     |  /   (1)     |  /
        //  | /            | /            | /
        //  |/             |/             |/
        //  1--------------2-------------12
        let mesh = join_two_meshes(&a, &b).unwrap();
        println!("{}", mesh);
        assert_eq!(mesh.ndim, 3);
        assert_eq!(mesh.points.len(), 18);
        assert_eq!(mesh.cells.len(), 4);
        mesh.check_ids_and_kind().unwrap();
        mesh.check_jacobian().unwrap();
        mesh.check_overlapping_points(0.01).unwrap();

        let sample = Samples::four_hex8();
        for i in 0..mesh.cells.len() {
            assert_eq!(mesh.cells[i].id, sample.cells[i].id);
            assert_eq!(mesh.cells[i].points, sample.cells[i].points);
        }
        for i in 0..mesh.points.len() {
            assert_eq!(mesh.points[i].id, sample.points[i].id);
            assert_eq!(mesh.points[i].coords, sample.points[i].coords);
        }

        assert_eq!(
            mesh.marked_edges,
            vec![
                (-4, 3, 7),
                (123, 7, 11),
                (-5, 8, 9),
                (-5, 10, 11),
                (-4, 13, 15),
                (123, 15, 17),
                (-5, 16, 17)
            ]
        );
        assert_eq!(
            mesh.marked_faces,
            vec![
                (-8, 2, 3, 6, 7),
                (-9, 8, 9, 10, 11),
                (-9, 10, 11, 16, 17),
                (-8, 12, 13, 14, 15)
            ]
        );
    }

    #[test]
    fn join_meshes_works() {
        //       -500        -600
        //   3-----------2-----------5
        //   |(-4)       |(-3)       |(-6)
        //   |    [0]    |    [1]    |
        //   |    (1)    |    (2)    |
        //   |(-1)       |(-2)       |(-5)
        //   0-----------1-----------4  → x
        //
        //             UNION
        //
        //       -300        -400
        //   3-----------2-----------5
        //   |(-4)       |(-3)       |(-6)
        //   |    [0]    |    [1]    |
        //   |    (1)    |    (2)    |
        //   |(-1)       |(-2)       |(-5)
        //   0-----------1-----------4
        //
        //             UNION
        //
        //       -100        -200
        //   3-----------2-----------5
        //   |(-4)       |(-3)       |(-6)
        //   |    [0]    |    [1]    |
        //   |    (1)    |    (2)    |
        //   |(-1)       |(-2)       |(-5)
        //   0-----------1-----------4
        let a = Samples::two_qua4();
        let mut b = Samples::two_qua4();
        let mut c = Samples::two_qua4();
        for m in 0..b.points.len() {
            b.points[m].coords[1] += 1.0;
            c.points[m].coords[1] += 2.0;
        }
        b.marked_edges[0].0 = -300;
        b.marked_edges[1].0 = -400;
        b.marked_edges.push((-123, 3, 2)); // will be overwritten by the first mesh
        b.marked_edges.push((-456, 2, 5)); // will be overwritten by the first mesh
        c.marked_edges[0].0 = -500;
        c.marked_edges[1].0 = -600;

        //        -500       -600
        //   10----------9----------11
        //   |(-4)       |(-3)       |(-6)
        //   |    [4]    |    [5]    |
        //   |    (1)    |    (2)    |
        //   |    -300   |   -400    |
        //   7-----------6-----------8
        //   |(-4)       |(-3)       |(-6)
        //   |    [2]    |    [3]    |
        //   |    (1)    |    (2)    |
        //   |    -100   |   -200    |
        //   3-----------2-----------5
        //   |(-4)       |(-3)       |(-6)
        //   |    [0]    |    [1]    |
        //   |    (1)    |    (2)    |
        //   |(-1)       |(-2)       |(-5)
        //   0-----------1-----------4
        let mesh = join_meshes(&[&a, &b, &c]).unwrap();
        mesh.check_overlapping_points(0.01).unwrap();
        assert_eq!(
            format!("{}", mesh),
            "# header\n\
             # ndim npoint ncell nmarked_edge nmarked_face\n\
             2 12 6 6 0\n\
             \n\
             # points\n\
             # id marker x y {z}\n\
             0 -1 0.0 0.0\n\
             1 -2 1.0 0.0\n\
             2 -3 1.0 1.0\n\
             3 -4 0.0 1.0\n\
             4 -5 2.0 0.0\n\
             5 -6 2.0 1.0\n\
             6 -3 1.0 2.0\n\
             7 -4 0.0 2.0\n\
             8 -6 2.0 2.0\n\
             9 -3 1.0 3.0\n\
             10 -4 0.0 3.0\n\
             11 -6 2.0 3.0\n\
             \n\
             # cells\n\
             # id attribute kind points\n\
             0 1 qua4 0 1 2 3\n\
             1 2 qua4 1 4 5 2\n\
             2 1 qua4 3 2 6 7\n\
             3 2 qua4 2 5 8 6\n\
             4 1 qua4 7 6 9 10\n\
             5 2 qua4 6 8 11 9\n\
             \n\
             # marked edges\n\
             # marker p1 p2\n\
             -100 2 3\n\
             -200 2 5\n\
             -300 6 7\n\
             -400 6 8\n\
             -500 9 10\n\
             -600 9 11\n"
        );
    }
}
