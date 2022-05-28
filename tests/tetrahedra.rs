use gemlab::mesh::{allocate_shapes, At, Boundary, Find, Mesh, NormalVector};
use gemlab::StrError;
use russell_chk::assert_vec_approx_eq;
use std::collections::HashSet;

fn check<T>(found: &HashSet<T>, correct: &[T])
where
    T: Copy + Ord + std::fmt::Debug,
{
    let mut ids: Vec<T> = found.iter().copied().collect();
    ids.sort();
    assert_eq!(ids, correct);
}

#[test]
fn five_tets_within_cube() -> Result<(), StrError> {
    // read mesh
    let mesh = Mesh::from_text_file("./data/meshes/five_tets_within_cube.msh")?;
    let npoint = mesh.points.len();
    assert_eq!(npoint, 8);
    assert_eq!(mesh.cells.len(), 5);

    // cells
    assert_eq!(mesh.cells[0].points, &[1, 2, 0, 5]);
    assert_eq!(mesh.cells[1].points, &[3, 0, 2, 7]);
    assert_eq!(mesh.cells[2].points, &[4, 7, 5, 0]);
    assert_eq!(mesh.cells[3].points, &[6, 5, 7, 2]);
    assert_eq!(mesh.cells[4].points, &[0, 2, 7, 5]);

    // boundary
    let shapes = allocate_shapes(&mesh)?;
    let boundary = Boundary::new(&mesh, &shapes)?;

    // x-min
    let face = boundary.faces.get(&(0, 3, 7, npoint)).unwrap();
    assert_eq!(face.points, &[3, 0, 7]);
    let face = boundary.faces.get(&(0, 4, 7, npoint)).unwrap();
    assert_eq!(face.points, &[4, 7, 0]);
    // x-max
    let face = boundary.faces.get(&(1, 2, 5, npoint)).unwrap();
    assert_eq!(face.points, &[1, 2, 5]);
    let face = boundary.faces.get(&(2, 5, 6, npoint)).unwrap();
    assert_eq!(face.points, &[6, 5, 2]);
    // y-min
    let face = boundary.faces.get(&(0, 1, 5, npoint)).unwrap();
    assert_eq!(face.points, &[1, 5, 0]);
    let face = boundary.faces.get(&(0, 4, 5, npoint)).unwrap();
    assert_eq!(face.points, &[4, 0, 5]);
    // y-max
    let face = boundary.faces.get(&(2, 3, 7, npoint)).unwrap();
    assert_eq!(face.points, &[3, 7, 2]);
    let face = boundary.faces.get(&(2, 6, 7, npoint)).unwrap();
    assert_eq!(face.points, &[6, 2, 7]);
    // z-min
    let face = boundary.faces.get(&(0, 1, 2, npoint)).unwrap();
    assert_eq!(face.points, &[1, 0, 2]);
    let face = boundary.faces.get(&(0, 2, 3, npoint)).unwrap();
    assert_eq!(face.points, &[3, 2, 0]);
    // z-max
    let face = boundary.faces.get(&(4, 5, 7, npoint)).unwrap();
    assert_eq!(face.points, &[4, 5, 7]);
    let face = boundary.faces.get(&(5, 6, 7, npoint)).unwrap();
    assert_eq!(face.points, &[6, 7, 5]);
    // internal
    assert!(boundary.faces.get(&(0, 2, 7, npoint)).is_none());

    // the norm of the normal vector should be equal to face_area / 0.5
    // where 0.5 corresponds to the face_area in the reference system
    let l = 2.0 / 0.5; // norm of normal vector

    // face keys and correct normal vectors (solutions)
    let face_keys_and_solutions = [
        // x-min
        (vec![(0, 3, 7, npoint), (0, 4, 7, npoint)], [-l, 0.0, 0.0]),
        // x-max
        (vec![(1, 2, 5, npoint), (2, 5, 6, npoint)], [l, 0.0, 0.0]),
        // y-min
        (vec![(0, 1, 5, npoint), (0, 4, 5, npoint)], [0.0, -l, 0.0]),
        // y-max
        (vec![(2, 3, 7, npoint), (2, 6, 7, npoint)], [0.0, l, 0.0]),
        // z-min
        (vec![(0, 1, 2, npoint), (0, 2, 3, npoint)], [0.0, 0.0, -l]),
        // z-max
        (vec![(4, 5, 7, npoint), (5, 6, 7, npoint)], [0.0, 0.0, l]),
    ];

    // check if the normal vectors at boundary are outward
    let ksi = &[0.0, 0.0, 0.0];
    for (face_keys, solution) in &face_keys_and_solutions {
        for face_key in face_keys {
            let mut normal = NormalVector::at_face(&mesh, &boundary, *face_key)?;
            normal.evaluate(ksi)?;
            assert_vec_approx_eq!(normal.value.as_data(), solution, 1e-15);
        }
    }

    // find
    let find = Find::new(&mesh, &boundary)?;

    // find points
    let points = find.points(At::X(0.0))?;
    check(&points, &[0, 3, 4, 7]);
    let points = find.points(At::Y(2.0))?;
    check(&points, &[2, 3, 6, 7]);

    // find edges
    let edges = find.edges(At::Z(0.0))?;
    check(&edges, &[(0, 1), (0, 2), (0, 3), (1, 2), (2, 3)]);
    let edges = find.edges(At::Y(2.0))?;
    check(&edges, &[(2, 3), (2, 6), (2, 7), (3, 7), (6, 7)]);

    // find faces
    let faces = find.faces(At::X(0.0))?;
    check(&faces, &[(0, 3, 7, npoint), (0, 4, 7, npoint)]);
    let faces = find.faces(At::X(2.0))?;
    check(&faces, &[(1, 2, 5, npoint), (2, 5, 6, npoint)]);
    let faces = find.faces(At::Y(0.0))?;
    check(&faces, &[(0, 1, 5, npoint), (0, 4, 5, npoint)]);
    let faces = find.faces(At::Y(2.0))?;
    check(&faces, &[(2, 3, 7, npoint), (2, 6, 7, npoint)]);
    let faces = find.faces(At::Z(0.0))?;
    check(&faces, &[(0, 1, 2, npoint), (0, 2, 3, npoint)]);
    let faces = find.faces(At::Z(2.0))?;
    check(&faces, &[(4, 5, 7, npoint), (5, 6, 7, npoint)]);
    Ok(())
}
