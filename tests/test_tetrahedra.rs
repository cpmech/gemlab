use gemlab::mesh::{check_face_normals, At, Extract, Features, Find, Mesh};
use gemlab::StrError;
use std::collections::HashMap;

fn any(_: &Vec<f64>) -> bool {
    true
}

#[test]
fn test_tetrahedra() -> Result<(), StrError> {
    // read mesh
    let mesh = Mesh::from_text_file("./data/meshes/five_tets_within_cube.msh")?;
    let features = Features::new(&mesh, Extract::Boundary);
    let npoint = mesh.points.len();
    assert_eq!(npoint, 8);
    assert_eq!(mesh.cells.len(), 5);

    // cells
    assert_eq!(mesh.cells[0].points, &[1, 2, 0, 5]);
    assert_eq!(mesh.cells[1].points, &[3, 0, 2, 7]);
    assert_eq!(mesh.cells[2].points, &[4, 7, 5, 0]);
    assert_eq!(mesh.cells[3].points, &[6, 5, 7, 2]);
    assert_eq!(mesh.cells[4].points, &[0, 2, 7, 5]);

    // x-min
    let face = features.faces.get(&(0, 3, 7, npoint)).unwrap();
    assert_eq!(face.points, &[3, 0, 7]);
    let face = features.faces.get(&(0, 4, 7, npoint)).unwrap();
    assert_eq!(face.points, &[4, 7, 0]);
    // x-max
    let face = features.faces.get(&(1, 2, 5, npoint)).unwrap();
    assert_eq!(face.points, &[1, 2, 5]);
    let face = features.faces.get(&(2, 5, 6, npoint)).unwrap();
    assert_eq!(face.points, &[6, 5, 2]);
    // y-min
    let face = features.faces.get(&(0, 1, 5, npoint)).unwrap();
    assert_eq!(face.points, &[1, 5, 0]);
    let face = features.faces.get(&(0, 4, 5, npoint)).unwrap();
    assert_eq!(face.points, &[4, 0, 5]);
    // y-max
    let face = features.faces.get(&(2, 3, 7, npoint)).unwrap();
    assert_eq!(face.points, &[3, 7, 2]);
    let face = features.faces.get(&(2, 6, 7, npoint)).unwrap();
    assert_eq!(face.points, &[6, 2, 7]);
    // z-min
    let face = features.faces.get(&(0, 1, 2, npoint)).unwrap();
    assert_eq!(face.points, &[1, 0, 2]);
    let face = features.faces.get(&(0, 2, 3, npoint)).unwrap();
    assert_eq!(face.points, &[3, 2, 0]);
    // z-max
    let face = features.faces.get(&(4, 5, 7, npoint)).unwrap();
    assert_eq!(face.points, &[4, 5, 7]);
    let face = features.faces.get(&(5, 6, 7, npoint)).unwrap();
    assert_eq!(face.points, &[6, 7, 5]);
    // internal
    assert!(features.faces.get(&(0, 2, 7, npoint)).is_none());

    // the norm of the normal vector should be equal to face_area / 0.5
    // where 0.5 corresponds to the face_area in the reference system
    let l = 2.0 / 0.5; // norm of normal vector

    // face keys and correct normal vectors (solutions)
    let solutions = HashMap::from([
        ((0, 3, 7, npoint), (l, [-1.0, 0.0, 0.0])), // x-min
        ((0, 4, 7, npoint), (l, [-1.0, 0.0, 0.0])), // x-min
        ((1, 2, 5, npoint), (l, [1.0, 0.0, 0.0])),  // x-max
        ((2, 5, 6, npoint), (l, [1.0, 0.0, 0.0])),  // x-max
        ((0, 1, 5, npoint), (l, [0.0, -1.0, 0.0])), // y-min
        ((0, 4, 5, npoint), (l, [0.0, -1.0, 0.0])), // y-min
        ((2, 3, 7, npoint), (l, [0.0, 1.0, 0.0])),  // y-max
        ((2, 6, 7, npoint), (l, [0.0, 1.0, 0.0])),  // y-max
        ((0, 1, 2, npoint), (l, [0.0, 0.0, -1.0])), // z-min
        ((0, 2, 3, npoint), (l, [0.0, 0.0, -1.0])), // z-min
        ((4, 5, 7, npoint), (l, [0.0, 0.0, 1.0])),  // z-max
        ((5, 6, 7, npoint), (l, [0.0, 0.0, 1.0])),  // z-max
    ]);

    // check if the normal vectors at boundary are outward
    check_face_normals(&mesh, &features.faces, &solutions, 1e-15).expect("ok");

    // find points
    let find = Find::new(&mesh, None);
    let points = find.point_ids(At::X(0.0), any)?;
    assert_eq!(&points, &[0, 3, 4, 7]);
    let points = find.point_ids(At::Y(2.0), any)?;
    assert_eq!(&points, &[2, 3, 6, 7]);

    // find edges
    let edges = find.edge_keys(At::Z(0.0), any)?;
    assert_eq!(&edges, &[(0, 1), (0, 2), (0, 3), (1, 2), (2, 3)]);
    let edges = find.edge_keys(At::Y(2.0), any)?;
    assert_eq!(&edges, &[(2, 3), (2, 6), (2, 7), (3, 7), (6, 7)]);

    // find faces
    let faces = find.face_keys(At::X(0.0), any)?;
    assert_eq!(&faces, &[(0, 3, 7, npoint), (0, 4, 7, npoint)]);
    let faces = find.face_keys(At::X(2.0), any)?;
    assert_eq!(&faces, &[(1, 2, 5, npoint), (2, 5, 6, npoint)]);
    let faces = find.face_keys(At::Y(0.0), any)?;
    assert_eq!(&faces, &[(0, 1, 5, npoint), (0, 4, 5, npoint)]);
    let faces = find.face_keys(At::Y(2.0), any)?;
    assert_eq!(&faces, &[(2, 3, 7, npoint), (2, 6, 7, npoint)]);
    let faces = find.face_keys(At::Z(0.0), any)?;
    assert_eq!(&faces, &[(0, 1, 2, npoint), (0, 2, 3, npoint)]);
    let faces = find.face_keys(At::Z(2.0), any)?;
    assert_eq!(&faces, &[(4, 5, 7, npoint), (5, 6, 7, npoint)]);
    Ok(())
}
