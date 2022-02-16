use gemlab::{mesh::Mesh, StrError};
use russell_chk::assert_vec_approx_eq;
use russell_lab::Vector;

#[test]
fn five_tets_within_cube() -> Result<(), StrError> {
    // read mesh
    let mesh = Mesh::from_text_file("./data/meshes/five_tets_within_cube.msh")?;
    let npoint = mesh.points.len();
    assert_eq!(npoint, 8);
    assert_eq!(mesh.cells.len(), 5);
    assert_eq!(mesh.boundary_points.len(), 8);
    assert_eq!(mesh.boundary_edges.len(), 18);
    assert_eq!(mesh.boundary_faces.len(), 12);

    // cells
    assert_eq!(mesh.cells[0].points, &[1, 2, 0, 5]);
    assert_eq!(mesh.cells[1].points, &[3, 0, 2, 7]);
    assert_eq!(mesh.cells[2].points, &[4, 7, 5, 0]);
    assert_eq!(mesh.cells[3].points, &[6, 5, 7, 2]);
    assert_eq!(mesh.cells[4].points, &[0, 2, 7, 5]);

    // x-min
    let face = mesh.boundary_faces.get(&(0, 3, 7, npoint)).unwrap();
    assert_eq!(face.points, &[3, 0, 7]);
    let face = mesh.boundary_faces.get(&(0, 4, 7, npoint)).unwrap();
    assert_eq!(face.points, &[4, 7, 0]);
    // x-max
    let face = mesh.boundary_faces.get(&(1, 2, 5, npoint)).unwrap();
    assert_eq!(face.points, &[1, 2, 5]);
    let face = mesh.boundary_faces.get(&(2, 5, 6, npoint)).unwrap();
    assert_eq!(face.points, &[6, 5, 2]);
    // y-min
    let face = mesh.boundary_faces.get(&(0, 1, 5, npoint)).unwrap();
    assert_eq!(face.points, &[1, 5, 0]);
    let face = mesh.boundary_faces.get(&(0, 4, 5, npoint)).unwrap();
    assert_eq!(face.points, &[4, 0, 5]);
    // y-max
    let face = mesh.boundary_faces.get(&(2, 3, 7, npoint)).unwrap();
    assert_eq!(face.points, &[3, 7, 2]);
    let face = mesh.boundary_faces.get(&(2, 6, 7, npoint)).unwrap();
    assert_eq!(face.points, &[6, 2, 7]);
    // z-min
    let face = mesh.boundary_faces.get(&(0, 1, 2, npoint)).unwrap();
    assert_eq!(face.points, &[1, 0, 2]);
    let face = mesh.boundary_faces.get(&(0, 2, 3, npoint)).unwrap();
    assert_eq!(face.points, &[3, 2, 0]);
    // z-max
    let face = mesh.boundary_faces.get(&(4, 5, 7, npoint)).unwrap();
    assert_eq!(face.points, &[4, 5, 7]);
    let face = mesh.boundary_faces.get(&(5, 6, 7, npoint)).unwrap();
    assert_eq!(face.points, &[6, 7, 5]);
    // internal
    assert!(mesh.boundary_faces.get(&(0, 2, 7, npoint)).is_none());

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
    let mut normal = Vector::new(mesh.space_ndim);
    let ksi = &[0.0, 0.0, 0.0];
    for (face_keys, solution) in &face_keys_and_solutions {
        for face_key in face_keys {
            let mut face_shape = mesh.alloc_shape_boundary_face(face_key)?;
            face_shape.calc_boundary_normal(&mut normal, ksi)?;
            assert_vec_approx_eq!(normal.as_data(), solution, 1e-15);
        }
    }
    Ok(())
}
