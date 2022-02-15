use gemlab::{mesh::Mesh, StrError};
use russell_chk::assert_vec_approx_eq;
use russell_lab::Vector;

#[test]
fn five_tets_within_cube() -> Result<(), StrError> {
    // read mesh
    let mesh = Mesh::from_text_file("./data/meshes/five_tets_within_cube.msh")?;
    println!("{}", mesh);
    let npoint = mesh.points.len();
    assert_eq!(npoint, 8);
    let face = mesh.boundary_faces.get(&(0, 1, 2, npoint)).unwrap();
    assert_eq!(face.points, &[1, 0, 2]);
    let face = mesh.boundary_faces.get(&(0, 2, 3, npoint)).unwrap();
    assert_eq!(face.points, &[3, 2, 0]);

    // the norm of the normal vector should be equal to face_area / 0.5
    // where 4.0 corresponds to the face_area in the reference system
    let l = 1.0; // norm of normal vector

    // face keys and correct normal vectors (solutions)
    let face_keys_and_solutions = [
        // bottom
        (vec![(0, 1, 2, npoint), (0, 2, 3, npoint)], [0.0, 0.0, -l]),
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
