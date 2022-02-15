use gemlab::{mesh::Mesh, StrError};

#[test]
fn column_distorted_tris_quads() -> Result<(), StrError> {
    // read mesh
    let mesh = Mesh::from_text_file("./data/meshes/column_distorted_tris_quads.msh")?;
    println!("{}", mesh);
    Ok(())
}

#[test]
fn rectangle_tris_quads() -> Result<(), StrError> {
    // read mesh
    let mesh = Mesh::from_text_file("./data/meshes/rectangle_tris_quads.msh")?;
    println!("{}", mesh);
    Ok(())
}
