/// Allocates a Shape instance for a given Cell and initializes the associated point coordinates
pub fn alloc_shape_cell(&self, cell_id: CellId) -> Result<Shape, StrError> {
    if !self.derived_props_computed {
        return Err("compute_derived_props must be called first");
    }
    if cell_id >= self.cells.len() {
        return Err("cell_id is out-of-bounds");
    }
    let cell = &self.cells[cell_id];
    let nnode = cell.points.len();
    let mut shape = Shape::new(self.space_ndim, cell.geo_ndim, nnode)?;
    for m in 0..nnode {
        let point_id = cell.points[m];
        for j in 0..self.space_ndim {
            shape.set_node(point_id, m, j, self.points[point_id].coords[j])?;
        }
    }
    Ok(shape)
}

/// Allocates a Shape instance for a given boundary Edge and initializes the associated point coordinates
pub fn alloc_shape_boundary_edge(&self, edge_key: &EdgeKey) -> Result<Shape, StrError> {
    if !self.derived_props_computed {
        return Err("compute_derived_props must be called first");
    }
    match self.boundary_edges.get(edge_key) {
        None => Err("edge_key is invalid"),
        Some(edge) => {
            let nnode = edge.points.len();
            let mut shape = Shape::new(self.space_ndim, 1, nnode)?;
            for m in 0..nnode {
                let point_id = edge.points[m];
                for j in 0..self.space_ndim {
                    shape.set_node(point_id, m, j, self.points[point_id].coords[j])?;
                }
            }
            Ok(shape)
        }
    }
}

/// Allocates a Shape instance for a given boundary Face and initializes the associated point coordinates
pub fn alloc_shape_boundary_face(&self, face_key: &FaceKey) -> Result<Shape, StrError> {
    if !self.derived_props_computed {
        return Err("compute_derived_props must be called first");
    }
    match self.boundary_faces.get(face_key) {
        None => Err("face_key is invalid"),
        Some(face) => {
            let nnode = face.points.len();
            let mut shape = Shape::new(self.space_ndim, 2, nnode)?;
            for m in 0..nnode {
                let point_id = face.points[m];
                for j in 0..self.space_ndim {
                    shape.set_node(point_id, m, j, self.points[point_id].coords[j])?;
                }
            }
            Ok(shape)
        }
    }
}

#[test]
fn alloc_shape_fns_are_correct_2d() -> Result<(), StrError> {
    //
    //  3--------2--------5
    //  |        |        |
    //  |        |        |
    //  |        |        |
    //  0--------1--------4
    //
    let mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;
    assert_eq!(mesh.alloc_shape_cell(2).err(), Some("cell_id is out-of-bounds"));
    assert_eq!(
        mesh.alloc_shape_boundary_edge(&(1, 0)).err(),
        Some("edge_key is invalid")
    );
    assert_eq!(
        mesh.alloc_shape_boundary_face(&(0, 1, 2, 3)).err(),
        Some("face_key is invalid")
    );

    let shape = mesh.alloc_shape_cell(0)?;
    assert_eq!(shape.kind, GeoKind::Qua4);
    assert_eq!(shape.node_to_point, &[0, 1, 2, 3]);
    assert_eq!(
        format!("{}", shape.coords_transp),
        "┌         ┐\n\
             │ 0 1 1 0 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
    );

    let shape = mesh.alloc_shape_cell(1)?;
    assert_eq!(shape.kind, GeoKind::Qua4);
    assert_eq!(shape.node_to_point, &[1, 4, 5, 2]);
    assert_eq!(
        format!("{}", shape.coords_transp),
        "┌         ┐\n\
             │ 1 2 2 1 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
    );

    let shape = mesh.alloc_shape_boundary_edge(&(0, 1))?;
    assert_eq!(shape.kind, GeoKind::Lin2);
    assert_eq!(shape.node_to_point, &[1, 0]);
    assert_eq!(
        format!("{}", shape.coords_transp),
        "┌     ┐\n\
             │ 1 0 │\n\
             │ 0 0 │\n\
             └     ┘"
    );

    let shape = mesh.alloc_shape_boundary_edge(&(2, 5))?;
    assert_eq!(shape.kind, GeoKind::Lin2);
    assert_eq!(shape.node_to_point, &[2, 5]);
    assert_eq!(
        format!("{}", shape.coords_transp),
        "┌     ┐\n\
             │ 1 2 │\n\
             │ 1 1 │\n\
             └     ┘"
    );
    Ok(())
}

#[test]
fn alloc_shape_fns_are_correct_3d() -> Result<(), StrError> {
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
    assert_eq!(mesh.alloc_shape_cell(2).err(), Some("cell_id is out-of-bounds"));
    assert_eq!(
        mesh.alloc_shape_boundary_edge(&(1, 0)).err(),
        Some("edge_key is invalid")
    );
    assert_eq!(
        mesh.alloc_shape_boundary_face(&(3, 2, 1, 0)).err(),
        Some("face_key is invalid")
    );

    let shape = mesh.alloc_shape_cell(0)?;
    assert_eq!(shape.kind, GeoKind::Hex8);
    assert_eq!(shape.node_to_point, &[0, 1, 2, 3, 4, 5, 6, 7]);
    assert_eq!(
        format!("{}", shape.coords_transp),
        "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 0 0 0 0 1 1 1 1 │\n\
             └                 ┘"
    );

    let shape = mesh.alloc_shape_cell(1)?;
    assert_eq!(shape.kind, GeoKind::Hex8);
    assert_eq!(shape.node_to_point, &[4, 5, 6, 7, 8, 9, 10, 11]);
    assert_eq!(
        format!("{}", shape.coords_transp),
        "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 1 1 1 1 2 2 2 2 │\n\
             └                 ┘"
    );

    let shape = mesh.alloc_shape_boundary_edge(&(0, 1))?;
    assert_eq!(shape.kind, GeoKind::Lin2);
    assert_eq!(shape.node_to_point, &[0, 1]);
    assert_eq!(
        format!("{}", shape.coords_transp),
        "┌     ┐\n\
             │ 0 1 │\n\
             │ 0 0 │\n\
             │ 0 0 │\n\
             └     ┘"
    );

    let shape = mesh.alloc_shape_boundary_edge(&(8, 11))?;
    assert_eq!(shape.kind, GeoKind::Lin2);
    assert_eq!(shape.node_to_point, &[11, 8]);
    assert_eq!(
        format!("{}", shape.coords_transp),
        "┌     ┐\n\
             │ 0 0 │\n\
             │ 1 0 │\n\
             │ 2 2 │\n\
             └     ┘"
    );

    let shape = mesh.alloc_shape_boundary_face(&(0, 1, 2, 3))?;
    assert_eq!(shape.kind, GeoKind::Qua4);
    assert_eq!(shape.node_to_point, &[0, 3, 2, 1]);
    assert_eq!(
        format!("{}", shape.coords_transp),
        "┌         ┐\n\
             │ 0 0 1 1 │\n\
             │ 0 1 1 0 │\n\
             │ 0 0 0 0 │\n\
             └         ┘"
    );

    let shape = mesh.alloc_shape_boundary_face(&(8, 9, 10, 11))?;
    assert_eq!(shape.kind, GeoKind::Qua4);
    assert_eq!(shape.node_to_point, &[8, 9, 10, 11]);
    assert_eq!(
        format!("{}", shape.coords_transp),
        "┌         ┐\n\
             │ 0 1 1 0 │\n\
             │ 0 0 1 1 │\n\
             │ 2 2 2 2 │\n\
             └         ┘"
    );
    Ok(())
}
