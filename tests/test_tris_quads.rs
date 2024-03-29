use gemlab::integ;
use gemlab::mesh::{At, Features, Mesh};
use gemlab::shapes::{GeoKind, Scratchpad};
use gemlab::util::any_x;
use gemlab::StrError;
use russell_lab::approx_eq;
use russell_lab::math::SQRT_2;
use std::collections::HashMap;

#[test]
fn test_column_distorted_tris_quads() -> Result<(), StrError> {
    // read mesh
    let mesh = Mesh::from_text_file("./data/meshes/column_distorted_tris_quads.msh")?;
    let features = Features::new(&mesh, false);

    // check sizes
    assert_eq!(mesh.points.len(), 13);
    assert_eq!(mesh.cells.len(), 7);

    // check cells
    assert_eq!(mesh.cells[0].points, &[0, 7, 8, 1]);
    assert_eq!(mesh.cells[1].points, &[1, 8, 9, 2]);
    assert_eq!(mesh.cells[2].points, &[2, 9, 10, 3]);
    assert_eq!(mesh.cells[3].points, &[3, 10, 4]);
    assert_eq!(mesh.cells[4].points, &[4, 10, 11]);
    assert_eq!(mesh.cells[5].points, &[4, 11, 5]);
    assert_eq!(mesh.cells[6].points, &[5, 11, 12, 6]);

    // check edges
    let edge = features.edges.get(&(0, 1)).unwrap();
    assert_eq!(edge.points, &[0, 1]);
    let edge = features.edges.get(&(9, 10)).unwrap();
    assert_eq!(edge.points, &[10, 9]);
    let edge = features.edges.get(&(3, 4)).unwrap();
    assert_eq!(edge.points, &[3, 4]);

    // the magnitude of the normal vector should be equal to edge_length / 2.0
    // for both tris or quas where 2.0 corresponds to the edge_length in the reference system
    // Note that the edge is mapped to [-1,+1] in both Lin or Qua

    // edge keys and correct normal vectors (solutions)
    let solutions = HashMap::from([
        // left
        ((0, 1), (0.25, [-1.0, 0.0])),
        ((1, 2), (0.25, [-1.0, 0.0])),
        ((2, 3), (0.50, [-1.0, 0.0])),
        ((3, 4), (0.25, [-1.0, 0.0])),
        ((4, 5), (0.25, [-1.0, 0.0])),
        ((5, 6), (0.05, [-1.0, 0.0])),
        // right
        ((7, 8), (0.1, [1.0, 0.0])),
        ((8, 9), (0.4, [1.0, 0.0])),
        ((9, 10), (0.4, [1.0, 0.0])),
        ((10, 11), (0.6, [1.0, 0.0])),
        ((11, 12), (0.05, [1.0, 0.0])),
        // bottom
        ((0, 7), (0.5, [0.0, -1.0])),
        // top
        ((6, 12), (0.5, [0.0, 1.0])),
    ]);

    // check if the normal vectors at boundary are outward
    mesh.check_2d_edge_normals(&features.edges, &solutions, 1e-15)
        .expect("ok");

    // search points
    let feat = Features::new(&mesh, false);
    let points = feat.search_point_ids(At::X(0.0), any_x)?;
    assert_eq!(&points, &[0, 1, 2, 3, 4, 5, 6]);
    let points = feat.search_point_ids(At::X(1.0), any_x)?;
    assert_eq!(&points, &[7, 8, 9, 10, 11, 12]);

    // find edges
    let edges = feat.search_edge_keys(At::X(0.0), any_x)?;
    assert_eq!(&edges, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]);
    let edges = feat.search_edge_keys(At::X(1.0), any_x)?;
    assert_eq!(&edges, &[(7, 8), (8, 9), (9, 10), (10, 11), (11, 12)]);

    // find faces
    assert_eq!(
        feat.search_face_keys(At::X(0.0), any_x).err(),
        Some("cannot find face keys in 2D")
    );
    Ok(())
}

#[test]
fn test_rectangle_tris_quads() -> Result<(), StrError> {
    // read mesh
    let mesh = Mesh::from_text_file("./data/meshes/rectangle_tris_quads.msh")?;
    let features = Features::new(&mesh, false);

    // the magnitude of the normal vector should be equal to edge_length / 2.0
    // for both tris or quas where 2.0 corresponds to the edge_length in the reference system
    // Note that the edge is mapped to [-1,+1] in both Lin or Qua

    // edge keys and correct normal vectors (solutions)
    let solutions = HashMap::from([
        ((0, 3), (0.5, [-1.0, 0.0])),  // left
        ((3, 7), (0.5, [-1.0, 0.0])),  // left
        ((7, 10), (0.5, [-1.0, 0.0])), // left
        ((2, 6), (0.5, [1.0, 0.0])),   // right
        ((6, 9), (0.5, [1.0, 0.0])),   // right
        ((9, 13), (0.5, [1.0, 0.0])),  // right
        ((0, 1), (1.0, [0.0, -1.0])),  // bottom
        ((1, 2), (1.0, [0.0, -1.0])),  // bottom
    ]);

    // check if the normal vectors at boundary are outward
    mesh.check_2d_edge_normals(&features.edges, &solutions, 1e-17)
        .expect("ok");

    // search edges
    let feat = Features::new(&mesh, false);
    let edges = feat.search_edge_keys(At::X(0.0), any_x)?;
    assert_eq!(&edges, &[(0, 3), (3, 7), (7, 10), (10, 14)]);
    let edges = feat.search_edge_keys(At::X(4.0), any_x)?;
    assert_eq!(&edges, &[(2, 6), (6, 9), (9, 13)]);

    // edge (7,11)
    let space_ndim = mesh.ndim;
    let p = &mesh.points;
    let mut pad_edge_7_11 = Scratchpad::new(space_ndim, GeoKind::Lin2)?;
    pad_edge_7_11.set_xx(0, 0, p[7].coords[0]);
    pad_edge_7_11.set_xx(0, 1, p[7].coords[1]);
    pad_edge_7_11.set_xx(1, 0, p[11].coords[0]);
    pad_edge_7_11.set_xx(1, 1, p[11].coords[1]);
    let ips = integ::default_points(pad_edge_7_11.kind);
    let mut length_numerical = 0.0;
    for index in 0..ips.len() {
        let iota = &ips[index];
        let weight = ips[index][3];
        let det_jac = pad_edge_7_11.calc_jacobian(iota)?;
        length_numerical += weight * det_jac;
    }
    approx_eq(length_numerical, SQRT_2, 1e-14);

    // TODO: numerical area of cell 5
    let cell = &mesh.cells[5];
    let mut pad_cell_5 = Scratchpad::new(2, cell.kind)?;
    for m in 0..cell.points.len() {
        for j in 0..space_ndim {
            pad_cell_5.set_xx(m, j, p[cell.points[m]].coords[j]);
        }
    }
    let ips = integ::default_points(pad_cell_5.kind);
    let mut area_numerical = 0.0;
    for p in 0..ips.len() {
        let iota = &ips[p];
        let weight = ips[p][3];
        let det_jac = pad_cell_5.calc_jacobian(iota)?;
        area_numerical += weight * det_jac;
    }
    approx_eq(area_numerical, 2.0, 1e-15);
    Ok(())
}
