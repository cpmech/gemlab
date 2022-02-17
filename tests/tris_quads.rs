use gemlab::mesh::{At, Mesh};
use gemlab::shapes::Shape;
use gemlab::util::SQRT_2;
use gemlab::StrError;
use russell_chk::{assert_approx_eq, assert_vec_approx_eq};
use russell_lab::Vector;

#[test]
fn column_distorted_tris_quads() -> Result<(), StrError> {
    // read mesh
    let mesh = Mesh::from_text_file("./data/meshes/column_distorted_tris_quads.msh")?;

    // check sizes
    assert_eq!(mesh.points.len(), 13);
    assert_eq!(mesh.cells.len(), 7);
    assert_eq!(mesh.boundary_points.len(), 13);
    assert_eq!(mesh.boundary_edges.len(), 13);
    assert_eq!(mesh.boundary_faces.len(), 0);

    // check cells
    assert_eq!(mesh.cells[0].points, &[0, 7, 8, 1]);
    assert_eq!(mesh.cells[1].points, &[1, 8, 9, 2]);
    assert_eq!(mesh.cells[2].points, &[2, 9, 10, 3]);
    assert_eq!(mesh.cells[3].points, &[3, 10, 4]);
    assert_eq!(mesh.cells[4].points, &[4, 10, 11]);
    assert_eq!(mesh.cells[5].points, &[4, 11, 5]);
    assert_eq!(mesh.cells[6].points, &[5, 11, 12, 6]);

    // check edges
    let edge = mesh.boundary_edges.get(&(0, 1)).unwrap();
    assert_eq!(edge.points, &[0, 1]);
    let edge = mesh.boundary_edges.get(&(9, 10)).unwrap();
    assert_eq!(edge.points, &[10, 9]);
    let edge = mesh.boundary_edges.get(&(3, 4)).unwrap();
    assert_eq!(edge.points, &[3, 4]);

    // the magnitude of the normal vector should be equal to edge_length / 2.0
    // for both tris or quas where 2.0 corresponds to the edge_length in the reference system
    // Note that the edge is mapped to [-1,+1] in both Lin or Qua

    // edge keys and correct normal vectors (solutions)
    let edge_keys_and_solutions = [
        // left
        (vec![(0, 1)], [-0.5 / 2.0, 0.0]),
        (vec![(1, 2)], [-0.5 / 2.0, 0.0]),
        (vec![(2, 3)], [-1.0 / 2.0, 0.0]),
        (vec![(3, 4)], [-0.5 / 2.0, 0.0]),
        (vec![(4, 5)], [-0.5 / 2.0, 0.0]),
        (vec![(5, 6)], [-0.1 / 2.0, 0.0]),
        // right
        (vec![(7, 8)], [0.2 / 2.0, 0.0]),
        (vec![(8, 9)], [0.8 / 2.0, 0.0]),
        (vec![(9, 10)], [0.8 / 2.0, 0.0]),
        (vec![(10, 11)], [1.2 / 2.0, 0.0]),
        (vec![(11, 12)], [0.1 / 2.0, 0.0]),
        // bottom
        (vec![(0, 7)], [0.0, -1.0 / 2.0]),
        // top
        (vec![(6, 12)], [0.0, 1.0 / 2.0]),
    ];

    // check if the normal vectors at boundary are outward
    let mut normal = Vector::new(mesh.space_ndim);
    let ksi = &[0.0, 0.0, 0.0];
    for (edge_keys, solution) in &edge_keys_and_solutions {
        for edge_key in edge_keys {
            let mut edge_shape = mesh.alloc_shape_boundary_edge(edge_key)?;
            assert_eq!(edge_shape.nnode, 2);
            edge_shape.calc_boundary_normal(&mut normal, ksi)?;
            assert_vec_approx_eq!(normal.as_data(), solution, 1e-14);
        }
    }

    // find points
    let points = mesh.find_boundary_points(At::X(0.0))?;
    assert_eq!(points, &[0, 1, 2, 3, 4, 5, 6]);
    let points = mesh.find_boundary_points(At::X(1.0))?;
    assert_eq!(points, &[7, 8, 9, 10, 11, 12]);

    // find edges
    let edges = mesh.find_boundary_edges(At::X(0.0))?;
    assert_eq!(edges, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]);
    let edges = mesh.find_boundary_edges(At::X(1.0))?;
    assert_eq!(edges, &[(7, 8), (8, 9), (9, 10), (10, 11), (11, 12)]);

    // find faces
    let faces = mesh.find_boundary_faces(At::X(0.0))?;
    assert_eq!(faces.len(), 0);
    Ok(())
}

#[test]
fn rectangle_tris_quads() -> Result<(), StrError> {
    // read mesh
    let mesh = Mesh::from_text_file("./data/meshes/rectangle_tris_quads.msh")?;

    // the magnitude of the normal vector should be equal to edge_length / 2.0
    // for both tris or quas where 2.0 corresponds to the edge_length in the reference system
    // Note that the edge is mapped to [-1,+1] in both Lin or Qua

    // edge keys and correct normal vectors (solutions)
    let edge_keys_and_solutions = [
        // left
        (vec![(0, 3), (3, 7), (7, 10)], [-1.0 / 2.0, 0.0]),
        // right
        (vec![(2, 6), (6, 9), (9, 13)], [1.0 / 2.0, 0.0]),
        // bottom
        (vec![(0, 1), (1, 2)], [0.0, -2.0 / 2.0]),
    ];

    // check if the normal vectors at boundary are outward
    let mut normal = Vector::new(mesh.space_ndim);
    let ksi = &[0.0, 0.0, 0.0];
    for (edge_keys, solution) in &edge_keys_and_solutions {
        for edge_key in edge_keys {
            let mut edge_shape = mesh.alloc_shape_boundary_edge(edge_key)?;
            assert_eq!(edge_shape.nnode, 2);
            edge_shape.calc_boundary_normal(&mut normal, ksi)?;
            assert_vec_approx_eq!(normal.as_data(), solution, 1e-14);
        }
    }

    // find edges
    let edges = mesh.find_boundary_edges(At::X(0.0))?;
    assert_eq!(edges, &[(0, 3), (3, 7), (7, 10), (10, 14)]);
    let edges = mesh.find_boundary_edges(At::X(4.0))?;
    assert_eq!(edges, &[(2, 6), (6, 9), (9, 13)]);

    // edge (4,7)
    let mut shape_edge_4_7 = Shape::new(2, 1, 2)?;
    shape_edge_4_7.set_node(7, 0, 0, mesh.points[7].coords[0])?;
    shape_edge_4_7.set_node(7, 0, 1, mesh.points[7].coords[1])?;
    shape_edge_4_7.set_node(11, 1, 0, mesh.points[11].coords[0])?;
    shape_edge_4_7.set_node(11, 1, 1, mesh.points[11].coords[1])?;
    let length_edge_4_7 = SQRT_2;
    let mut length_numerical = 0.0;
    for index in 0..shape_edge_4_7.integ_points.len() {
        let iota = &shape_edge_4_7.integ_points[index];
        let weight = shape_edge_4_7.integ_points[index][3];
        let det_jac = shape_edge_4_7.calc_jacobian(iota)?;
        length_numerical += weight * det_jac;
    }
    assert_approx_eq!(length_edge_4_7, length_numerical, 1e-14);
    let mut normal = Vector::new(2);
    shape_edge_4_7.calc_boundary_normal(&mut normal, &[0.0, 0.0])?;
    let l = length_edge_4_7 / 2.0;
    assert_vec_approx_eq!(normal.as_data(), &[-l / SQRT_2, l / SQRT_2], 1e-14);

    // cell 5
    let mut shape_cell_5 = mesh.alloc_shape_cell(5)?;
    let mut area_numerical = 0.0;
    for index in 0..shape_cell_5.integ_points.len() {
        let iota = &shape_cell_5.integ_points[index];
        let weight = shape_cell_5.integ_points[index][3];
        let det_jac = shape_cell_5.calc_jacobian(iota)?;
        area_numerical += weight * det_jac;
    }
    assert_approx_eq!(area_numerical, length_edge_4_7 * length_edge_4_7, 1e-15);
    Ok(())
}
