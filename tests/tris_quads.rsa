use gemlab::integ::default_integ_points;
use gemlab::mesh::{all_edges_2d, allocate_shapes, allocate_states, At, Extract, Features, Find, Mesh, NormalVector};
use gemlab::shapes::{Shape, StateOfShape};
use gemlab::util::SQRT_2;
use gemlab::StrError;
use russell_chk::{assert_approx_eq, assert_vec_approx_eq};
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
fn column_distorted_tris_quads() -> Result<(), StrError> {
    // read mesh
    let mesh = Mesh::from_text_file("./data/meshes/column_distorted_tris_quads.msh")?;

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

    // boundary
    let shapes = allocate_shapes(&mesh)?;
    let edges = all_edges_2d(&mesh, &shapes)?;
    let boundary = Features::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;

    // check edges
    let edge = boundary.edges.get(&(0, 1)).unwrap();
    assert_eq!(edge.points, &[0, 1]);
    let edge = boundary.edges.get(&(9, 10)).unwrap();
    assert_eq!(edge.points, &[10, 9]);
    let edge = boundary.edges.get(&(3, 4)).unwrap();
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
    let ksi = &[0.0, 0.0, 0.0];
    for (edge_keys, solution) in &edge_keys_and_solutions {
        for edge_key in edge_keys {
            let mut normal = NormalVector::at_edge(&mesh, &boundary, *edge_key)?;
            normal.evaluate(ksi)?;
            assert_vec_approx_eq!(normal.value.as_data(), solution, 1e-14);
        }
    }

    // find
    let find = Find::new(&mesh, &boundary)?;

    // find points
    let points = find.points(At::X(0.0))?;
    check(&points, &[0, 1, 2, 3, 4, 5, 6]);
    let points = find.points(At::X(1.0))?;
    check(&points, &[7, 8, 9, 10, 11, 12]);

    // find edges
    let edges = find.edges(At::X(0.0))?;
    check(&edges, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]);
    let edges = find.edges(At::X(1.0))?;
    check(&edges, &[(7, 8), (8, 9), (9, 10), (10, 11), (11, 12)]);

    // find faces
    let faces = find.faces(At::X(0.0))?;
    assert_eq!(faces.len(), 0);
    Ok(())
}

#[test]
fn rectangle_tris_quads() -> Result<(), StrError> {
    // read mesh
    let mesh = Mesh::from_text_file("./data/meshes/rectangle_tris_quads.msh")?;
    let shapes = allocate_shapes(&mesh)?;
    let edges = all_edges_2d(&mesh, &shapes)?;
    let boundary = Features::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;

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
    let ksi = &[0.0, 0.0, 0.0];
    for (edge_keys, solution) in &edge_keys_and_solutions {
        for edge_key in edge_keys {
            let mut normal = NormalVector::at_edge(&mesh, &boundary, *edge_key)?;
            normal.evaluate(ksi)?;
            assert_vec_approx_eq!(normal.value.as_data(), solution, 1e-14);
        }
    }

    // find
    let find = Find::new(&mesh, &boundary)?;

    // find edges
    let edges = find.edges(At::X(0.0))?;
    check(&edges, &[(0, 3), (3, 7), (7, 10), (10, 14)]);
    let edges = find.edges(At::X(4.0))?;
    check(&edges, &[(2, 6), (6, 9), (9, 13)]);

    // edge (7,11)
    let shape_edge_7_11 = Shape::new(2, 1, 2)?;
    let mut state_edge_7_11 = StateOfShape::new(
        1,
        &[
            [mesh.points[7].coords[0], mesh.points[7].coords[1]],
            [mesh.points[11].coords[0], mesh.points[11].coords[1]],
        ],
    )?;
    let ips = default_integ_points(shape_edge_7_11.kind);
    let mut length_numerical = 0.0;
    for index in 0..ips.len() {
        let iota = &ips[index];
        let weight = ips[index][3];
        let det_jac = shape_edge_7_11.calc_jacobian(&mut state_edge_7_11, iota)?;
        length_numerical += weight * det_jac;
    }
    assert_approx_eq!(length_numerical, SQRT_2, 1e-14);

    // states
    let mut states = allocate_states(&mesh)?;

    // cell 5
    let shape_cell_5 = Shape::new(2, 2, 4)?;
    let state_cell_5 = &mut states[5];
    let ips = default_integ_points(shape_cell_5.kind);
    let mut area_numerical = 0.0;
    for p in 0..ips.len() {
        let iota = &ips[p];
        let weight = ips[p][3];
        let det_jac = shape_cell_5.calc_jacobian(state_cell_5, iota)?;
        area_numerical += weight * det_jac;
    }
    assert_approx_eq!(area_numerical, 2.0, 1e-15);
    Ok(())
}
