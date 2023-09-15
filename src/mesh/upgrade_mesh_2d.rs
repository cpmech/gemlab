use super::{Extract, Features, Mesh, Point};
use crate::shapes::{GeoClass, GeoKind, Scratchpad};
use crate::StrError;
use russell_lab::Vector;

/// Upgrades a mesh with triangles or quadrilaterals to a higher order
///
/// # Notes
///
/// 1. All cells must have the same GeoKind
/// 2. Only [GeoClass::Tri] and [GeoClass::Qua] are allowed
/// 3. The target GeoKind must have more nodes that the current GeoKind
/// 4. The current cells must not have interior nodes, i.e., the current cells must **not** be Tri10, Qua9, or Qua16
pub fn upgrade_mesh_2d(mesh: &mut Mesh, target: GeoKind) -> Result<(), StrError> {
    //        2,
    //  s     | ',
    //        |   ',
    //  i    10     9,
    //        |       ',
    //  d     |         ',   side = 1
    //        5    14     4,
    //  e     |             ',
    //        |               ',         NOTE that two neighboring triangles
    //  =    11    12    13     8,       will have reversed pairs of local
    //        |                   ',     points, e.g. (6,7) <-> (7,6)
    //  2     |      side = 0       ',   Thus, the first point along the
    //        0-----6-----3-----7-----1  edge of one triangle will be the
    //        1-----7-----3-----6-----0  second point along the edge of the
    //  2     |      side = 0       ,'   neighbor triangle, and vice-versa
    //        |                   ,'     Solution: Flip Normals
    //  =     8    13    12    11'
    //        |               ,'
    //  e     |             ,'
    //        4    14     5'
    //  d     |         ,'   side = 1
    //        |       ,'
    //  i     9    10'
    //        |   ,'
    //  s     | ,'
    //        2'
    if mesh.ndim != 2 {
        return Err("mesh ndim must be equal to 2");
    }

    // current number of cells
    let ncell = mesh.cells.len();
    if ncell < 1 {
        return Err("the conversion requires at least one cell");
    }

    // get info about first cell in the mesh
    let source = mesh.cells[0].kind;
    let source_nnode = source.nnode();
    let target_nnode = target.nnode();
    if target.class() != source.class() {
        return Err("target class must equal the GeoClass of current cells");
    }
    if target_nnode <= source_nnode {
        return Err("target GeoKind must have more nodes than the current GeoKind");
    }
    if !(target.class() == GeoClass::Tri || target.class() == GeoClass::Qua) {
        return Err("target GeoClass must be Tri or Qua");
    }

    // note: we cannot allow interior nodes in the source type because their position
    // in the coordinates array of the cell are inconsistent with the position in the
    // target array. Thus, we would have to swap point ids and do other checks...
    if source.n_interior_nodes() > 0 {
        return Err("source class must not have interior nodes");
    }

    // info about the edge nodes
    let nedge = source.nedge();
    let source_edge_nnode = source.edge_nnode();
    let target_edge_nnode = target.edge_nnode();
    assert!(target_edge_nnode >= source_edge_nnode);

    // find neighbors (must be done before expanding the connectivity arrays)
    let features = Features::new(mesh, Extract::All);
    let edges = features.all_2d_edges.unwrap();

    // calculate not_set_yet, a number that indicates that a point has not been set yet
    let not_set_yet = ncell * target_nnode;
    let needs_to_move = 1 + not_set_yet;

    // expand connectivity arrays and mark points that aren't set yet or need to move
    for i in 0..ncell {
        if mesh.cells[i].kind != source {
            return Err("all cells must have the same GeoKind");
        }
        mesh.cells[i].points.resize(target_nnode, not_set_yet);
        // mark points that need to be moved
        if source_edge_nnode == 3 && target_edge_nnode == 4 {
            for e in 0..nedge {
                let m = source.edge_node_id(e, 2);
                mesh.cells[i].points[m] += needs_to_move; // we use += here because we need to recover the original id
            }
        }
        if source_edge_nnode == 4 && target_edge_nnode == 5 {
            for e in 0..nedge {
                let m = source.edge_node_id(e, 2);
                mesh.cells[i].points[m] += needs_to_move; // we use += here because we need to recover the original id
                let m = source.edge_node_id(e, 3);
                mesh.cells[i].points[m] += needs_to_move; // we use += here because we need to recover the original id
            }
        }
    }

    // define function to get the original point id, in case it has been marked as need_to_move
    let original_point_id = |current_point_id: usize| {
        if current_point_id > not_set_yet {
            current_point_id - needs_to_move
        } else {
            current_point_id
        }
    };

    // scratchpad for point interpolation
    const NDIM: usize = 2;
    let mut pad = Scratchpad::new(NDIM, source).unwrap();

    // real coordinates of new points
    let mut x = Vector::new(NDIM);

    // get info about target
    let target_first_interior_node = target_nnode - target.n_interior_nodes();

    // add interior points
    // this must be done before (eventually) moving the edge nodes (if any),
    // because the coordinates matrix will change when the nodes are moved
    if target.n_interior_nodes() > 0 {
        for this_cell_id in 0..ncell {
            // set coordinates matrix of the source
            for m in 0..source_nnode {
                let p = original_point_id(mesh.cells[this_cell_id].points[m]);
                let x = mesh.points[p].coords[0];
                let y = mesh.points[p].coords[1];
                pad.set_xx(m, 0, x);
                pad.set_xx(m, 1, y);
            }

            // add the new interior points
            for m in target_first_interior_node..target_nnode {
                // coordinates of new point
                pad.calc_coords(&mut x, target.reference_coords(m)).unwrap();

                // new point
                let new_point_id = mesh.points.len();
                mesh.points.push(Point {
                    id: new_point_id,
                    marker: 0,
                    coords: x.as_data().clone(),
                });

                // set cell point
                mesh.cells[this_cell_id].points[m] = new_point_id;
            }
        }
    }

    // add or move edge points
    // a move occurs, for example, when the middle node of a Lin3 becomes a side node of a Lin4
    for this_cell_id in 0..ncell {
        // set coordinates matrix of the source
        for m in 0..source_nnode {
            let p = original_point_id(mesh.cells[this_cell_id].points[m]);
            let x = mesh.points[p].coords[0];
            let y = mesh.points[p].coords[1];
            pad.set_xx(m, 0, x);
            pad.set_xx(m, 1, y);
        }

        // new points on the edge of the target
        if target_edge_nnode > source_edge_nnode {
            for e in 0..nedge {
                // find the first non-zero marker of a node on the source edge
                let mut marker = 0;
                for j in 2..source_edge_nnode {
                    let n = source.edge_node_id(e, j);
                    let p = original_point_id(mesh.cells[this_cell_id].points[n]);
                    marker = mesh.points[p].marker;
                    if marker != 0 {
                        break;
                    }
                }

                for i in 2..target_edge_nnode {
                    // local point on the edge of target
                    let m = target.edge_node_id_inward(e, i); // order => inward normals

                    // do this next only if the point is marked as not_set_yet or
                    // it is marked as need_to_move. these cases have current_new_point >= not_set_yet
                    let current_new_point = mesh.cells[this_cell_id].points[m];
                    if current_new_point >= not_set_yet {
                        // coordinates of new point
                        pad.calc_coords(&mut x, target.reference_coords(m)).unwrap();

                        // new point
                        let point_id = if current_new_point == not_set_yet {
                            // needs to be added
                            let new_point_id = mesh.points.len();
                            mesh.points.push(Point {
                                id: new_point_id,
                                marker,
                                coords: x.as_data().clone(),
                            });
                            new_point_id
                        } else {
                            // needs to move
                            let original_point_id = current_new_point - needs_to_move;

                            // since the matrix of coordinates (in pad) is already set,
                            // we can change the coordinates now
                            mesh.points[original_point_id].coords[0] = x[0];
                            mesh.points[original_point_id].coords[1] = x[1];
                            original_point_id
                        };

                        // set cell point
                        mesh.cells[this_cell_id].points[m] = point_id;

                        // search cells sharing this edge
                        let local_a = source.edge_node_id(e, 0);
                        let local_b = source.edge_node_id(e, 1);
                        let mut this_a = mesh.cells[this_cell_id].points[local_a];
                        let mut this_b = mesh.cells[this_cell_id].points[local_b];
                        if this_b < this_a {
                            let temp = this_a;
                            this_a = this_b;
                            this_b = temp;
                        }
                        let edge_shares = edges.get(&(this_a, this_b)).unwrap();
                        assert!(edge_shares.len() == 1 || edge_shares.len() == 2);
                        let shares: Vec<_> = edge_shares.iter().filter(|v| v.0 != this_cell_id).collect();

                        // found neighbor
                        if shares.len() == 1 {
                            let neigh_cell_id = shares[0].0;
                            let neigh_edge_id = shares[0].1;
                            let n = target.edge_node_id(neigh_edge_id, i); // order => outward normals

                            // set new point of neighbor
                            mesh.cells[neigh_cell_id].points[n] = point_id;
                        }
                    }
                }
            }
        }

        // set new geo code
        mesh.cells[this_cell_id].kind = target;
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::upgrade_mesh_2d;
    use crate::mesh::{check_all, check_overlapping_points, Cell, Mesh, Point, Samples};
    use crate::shapes::GeoKind;
    use russell_chk::vec_approx_eq;

    #[allow(unused_imports)]
    use crate::mesh::draw_mesh;

    const SAVE_FIGURE: bool = false;

    #[test]
    fn upgrade_mesh_2d_captures_errors() {
        #[rustfmt::skip]
        let mut mesh = Mesh {
            ndim: 3,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0, 0.0, 0.0] },
                Point { id: 1, marker: 0, coords: vec![1.0, 0.0, 0.0] },
                Point { id: 2, marker: 0, coords: vec![1.0, 1.0, 0.0] },
                Point { id: 3, marker: 0, coords: vec![0.0, 1.0, 0.0] },
                Point { id: 4, marker: 0, coords: vec![0.0, 0.0, 1.0] },
                Point { id: 5, marker: 0, coords: vec![1.0, 0.0, 1.0] },
                Point { id: 6, marker: 0, coords: vec![1.0, 1.0, 1.0] },
                Point { id: 7, marker: 0, coords: vec![0.0, 1.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Hex8, points: vec![0,1,2,3, 4,5,6,7] },
            ],
        };
        assert_eq!(
            upgrade_mesh_2d(&mut mesh, GeoKind::Tri15).err(),
            Some("mesh ndim must be equal to 2")
        );
        let mut mesh = Mesh {
            ndim: 2,
            points: vec![],
            cells: vec![],
        };
        assert_eq!(
            upgrade_mesh_2d(&mut mesh, GeoKind::Tri15).err(),
            Some("the conversion requires at least one cell")
        );
        #[rustfmt::skip]
        let mut mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0, 0.0 ] },
                Point { id: 1, marker: 0, coords: vec![1.0, 0.0 ] },
                Point { id: 2, marker: 0, coords: vec![0.5, 0.85] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Tri3, points: vec![0, 1, 2] },
            ],
        };
        assert_eq!(
            upgrade_mesh_2d(&mut mesh, GeoKind::Qua8).err(),
            Some("target class must equal the GeoClass of current cells")
        );
        assert_eq!(
            upgrade_mesh_2d(&mut mesh, GeoKind::Tri3).err(),
            Some("target GeoKind must have more nodes than the current GeoKind")
        );
        #[rustfmt::skip]
        let mut mesh =Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, marker: 0, coords: vec![1.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
            ],
        };
        assert_eq!(
            upgrade_mesh_2d(&mut mesh, GeoKind::Lin3).err(),
            Some("target GeoClass must be Tri or Qua")
        );
        #[rustfmt::skip]
        let mut mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0, 0.0 ] }, // 0
                Point { id: 1, marker: 0, coords: vec![0.8, 0.0 ] }, // 1
                Point { id: 2, marker: 0, coords: vec![0.8, 0.8 ] }, // 2
                Point { id: 3, marker: 0, coords: vec![0.0, 0.8 ] }, // 3
                Point { id: 4, marker: 0, coords: vec![0.4, 0.05] }, // 4
                Point { id: 5, marker: 0, coords: vec![0.8, 0.4 ] }, // 5
                Point { id: 6, marker: 0, coords: vec![0.4, 0.85] }, // 6
                Point { id: 7, marker: 0, coords: vec![0.0, 0.4 ] }, // 7
                Point { id: 8, marker: 0, coords: vec![0.4, 0.4 ] }, // 8
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Qua9, points: vec![0, 1, 2, 3, 4, 5, 6, 7, 8] },
            ],
        };
        assert_eq!(
            upgrade_mesh_2d(&mut mesh, GeoKind::Qua12).err(),
            Some("source class must not have interior nodes")
        );
        #[rustfmt::skip]
        let mut mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0, 0.0 ] },
                Point { id: 1, marker: 0, coords: vec![1.0, 0.0 ] },
                Point { id: 2, marker: 0, coords: vec![0.5, 0.85] },
                Point { id: 3, marker: 0, coords: vec![1.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Tri3, points: vec![0, 1, 2] },
                Cell { id: 1, attribute: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
            ],
        };
        assert_eq!(
            upgrade_mesh_2d(&mut mesh, GeoKind::Tri6).err(),
            Some("all cells must have the same GeoKind")
        );
    }

    #[test]
    fn upgrade_tri6_to_tri15_works_1() {
        #[rustfmt::skip]
        let mut mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords:   vec![0.0, 0.0] },
                Point { id: 1, marker: 0, coords:   vec![4.0, 0.0] },
                Point { id: 2, marker: 0, coords:   vec![0.0, 4.0] },
                Point { id: 3, marker: 0, coords:   vec![2.0, 0.0] },
                Point { id: 4, marker: -20, coords: vec![2.0, 2.0] },
                Point { id: 5, marker: -10, coords: vec![0.0, 2.0] },
                Point { id: 6, marker: 0, coords:   vec![0.0,-4.0] },
                Point { id: 7, marker: -30, coords: vec![2.0,-2.0] },
                Point { id: 8, marker: -10, coords: vec![0.0,-2.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Tri6, points: vec![1, 2, 0, 4, 5, 3] },
                Cell { id: 1, attribute: 2, kind: GeoKind::Tri6, points: vec![0, 6, 1, 8, 7, 3] },
            ],
        };
        check_all(&mesh).unwrap();

        if SAVE_FIGURE {
            draw_mesh(&mesh, true, true, false, "/tmp/gemlab/test_tri6_to_tri15_1_before.svg").unwrap();
        }

        upgrade_mesh_2d(&mut mesh, GeoKind::Tri15).unwrap();

        if SAVE_FIGURE {
            draw_mesh(&mesh, true, true, false, "/tmp/gemlab/test_tri6_to_tri15_1_after.svg").unwrap();
        }

        check_all(&mesh).unwrap();
        check_overlapping_points(&mesh, 0.2).unwrap();

        assert_eq!(mesh.points.len(), 25);
        assert_eq!(
            mesh.cells[0].points,
            &[1, 2, 0, 4, 5, 3, 15, 16, 17, 18, 19, 20, 9, 10, 11]
        );
        assert_eq!(
            mesh.cells[1].points,
            &[0, 6, 1, 8, 7, 3, 21, 22, 23, 24, 20, 19, 12, 13, 14]
        );

        // original
        vec_approx_eq(&mesh.points[0].coords, &[0.0, 0.0], 1e-15);
        vec_approx_eq(&mesh.points[1].coords, &[4.0, 0.0], 1e-15);
        vec_approx_eq(&mesh.points[2].coords, &[0.0, 4.0], 1e-15);
        vec_approx_eq(&mesh.points[3].coords, &[2.0, 0.0], 1e-15);
        vec_approx_eq(&mesh.points[4].coords, &[2.0, 2.0], 1e-15);
        vec_approx_eq(&mesh.points[5].coords, &[0.0, 2.0], 1e-15);
        vec_approx_eq(&mesh.points[6].coords, &[0.0, -4.0], 1e-15);
        vec_approx_eq(&mesh.points[7].coords, &[2.0, -2.0], 1e-15);
        vec_approx_eq(&mesh.points[8].coords, &[0.0, -2.0], 1e-15);

        // interior
        vec_approx_eq(&mesh.points[9].coords, &[2.0, 1.0], 1e-15);
        vec_approx_eq(&mesh.points[10].coords, &[1.0, 2.0], 1e-15);
        vec_approx_eq(&mesh.points[11].coords, &[1.0, 1.0], 1e-15);
        vec_approx_eq(&mesh.points[12].coords, &[1.0, -1.0], 1e-15);
        vec_approx_eq(&mesh.points[13].coords, &[1.0, -2.0], 1e-15);
        vec_approx_eq(&mesh.points[14].coords, &[2.0, -1.0], 1e-15);

        // edges of cell # 0 and cell # 1
        vec_approx_eq(&mesh.points[15].coords, &[3.0, 1.0], 1e-15);
        vec_approx_eq(&mesh.points[16].coords, &[1.0, 3.0], 1e-15);
        vec_approx_eq(&mesh.points[17].coords, &[0.0, 3.0], 1e-15);
        vec_approx_eq(&mesh.points[18].coords, &[0.0, 1.0], 1e-15);
        vec_approx_eq(&mesh.points[19].coords, &[1.0, 0.0], 1e-15);
        vec_approx_eq(&mesh.points[20].coords, &[3.0, 0.0], 1e-15);

        // edges of cell # 1
        vec_approx_eq(&mesh.points[21].coords, &[0.0, -1.0], 1e-15);
        vec_approx_eq(&mesh.points[22].coords, &[0.0, -3.0], 1e-15);
        vec_approx_eq(&mesh.points[23].coords, &[1.0, -3.0], 1e-15);
        vec_approx_eq(&mesh.points[24].coords, &[3.0, -1.0], 1e-15);

        // original
        assert_eq!(mesh.points[0].marker, 0);
        assert_eq!(mesh.points[1].marker, 0);
        assert_eq!(mesh.points[2].marker, 0);
        assert_eq!(mesh.points[3].marker, 0);
        assert_eq!(mesh.points[4].marker, -20);
        assert_eq!(mesh.points[5].marker, -10);
        assert_eq!(mesh.points[6].marker, 0);
        assert_eq!(mesh.points[7].marker, -30);
        assert_eq!(mesh.points[8].marker, -10);

        // interior
        assert_eq!(mesh.points[9].marker, 0);
        assert_eq!(mesh.points[10].marker, 0);
        assert_eq!(mesh.points[11].marker, 0);
        assert_eq!(mesh.points[12].marker, 0);
        assert_eq!(mesh.points[13].marker, 0);
        assert_eq!(mesh.points[14].marker, 0);

        // edges of cell # 0 and cell # 1
        assert_eq!(mesh.points[15].marker, -20);
        assert_eq!(mesh.points[16].marker, -20);
        assert_eq!(mesh.points[17].marker, -10);
        assert_eq!(mesh.points[18].marker, -10);
        assert_eq!(mesh.points[19].marker, 0);
        assert_eq!(mesh.points[20].marker, 0);

        // edges of cell # 1
        assert_eq!(mesh.points[21].marker, -10);
        assert_eq!(mesh.points[22].marker, -10);
        assert_eq!(mesh.points[23].marker, -30);
        assert_eq!(mesh.points[24].marker, -30);
    }

    #[test]
    fn upgrade_tri6_to_tri10_works() {
        #[rustfmt::skip]
        let mut mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, marker: 0, coords: vec![4.0, 0.0] },
                Point { id: 2, marker: 0, coords: vec![0.0, 4.0] },
                Point { id: 3, marker: 0, coords: vec![2.0, 0.0] },
                Point { id: 4, marker: -20, coords: vec![2.0, 2.0] },
                Point { id: 5, marker: -10, coords: vec![0.0, 2.0] },
                Point { id: 6, marker: 0, coords: vec![0.0,-4.0] },
                Point { id: 7, marker: -30, coords: vec![2.0,-2.0] },
                Point { id: 8, marker: -10, coords: vec![0.0,-2.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Tri6, points: vec![1, 2, 0, 4, 5, 3] },
                Cell { id: 1, attribute: 2, kind: GeoKind::Tri6, points: vec![0, 6, 1, 8, 7, 3] },
            ],
        };
        check_all(&mesh).unwrap();

        if SAVE_FIGURE {
            draw_mesh(&mesh, true, true, false, "/tmp/gemlab/test_tri6_to_tri10_before.svg").unwrap();
        }

        upgrade_mesh_2d(&mut mesh, GeoKind::Tri10).unwrap();

        if SAVE_FIGURE {
            draw_mesh(&mesh, true, true, false, "/tmp/gemlab/test_tri6_to_tri10_after.svg").unwrap();
        }

        check_all(&mesh).unwrap();
        check_overlapping_points(&mesh, 0.2).unwrap();

        assert_eq!(mesh.points.len(), 16);

        // original
        let ft = 4.0 / 3.0;
        let et = 8.0 / 3.0;
        vec_approx_eq(&mesh.points[0].coords, &[0.0, 0.0], 1e-15);
        vec_approx_eq(&mesh.points[1].coords, &[4.0, 0.0], 1e-15);
        vec_approx_eq(&mesh.points[2].coords, &[0.0, 4.0], 1e-15);
        vec_approx_eq(&mesh.points[3].coords, &[ft, 0.0], 1e-15); // moved
        vec_approx_eq(&mesh.points[4].coords, &[et, ft], 1e-15); // moved
        vec_approx_eq(&mesh.points[5].coords, &[0.0, et], 1e-15); // moved
        vec_approx_eq(&mesh.points[6].coords, &[0.0, -4.0], 1e-15);
        vec_approx_eq(&mesh.points[7].coords, &[ft, -et], 1e-15); // moved
        vec_approx_eq(&mesh.points[8].coords, &[0.0, -ft], 1e-15); // moved

        // interior
        vec_approx_eq(&mesh.points[9].coords, &[4.0 / 3.0, 4.0 / 3.0], 1e-15);
        vec_approx_eq(&mesh.points[10].coords, &[4.0 / 3.0, -4.0 / 3.0], 1e-15);

        // edges of cell # 0 and cell # 1
        vec_approx_eq(&mesh.points[11].coords, &[ft, et], 1e-15);
        vec_approx_eq(&mesh.points[12].coords, &[0.0, ft], 1e-15);
        vec_approx_eq(&mesh.points[13].coords, &[et, 0.0], 1e-15);

        // edges of cell # 1
        vec_approx_eq(&mesh.points[14].coords, &[0.0, -et], 1e-15);
        vec_approx_eq(&mesh.points[15].coords, &[et, -ft], 1e-15);

        // original
        assert_eq!(mesh.points[0].marker, 0);
        assert_eq!(mesh.points[1].marker, 0);
        assert_eq!(mesh.points[2].marker, 0);
        assert_eq!(mesh.points[3].marker, 0);
        assert_eq!(mesh.points[4].marker, -20);
        assert_eq!(mesh.points[5].marker, -10);
        assert_eq!(mesh.points[6].marker, 0);
        assert_eq!(mesh.points[7].marker, -30);
        assert_eq!(mesh.points[8].marker, -10);

        // interior
        assert_eq!(mesh.points[9].marker, 0);
        assert_eq!(mesh.points[10].marker, 0);

        // edges of cell # 0 and cell # 1
        assert_eq!(mesh.points[11].marker, -20);
        assert_eq!(mesh.points[12].marker, -10);
        assert_eq!(mesh.points[13].marker, 0);

        // edges of cell # 1
        assert_eq!(mesh.points[14].marker, -10);
        assert_eq!(mesh.points[15].marker, -30);
    }

    #[test]
    fn upgrade_tri3_to_tri6_works() {
        let mut mesh = Samples::two_tri3().clone();
        upgrade_mesh_2d(&mut mesh, GeoKind::Tri6).unwrap();
        if SAVE_FIGURE {
            draw_mesh(&mesh, true, true, false, "/tmp/gemlab/test_tri3_to_tri6_after.svg").unwrap();
        }
        check_all(&mesh).unwrap();
        check_overlapping_points(&mesh, 0.1).unwrap();
        assert_eq!(mesh.points.len(), 9);
        assert_eq!(mesh.cells[0].points, &[0, 1, 3, 4, 5, 6]);
        assert_eq!(mesh.cells[1].points, &[2, 3, 1, 7, 5, 8]);
    }

    #[test]
    fn upgrade_tri3_to_tri10_works() {
        let mut mesh = Samples::two_tri3().clone();
        upgrade_mesh_2d(&mut mesh, GeoKind::Tri10).unwrap();
        if SAVE_FIGURE {
            draw_mesh(&mesh, true, true, false, "/tmp/gemlab/test_tri3_to_tri10_after.svg").unwrap();
        }
        check_all(&mesh).unwrap();
        check_overlapping_points(&mesh, 0.1).unwrap();
        assert_eq!(mesh.points.len(), 16);
        assert_eq!(mesh.cells[0].points, &[0, 1, 3, 6, 8, 10, 7, 9, 11, 4]);
        assert_eq!(mesh.cells[1].points, &[2, 3, 1, 12, 9, 14, 13, 8, 15, 5]);
    }

    #[test]
    fn upgrade_qua4_to_qua8_works() {
        let mut mesh = Samples::two_qua4().clone();
        upgrade_mesh_2d(&mut mesh, GeoKind::Qua8).unwrap();
        if SAVE_FIGURE {
            draw_mesh(&mesh, true, true, false, "/tmp/gemlab/test_qua4_to_qua8_after.svg").unwrap();
        }
        check_all(&mesh).unwrap();
        check_overlapping_points(&mesh, 0.2).unwrap();
        assert_eq!(mesh.points.len(), 13);
        assert_eq!(mesh.cells[0].points, &[0, 1, 2, 3, 6, 7, 8, 9]);
        assert_eq!(mesh.cells[1].points, &[1, 4, 5, 2, 10, 11, 12, 7]);
    }

    #[test]
    fn upgrade_qua12_to_qua16_works() {
        let mut mesh = Samples::block_2d_four_qua12().clone();
        upgrade_mesh_2d(&mut mesh, GeoKind::Qua16).unwrap();
        if SAVE_FIGURE {
            draw_mesh(&mesh, true, true, false, "/tmp/gemlab/test_qua12_to_qua16_after.svg").unwrap();
        }
        check_all(&mesh).unwrap();
        check_overlapping_points(&mesh, 0.2).unwrap();
        assert_eq!(mesh.points.len(), 33 + 4 * 4);
        assert_eq!(
            mesh.cells[0].points,
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 33, 34, 35, 36]
        );
        assert_eq!(
            mesh.cells[1].points,
            [1, 12, 13, 2, 14, 15, 16, 9, 17, 18, 19, 5, 37, 38, 39, 40]
        );
        assert_eq!(
            mesh.cells[2].points,
            [3, 2, 20, 21, 10, 22, 23, 24, 6, 25, 26, 27, 41, 42, 43, 44]
        );
        assert_eq!(
            mesh.cells[3].points,
            [2, 13, 28, 20, 19, 29, 30, 25, 16, 31, 32, 22, 45, 46, 47, 48]
        );
    }
}
