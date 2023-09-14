use super::{Extract, Features, Mesh, Point};
use crate::shapes::{GeoClass, GeoKind, Scratchpad};
use crate::StrError;
use russell_lab::Vector;

/// Upgrades a mesh with triangles or quadrilaterals to a higher order
pub fn upgrade_tri_mesh(mesh: &mut Mesh, target: GeoKind) -> Result<(), StrError> {
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

    // current number of points and cells
    let npoint = mesh.points.len();
    let ncell = mesh.cells.len();
    if ncell < 1 {
        return Err("the conversion requires at least one cell");
    }

    // get info about first cell in the mesh
    let source = mesh.cells[0].kind;
    let source_nnode = source.nnode();
    let target_nnode = target.nnode();
    let delta_nnode = target_nnode - source_nnode;
    if target.class() != source.class() {
        return Err("target class must equal the GeoClass of current cells");
    }
    if !(target.class() == GeoClass::Tri || target.class() == GeoClass::Qua) {
        return Err("target GeoClass must be Tri or Qua");
    }
    if delta_nnode < 1 {
        return Err("target GeoKind must have more nodes than the current GeoKind");
    }
    let nedge = source.nedge();
    let source_edge_nnode = source.edge_nnode();
    let target_edge_nnode = target.edge_nnode();
    assert!(target_edge_nnode >= source_edge_nnode);

    // not_set_yet: number that indicates that a point has not been set yet
    // here we use 1 + max number of points possible, ignoring edge overlaps
    // for example, in a mesh with two Tri6 cells (ncell = 2, npoint = 9), if target = Tri15,
    // then delta_nnode = 9, and not_set_yet will be equal to 1 + 9 + 2 * 9 = 28,
    // however, the final total number of points will be 25
    let not_set_yet = 1 + npoint + ncell * delta_nnode;

    // expand connectivity array
    for i in 0..ncell {
        if mesh.cells[i].kind != source {
            return Err("all current cells must have the same GeoKind");
        }
        mesh.cells[i].points.resize(target_nnode, not_set_yet);
    }

    // scratchpad for point interpolation
    const NDIM: usize = 2;
    let mut pad = Scratchpad::new(NDIM, source).unwrap();

    // real coordinates of new points
    let mut x = Vector::new(NDIM);

    // find neighbors
    let features = Features::new(mesh, Extract::All);
    let edges = features.all_2d_edges.unwrap();

    // get info about target
    let target_first_interior_node = target_nnode - target.n_interior_nodes();

    // perform the conversion
    for this_cell_id in 0..ncell {
        // set coordinates matrix of the source kind
        for m in 0..source_nnode {
            let p = mesh.cells[this_cell_id].points[m];
            let x = mesh.points[p].coords[0];
            let y = mesh.points[p].coords[1];
            pad.set_xx(m, 0, x);
            pad.set_xx(m, 1, y);
        }

        // new points on the edge of the target kind
        if target_edge_nnode > source_edge_nnode {
            for e in 0..nedge {
                for i in 2..target_edge_nnode {
                    // local point on the edge of target
                    let m = target.edge_node_id_inward(e, i); // order => inward normals

                    let current_new_point = mesh.cells[this_cell_id].points[m];
                    if current_new_point == not_set_yet {
                        // use marker of the middle point, if any
                        let marker = if source_edge_nnode > 2 {
                            let local_mid = source.edge_node_id_inward(e, 2);
                            let mid_point = mesh.cells[this_cell_id].points[local_mid];
                            mesh.points[mid_point].marker
                        } else {
                            0
                        };

                        // coordinates of new point
                        pad.calc_coords(&mut x, target.reference_coords(m)).unwrap();

                        // new point
                        let new_point_id = mesh.points.len();
                        mesh.points.push(Point {
                            id: new_point_id,
                            marker,
                            coords: x.as_data().clone(),
                        });

                        // set cell point
                        mesh.cells[this_cell_id].points[m] = new_point_id;

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
                        let shares: Vec<_> = edge_shares.iter().filter(|v| v.0 != this_cell_id).collect();

                        // found neighbor
                        if shares.len() == 1 {
                            let neigh_cell_id = shares[0].0;
                            let neigh_edge_id = shares[0].1;
                            let n = target.edge_node_id(neigh_edge_id, i); // order => outward normals

                            // set new point of neighbor
                            mesh.cells[neigh_cell_id].points[n] = new_point_id;
                        }
                    }
                }
            }
        }

        // new points in the interior
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

        // set new geo code
        mesh.cells[this_cell_id].kind = target;
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::upgrade_tri_mesh;
    use crate::mesh::{check_all, Cell, Mesh, Point};
    use crate::shapes::GeoKind;
    use russell_chk::vec_approx_eq;

    #[allow(unused_imports)]
    use crate::mesh::draw_mesh;

    const SAVE_FIGURE: bool = false;

    #[test]
    fn upgrade_tri6_to_tri15_captures_errors() {
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
            upgrade_tri_mesh(&mut mesh, GeoKind::Tri15).err(),
            Some("mesh ndim must be equal to 2")
        );
    }

    #[test]
    fn upgrade_tri6_to_tri15_works_1() {
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
            draw_mesh(&mesh, true, true, false, "/tmp/gemlab/test_tri6_to_tri15_1_before.svg").unwrap();
        }

        upgrade_tri_mesh(&mut mesh, GeoKind::Tri15).unwrap();

        if SAVE_FIGURE {
            draw_mesh(&mesh, true, true, false, "/tmp/gemlab/test_tri6_to_tri15_1_after.svg").unwrap();
        }

        check_all(&mesh).unwrap();
        assert_eq!(mesh.points.len(), 25);
        assert_eq!(
            mesh.cells[0].points,
            &[1, 2, 0, 4, 5, 3, 9, 10, 11, 12, 13, 14, 15, 16, 17]
        );
        assert_eq!(
            mesh.cells[1].points,
            &[0, 6, 1, 8, 7, 3, 18, 19, 20, 21, 14, 13, 22, 23, 24]
        );
        vec_approx_eq(&mesh.points[9].coords, &[3.0, 1.0], 1e-15);
        vec_approx_eq(&mesh.points[10].coords, &[1.0, 3.0], 1e-15);
        vec_approx_eq(&mesh.points[11].coords, &[0.0, 3.0], 1e-15);
        vec_approx_eq(&mesh.points[12].coords, &[0.0, 1.0], 1e-15);
        vec_approx_eq(&mesh.points[13].coords, &[1.0, 0.0], 1e-15);
        vec_approx_eq(&mesh.points[14].coords, &[3.0, 0.0], 1e-15);
        vec_approx_eq(&mesh.points[15].coords, &[2.0, 1.0], 1e-15);
        vec_approx_eq(&mesh.points[16].coords, &[1.0, 2.0], 1e-15);
        vec_approx_eq(&mesh.points[17].coords, &[1.0, 1.0], 1e-15);
        vec_approx_eq(&mesh.points[18].coords, &[0.0, -1.0], 1e-15);
        vec_approx_eq(&mesh.points[19].coords, &[0.0, -3.0], 1e-15);
        vec_approx_eq(&mesh.points[20].coords, &[1.0, -3.0], 1e-15);
        vec_approx_eq(&mesh.points[21].coords, &[3.0, -1.0], 1e-15);
        vec_approx_eq(&mesh.points[22].coords, &[1.0, -1.0], 1e-15);
        vec_approx_eq(&mesh.points[23].coords, &[1.0, -2.0], 1e-15);
        vec_approx_eq(&mesh.points[24].coords, &[2.0, -1.0], 1e-15);
        assert_eq!(mesh.points[9].marker, -20);
        assert_eq!(mesh.points[10].marker, -20);
        assert_eq!(mesh.points[11].marker, -10);
        assert_eq!(mesh.points[12].marker, -10);
        assert_eq!(mesh.points[13].marker, 0);
        assert_eq!(mesh.points[14].marker, 0);
        assert_eq!(mesh.points[15].marker, 0);
        assert_eq!(mesh.points[16].marker, 0);
        assert_eq!(mesh.points[17].marker, 0);
        assert_eq!(mesh.points[18].marker, -10);
        assert_eq!(mesh.points[19].marker, -10);
        assert_eq!(mesh.points[20].marker, -30);
        assert_eq!(mesh.points[21].marker, -30);
        assert_eq!(mesh.points[22].marker, 0);
        assert_eq!(mesh.points[23].marker, 0);
        assert_eq!(mesh.points[24].marker, 0);
    }
}
