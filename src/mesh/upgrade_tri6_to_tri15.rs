use super::{Extract, Features, Mesh, Point};
use crate::shapes::{GeoKind, Scratchpad, Tri15, Tri6};
use crate::StrError;
use russell_lab::Vector;

/// Converts a mesh with only Tri6 cells to a mesh with Tri15 cells
pub fn upgrade_tri6_to_tri15(mesh: &mut Mesh) -> Result<(), StrError> {
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
    //        |                   ,'
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

    // current number of points and cells
    let npoint = mesh.points.len();
    let ncell = mesh.cells.len();

    // 1 + max possible number of points (neglecting overlap)
    // this number will indicate that a point has not been set yet
    // if there is only one Tri6 cell, not_set_yet will be equal to 16
    let not_set_yet = 1 + npoint + ncell * 9;

    // expand connectivity array
    for i in 0..ncell {
        if mesh.cells[i].kind != GeoKind::Tri6 {
            return Err("only Tri6 cells are allowed in the conversion");
        }
        mesh.cells[i].points.resize(15, not_set_yet);
    }

    // scratchpad for point interpolation
    const NDIM: usize = 2;
    let mut pad = Scratchpad::new(NDIM, GeoKind::Tri6)?;

    // real coordinates of new points
    let mut x = Vector::new(NDIM);

    // find neighbors
    let features = Features::new(mesh, Extract::All);
    let edges = features.all_2d_edges.unwrap();

    // perform the conversion
    for this_cell_id in 0..ncell {
        // access this cell (opposite to neighbor cell)
        // let this_cell = &mut mesh.cells[this_cell_id];

        // set coordinates matrix of Tri6
        for m in 0..6 {
            let p = mesh.cells[this_cell_id].points[m];
            let x = mesh.points[p].coords[0];
            let y = mesh.points[p].coords[1];
            pad.set_xx(m, 0, x);
            pad.set_xx(m, 1, y);
        }

        // new points on the edge of the Tri15
        for m in 6..12 {
            let current_new_point = mesh.cells[this_cell_id].points[m];
            if current_new_point == not_set_yet {
                // id and marker of middle point
                let local_mid = m / 2; // 3 3  4 4  5 5
                let mid_point = mesh.cells[this_cell_id].points[local_mid];
                let marker = mesh.points[mid_point].marker;

                // coordinates of new point
                pad.calc_coords(&mut x, &Tri15::NODE_REFERENCE_COORDS[m])?;

                // new point
                let new_point_id = mesh.points.len();
                mesh.points.push(Point {
                    id: new_point_id,
                    marker,
                    coords: x.as_data().clone(),
                });

                // set cell point
                mesh.cells[this_cell_id].points[m] = new_point_id;

                // parameters
                // m         = 6 7  8 9 10 11
                // m - 6     = 0 1  2 3  4 5
                // this_side =  0    1    2   = (m - 6) / 2
                // this_pos  = 0 1  0 1  0 1  = (m - 6) % 2
                // neigh_pos = 1 0  1 0  1 0  = 1 - this_pos
                let this_side = (m - 6) / 2; // side of the new vertex on this triangle
                let this_pos = (m - 6) % 2; // 0 or 1 (first or second node on the side of this triangle)
                let neigh_pos = 1 - this_pos; // 1 or 0 (second or first node on neighbor's side)

                // find neighboring triangle and the side of edge on the neighbor triangle
                let local_a = Tri6::EDGE_NODE_IDS[this_side][0];
                let local_b = Tri6::EDGE_NODE_IDS[this_side][1];
                let mut this_a = mesh.cells[this_cell_id].points[local_a];
                let mut this_b = mesh.cells[this_cell_id].points[local_b];
                if this_b < this_a {
                    let temp = this_a;
                    this_a = this_b;
                    this_b = temp;
                }
                let edge_shares = edges.get(&(this_a, this_b)).unwrap();
                let shares: Vec<_> = edge_shares.iter().filter(|v| v.0 != this_cell_id).collect();
                if shares.len() == 1 {
                    // n is the index of new point on neighbor's edge
                    let neigh_cell_id = shares[0].0;
                    let neigh_side = shares[0].1;
                    let n = 6 + neigh_side * 2 + neigh_pos;

                    // set new point of neighbor
                    mesh.cells[neigh_cell_id].points[n] = new_point_id;
                } // else => the edge is on the boundary (shared by one triangle only)
            }
        }

        // new points in the interior of the Tri15
        for m in 12..15 {
            // coordinates of new point
            pad.calc_coords(&mut x, &Tri15::NODE_REFERENCE_COORDS[m])?;

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
        mesh.cells[this_cell_id].kind = GeoKind::Tri15;
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use russell_chk::vec_approx_eq;

    use super::upgrade_tri6_to_tri15;
    use crate::mesh::{check_all, Cell, Mesh, Point};
    use crate::shapes::GeoKind;

    #[allow(unused_imports)]
    use crate::mesh::draw_mesh;

    const SAVE_FIGURE: bool = false;

    #[test]
    fn tri6_to_tri15_works_1() {
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

        upgrade_tri6_to_tri15(&mut mesh).unwrap();

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
