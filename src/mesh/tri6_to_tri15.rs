#![allow(unused)]

use super::{Extract, Features, Mesh, Point};
use crate::shapes::{GeoKind, Scratchpad, Tri15, Tri6};
use crate::StrError;
use russell_lab::Vector;

/// Converts a mesh with only Tri6 cells to a mesh with Tri15 cells
pub fn tri6_to_tri15(mesh: &mut Mesh) -> Result<(), StrError> {
    // current number of points and cells
    let npoint = mesh.points.len();
    let ncell = mesh.cells.len();

    // 1 + max possible number of points (neglecting overlap)
    // this number will indicate that a point has not been set yet
    // if there is only one Tri6 cell, not_set will be equal to 16
    let not_set = 1 + npoint + ncell * 9;

    // counter for the next point id
    // let mut next_point_id = npoint;

    // expand connectivity array
    for i in 0..ncell {
        if mesh.cells[i].kind != GeoKind::Tri6 {
            return Err("only Tri6 cells are allowed in the conversion");
        }
        mesh.cells[i].points.resize(15, not_set);
    }

    // scratchpad for point interpolation
    const NDIM: usize = 2;
    let mut pad = Scratchpad::new(NDIM, GeoKind::Tri6)?;

    // natural coordinates of new points
    // let mut ksi = vec![0.0; NDIM];

    // real coordinates of new points
    let mut x = Vector::new(NDIM);

    // natural r coord of new nodes (6->14)
    // Vec_t R(9);
    // natural s coord of new nodes (6->14)
    // Vec_t S(9);
    // let rr =vec![1.0/4.0 , 3.0/4.0 , 3.0/4.0 , 1.0/4.0 , 0.0     , 0.0     , 1.0/4.0 , 1.0/2.0 , 1.0/4.0];
    // let ss =vec![0.0     , 0.0     , 1.0/4.0 , 3.0/4.0 , 3.0/4.0 , 1.0/4.0 , 1.0/4.0 , 1.0/4.0 , 1.0/2.0];

    const NEW_NODE_TO_EDGE: [usize; 15] = [
        0, //  0 (old node => ignored)
        0, //  1 (old node => ignored)
        0, //  2 (old node => ignored)
        0, //  3 (old node => ignored)
        0, //  4 (old node => ignored)
        0, //  5 (old node => ignored)
        0, //  6
        0, //  7
        1, //  8
        1, //  9
        2, // 10
        2, // 11
        0, // 12 (interior => ignored)
        0, // 13 (interior => ignored)
        0, // 14 (interior => ignored)
    ];

    // find neighbors
    let features = Features::new(mesh, Extract::All);
    let edges = match features.all_2d_edges {
        Some(e) => e,
        None => return Err("cannot extract 2d edges"),
    };

    // perform the conversion
    for i in 0..ncell {
        // set transformation matrix
        for m in 0..6 {
            let p = mesh.cells[i].points[m];
            let x = mesh.points[p].coords[0];
            let y = mesh.points[p].coords[1];
            pad.set_xx(m, 0, x);
            pad.set_xx(m, 1, y);
        }

        // points on the edge of the Tri15
        for m in 6..12 {
            let p = mesh.cells[i].points[m];
            if p == not_set {
                // new point
                pad.calc_coords(&mut x, &Tri15::NODE_REFERENCE_COORDS[m])?;
                let new_point_id = mesh.points.len();
                mesh.points.push(Point {
                    id: new_point_id,
                    marker: 0,
                    coords: x.as_data().clone(),
                });

                // set cell vertices
                mesh.cells[i].points[m] = new_point_id;

                // find neighboring triangle
                let local_edge_id = NEW_NODE_TO_EDGE[m];
                let tri6_local_node_a = Tri6::EDGE_NODE_IDS[local_edge_id][0];
                let tri6_local_node_b = Tri6::EDGE_NODE_IDS[local_edge_id][1];
                let a = mesh.cells[i].points[tri6_local_node_a];
                let b = mesh.cells[i].points[tri6_local_node_b];
                // let edge = features.get_edge(a, b);
                let pairs = edges.get(&(a, b)).unwrap();
                for (cell_id, e) in pairs {
                    if *cell_id != i {
                        // it's not me
                        println!("{} {}", cell_id, e)
                    }
                }

                // set neighbor cell vertices
                // TODO
            }
        }
    }

    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::tri6_to_tri15;
    use crate::mesh::{check_all, Cell, Mesh, Point};
    use crate::shapes::GeoKind;

    #[allow(unused_imports)]
    use crate::mesh::draw_mesh;

    const SAVE_FIGURE: bool = true;

    #[test]
    fn tri6_to_tri15_works_single_cell() {
        #[rustfmt::skip]
        let mut mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, marker: 0, coords: vec![4.0, 0.0] },
                Point { id: 2, marker: 0, coords: vec![0.0, 4.0] },
                Point { id: 3, marker: 0, coords: vec![2.0, 0.0] },
                Point { id: 4, marker: 0, coords: vec![2.0, 2.0] },
                Point { id: 5, marker: 0, coords: vec![0.0, 2.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Tri6, points: vec![0, 1, 2, 3, 4, 5] },
            ],
        };
        #[rustfmt::skip]
        let mut target = Mesh {
            ndim: 2,
            points: vec![
                Point { id:  0, marker: 0, coords: vec![0.0, 0.0] },
                Point { id:  1, marker: 0, coords: vec![4.0, 0.0] },
                Point { id:  2, marker: 0, coords: vec![0.0, 4.0] },
                Point { id:  3, marker: 0, coords: vec![2.0, 0.0] },
                Point { id:  4, marker: 0, coords: vec![2.0, 2.0] },
                Point { id:  5, marker: 0, coords: vec![0.0, 2.0] },
                Point { id:  6, marker: 0, coords: vec![1.0, 0.0] },
                Point { id:  7, marker: 0, coords: vec![3.0, 0.0] },
                Point { id:  8, marker: 0, coords: vec![3.0, 1.0] },
                Point { id:  9, marker: 0, coords: vec![1.0, 3.0] },
                Point { id: 10, marker: 0, coords: vec![0.0, 3.0] },
                Point { id: 11, marker: 0, coords: vec![0.0, 1.0] },
                Point { id: 12, marker: 0, coords: vec![1.0, 1.0] },
                Point { id: 13, marker: 0, coords: vec![2.0, 1.0] },
                Point { id: 14, marker: 0, coords: vec![1.0, 2.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Tri15, points: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14] },
            ],
        };
        check_all(&mesh).unwrap();
        check_all(&target).unwrap();

        if SAVE_FIGURE {
            draw_mesh(&mesh, true, true, false, "/tmp/gemlab/test_tri6_to_tri15.svg").unwrap();
        }

        tri6_to_tri15(&mut mesh).unwrap();
    }
}
