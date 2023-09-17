use super::{get_neighbors_2d, Cell, Extract, Features, Mesh, Point};
use crate::prelude::PointId;
use crate::shapes::{GeoClass, GeoKind, Scratchpad};
use crate::StrError;
use russell_lab::Vector;
use std::collections::HashMap;

impl Mesh {
    /// Upgrades or downgrades a mesh with triangles or quadrilaterals
    ///
    /// # Notes
    ///
    /// 1. All cells must have the same GeoKind
    /// 2. Only [GeoClass::Tri] and [GeoClass::Qua] are allowed
    /// 3. The points will be completely renumbered
    /// 4. The corner tags will be replicated into the new mesh
    /// 5. The points at the middle of edges will inherit the tag of ONE corresponding middle point
    pub fn convert_2d(&self, target: GeoKind) -> Result<Mesh, StrError> {
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
        if self.ndim != 2 {
            return Err("ndim must be equal to 2");
        }

        // current number of cells
        let ncell = self.cells.len();
        if ncell < 1 {
            return Err("the conversion requires at least one cell");
        }

        // get info about first cell in the mesh
        let source = self.cells[0].kind;
        let class = source.class();
        let source_nnode = source.nnode();
        let target_nnode = target.nnode();
        if target.class() != class {
            return Err("target class must equal the GeoClass of current cells");
        }
        if target.class() != GeoClass::Tri && target.class() != GeoClass::Qua {
            return Err("target GeoClass must be Tri or Qua");
        }

        // info about the edges and corners
        let source_edge_nnode = source.edge_nnode();
        let target_edge_nnode = target.edge_nnode();
        let nedge = target.nedge();
        let ncorner = if class == GeoClass::Tri { 3 } else { 4 };
        let mut edge_point_markers = vec![0; nedge];

        // maps target local node id to edge number (to set markers of points on edges)
        const INTERIOR_OR_CORNER: usize = usize::MAX;
        let mut target_node_to_edge = vec![INTERIOR_OR_CORNER; target_nnode];
        for e in 0..nedge {
            for i in 2..target_edge_nnode {
                let m = target.edge_node_id(e, i);
                target_node_to_edge[m] = e;
            }
        }

        // find neighbors in the original mesh
        let source_features = Features::new(self, Extract::All);
        let source_edges = source_features.all_2d_edges.unwrap();

        // flag indicating that a connectivity point has not been set
        const UNSET: usize = usize::MAX;
        assert!(ncell * target_nnode < UNSET);

        // zeroed new cell
        let zero_cell = Cell {
            id: 0,
            attribute: 0,
            kind: target,
            points: vec![UNSET; target_nnode],
        };

        // allocate destination mesh
        let mut dest = Mesh {
            ndim: self.ndim,
            points: Vec::new(),
            cells: vec![zero_cell; ncell],
        };

        // scratchpad for point interpolation (based on the original mesh)
        let mut pad = Scratchpad::new(self.ndim, source).unwrap();

        // coordinates of new points
        let mut x = Vector::new(self.ndim);

        // maps old point id to new point id
        let mut corners: HashMap<PointId, PointId> = HashMap::new();

        // upgrade or downgrade mesh
        for cell_id in 0..ncell {
            // check GeoKind
            if self.cells[cell_id].kind != source {
                return Err("all cells must have the same GeoKind");
            }

            // set the coordinates matrix of the original mesh
            for m in 0..source_nnode {
                let p = self.cells[cell_id].points[m];
                for j in 0..self.ndim {
                    pad.set_xx(m, j, self.points[p].coords[j]);
                }
            }

            // find markers of points on edges
            if source_edge_nnode > 2 {
                for e in 0..nedge {
                    for i in 2..source_edge_nnode {
                        let m = source.edge_node_id(e, i);
                        let p = self.cells[cell_id].points[m];
                        edge_point_markers[e] = self.points[p].marker;
                    }
                }
            }

            // set the new cell data
            dest.cells[cell_id].id = cell_id;
            dest.cells[cell_id].attribute = self.cells[cell_id].attribute;

            // handle corner nodes
            for m in 0..ncorner {
                let old_point_id = self.cells[cell_id].points[m];
                if let Some(new_point_id) = corners.get(&old_point_id) {
                    dest.cells[cell_id].points[m] = *new_point_id;
                } else {
                    let new_point_id = dest.points.len();
                    dest.points.push(Point {
                        id: new_point_id,
                        marker: self.points[old_point_id].marker,
                        coords: self.points[old_point_id].coords.clone(),
                    });
                    dest.cells[cell_id].points[m] = new_point_id;
                    corners.insert(old_point_id, new_point_id);
                }
            }

            // consult neighbors to see if there are points (in the middle of the edge) set already
            if target_edge_nnode > 2 {
                let neighbors = get_neighbors_2d(self, &source_edges, cell_id);
                for (e, neigh_cell_id, neigh_e) in neighbors {
                    // only deal with the centre edge points, not the corner ones (start at 2)
                    for i in 2..target_edge_nnode {
                        let n = target.edge_node_id(neigh_e, i); // outward normal
                        let p = dest.cells[neigh_cell_id].points[n];
                        if p != UNSET {
                            let m = target.edge_node_id_inward(e, i); // inward normal
                            dest.cells[cell_id].points[m] = p;
                        }
                    }
                }
            }

            // add new points (skip corner nodes)
            for m in ncorner..target_nnode {
                if dest.cells[cell_id].points[m] == UNSET {
                    let e = target_node_to_edge[m];
                    let marker = if e == INTERIOR_OR_CORNER {
                        0
                    } else {
                        edge_point_markers[e]
                    };
                    let new_point_id = dest.points.len();
                    pad.calc_coords(&mut x, target.reference_coords(m)).unwrap();
                    dest.points.push(Point {
                        id: new_point_id,
                        marker,
                        coords: x.as_data().clone(),
                    });
                    dest.cells[cell_id].points[m] = new_point_id;
                }
            }
        }
        Ok(dest)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::mesh::{Cell, Figure, Mesh, Point, Samples};
    use crate::shapes::GeoKind;
    use russell_chk::vec_approx_eq;

    const SAVE_FIGURE: bool = false;

    fn draw(mesh: &Mesh, larger: bool, filename: &str) {
        let mut fig = Figure::new();
        fig.param_cell_ids = true;
        fig.param_point_ids = true;
        if larger {
            fig.param_figure_size = Some((600.0, 600.0));
        }
        mesh.draw(Some(fig), filename).unwrap();
    }

    #[test]
    fn convert_mesh_2d_captures_errors() {
        #[rustfmt::skip]
        let mesh = Mesh {
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
        assert_eq!(mesh.convert_2d(GeoKind::Tri15).err(), Some("ndim must be equal to 2"));

        let mesh = Mesh {
            ndim: 2,
            points: vec![],
            cells: vec![],
        };
        assert_eq!(
            mesh.convert_2d(GeoKind::Tri15).err(),
            Some("the conversion requires at least one cell")
        );

        #[rustfmt::skip]
        let mesh = Mesh {
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
            mesh.convert_2d(GeoKind::Qua8).err(),
            Some("target class must equal the GeoClass of current cells")
        );

        #[rustfmt::skip]
        let mesh =Mesh {
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
            mesh.convert_2d(GeoKind::Lin3).err(),
            Some("target GeoClass must be Tri or Qua")
        );

        #[rustfmt::skip]
        let mesh = Mesh {
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
            mesh.convert_2d(GeoKind::Tri6).err(),
            Some("all cells must have the same GeoKind")
        );
    }

    #[test]
    fn convert_tri6_to_tri15_works() {
        #[rustfmt::skip]
        let mesh = Mesh {
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
        mesh.check_all().unwrap();

        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_convert_2d_tri6_to_tri15_before.svg");
        }

        let res = mesh.convert_2d(GeoKind::Tri15).unwrap();

        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_convert_2d_tri6_to_tri15_after.svg");
        }

        res.check_all().unwrap();
        res.check_overlapping_points(0.2).unwrap();

        assert_eq!(res.points.len(), 2 * 15 - 5);
        assert_eq!(res.cells[0].points, (0..15).collect::<Vec<_>>());
        assert_eq!(
            res.cells[1].points,
            &[2, 15, 0, 16, 17, 5, 18, 19, 20, 21, 11, 10, 22, 23, 24],
        );

        vec_approx_eq(&res.points[0].coords, &[4.0, 0.0], 1e-15);
        vec_approx_eq(&res.points[1].coords, &[0.0, 4.0], 1e-15);
        vec_approx_eq(&res.points[2].coords, &[0.0, 0.0], 1e-15);
        vec_approx_eq(&res.points[3].coords, &[2.0, 2.0], 1e-15);
        vec_approx_eq(&res.points[4].coords, &[0.0, 2.0], 1e-15);
        vec_approx_eq(&res.points[5].coords, &[2.0, 0.0], 1e-15);
        vec_approx_eq(&res.points[6].coords, &[3.0, 1.0], 1e-15);
        vec_approx_eq(&res.points[7].coords, &[1.0, 3.0], 1e-15);
        vec_approx_eq(&res.points[8].coords, &[0.0, 3.0], 1e-15);
        vec_approx_eq(&res.points[9].coords, &[0.0, 1.0], 1e-15);
        vec_approx_eq(&res.points[10].coords, &[1.0, 0.0], 1e-15);
        vec_approx_eq(&res.points[11].coords, &[3.0, 0.0], 1e-15);
        vec_approx_eq(&res.points[12].coords, &[2.0, 1.0], 1e-15);
        vec_approx_eq(&res.points[13].coords, &[1.0, 2.0], 1e-15);
        vec_approx_eq(&res.points[14].coords, &[1.0, 1.0], 1e-15);
        vec_approx_eq(&res.points[15].coords, &[0.0, -4.0], 1e-15);
        vec_approx_eq(&res.points[16].coords, &[0.0, -2.0], 1e-15);
        vec_approx_eq(&res.points[17].coords, &[2.0, -2.0], 1e-15);
        vec_approx_eq(&res.points[18].coords, &[0.0, -1.0], 1e-15);
        vec_approx_eq(&res.points[19].coords, &[0.0, -3.0], 1e-15);
        vec_approx_eq(&res.points[20].coords, &[1.0, -3.0], 1e-15);
        vec_approx_eq(&res.points[21].coords, &[3.0, -1.0], 1e-15);
        vec_approx_eq(&res.points[22].coords, &[1.0, -1.0], 1e-15);
        vec_approx_eq(&res.points[23].coords, &[1.0, -2.0], 1e-15);
        vec_approx_eq(&res.points[24].coords, &[2.0, -1.0], 1e-15);

        assert_eq!(res.points[0].marker, 0);
        assert_eq!(res.points[1].marker, 0);
        assert_eq!(res.points[2].marker, 0);
        assert_eq!(res.points[3].marker, -20);
        assert_eq!(res.points[4].marker, -10);
        assert_eq!(res.points[5].marker, 0);
        assert_eq!(res.points[6].marker, -20);
        assert_eq!(res.points[7].marker, -20);
        assert_eq!(res.points[8].marker, -10);
        assert_eq!(res.points[9].marker, -10);
        assert_eq!(res.points[10].marker, 0);
        assert_eq!(res.points[11].marker, 0);
        assert_eq!(res.points[12].marker, 0);
        assert_eq!(res.points[13].marker, 0);
        assert_eq!(res.points[14].marker, 0);
        assert_eq!(res.points[15].marker, 0);
        assert_eq!(res.points[16].marker, -10);
        assert_eq!(res.points[17].marker, -30);
        assert_eq!(res.points[18].marker, -10);
        assert_eq!(res.points[19].marker, -10);
        assert_eq!(res.points[20].marker, -30);
        assert_eq!(res.points[21].marker, -30);
        assert_eq!(res.points[22].marker, 0);
        assert_eq!(res.points[23].marker, 0);
        assert_eq!(res.points[24].marker, 0);
    }

    #[test]
    fn convert_tri6_to_tri10_works() {
        #[rustfmt::skip]
        let mesh = Mesh {
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
        mesh.check_all().unwrap();

        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_convert_2d_tri6_to_tri10_before.svg");
        }

        let res = mesh.convert_2d(GeoKind::Tri10).unwrap();

        if SAVE_FIGURE {
            draw(&res, false, "/tmp/gemlab/test_convert_2d_tri6_to_tri10_after.svg");
        }

        res.check_all().unwrap();
        res.check_overlapping_points(0.2).unwrap();

        assert_eq!(res.points.len(), 16);

        let ft = 4.0 / 3.0;
        let et = 8.0 / 3.0;
        vec_approx_eq(&res.points[0].coords, &[4.0, 0.0], 1e-15);
        vec_approx_eq(&res.points[1].coords, &[0.0, 4.0], 1e-15);
        vec_approx_eq(&res.points[2].coords, &[0.0, 0.0], 1e-15);
        vec_approx_eq(&res.points[3].coords, &[et, ft], 1e-15);
        vec_approx_eq(&res.points[4].coords, &[0.0, et], 1e-15);
        vec_approx_eq(&res.points[5].coords, &[ft, 0.0], 1e-15);
        vec_approx_eq(&res.points[6].coords, &[ft, et], 1e-15);
        vec_approx_eq(&res.points[7].coords, &[0.0, ft], 1e-15);
        vec_approx_eq(&res.points[8].coords, &[et, 0.0], 1e-15);
        vec_approx_eq(&res.points[9].coords, &[ft, ft], 1e-15);
        vec_approx_eq(&res.points[10].coords, &[0.0, -4.0], 1e-15);
        vec_approx_eq(&res.points[11].coords, &[0.0, -ft], 1e-15);
        vec_approx_eq(&res.points[12].coords, &[ft, -et], 1e-15);
        vec_approx_eq(&res.points[13].coords, &[0.0, -et], 1e-15);
        vec_approx_eq(&res.points[14].coords, &[et, -ft], 1e-15);
        vec_approx_eq(&res.points[15].coords, &[ft, -ft], 1e-15);

        assert_eq!(res.points[0].marker, 0);
        assert_eq!(res.points[1].marker, 0);
        assert_eq!(res.points[2].marker, 0);
        assert_eq!(res.points[3].marker, -20);
        assert_eq!(res.points[4].marker, -10);
        assert_eq!(res.points[5].marker, 0);
        assert_eq!(res.points[6].marker, -20);
        assert_eq!(res.points[7].marker, -10);
        assert_eq!(res.points[8].marker, 0);
        assert_eq!(res.points[9].marker, 0);
        assert_eq!(res.points[10].marker, 0);
        assert_eq!(res.points[11].marker, -10);
        assert_eq!(res.points[12].marker, -30);
        assert_eq!(res.points[13].marker, -10);
        assert_eq!(res.points[14].marker, -30);
        assert_eq!(res.points[15].marker, 0);
    }

    #[test]
    fn convert_tri6_to_tri3_works() {
        #[rustfmt::skip]
        let mesh = Mesh {
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
        mesh.check_all().unwrap();

        let res = mesh.convert_2d(GeoKind::Tri3).unwrap();

        if SAVE_FIGURE {
            draw(&res, false, "/tmp/gemlab/test_convert_2d_tri6_to_tri3_after.svg");
        }

        res.check_all().unwrap();
        res.check_overlapping_points(0.1).unwrap();

        assert_eq!(res.points.len(), 4);
        assert_eq!(res.cells[0].points, &[0, 1, 2]);
        assert_eq!(res.cells[1].points, &[2, 3, 0]);
    }

    #[test]
    fn convert_tri3_to_tri6_works() {
        let mesh = Samples::two_tri3().clone();
        let res = mesh.convert_2d(GeoKind::Tri6).unwrap();
        if SAVE_FIGURE {
            draw(&res, false, "/tmp/gemlab/test_convert_2d_tri3_to_tri6_after.svg");
        }
        res.check_all().unwrap();
        res.check_overlapping_points(0.1).unwrap();
        assert_eq!(res.points.len(), 9);
        assert_eq!(res.cells[0].points, (0..6).collect::<Vec<_>>());
        assert_eq!(res.cells[1].points, &[6, 2, 1, 7, 4, 8]);
    }

    #[test]
    fn convert_four_tri3_to_tri6_works() {
        let mesh = Samples::four_tri3().clone();
        let res = mesh.convert_2d(GeoKind::Tri6).unwrap();
        if SAVE_FIGURE {
            draw(&res, false, "/tmp/gemlab/test_convert_2d_four_tri3_to_tri6_after.svg");
        }
        res.check_all().unwrap();
        res.check_overlapping_points(0.1).unwrap();
        assert_eq!(res.points.len(), 13);
        assert_eq!(res.cells[0].points, (0..6).collect::<Vec<_>>());
        assert_eq!(res.cells[1].points, &[0, 6, 7, 8, 9, 10]);
        assert_eq!(res.cells[2].points, &[0, 2, 6, 5, 11, 8]);
        assert_eq!(res.cells[3].points, &[0, 7, 1, 10, 12, 3]);
    }

    #[test]
    fn convert_tri3_to_tri10_works() {
        let mesh = Samples::two_tri3().clone();
        let res = mesh.convert_2d(GeoKind::Tri10).unwrap();
        if SAVE_FIGURE {
            draw(&res, false, "/tmp/gemlab/test_convert_2d_tri3_to_tri10_after.svg");
        }
        res.check_all().unwrap();
        res.check_overlapping_points(0.1).unwrap();
        assert_eq!(res.points.len(), 16);
        assert_eq!(res.cells[0].points, (0..10).collect::<Vec<_>>());
        assert_eq!(res.cells[1].points, &[10, 2, 1, 11, 7, 12, 13, 4, 14, 15]);
    }

    #[test]
    fn convert_tri6_arrow_to_tri10_works() {
        let mesh = Samples::three_tri6_arrow().clone();
        let res = mesh.convert_2d(GeoKind::Tri10).unwrap();
        if SAVE_FIGURE {
            draw(&res, false, "/tmp/gemlab/test_convert_2d_tri6_arrow_to_tri10_after.svg");
        }
        res.check_all().unwrap();
        res.check_overlapping_points(0.1).unwrap();
        assert_eq!(res.points.len(), 19);
        assert_eq!(res.cells[0].points, (0..10).collect::<Vec<_>>());
        assert_eq!(res.cells[1].points, &[2, 1, 10, 7, 11, 12, 4, 13, 14, 15]);
        assert_eq!(res.cells[2].points, &[0, 2, 10, 8, 14, 16, 5, 12, 17, 18]);
        assert_eq!(res.points[0].marker, -1);
        assert_eq!(res.points[1].marker, -2);
        assert_eq!(res.points[2].marker, -3);
        assert_eq!(res.points[3].marker, -10);
        assert_eq!(res.points[4].marker, -5);
        assert_eq!(res.points[5].marker, -6);
        assert_eq!(res.points[6].marker, -10);
        assert_eq!(res.points[7].marker, -5);
        assert_eq!(res.points[8].marker, -6);
        assert_eq!(res.points[9].marker, 0);
        assert_eq!(res.points[10].marker, -7);
        assert_eq!(res.points[11].marker, -20);
        assert_eq!(res.points[12].marker, 0);
        assert_eq!(res.points[13].marker, -20);
        assert_eq!(res.points[14].marker, 0);
        assert_eq!(res.points[15].marker, 0);
        assert_eq!(res.points[16].marker, -30);
        assert_eq!(res.points[17].marker, -30);
        assert_eq!(res.points[18].marker, 0);
    }

    #[test]
    fn convert_qua4_to_qua8_works() {
        let mesh = Samples::two_qua4().clone();
        let res = mesh.convert_2d(GeoKind::Qua8).unwrap();
        if SAVE_FIGURE {
            draw(&res, false, "/tmp/gemlab/test_convert_2d_qua4_to_qua8_after.svg");
        }
        res.check_all().unwrap();
        res.check_overlapping_points(0.2).unwrap();
        assert_eq!(res.points.len(), 13);
        assert_eq!(res.cells[0].points, (0..8).collect::<Vec<_>>());
        assert_eq!(res.cells[1].points, &[1, 8, 9, 2, 10, 11, 12, 5]);
    }

    #[test]
    fn convert_qua12_to_qua16_works() {
        let mesh = Samples::block_2d_four_qua12().clone();
        let res = mesh.convert_2d(GeoKind::Qua16).unwrap();
        if SAVE_FIGURE {
            draw(&res, false, "/tmp/gemlab/test_convert_2d_qua12_to_qua16_after.svg");
        }
        res.check_all().unwrap();
        res.check_overlapping_points(0.2).unwrap();
        assert_eq!(res.points.len(), 33 + 4 * 4);
        assert_eq!(res.cells[0].points, (0..16).collect::<Vec<_>>());
        assert_eq!(
            res.cells[1].points,
            [1, 16, 17, 2, 18, 19, 20, 9, 21, 22, 23, 5, 24, 25, 26, 27]
        );
        assert_eq!(
            res.cells[2].points,
            [3, 2, 28, 29, 10, 30, 31, 32, 6, 33, 34, 35, 36, 37, 38, 39]
        );
        assert_eq!(
            res.cells[3].points,
            [2, 17, 40, 28, 23, 41, 42, 33, 20, 43, 44, 30, 45, 46, 47, 48]
        );
    }

    #[test]
    fn convert_qua12_to_qua17_works() {
        let mesh = Samples::block_2d_four_qua12().clone();
        let res = mesh.convert_2d(GeoKind::Qua17).unwrap();
        if SAVE_FIGURE {
            draw(&res, false, "/tmp/gemlab/test_convert_2d_qua12_to_qua17_after.svg");
        }
        res.check_all().unwrap();
        res.check_overlapping_points(0.2).unwrap();
        assert_eq!(res.points.len(), 3 * 9 + 6 * 3 + 4);
        assert_eq!(res.cells[0].points, (0..17).collect::<Vec<_>>());
        assert_eq!(
            res.cells[1].points,
            [1, 17, 18, 2, 19, 20, 21, 5, 22, 23, 24, 25, 26, 27, 28, 12, 11]
        );
        assert_eq!(
            res.cells[2].points,
            [3, 2, 29, 30, 6, 31, 32, 33, 34, 14, 13, 35, 36, 37, 38, 39, 40]
        );
        assert_eq!(
            res.cells[3].points,
            [2, 18, 41, 29, 21, 42, 43, 31, 44, 28, 27, 45, 46, 47, 48, 36, 35]
        );
        let correct = Samples::block_2d_four_qua17();
        for i in 0..correct.cells.len() {
            assert_eq!(&res.cells[i].points, &correct.cells[i].points);
        }
        let sf = 3.0 / 4.0;
        for i in 0..correct.points.len() {
            let scaled = &[sf * correct.points[i].coords[0], sf * correct.points[i].coords[1]];
            vec_approx_eq(&res.points[i].coords, scaled, 1e-15);
        }
    }

    #[test]
    fn convert_qua17_to_qua4_works() {
        let mesh = Samples::block_2d_four_qua17().clone();
        let res = mesh.convert_2d(GeoKind::Qua4).unwrap();
        if SAVE_FIGURE {
            draw(&res, false, "/tmp/gemlab/test_convert_2d_qua17_to_qua4_after.svg");
        }
        res.check_all().unwrap();
        res.check_overlapping_points(0.2).unwrap();
        assert_eq!(res.points.len(), 9);
        let correct = Samples::block_2d_four_qua4();
        for i in 0..correct.cells.len() {
            assert_eq!(&res.cells[i].points, &correct.cells[i].points);
        }
        let sf = 4.0 / 2.0;
        for i in 0..correct.points.len() {
            let scaled = &[sf * correct.points[i].coords[0], sf * correct.points[i].coords[1]];
            vec_approx_eq(&res.points[i].coords, scaled, 1e-15);
        }
    }
}
