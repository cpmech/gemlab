#![allow(dead_code)]
#![allow(unused_variables)]

use crate::{GridSearch, Mesh};

pub struct Generator {
    mesh: Mesh,
}

impl Generator {
    pub fn from_mesh(mesh: Mesh) -> Generator {
        Generator { mesh }
    }

    pub fn merge(&mut self, right: &Mesh) -> Mesh {
        let ndim = self.mesh.ndim;
        let mesh = Mesh::new(ndim).unwrap();

        const NDIV: usize = 2;
        const GRID_NDIV: usize = 20;
        const GRID_MIN: f64 = -1.0;
        const GRID_MAX: f64 = 1.0;
        const TOLERANCE: f64 = 1e-4;

        let mut grid_bry = GridSearch::new(
            &vec![GRID_NDIV; ndim],
            &vec![GRID_MIN; ndim],
            &vec![GRID_MAX; ndim],
            &vec![TOLERANCE; ndim],
        )
        .unwrap();

        for (edge_key, edge) in &self.mesh.boundary_edges {
            for i in 0..2 {
                let point_id = edge.point_ids[i];
                let x = &self.mesh.points[point_id].coords;
                let existing_point_id = grid_bry.maybe_insert(point_id, x).unwrap();
                if point_id != existing_point_id {
                    panic!("problem with duplicate points");
                }
            }
        }

        println!("{}", grid_bry);

        for (edge_key, edge) in &right.boundary_edges {
            let mut npoint_found = 0;
            for i in 0..2 {
                let point_id = edge.point_ids[i];
                let x = &right.points[point_id].coords;
                if let Some(left_point_id) = grid_bry.find(x).unwrap() {
                    npoint_found += 1;
                }
            }
            if npoint_found == 2 {
                // existing edge
            }
        }

        mesh
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Cell, Edge, Point};

    #[test]
    fn merge_works() {
        let left = gen_mesh_left();
        // let right = gen_mesh_right();
        // let mut generator = Generator::from_mesh(left);
        // generator.merge(&right);
        left.write_serialized("/tmp/gemlab/mesh.out").unwrap();
    }

    #[rustfmt::skip]
    fn gen_mesh_left() -> Mesh {
        let mut mesh = Mesh::new(2).unwrap();

        mesh.points = vec![
            Point { id: 0, group: 1, coords: vec![0.0,0.0], shared_by_cell_ids: vec![0] },
            Point { id: 1, group: 1, coords: vec![1.0,0.0], shared_by_cell_ids: vec![0] },
            Point { id: 2, group: 1, coords: vec![1.0,1.0], shared_by_cell_ids: vec![0] },
            Point { id: 3, group: 1, coords: vec![0.0,1.0], shared_by_cell_ids: vec![0] },
        ];

        mesh.cells = vec![Cell { id: 0, group: 1, point_ids: vec![0,1,2,3], boundary_edge_ids: vec![0,1,2,3], boundary_face_ids: Vec::new() }];

        mesh.boundary_points.insert(0, true);
        mesh.boundary_points.insert(1, true);
        mesh.boundary_points.insert(2, true);
        mesh.boundary_points.insert(3, true);

        mesh.boundary_edges.insert((0,1), Edge { id: 0, group: 1, point_ids: vec![0,1], shared_by_cell_ids: vec![0] });
        mesh.boundary_edges.insert((1,2), Edge { id: 1, group: 1, point_ids: vec![1,2], shared_by_cell_ids: vec![0] });
        mesh.boundary_edges.insert((2,3), Edge { id: 2, group: 1, point_ids: vec![2,3], shared_by_cell_ids: vec![0] });
        mesh.boundary_edges.insert((0,3), Edge { id: 3, group: 1, point_ids: vec![3,0], shared_by_cell_ids: vec![0] });

        mesh
    }

    #[rustfmt::skip]
    fn gen_mesh_right() -> Mesh {
        let mut mesh = Mesh::new(2).unwrap();

        mesh.points = vec![
            Point { id: 0, group: 2, coords: vec![1.0,0.0], shared_by_cell_ids: vec![0] },
            Point { id: 1, group: 2, coords: vec![2.0,0.0], shared_by_cell_ids: vec![0] },
            Point { id: 2, group: 2, coords: vec![2.0,1.0], shared_by_cell_ids: vec![0] },
            Point { id: 3, group: 2, coords: vec![1.0,1.0], shared_by_cell_ids: vec![0] },
        ];

        mesh.cells = vec![Cell { id: 0, group: 1, point_ids: vec![0,1,2,3], boundary_edge_ids: vec![0,1,2,3], boundary_face_ids: Vec::new() }];

        mesh.boundary_points.insert(0, true);
        mesh.boundary_points.insert(1, true);
        mesh.boundary_points.insert(2, true);
        mesh.boundary_points.insert(3, true);

        mesh.boundary_edges.insert((0,1), Edge { id: 0, group: 1, point_ids: vec![0,1], shared_by_cell_ids: vec![0] });
        mesh.boundary_edges.insert((1,2), Edge { id: 1, group: 1, point_ids: vec![1,2], shared_by_cell_ids: vec![0] });
        mesh.boundary_edges.insert((2,3), Edge { id: 2, group: 1, point_ids: vec![2,3], shared_by_cell_ids: vec![0] });
        mesh.boundary_edges.insert((0,3), Edge { id: 3, group: 1, point_ids: vec![3,0], shared_by_cell_ids: vec![0] });

        mesh
    }
}
