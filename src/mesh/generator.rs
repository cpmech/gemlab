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
                let point_id = edge.points[i];
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
                let point_id = edge.points[i];
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
    // use super::*;

    #[test]
    fn merge_works() {
        // todo
    }
}
