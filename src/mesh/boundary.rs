use super::{all_edges_2d, all_faces_3d, alloc_cell_shapes, CellId, Edge, EdgeKey, Face, FaceKey, Mesh, PointId};
use crate::{shapes::Shape, StrError};
use russell_lab::sort2;
use std::collections::HashMap;

/// Holds points, edges and faces on the boundaries of a mesh
pub struct Boundary {
    /// Set of edges on the boundaries
    ///
    /// Note:
    ///
    /// * In 2D, a boundary edge is such that it is shared by one 2D cell only (1D cells are ignored)
    /// * In 3D, a boundary edge belongs to a boundary face
    pub edges: HashMap<EdgeKey, Edge>,

    /// Set of faces on the boundaries
    ///
    /// Note: A boundary face is such that it is shared by one 3D cell only
    pub faces: HashMap<FaceKey, Face>,
}

impl Boundary {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh) -> Result<Self, StrError> {
        let shapes = alloc_cell_shapes(mesh)?;
        let (boundary_edges, boundary_faces) = match mesh.space_ndim {
            2 => {
                let edges = all_edges_2d(mesh, &shapes)?;
                let boundary_edges = Boundary::two_dim(mesh, &shapes, &edges)?;
                (boundary_edges, HashMap::new())
            }
            3 => {
                let faces = all_faces_3d(mesh, &shapes)?;
                Boundary::three_dim(mesh, &shapes, &faces)?
            }
            _ => panic!("space_ndim must be 2 or 3"),
        };
        Ok(Boundary {
            edges: boundary_edges,
            faces: boundary_faces,
        })
    }

    /// Finds boundary entities in 2D
    ///
    /// **Note:** Call this function after `all_edges_2d`.
    ///
    /// # Input
    ///
    /// * `mesh` -- the Mesh
    /// * `shapes` -- the shapes of cells (len == cells.len())
    /// * `edges` -- all edges (internal and boundary)
    ///
    /// # Output
    ///
    /// * Returns a map relating edge keys to Edge
    pub fn two_dim(
        mesh: &Mesh,
        shapes: &Vec<Shape>,
        edges: &HashMap<EdgeKey, Vec<(CellId, usize)>>,
    ) -> Result<HashMap<EdgeKey, Edge>, StrError> {
        if mesh.space_ndim != 2 {
            return Err("this function works in 2D only");
        }
        let mut boundary_edges: HashMap<EdgeKey, Edge> = HashMap::new();
        for (edge_key, shared_by) in edges {
            if shared_by.len() != 1 {
                continue; // skip internal edges (those shared by multiple cells)
            }
            let (cell_id, e) = shared_by[0];
            let cell = &mesh.cells[cell_id];
            let shape = &shapes[cell_id];
            boundary_edges.entry(*edge_key).or_insert(Edge {
                points: (0..shape.edge_nnode)
                    .into_iter()
                    .map(|i| cell.points[shape.edge_node_id(e, i)])
                    .collect(),
            });
        }
        Ok(boundary_edges)
    }

    /// Finds boundary entities in 3D
    ///
    /// **Note:** Call this function after `all_faces`.
    ///
    /// # Input
    ///
    /// * `mesh` -- the Mesh
    /// * `shapes` -- the shapes of cells (len == cells.len())
    /// * `faces` -- all faces (internal and boundary)
    ///
    /// # Output
    ///
    /// * Returns:
    ///     - a map relating edge keys to Edge
    ///     - a map relating face keys to Face
    pub fn three_dim(
        mesh: &Mesh,
        shapes: &Vec<Shape>,
        faces: &HashMap<FaceKey, Vec<(CellId, usize)>>,
    ) -> Result<(HashMap<EdgeKey, Edge>, HashMap<FaceKey, Face>), StrError> {
        if mesh.space_ndim != 3 {
            return Err("this function works in 3D only");
        }

        // sort face keys just so the next loop is deterministic
        let mut face_keys: Vec<_> = faces.keys().collect();
        face_keys.sort();

        // loop over all faces
        let mut boundary_edges: HashMap<EdgeKey, Edge> = HashMap::new();
        let mut boundary_faces: HashMap<FaceKey, Face> = HashMap::new();
        for face_key in face_keys {
            // skip internal faces (those shared by multiple cells)
            let shared_by = faces.get(face_key).unwrap();
            if shared_by.len() != 1 {
                continue;
            }

            // cell and face
            let (cell_id, f) = shared_by[0];
            let cell = &mesh.cells[cell_id];
            let shape = &shapes[cell_id];
            let face_nnode = shape.face_nnode;
            let face_shape = Shape::new(mesh.space_ndim, 2, face_nnode)?;

            // ids of points on faces
            let face_points: Vec<PointId> = (0..face_nnode)
                .into_iter()
                .map(|i| cell.points[shape.face_node_id(f, i)])
                .collect();

            // loop over all face edges
            for e in 0..face_shape.nedge {
                // define edge key (sorted point ids)
                let mut edge_key: EdgeKey = (
                    face_points[face_shape.edge_node_id(e, 0)],
                    face_points[face_shape.edge_node_id(e, 1)],
                );
                sort2(&mut edge_key);

                // insert boundary edge
                boundary_edges.entry(edge_key).or_insert(Edge {
                    points: (0..face_shape.edge_nnode)
                        .into_iter()
                        .map(|i| face_points[face_shape.edge_node_id(e, i)])
                        .collect(),
                });
            }

            // insert boundary face
            boundary_faces.entry(*face_key).or_insert(Face { points: face_points });
        }
        Ok((boundary_edges, boundary_faces))
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Boundary;
    use crate::mesh::Samples;
    use crate::StrError;

    #[test]
    fn boundary_2d_works() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let boundary = Boundary::new(&mesh)?;
        let mut edge_keys: Vec<_> = boundary.edges.keys().collect();
        edge_keys.sort();
        assert_eq!(edge_keys, [&(0, 1), &(0, 3), &(1, 4), &(2, 3), &(2, 5), &(4, 5)]);
        assert_eq!(boundary.edges.get(&(0, 1)).unwrap().points, &[1, 0]);
        assert_eq!(boundary.edges.get(&(0, 3)).unwrap().points, &[0, 3]);
        assert_eq!(boundary.edges.get(&(1, 4)).unwrap().points, &[4, 1]);
        assert_eq!(boundary.edges.get(&(2, 3)).unwrap().points, &[3, 2]);
        assert_eq!(boundary.edges.get(&(2, 5)).unwrap().points, &[2, 5]);
        assert_eq!(boundary.edges.get(&(4, 5)).unwrap().points, &[5, 4]);

        //           4---------3
        //           |         |
        //           |   [1]   |
        //           |         |
        //  0--------1---------2
        let mesh = Samples::mixed_shapes_2d();
        let boundary = Boundary::new(&mesh)?;
        let mut edge_keys: Vec<_> = boundary.edges.keys().collect();
        edge_keys.sort();
        assert_eq!(edge_keys, [&(1, 2), &(1, 4), &(2, 3), &(3, 4)]);
        assert_eq!(boundary.edges.get(&(1, 2)).unwrap().points, &[2, 1]);
        assert_eq!(boundary.edges.get(&(1, 4)).unwrap().points, &[1, 4]);
        assert_eq!(boundary.edges.get(&(2, 3)).unwrap().points, &[3, 2]);
        assert_eq!(boundary.edges.get(&(3, 4)).unwrap().points, &[4, 3]);
        Ok(())
    }

    #[test]
    fn boundary_3d_works() -> Result<(), StrError> {
        //      8-------------11
        //     /.             /|
        //    / .            / |
        //   /  .           /  |
        //  /   .          /   |       id = 1
        // 9-------------10    |       attribute_id = 2
        // |    .         |    |
        // |    4---------|----7
        // |   /.         |   /|
        // |  / .         |  / |
        // | /  .         | /  |
        // |/   .         |/   |
        // 5--------------6    |       id = 0
        // |    .         |    |       attribute_id = 1
        // |    0---------|----3
        // |   /          |   /
        // |  /           |  /
        // | /            | /
        // |/             |/
        // 1--------------2
        let mesh = Samples::two_cubes_vertical();
        let boundary = Boundary::new(&mesh)?;
        let mut edge_keys: Vec<_> = boundary.edges.keys().collect();
        let mut face_keys: Vec<_> = boundary.faces.keys().collect();
        edge_keys.sort();
        face_keys.sort();
        assert_eq!(
            edge_keys,
            [
                &(0, 1),
                &(0, 3),
                &(0, 4),
                &(1, 2),
                &(1, 5),
                &(2, 3),
                &(2, 6),
                &(3, 7),
                &(4, 5),
                &(4, 7),
                &(4, 8),
                &(5, 6),
                &(5, 9),
                &(6, 7),
                &(6, 10),
                &(7, 11),
                &(8, 9),
                &(8, 11),
                &(9, 10),
                &(10, 11),
            ]
        );
        assert_eq!(
            face_keys,
            [
                &(0, 1, 2, 3),
                &(0, 1, 4, 5),
                &(0, 3, 4, 7),
                &(1, 2, 5, 6),
                &(2, 3, 6, 7),
                &(4, 5, 8, 9),
                &(4, 7, 8, 11),
                &(5, 6, 9, 10),
                &(6, 7, 10, 11),
                &(8, 9, 10, 11),
            ]
        );
        assert_eq!(boundary.edges.get(&(0, 1)).unwrap().points, &[0, 1]);
        assert_eq!(boundary.edges.get(&(0, 3)).unwrap().points, &[3, 0]);
        assert_eq!(boundary.edges.get(&(0, 4)).unwrap().points, &[0, 4]);
        assert_eq!(boundary.edges.get(&(1, 2)).unwrap().points, &[1, 2]);
        assert_eq!(boundary.edges.get(&(1, 5)).unwrap().points, &[5, 1]);
        assert_eq!(boundary.edges.get(&(2, 3)).unwrap().points, &[2, 3]);
        assert_eq!(boundary.edges.get(&(2, 6)).unwrap().points, &[6, 2]);
        assert_eq!(boundary.edges.get(&(3, 7)).unwrap().points, &[3, 7]);
        assert_eq!(boundary.edges.get(&(4, 5)).unwrap().points, &[4, 5]);
        assert_eq!(boundary.edges.get(&(4, 7)).unwrap().points, &[7, 4]);
        assert_eq!(boundary.edges.get(&(4, 8)).unwrap().points, &[4, 8]);
        assert_eq!(boundary.edges.get(&(5, 6)).unwrap().points, &[5, 6]);
        assert_eq!(boundary.edges.get(&(5, 9)).unwrap().points, &[9, 5]);
        assert_eq!(boundary.edges.get(&(6, 7)).unwrap().points, &[6, 7]);
        assert_eq!(boundary.edges.get(&(6, 10)).unwrap().points, &[10, 6]);
        assert_eq!(boundary.edges.get(&(7, 11)).unwrap().points, &[7, 11]);
        assert_eq!(boundary.edges.get(&(8, 9)).unwrap().points, &[8, 9]);
        assert_eq!(boundary.edges.get(&(8, 11)).unwrap().points, &[11, 8]);
        assert_eq!(boundary.edges.get(&(9, 10)).unwrap().points, &[9, 10]);
        assert_eq!(boundary.edges.get(&(10, 11)).unwrap().points, &[10, 11]);
        assert_eq!(boundary.faces.get(&(0, 1, 2, 3)).unwrap().points, &[0, 3, 2, 1]);
        assert_eq!(boundary.faces.get(&(0, 1, 4, 5)).unwrap().points, &[0, 1, 5, 4]);
        assert_eq!(boundary.faces.get(&(0, 3, 4, 7)).unwrap().points, &[0, 4, 7, 3]);
        assert_eq!(boundary.faces.get(&(1, 2, 5, 6)).unwrap().points, &[1, 2, 6, 5]);
        assert_eq!(boundary.faces.get(&(2, 3, 6, 7)).unwrap().points, &[2, 3, 7, 6]);
        assert_eq!(boundary.faces.get(&(4, 5, 8, 9)).unwrap().points, &[4, 5, 9, 8]);
        assert_eq!(boundary.faces.get(&(4, 7, 8, 11)).unwrap().points, &[4, 8, 11, 7]);
        assert_eq!(boundary.faces.get(&(5, 6, 9, 10)).unwrap().points, &[5, 6, 10, 9]);
        assert_eq!(boundary.faces.get(&(6, 7, 10, 11)).unwrap().points, &[6, 7, 11, 10]);
        assert_eq!(boundary.faces.get(&(8, 9, 10, 11)).unwrap().points, &[8, 9, 10, 11]);

        //                       4------------7-----------10
        //                      /.           /|            |
        //                     / .          / |            |
        //                    /  .         /  |            |
        //                   /   .        /   |            |
        //                  5------------6    |            |
        //                  |    .       |`.  |            |
        //                  |    0-------|--`.3------------9
        //                  |   /        |   /`.          /
        //                  |  /         |  /   `.       /
        //                  | /          | /      `.    /
        //                  |/           |/         `. /
        //  12-----11-------1------------2------------8
        //
        let mesh = Samples::mixed_shapes_3d();
        let boundary = Boundary::new(&mesh)?;
        let mut edge_keys: Vec<_> = boundary.edges.keys().collect();
        let mut face_keys: Vec<_> = boundary.faces.keys().collect();
        edge_keys.sort();
        face_keys.sort();
        assert_eq!(
            edge_keys,
            [
                &(0, 1),
                &(0, 3),
                &(0, 4),
                &(1, 2),
                &(1, 5),
                &(2, 3),
                &(2, 6),
                &(2, 8),
                &(3, 6),
                &(3, 7),
                &(3, 8),
                &(4, 5),
                &(4, 7),
                &(5, 6),
                &(6, 7),
                &(6, 8),
            ]
        );
        assert_eq!(
            face_keys,
            [
                &(0, 1, 2, 3),
                &(0, 1, 4, 5),
                &(0, 3, 4, 7),
                &(1, 2, 5, 6),
                &(2, 3, 6, 7),
                &(2, 3, 6, 13),
                &(2, 3, 8, 13),
                &(2, 6, 8, 13),
                &(3, 6, 8, 13),
                &(4, 5, 6, 7),
            ]
        );
        assert_eq!(boundary.edges.get(&(0, 1)).unwrap().points, &[0, 1]);
        assert_eq!(boundary.edges.get(&(0, 3)).unwrap().points, &[3, 0]);
        assert_eq!(boundary.edges.get(&(0, 4)).unwrap().points, &[0, 4]);
        assert_eq!(boundary.edges.get(&(1, 2)).unwrap().points, &[1, 2]);
        assert_eq!(boundary.edges.get(&(1, 5)).unwrap().points, &[5, 1]);
        assert_eq!(boundary.edges.get(&(2, 3)).unwrap().points, &[2, 3]);
        assert_eq!(boundary.edges.get(&(2, 6)).unwrap().points, &[6, 2]);
        assert_eq!(boundary.edges.get(&(2, 8)).unwrap().points, &[2, 8]);
        assert_eq!(boundary.edges.get(&(3, 6)).unwrap().points, &[3, 6]);
        assert_eq!(boundary.edges.get(&(3, 7)).unwrap().points, &[3, 7]);
        assert_eq!(boundary.edges.get(&(3, 8)).unwrap().points, &[8, 3]);
        assert_eq!(boundary.edges.get(&(4, 5)).unwrap().points, &[4, 5]);
        assert_eq!(boundary.edges.get(&(4, 7)).unwrap().points, &[7, 4]);
        assert_eq!(boundary.edges.get(&(5, 6)).unwrap().points, &[5, 6]);
        assert_eq!(boundary.edges.get(&(6, 7)).unwrap().points, &[6, 7]);
        assert_eq!(boundary.edges.get(&(6, 8)).unwrap().points, &[6, 8]);
        assert_eq!(boundary.faces.get(&(0, 1, 2, 3)).unwrap().points, &[0, 3, 2, 1]);
        assert_eq!(boundary.faces.get(&(0, 1, 4, 5)).unwrap().points, &[0, 1, 5, 4]);
        assert_eq!(boundary.faces.get(&(0, 3, 4, 7)).unwrap().points, &[0, 4, 7, 3]);
        assert_eq!(boundary.faces.get(&(1, 2, 5, 6)).unwrap().points, &[1, 2, 6, 5]);
        assert_eq!(boundary.faces.get(&(2, 3, 6, 7)).unwrap().points, &[2, 3, 7, 6]);
        assert_eq!(boundary.faces.get(&(2, 3, 6, 13)).unwrap().points, &[2, 6, 3]);
        assert_eq!(boundary.faces.get(&(2, 3, 8, 13)).unwrap().points, &[2, 3, 8]);
        assert_eq!(boundary.faces.get(&(2, 6, 8, 13)).unwrap().points, &[2, 8, 6]);
        assert_eq!(boundary.faces.get(&(3, 6, 8, 13)).unwrap().points, &[8, 3, 6]);
        assert_eq!(boundary.faces.get(&(4, 5, 6, 7)).unwrap().points, &[4, 5, 6, 7]);
        Ok(())
    }
}
