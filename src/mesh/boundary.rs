use super::{all_edges_2d, all_faces_3d, alloc_cell_shapes, CellId, Edge, EdgeKey, Face, FaceKey, Mesh, PointId};
use crate::{shapes::Shape, StrError};
use russell_lab::sort2;
use std::collections::{HashMap, HashSet};

/// Holds points, edges and faces on the boundaries of a mesh
pub struct Boundary {
    /// Set of points on the boundaries
    ///
    /// Note: a boundary point belongs to a boundary edge or a boundary face
    pub points: HashSet<PointId>,

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

    /// The minimum coordinates; len = space_ndim
    pub min: Vec<f64>,

    /// The maximum coordinates; len = space_ndim
    pub max: Vec<f64>,
}

impl Boundary {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh) -> Result<Self, StrError> {
        let shapes = alloc_cell_shapes(mesh)?;
        let mut boundary = match mesh.space_ndim {
            2 => {
                let edges = all_edges_2d(mesh, &shapes)?;
                Boundary::two_dim(mesh, &shapes, &edges)?
            }
            3 => {
                let faces = all_faces_3d(mesh, &shapes)?;
                Boundary::three_dim(mesh, &shapes, &faces)?
            }
            _ => panic!("space_ndim must be 2 or 3"),
        };
        // handle points of (rods in 2D or 3D) or (shells in 3D)
        mesh.cells.iter().for_each(|cell| {
            if cell.geo_ndim == 1 || (cell.geo_ndim == 2 && mesh.space_ndim == 3) {
                cell.points.iter().for_each(|id| {
                    boundary.points.insert(*id);
                    for j in 0..mesh.space_ndim {
                        if mesh.points[*id].coords[j] < boundary.min[j] {
                            boundary.min[j] = mesh.points[*id].coords[j];
                        }
                        if mesh.points[*id].coords[j] > boundary.max[j] {
                            boundary.max[j] = mesh.points[*id].coords[j];
                        }
                    }
                });
            }
        });
        Ok(boundary)
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
    pub fn two_dim(
        mesh: &Mesh,
        shapes: &Vec<Shape>,
        edges: &HashMap<EdgeKey, Vec<(CellId, usize)>>,
    ) -> Result<Boundary, StrError> {
        // check
        if mesh.space_ndim != 2 {
            return Err("this function works in 2D only");
        }

        // output
        let mut boundary = Boundary {
            points: HashSet::new(),
            edges: HashMap::new(),
            faces: HashMap::new(),
            min: vec![f64::MAX; mesh.space_ndim],
            max: vec![f64::MIN; mesh.space_ndim],
        };

        // loop over all edges
        for (edge_key, shared_by) in edges {
            // skip internal edges (those shared by multiple cells)
            if shared_by.len() != 1 {
                continue;
            }

            // skip already handled edges
            if boundary.edges.contains_key(edge_key) {
                continue;
            }

            // cell and edge
            let (cell_id, e) = shared_by[0];
            let cell = &mesh.cells[cell_id];
            let shape = &shapes[cell_id];
            let mut edge = Edge {
                points: vec![0; shape.edge_nnode],
            };

            // process points on edge
            for i in 0..shape.edge_nnode {
                edge.points[i] = cell.points[shape.edge_node_id(e, i)];
                boundary.points.insert(edge.points[i]);
                for j in 0..mesh.space_ndim {
                    if mesh.points[edge.points[i]].coords[j] < boundary.min[j] {
                        boundary.min[j] = mesh.points[edge.points[i]].coords[j];
                    }
                    if mesh.points[edge.points[i]].coords[j] > boundary.max[j] {
                        boundary.max[j] = mesh.points[edge.points[i]].coords[j];
                    }
                }
            }

            // new edge
            boundary.edges.insert(*edge_key, edge);
        }
        Ok(boundary)
    }

    /// Finds boundary entities in 3D
    ///
    /// **Note:** Call this function after `all_faces_3d`.
    ///
    /// # Input
    ///
    /// * `mesh` -- the Mesh
    /// * `shapes` -- the shapes of cells (len == cells.len())
    /// * `faces` -- all faces (internal and boundary)
    pub fn three_dim(
        mesh: &Mesh,
        shapes: &Vec<Shape>,
        faces: &HashMap<FaceKey, Vec<(CellId, usize)>>,
    ) -> Result<Boundary, StrError> {
        // check
        if mesh.space_ndim != 3 {
            return Err("this function works in 3D only");
        }

        // output
        let mut boundary = Boundary {
            points: HashSet::new(),
            edges: HashMap::new(),
            faces: HashMap::new(),
            min: vec![f64::MAX; mesh.space_ndim],
            max: vec![f64::MIN; mesh.space_ndim],
        };

        // sort face keys just so the next loop is deterministic
        let mut face_keys: Vec<_> = faces.keys().collect();
        face_keys.sort();

        // loop over all faces
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
            let mut face = Face {
                points: vec![0; shape.face_nnode],
            };

            // process points on face
            for i in 0..shape.face_nnode {
                face.points[i] = cell.points[shape.face_node_id(f, i)];
                boundary.points.insert(face.points[i]);
                for j in 0..mesh.space_ndim {
                    if mesh.points[face.points[i]].coords[j] < boundary.min[j] {
                        boundary.min[j] = mesh.points[face.points[i]].coords[j];
                    }
                    if mesh.points[face.points[i]].coords[j] > boundary.max[j] {
                        boundary.max[j] = mesh.points[face.points[i]].coords[j];
                    }
                }
            }

            // loop over all edges on face
            let face_shape = Shape::new(mesh.space_ndim, 2, shape.face_nnode)?;
            for e in 0..face_shape.nedge {
                // define edge key (sorted point ids)
                let mut edge_key: EdgeKey = (
                    face.points[face_shape.edge_node_id(e, 0)],
                    face.points[face_shape.edge_node_id(e, 1)],
                );
                sort2(&mut edge_key);

                // skip already handled edge
                if boundary.edges.contains_key(&edge_key) {
                    continue;
                }

                // new edge
                let mut edge = Edge {
                    points: vec![0; face_shape.edge_nnode],
                };
                for i in 0..face_shape.edge_nnode {
                    edge.points[i] = face.points[face_shape.edge_node_id(e, i)];
                }
                boundary.edges.insert(edge_key, edge);
            }

            // new face
            boundary.faces.insert(*face_key, face);
        }
        Ok(boundary)
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
        assert_eq!(boundary.min, &[0.0, 0.0]);
        assert_eq!(boundary.max, &[2.0, 1.0]);
        let mut points: Vec<_> = boundary.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[0, 1, 2, 3, 4, 5]);
        Ok(())
    }

    #[test]
    fn boundary_2d_mixed_works() -> Result<(), StrError> {
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
        assert_eq!(boundary.min, &[0.0, 0.0]);
        assert_eq!(boundary.max, &[2.0, 1.0]);
        let mut points: Vec<_> = boundary.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[0, 1, 2, 3, 4]);
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
        assert_eq!(boundary.min, &[0.0, 0.0, 0.0]);
        assert_eq!(boundary.max, &[1.0, 1.0, 2.0]);
        let mut points: Vec<_> = boundary.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
        Ok(())
    }

    #[test]
    fn boundary_3d_mixed_works() -> Result<(), StrError> {
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
        assert_eq!(boundary.min, &[0.0, -1.0, 0.0]);
        assert_eq!(boundary.max, &[1.0, 2.0, 1.0]);
        let mut points: Vec<_> = boundary.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
        Ok(())
    }
}
