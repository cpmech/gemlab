use super::{CellId, Mesh, PointId};
use crate::shapes::Shape;
use crate::StrError;
use russell_lab::{sort2, sort4};
use std::collections::HashMap;

/// Aliases (usize,usize) as the key of edges
///
/// # Note
///
/// Since the local numbering scheme runs over "corners" first, we can compare
/// edges using only two points; i.e., the middle points don't matter.
pub type EdgeKey = (usize, usize);

/// Aliases (usize,usize,usize,usize) as the key of faces
///
/// # Note
///
/// If a face has at most 3 points, the fourth entry in the key will be set to the total number of points.
/// In this way, we can compare 4-node (or more nodes) faces with each other. Since the local numbering
/// scheme runs over the "corners" first, the middle points don't matter.
pub type FaceKey = (usize, usize, usize, usize);

/// Holds the point ids of an edge (an entity belonging to a solid cell in 2D or a face in 3D)
#[derive(Clone, Debug)]
pub struct Edge {
    /// List of points defining this edge; in the right order (unsorted)
    pub points: Vec<PointId>,
}

/// Holds the point ids of a face (an entity belonging to a solid cell in 3D)
#[derive(Clone, Debug)]
pub struct Face {
    /// List of points defining this face; in the right order (unsorted)
    pub points: Vec<PointId>,
}

/// Extracts all edges (internal and boundary) in 2D
///
/// # Input
///
/// * `mesh` -- the Mesh
/// * `shapes` -- the shapes of cells (len == cells.len())
///
/// # Output
///
/// * Returns a map relating edge keys to `Vec<(cell_id, e)>` where:
///     - `cell_id` -- the id of the cell sharing the edge
///     - `e` -- is the cell's local edge index
pub fn all_edges_2d(mesh: &Mesh, shapes: &Vec<Shape>) -> Result<HashMap<EdgeKey, Vec<(CellId, usize)>>, StrError> {
    if mesh.space_ndim != 2 {
        return Err("this function works in 2D only");
    }
    let mut edges: HashMap<EdgeKey, Vec<(CellId, usize)>> = HashMap::new();
    mesh.cells.iter().zip(shapes).for_each(|(cell, shape)| {
        if shape.geo_ndim == 2 {
            for e in 0..shape.nedge {
                let mut edge_key: EdgeKey = (
                    cell.points[shape.edge_node_id(e, 0)],
                    cell.points[shape.edge_node_id(e, 1)],
                );
                sort2(&mut edge_key);
                let data = edges.entry(edge_key).or_insert(Vec::new());
                data.push((cell.id, e));
            }
        }
    });
    Ok(edges)
}

/// Extracts all faces (internal and boundary) in 3D
///
/// # Input
///
/// * `mesh` -- the Mesh
/// * `shapes` -- the shapes of cells (len == cells.len())
///
/// # Output
///
/// * Returns a map relating face keys to `Vec<(cell_id, f)>` where:
///     - `cell_id` -- the id of the cell sharing the face
///     - `f` -- is the cell's local face index
pub fn all_faces_3d(mesh: &Mesh, shapes: &Vec<Shape>) -> Result<HashMap<FaceKey, Vec<(CellId, usize)>>, StrError> {
    if mesh.space_ndim != 3 {
        return Err("this function works in 3D only");
    }
    let mut faces: HashMap<FaceKey, Vec<(CellId, usize)>> = HashMap::new();
    mesh.cells.iter().zip(shapes).for_each(|(cell, shape)| {
        if shape.geo_ndim == 3 {
            for f in 0..shape.nface {
                let mut face_key: FaceKey = if shape.face_nnode > 3 {
                    (
                        cell.points[shape.face_node_id(f, 0)],
                        cell.points[shape.face_node_id(f, 1)],
                        cell.points[shape.face_node_id(f, 2)],
                        cell.points[shape.face_node_id(f, 3)],
                    )
                } else {
                    (
                        cell.points[shape.face_node_id(f, 0)],
                        cell.points[shape.face_node_id(f, 1)],
                        cell.points[shape.face_node_id(f, 2)],
                        mesh.points.len(),
                    )
                };
                sort4(&mut face_key);
                let data = faces.entry(face_key).or_insert(Vec::new());
                data.push((cell.id, f));
            }
        }
    });
    Ok(faces)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{all_edges_2d, all_faces_3d};
    use crate::mesh::{Samples, Shapes};
    use crate::StrError;

    #[test]
    fn all_edges_2d_works() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let shapes = Shapes::new(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let mut keys: Vec<_> = edges.keys().collect();
        keys.sort();
        assert_eq!(keys, [&(0, 1), &(0, 3), &(1, 2), &(1, 4), &(2, 3), &(2, 5), &(4, 5)]);
        assert_eq!(edges.get(&(0, 1)).unwrap(), &[(0, 0)]);
        assert_eq!(edges.get(&(0, 3)).unwrap(), &[(0, 3)]);
        assert_eq!(edges.get(&(1, 2)).unwrap(), &[(0, 1), (1, 3)]);
        assert_eq!(edges.get(&(1, 4)).unwrap(), &[(1, 0)]);
        assert_eq!(edges.get(&(2, 3)).unwrap(), &[(0, 2)]);
        assert_eq!(edges.get(&(2, 5)).unwrap(), &[(1, 2)]);
        assert_eq!(edges.get(&(4, 5)).unwrap(), &[(1, 1)]);
        Ok(())
    }

    #[test]
    fn all_edges_2d_mixed_works() -> Result<(), StrError> {
        //           4---------3
        //           |         |
        //           |   [1]   |
        //           |         |
        //  0--------1---------2
        let mesh = Samples::mixed_shapes_2d();
        let shapes = Shapes::new(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let mut keys: Vec<_> = edges.keys().collect();
        keys.sort();
        assert_eq!(keys, [&(1, 2), &(1, 4), &(2, 3), &(3, 4)]);
        assert_eq!(edges.get(&(1, 2)).unwrap(), &[(1, 0)]);
        assert_eq!(edges.get(&(1, 4)).unwrap(), &[(1, 3)]);
        assert_eq!(edges.get(&(2, 3)).unwrap(), &[(1, 1)]);
        assert_eq!(edges.get(&(3, 4)).unwrap(), &[(1, 2)]);
        Ok(())
    }

    #[test]
    fn all_faces_3d_works() -> Result<(), StrError> {
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
        let shapes = Shapes::new(&mesh)?;
        let faces = all_faces_3d(&mesh, &shapes)?;
        let mut keys: Vec<_> = faces.keys().collect();
        keys.sort();
        assert_eq!(
            keys,
            [
                &(0, 1, 2, 3),
                &(0, 1, 4, 5),
                &(0, 3, 4, 7),
                &(1, 2, 5, 6),
                &(2, 3, 6, 7),
                &(4, 5, 6, 7),
                &(4, 5, 8, 9),
                &(4, 7, 8, 11),
                &(5, 6, 9, 10),
                &(6, 7, 10, 11),
                &(8, 9, 10, 11),
            ]
        );
        assert_eq!(faces.get(&(0, 1, 2, 3)).unwrap(), &[(0, 4)]);
        assert_eq!(faces.get(&(0, 1, 4, 5)).unwrap(), &[(0, 2)]);
        assert_eq!(faces.get(&(0, 3, 4, 7)).unwrap(), &[(0, 0)]);
        assert_eq!(faces.get(&(1, 2, 5, 6)).unwrap(), &[(0, 1)]);
        assert_eq!(faces.get(&(2, 3, 6, 7)).unwrap(), &[(0, 3)]);
        assert_eq!(faces.get(&(4, 5, 6, 7)).unwrap(), &[(0, 5), (1, 4)]);
        assert_eq!(faces.get(&(4, 5, 8, 9)).unwrap(), &[(1, 2)]);
        assert_eq!(faces.get(&(4, 7, 8, 11)).unwrap(), &[(1, 0)]);
        assert_eq!(faces.get(&(5, 6, 9, 10)).unwrap(), &[(1, 1)]);
        assert_eq!(faces.get(&(6, 7, 10, 11)).unwrap(), &[(1, 3)]);
        assert_eq!(faces.get(&(8, 9, 10, 11)).unwrap(), &[(1, 5)]);
        Ok(())
    }

    #[test]
    fn all_faces_3d_mixed_works() -> Result<(), StrError> {
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
        let shapes = Shapes::new(&mesh)?;
        let faces = all_faces_3d(&mesh, &shapes)?;
        let mut keys: Vec<_> = faces.keys().collect();
        keys.sort();
        assert_eq!(
            keys,
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
        assert_eq!(faces.get(&(0, 1, 2, 3)).unwrap(), &[(0, 4)]);
        assert_eq!(faces.get(&(0, 1, 4, 5)).unwrap(), &[(0, 2)]);
        assert_eq!(faces.get(&(0, 3, 4, 7)).unwrap(), &[(0, 0)]);
        assert_eq!(faces.get(&(1, 2, 5, 6)).unwrap(), &[(0, 1)]);
        assert_eq!(faces.get(&(2, 3, 6, 7)).unwrap(), &[(0, 3)]);
        assert_eq!(faces.get(&(2, 3, 6, 13)).unwrap(), &[(1, 0)]);
        assert_eq!(faces.get(&(2, 3, 8, 13)).unwrap(), &[(1, 2)]);
        assert_eq!(faces.get(&(2, 6, 8, 13)).unwrap(), &[(1, 1)]);
        assert_eq!(faces.get(&(3, 6, 8, 13)).unwrap(), &[(1, 3)]);
        assert_eq!(faces.get(&(4, 5, 6, 7)).unwrap(), &[(0, 5)]);
        Ok(())
    }
}
