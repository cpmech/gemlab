use crate::mesh::{FaceKey, MapFaceToCells, Mesh};
use crate::shapes::Shape;
use russell_lab::sort4;
use std::collections::HashMap;

/// Extracts all faces (internal and boundary)
///
/// # Input
///
/// * `mesh` -- the Mesh
/// * `shapes` -- the shapes of cells (len == cells.len())
///
/// # Output
///
/// * Returns a map relating face keys to `Vec<(cell_id, f)>` where:
///     - `cell_id` -- is the id of the cell sharing the face
///     - `f` -- is the cell's local face index
///
/// # Panics
///
/// 1. It panics if `shapes.len() != mesh.cells.len()` (i.e., the shapes vector and the mesh must be compatible)
/// 2. It panics if `mesh.space_ndim != 3` (i.e., this function works in 3D only)
pub(crate) fn extract_all_faces(mesh: &Mesh, shapes: &Vec<Shape>) -> MapFaceToCells {
    assert_eq!(mesh.space_ndim, 3);
    let mut faces: MapFaceToCells = HashMap::new();
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
    faces
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::extract_all_faces;
    use crate::mesh::{allocate_shapes, Samples};

    #[test]
    fn extract_all_faces_works() {
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
        let shapes = allocate_shapes(&mesh).unwrap();
        let faces = extract_all_faces(&mesh, &shapes);
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
    }

    #[test]
    fn extract_all_faces_mixed_works() {
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
        let shapes = allocate_shapes(&mesh).unwrap();
        let faces = extract_all_faces(&mesh, &shapes);
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
    }
}
