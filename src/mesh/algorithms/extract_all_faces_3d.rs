use crate::mesh::{FaceKey, MapFaceToCells, Mesh};
use crate::shapes::Shape;
use crate::StrError;
use russell_lab::sort4;
use std::collections::HashMap;

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
pub(crate) fn extract_all_faces_3d(mesh: &Mesh, shapes: &Vec<Shape>) -> Result<MapFaceToCells, StrError> {
    if mesh.space_ndim != 3 {
        return Err("this function works in 3D only");
    }
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
    Ok(faces)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::extract_all_faces_3d;
    use crate::mesh::{allocate_shapes, Mesh, Samples};
    use crate::StrError;

    #[test]
    fn capture_some_wrong_input() {
        let mesh = Mesh {
            space_ndim: 1,
            points: Vec::new(),
            cells: Vec::new(),
        };
        assert_eq!(
            extract_all_faces_3d(&mesh, &Vec::new()).err(),
            Some("this function works in 3D only")
        );
    }

    #[test]
    fn extract_all_faces_3d_works() -> Result<(), StrError> {
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
        let shapes = allocate_shapes(&mesh)?;
        let faces = extract_all_faces_3d(&mesh, &shapes)?;
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
    fn extract_all_faces_3d_mixed_works() -> Result<(), StrError> {
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
        let shapes = allocate_shapes(&mesh)?;
        let faces = extract_all_faces_3d(&mesh, &shapes)?;
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
