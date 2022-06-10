use crate::mesh::{Edge, EdgeKey, Extract, Face, Features, MapFaceToCells, Mesh};
use russell_lab::sort2;
use std::collections::{HashMap, HashSet};

/// Extracts mesh features in 3D
///
/// # Panics
///
/// 1. It panics if `mesh.space != 3` (i.e., this function works in 3D only)
/// 2. It panics if the mesh data is invalid, e.g., the cell points array doesn't contain enough points
pub(crate) fn extract_features_3d(mesh: &Mesh, faces: &MapFaceToCells, extract: Extract) -> Features {
    assert_eq!(mesh.ndim, 3);

    // output
    let mut features = Features {
        points: HashSet::new(),
        edges: HashMap::new(),
        faces: HashMap::new(),
        min: vec![f64::MAX; mesh.ndim],
        max: vec![f64::MIN; mesh.ndim],
    };

    // sort face keys just so the next loop is deterministic
    let mut face_keys: Vec<_> = faces.keys().collect();
    face_keys.sort();

    // loop over all faces
    for face_key in face_keys {
        let shared_by = faces.get(face_key).unwrap();

        // accept feature?
        let accept = match extract {
            Extract::All => true,
            Extract::Boundary => shared_by.len() == 1, // boundary face; with only one shared cell
            Extract::Interior => shared_by.len() != 1, // interior face; shared by multiple cells
        };
        if !accept {
            continue;
        }

        // cell and face
        let (cell_id, f) = shared_by[0];
        let cell = &mesh.cells[cell_id];
        let mut face = Face {
            kind: cell.kind.face_kind().unwrap(),
            points: vec![0; cell.kind.face_nnode()],
        };

        // process points on face
        for i in 0..face.points.len() {
            face.points[i] = cell.points[cell.kind.face_node_id(f, i)];
            features.points.insert(face.points[i]);
            for j in 0..mesh.ndim {
                features.min[j] = f64::min(features.min[j], mesh.points[face.points[i]].coords[j]);
                features.max[j] = f64::max(features.max[j], mesh.points[face.points[i]].coords[j]);
            }
        }

        // loop over all edges on face
        for e in 0..face.kind.nedge() {
            // define edge key (sorted point ids)
            let mut edge_key: EdgeKey = (
                face.points[face.kind.edge_node_id(e, 0)],
                face.points[face.kind.edge_node_id(e, 1)],
            );
            sort2(&mut edge_key);

            // skip already handled edge
            if features.edges.contains_key(&edge_key) {
                continue;
            }

            // new edge
            let mut edge = Edge {
                kind: face.kind.edge_kind().unwrap(),
                points: vec![0; face.kind.edge_nnode()],
            };
            for i in 0..edge.points.len() {
                edge.points[i] = face.points[face.kind.edge_node_id(e, i)];
            }
            features.edges.insert(edge_key, edge);
        }

        // new face
        features.faces.insert(*face_key, face);
    }
    features
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::extract_features_3d;
    use crate::mesh::algorithms::extract_all_faces;
    use crate::mesh::{EdgeKey, Extract, FaceKey, Features, PointId, Samples};
    use crate::util::AsArray2D;
    use crate::StrError;

    fn validate_edges<'a, T>(
        features: &Features,
        correct_keys: &[EdgeKey], // sorted
        correct_points: &'a T,
    ) where
        T: AsArray2D<'a, PointId>,
    {
        let mut keys: Vec<_> = features.edges.keys().map(|k| *k).collect();
        keys.sort();
        assert_eq!(keys, correct_keys);
        for i in 0..keys.len() {
            assert_eq!(features.edges.get(&keys[i]).unwrap().points, correct_points.row(i));
        }
    }

    fn validate_faces<'a, T>(
        features: &Features,
        correct_keys: &[FaceKey], // sorted
        correct_points: &'a T,
    ) where
        T: AsArray2D<'a, PointId>,
    {
        let mut keys: Vec<_> = features.faces.keys().map(|k| *k).collect();
        keys.sort();
        assert_eq!(keys, correct_keys);
        for i in 0..keys.len() {
            assert_eq!(features.faces.get(&keys[i]).unwrap().points, correct_points.row(i));
        }
    }

    #[test]
    fn extract_features_3d_works() -> Result<(), StrError> {
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
        let faces = extract_all_faces(&mesh);
        let features = extract_features_3d(&mesh, &faces, Extract::Boundary);
        let correct_edge_keys = [
            (0, 1),
            (0, 3),
            (0, 4),
            (1, 2),
            (1, 5),
            (2, 3),
            (2, 6),
            (3, 7),
            (4, 5),
            (4, 7),
            (4, 8),
            (5, 6),
            (5, 9),
            (6, 7),
            (6, 10),
            (7, 11),
            (8, 9),
            (8, 11),
            (9, 10),
            (10, 11),
        ];
        let correct_edge_points = [
            [0, 1],
            [3, 0],
            [0, 4],
            [1, 2],
            [5, 1],
            [2, 3],
            [6, 2],
            [3, 7],
            [4, 5],
            [7, 4],
            [4, 8],
            [5, 6],
            [9, 5],
            [6, 7],
            [10, 6],
            [7, 11],
            [8, 9],
            [11, 8],
            [9, 10],
            [10, 11],
        ];
        validate_edges(&features, &correct_edge_keys, &correct_edge_points);
        let correct_face_keys = [
            (0, 1, 2, 3),
            (0, 1, 4, 5),
            (0, 3, 4, 7),
            (1, 2, 5, 6),
            (2, 3, 6, 7),
            (4, 5, 8, 9),
            (4, 7, 8, 11),
            (5, 6, 9, 10),
            (6, 7, 10, 11),
            (8, 9, 10, 11),
        ];
        let correct_face_points = [
            [0, 3, 2, 1],
            [0, 1, 5, 4],
            [0, 4, 7, 3],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [4, 5, 9, 8],
            [4, 8, 11, 7],
            [5, 6, 10, 9],
            [6, 7, 11, 10],
            [8, 9, 10, 11],
        ];
        validate_faces(&features, &correct_face_keys, &correct_face_points);
        assert_eq!(features.min, &[0.0, 0.0, 0.0]);
        assert_eq!(features.max, &[1.0, 1.0, 2.0]);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
        Ok(())
    }

    #[test]
    fn extract_features_3d_all_works() -> Result<(), StrError> {
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
        let faces = extract_all_faces(&mesh);
        let features = extract_features_3d(&mesh, &faces, Extract::All);
        let correct_edge_keys = [
            (0, 1),
            (0, 3),
            (0, 4),
            (1, 2),
            (1, 5),
            (2, 3),
            (2, 6),
            (3, 7),
            (4, 5),
            (4, 7),
            (4, 8),
            (5, 6),
            (5, 9),
            (6, 7),
            (6, 10),
            (7, 11),
            (8, 9),
            (8, 11),
            (9, 10),
            (10, 11),
        ];
        let correct_edge_points = [
            [0, 1],
            [3, 0],
            [0, 4],
            [1, 2],
            [5, 1],
            [2, 3],
            [6, 2],
            [3, 7],
            [4, 5],
            [7, 4],
            [4, 8],
            [5, 6],
            [9, 5],
            [6, 7],
            [10, 6],
            [7, 11],
            [8, 9],
            [11, 8],
            [9, 10],
            [10, 11],
        ];
        validate_edges(&features, &correct_edge_keys, &correct_edge_points);
        let correct_face_keys = [
            (0, 1, 2, 3),
            (0, 1, 4, 5),
            (0, 3, 4, 7),
            (1, 2, 5, 6),
            (2, 3, 6, 7),
            (4, 5, 6, 7),
            (4, 5, 8, 9),
            (4, 7, 8, 11),
            (5, 6, 9, 10),
            (6, 7, 10, 11),
            (8, 9, 10, 11),
        ];
        let correct_face_points = [
            [0, 3, 2, 1],
            [0, 1, 5, 4],
            [0, 4, 7, 3],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [4, 5, 6, 7],
            [4, 5, 9, 8],
            [4, 8, 11, 7],
            [5, 6, 10, 9],
            [6, 7, 11, 10],
            [8, 9, 10, 11],
        ];
        validate_faces(&features, &correct_face_keys, &correct_face_points);
        assert_eq!(features.min, &[0.0, 0.0, 0.0]);
        assert_eq!(features.max, &[1.0, 1.0, 2.0]);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
        Ok(())
    }

    #[test]
    fn extract_features_3d_interior_works() -> Result<(), StrError> {
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
        let faces = extract_all_faces(&mesh);
        let features = extract_features_3d(&mesh, &faces, Extract::Interior);
        let correct_edge_keys = [(4, 5), (4, 7), (5, 6), (6, 7)];
        let correct_edge_points = [[5, 4], [4, 7], [6, 5], [7, 6]];
        validate_edges(&features, &correct_edge_keys, &correct_edge_points);
        let correct_face_keys = [(4, 5, 6, 7)];
        let correct_face_points = [[4, 5, 6, 7]];
        validate_faces(&features, &correct_face_keys, &correct_face_points);
        assert_eq!(features.min, &[0.0, 0.0, 1.0]);
        assert_eq!(features.max, &[1.0, 1.0, 1.0]);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[4, 5, 6, 7]);
        Ok(())
    }

    #[test]
    fn extract_features_3d_mixed_works() -> Result<(), StrError> {
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
        let faces = extract_all_faces(&mesh);
        let features = extract_features_3d(&mesh, &faces, Extract::Boundary);
        let correct_edge_keys = [
            (0, 1),
            (0, 3),
            (0, 4),
            (1, 2),
            (1, 5),
            (2, 3),
            (2, 6),
            (2, 8),
            (3, 6),
            (3, 7),
            (3, 8),
            (4, 5),
            (4, 7),
            (5, 6),
            (6, 7),
            (6, 8),
        ];
        let correct_edge_points = [
            [0, 1],
            [3, 0],
            [0, 4],
            [1, 2],
            [5, 1],
            [2, 3],
            [6, 2],
            [2, 8],
            [3, 6],
            [3, 7],
            [8, 3],
            [4, 5],
            [7, 4],
            [5, 6],
            [6, 7],
            [6, 8],
        ];
        validate_edges(&features, &correct_edge_keys, &correct_edge_points);
        let correct_face_keys = [
            (0, 1, 2, 3),
            (0, 1, 4, 5),
            (0, 3, 4, 7),
            (1, 2, 5, 6),
            (2, 3, 6, 7),
            (2, 3, 6, 13),
            (2, 3, 8, 13),
            (2, 6, 8, 13),
            (3, 6, 8, 13),
            (4, 5, 6, 7),
        ];
        let correct_face_points: &[&[PointId]] = &[
            &[0, 3, 2, 1],
            &[0, 1, 5, 4],
            &[0, 4, 7, 3],
            &[1, 2, 6, 5],
            &[2, 3, 7, 6],
            &[2, 6, 3],
            &[2, 3, 8],
            &[2, 8, 6],
            &[8, 3, 6],
            &[4, 5, 6, 7],
        ];
        validate_faces(&features, &correct_face_keys, &correct_face_points);
        assert_eq!(features.min, &[0.0, 0.0, 0.0]);
        assert_eq!(features.max, &[1.0, 2.0, 1.0]);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[0, 1, 2, 3, 4, 5, 6, 7, 8]);
        Ok(())
    }

    #[test]
    fn extract_features_3d_block_works() -> Result<(), StrError> {
        //              51--------58--------54--------74--------71
        //              /.                  /.                  /|
        //             / .                 / .                 / |
        //           55  .               57  .               73  |
        //           /   .               /   .               /   |
        //          /    .              /    .              /    |
        //        52--------56--------53--------72--------70     |
        //        /.     .            /.     .            /|     |
        //       / .    59           / .    62           / |    76
        //     65  .     .         67  .     .         79  |     |
        //     /   .     .         /   .     .         /   |     |
        //    /    .     .        /    .     .        /    |     |
        //  63========66========64========78========77     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |    60     .       |    61     .       |    75     |
        //   |     .     4 - - - |15 - . - - 7 - - - |41 - | - -35
        //   |     .    /.       |     .    /.       |     |    /|
        //   |     .   / .       |     .   / .       |     |   / |
        //   |     . 12  .       |     . 14  .       |     | 40  |
        //   |     . /   .       |     . /   .       |     | /   |
        //  68     ./    .      69     ./    .      80     |/    |
        //   |     5 - - - -13 - | - - 6 - - - -39 - | - -34     |
        //   |    /.     .       |    /.     .       |    /|     |
        //   |   / .    16       |   / .    19       |   / |    43
        //   | 27  .     .       | 29  .     .       | 49  |     |
        //   | /   .     .       | /   .     .       | /   |     |
        //   |/    .     .       |/    .     .       |/    |     |
        //  22========28========23========48========45     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |    17     .       |    18     .       |    42     |
        //   |     .     0 - - - |11 - . - - 3 - - - |38 - | - -33
        //   |     .    /        |     .    /        |     |    /
        //   |     .   /         |     .   /         |     |   /
        //   |     .  8          |     . 10          |     | 37
        //   |     . /           |     . /           |     | /
        //  30     ./           31     ./           50     |/
        //   |     1 - - - - 9 - | - - 2 - - - -36 - | - -32
        //   |    /              |    /              |    /
        //   |   /               |   /               |   /
        //   | 24                | 26                | 47
        //   | /                 | /                 | /
        //   |/                  |/                  |/
        //  20========25========21========46========44
        let mesh = Samples::block_3d_eight_hex20();
        let faces = extract_all_faces(&mesh);
        let features = extract_features_3d(&mesh, &faces, Extract::Boundary);
        let correct_edge_keys = [
            (0, 1),
            (0, 3),
            (0, 4),
            (1, 2),
            (1, 5),
            (1, 20),
            (2, 3),
            (2, 21),
            (2, 32),
            (3, 7),
            (3, 33),
            (4, 5),
            (4, 7),
            (4, 51),
            (5, 22),
            (5, 52),
            (7, 35),
            (7, 54),
            (20, 21),
            (20, 22),
            (21, 23),
            (21, 44),
            (22, 23),
            (22, 63),
            (23, 45),
            (23, 64),
            (32, 33),
            (32, 34),
            (32, 44),
            (33, 35),
            (34, 35),
            (34, 45),
            (34, 70),
            (35, 71),
            (44, 45),
            (45, 77),
            (51, 52),
            (51, 54),
            (52, 53),
            (52, 63),
            (53, 54),
            (53, 64),
            (53, 70),
            (54, 71),
            (63, 64),
            (64, 77),
            (70, 71),
            (70, 77),
        ];
        let correct_edge_points = [
            [0, 1, 8],
            [3, 0, 11],
            [0, 4, 16],
            [1, 2, 9],
            [5, 1, 17],
            [1, 20, 24],
            [2, 3, 10],
            [21, 2, 26],
            [2, 32, 36],
            [3, 7, 19],
            [33, 3, 38],
            [4, 5, 12],
            [7, 4, 15],
            [4, 51, 59],
            [5, 22, 27],
            [52, 5, 60],
            [35, 7, 41],
            [7, 54, 62],
            [20, 21, 25],
            [22, 20, 30],
            [23, 21, 31],
            [21, 44, 46],
            [22, 23, 28],
            [63, 22, 68],
            [23, 45, 48],
            [64, 23, 69],
            [32, 33, 37],
            [32, 34, 42],
            [44, 32, 47],
            [33, 35, 43],
            [34, 35, 40],
            [45, 34, 49],
            [34, 70, 75],
            [35, 71, 76],
            [45, 44, 50],
            [77, 45, 80],
            [51, 52, 55],
            [54, 51, 58],
            [53, 52, 56],
            [52, 63, 65],
            [54, 53, 57],
            [53, 64, 67],
            [70, 53, 72],
            [71, 54, 74],
            [63, 64, 66],
            [64, 77, 78],
            [70, 71, 73],
            [77, 70, 79],
        ];
        validate_edges(&features, &correct_edge_keys, &correct_edge_points);
        let correct_face_keys = [
            (0, 1, 2, 3),
            (0, 1, 4, 5),
            (0, 3, 4, 7),
            (1, 2, 20, 21),
            (1, 5, 20, 22),
            (2, 3, 32, 33),
            (2, 21, 32, 44),
            (3, 7, 33, 35),
            (4, 5, 51, 52),
            (4, 7, 51, 54),
            (5, 22, 52, 63),
            (7, 35, 54, 71),
            (20, 21, 22, 23),
            (21, 23, 44, 45),
            (22, 23, 63, 64),
            (23, 45, 64, 77),
            (32, 33, 34, 35),
            (32, 34, 44, 45),
            (34, 35, 70, 71),
            (34, 45, 70, 77),
            (51, 52, 53, 54),
            (52, 53, 63, 64),
            (53, 54, 70, 71),
            (53, 64, 70, 77),
        ];
        let correct_face_points = [
            [0, 3, 2, 1, 11, 10, 9, 8],
            [0, 1, 5, 4, 8, 17, 12, 16],
            [0, 4, 7, 3, 16, 15, 19, 11],
            [1, 2, 21, 20, 9, 26, 25, 24],
            [1, 20, 22, 5, 24, 30, 27, 17],
            [3, 33, 32, 2, 38, 37, 36, 10],
            [2, 32, 44, 21, 36, 47, 46, 26],
            [3, 7, 35, 33, 19, 41, 43, 38],
            [4, 5, 52, 51, 12, 60, 55, 59],
            [4, 51, 54, 7, 59, 58, 62, 15],
            [5, 22, 63, 52, 27, 68, 65, 60],
            [7, 54, 71, 35, 62, 74, 76, 41],
            [20, 21, 23, 22, 25, 31, 28, 30],
            [21, 44, 45, 23, 46, 50, 48, 31],
            [22, 23, 64, 63, 28, 69, 66, 68],
            [23, 45, 77, 64, 48, 80, 78, 69],
            [32, 33, 35, 34, 37, 43, 40, 42],
            [44, 32, 34, 45, 47, 42, 49, 50],
            [34, 35, 71, 70, 40, 76, 73, 75],
            [45, 34, 70, 77, 49, 75, 79, 80],
            [51, 52, 53, 54, 55, 56, 57, 58],
            [52, 63, 64, 53, 65, 66, 67, 56],
            [54, 53, 70, 71, 57, 72, 73, 74],
            [53, 64, 77, 70, 67, 78, 79, 72],
        ];
        validate_faces(&features, &correct_face_keys, &correct_face_points);
        assert_eq!(features.min, &[0.0, 0.0, 0.0]);
        assert_eq!(features.max, &[2.0, 2.0, 4.0]);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(
            points,
            &[
                0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 31, 32,
                33, 34, 35, 36, 37, 38, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                60, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80
            ]
        );
        Ok(())
    }

    #[test]
    fn extract_features_3d_all_block_works() -> Result<(), StrError> {
        //              51--------58--------54--------74--------71
        //              /.                  /.                  /|
        //             / .                 / .                 / |
        //           55  .               57  .               73  |
        //           /   .               /   .               /   |
        //          /    .              /    .              /    |
        //        52--------56--------53--------72--------70     |
        //        /.     .            /.     .            /|     |
        //       / .    59           / .    62           / |    76
        //     65  .     .         67  .     .         79  |     |
        //     /   .     .         /   .     .         /   |     |
        //    /    .     .        /    .     .        /    |     |
        //  63========66========64========78========77     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |    60     .       |    61     .       |    75     |
        //   |     .     4 - - - |15 - . - - 7 - - - |41 - | - -35
        //   |     .    /.       |     .    /.       |     |    /|
        //   |     .   / .       |     .   / .       |     |   / |
        //   |     . 12  .       |     . 14  .       |     | 40  |
        //   |     . /   .       |     . /   .       |     | /   |
        //  68     ./    .      69     ./    .      80     |/    |
        //   |     5 - - - -13 - | - - 6 - - - -39 - | - -34     |
        //   |    /.     .       |    /.     .       |    /|     |
        //   |   / .    16       |   / .    19       |   / |    43
        //   | 27  .     .       | 29  .     .       | 49  |     |
        //   | /   .     .       | /   .     .       | /   |     |
        //   |/    .     .       |/    .     .       |/    |     |
        //  22========28========23========48========45     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |    17     .       |    18     .       |    42     |
        //   |     .     0 - - - |11 - . - - 3 - - - |38 - | - -33
        //   |     .    /        |     .    /        |     |    /
        //   |     .   /         |     .   /         |     |   /
        //   |     .  8          |     . 10          |     | 37
        //   |     . /           |     . /           |     | /
        //  30     ./           31     ./           50     |/
        //   |     1 - - - - 9 - | - - 2 - - - -36 - | - -32
        //   |    /              |    /              |    /
        //   |   /               |   /               |   /
        //   | 24                | 26                | 47
        //   | /                 | /                 | /
        //   |/                  |/                  |/
        //  20========25========21========46========44
        let mesh = Samples::block_3d_eight_hex20();
        let faces = extract_all_faces(&mesh);
        let features = extract_features_3d(&mesh, &faces, Extract::All);
        let correct_edge_keys = [
            (0, 1),
            (0, 3),
            (0, 4),
            (1, 2),
            (1, 5),
            (1, 20),
            (2, 3),
            (2, 6), // interior
            (2, 21),
            (2, 32),
            (3, 7),
            (3, 33),
            (4, 5),
            (4, 7),
            (4, 51),
            (5, 6), // interior
            (5, 22),
            (5, 52),
            (6, 7),  // interior
            (6, 23), // interior
            (6, 34), // interior
            (6, 53), // interior
            (7, 35),
            (7, 54),
            (20, 21),
            (20, 22),
            (21, 23),
            (21, 44),
            (22, 23),
            (22, 63),
            (23, 45),
            (23, 64),
            (32, 33),
            (32, 34),
            (32, 44),
            (33, 35),
            (34, 35),
            (34, 45),
            (34, 70),
            (35, 71),
            (44, 45),
            (45, 77),
            (51, 52),
            (51, 54),
            (52, 53),
            (52, 63),
            (53, 54),
            (53, 64),
            (53, 70),
            (54, 71),
            (63, 64),
            (64, 77),
            (70, 71),
            (70, 77),
        ];
        let correct_edge_points = [
            [0, 1, 8],
            [3, 0, 11],
            [0, 4, 16],
            [1, 2, 9],
            [5, 1, 17],
            [1, 20, 24],
            [2, 3, 10],
            [6, 2, 18], // interior
            [21, 2, 26],
            [2, 32, 36],
            [3, 7, 19],
            [33, 3, 38],
            [4, 5, 12],
            [7, 4, 15],
            [4, 51, 59],
            [5, 6, 13], // interior
            [5, 22, 27],
            [52, 5, 60],
            [6, 7, 14],  // interior
            [23, 6, 29], // interior
            [6, 34, 39], // interior
            [53, 6, 61], // interior
            [35, 7, 41],
            [7, 54, 62],
            [20, 21, 25],
            [22, 20, 30],
            [21, 23, 31], // [23, 21, 31] => the order changed because an interior face now can set this edge
            [21, 44, 46],
            [23, 22, 28], // [22, 23, 28] => the order changed because an interior face now can set this edge
            [63, 22, 68],
            [45, 23, 48], // [23, 45, 48] => the order changed because an interior face now can set this edge
            [23, 64, 69], // [64, 23, 69] => the order changed because an interior face now can set this edge
            [32, 33, 37],
            [34, 32, 42], // [32, 34, 42] => the order changed because an interior face now can set this edge
            [44, 32, 47],
            [33, 35, 43],
            [35, 34, 40], // [34, 35, 40] => the order changed because an interior face now can set this edge
            [34, 45, 49], // [45, 34, 49] => the order changed because an interior face now can set this edge
            [70, 34, 75], // [34, 70, 75] => the order changed because an interior face now can set this edge
            [35, 71, 76],
            [45, 44, 50],
            [77, 45, 80],
            [51, 52, 55],
            [54, 51, 58],
            [52, 53, 56], // [53, 52, 56] => the order changed because an interior face now can set this edge
            [52, 63, 65],
            [53, 54, 57], // [54, 53, 57] => the order changed because an interior face now can set this edge
            [64, 53, 67], // [53, 64, 67] => the order changed because an interior face now can set this edge
            [53, 70, 72], // [70, 53, 72] => the order changed because an interior face now can set this edge
            [71, 54, 74],
            [63, 64, 66],
            [64, 77, 78],
            [70, 71, 73],
            [77, 70, 79],
        ];
        validate_edges(&features, &correct_edge_keys, &correct_edge_points);
        let correct_face_keys = [
            (0, 1, 2, 3),
            (0, 1, 4, 5),
            (0, 3, 4, 7),
            (1, 2, 5, 6), // interior
            (1, 2, 20, 21),
            (1, 5, 20, 22),
            (2, 3, 6, 7), // interior
            (2, 3, 32, 33),
            (2, 6, 21, 23), // interior
            (2, 6, 32, 34), // interior
            (2, 21, 32, 44),
            (3, 7, 33, 35),
            (4, 5, 6, 7), // interior
            (4, 5, 51, 52),
            (4, 7, 51, 54),
            (5, 6, 22, 23), // interior
            (5, 6, 52, 53), // interior
            (5, 22, 52, 63),
            (6, 7, 34, 35),  // interior
            (6, 7, 53, 54),  // interior
            (6, 23, 34, 45), // interior
            (6, 23, 53, 64), // interior
            (6, 34, 53, 70), // interior
            (7, 35, 54, 71),
            (20, 21, 22, 23),
            (21, 23, 44, 45),
            (22, 23, 63, 64),
            (23, 45, 64, 77),
            (32, 33, 34, 35),
            (32, 34, 44, 45),
            (34, 35, 70, 71),
            (34, 45, 70, 77),
            (51, 52, 53, 54),
            (52, 53, 63, 64),
            (53, 54, 70, 71),
            (53, 64, 70, 77),
        ];
        let correct_face_points = [
            [0, 3, 2, 1, 11, 10, 9, 8],
            [0, 1, 5, 4, 8, 17, 12, 16],
            [0, 4, 7, 3, 16, 15, 19, 11],
            [1, 2, 6, 5, 9, 18, 13, 17], // interior
            [1, 2, 21, 20, 9, 26, 25, 24],
            [1, 20, 22, 5, 24, 30, 27, 17],
            [2, 3, 7, 6, 10, 19, 14, 18], // interior
            [3, 33, 32, 2, 38, 37, 36, 10],
            [21, 2, 6, 23, 26, 18, 29, 31], // interior
            [2, 32, 34, 6, 36, 42, 39, 18], // interior
            [2, 32, 44, 21, 36, 47, 46, 26],
            [3, 7, 35, 33, 19, 41, 43, 38],
            [4, 5, 6, 7, 12, 13, 14, 15], // interior
            [4, 5, 52, 51, 12, 60, 55, 59],
            [4, 51, 54, 7, 59, 58, 62, 15],
            [5, 22, 23, 6, 27, 28, 29, 13], // interior
            [5, 6, 53, 52, 13, 61, 56, 60], // interior
            [5, 22, 63, 52, 27, 68, 65, 60],
            [7, 6, 34, 35, 14, 39, 40, 41],  // interior
            [6, 7, 54, 53, 14, 62, 57, 61],  // interior
            [6, 23, 45, 34, 29, 48, 49, 39], // interior
            [23, 6, 53, 64, 29, 61, 67, 69], // interior
            [6, 34, 70, 53, 39, 75, 72, 61], // interior
            [7, 54, 71, 35, 62, 74, 76, 41],
            [20, 21, 23, 22, 25, 31, 28, 30],
            [21, 44, 45, 23, 46, 50, 48, 31],
            [22, 23, 64, 63, 28, 69, 66, 68],
            [23, 45, 77, 64, 48, 80, 78, 69],
            [32, 33, 35, 34, 37, 43, 40, 42],
            [44, 32, 34, 45, 47, 42, 49, 50],
            [34, 35, 71, 70, 40, 76, 73, 75],
            [45, 34, 70, 77, 49, 75, 79, 80],
            [51, 52, 53, 54, 55, 56, 57, 58],
            [52, 63, 64, 53, 65, 66, 67, 56],
            [54, 53, 70, 71, 57, 72, 73, 74],
            [53, 64, 77, 70, 67, 78, 79, 72],
        ];
        validate_faces(&features, &correct_face_keys, &correct_face_points);
        assert_eq!(features.min, &[0.0, 0.0, 0.0]);
        assert_eq!(features.max, &[2.0, 2.0, 4.0]);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(
            points,
            &[
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
                54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                80
            ]
        );
        Ok(())
    }
}
