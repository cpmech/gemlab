use crate::mesh::{Edge, Extract, Features, MapEdge2dToCells, Mesh};
use std::collections::{HashMap, HashSet};

/// Extracts mesh features in 2D
///
/// # Panics
///
/// 1. It panics if `mesh.ndim != 2` (i.e., this function works in 2D only)
pub(crate) fn extract_features_2d(mesh: &Mesh, edges: &MapEdge2dToCells, extract: Extract) -> Features {
    assert_eq!(mesh.ndim, 2);

    // output
    let mut features = Features {
        points: HashSet::new(),
        edges: HashMap::new(),
        faces: HashMap::new(),
        min: vec![f64::MAX; mesh.ndim],
        max: vec![f64::MIN; mesh.ndim],
    };

    // loop over all edges
    for (edge_key, shared_by) in edges {
        // accept feature?
        let accept = match extract {
            Extract::All => true,
            Extract::Boundary => shared_by.len() == 1, // boundary edge; with only one shared cell
            Extract::Interior => shared_by.len() != 1, // interior edge; shared by multiple cells
        };
        if !accept {
            continue;
        }

        // cell and edge
        let (cell_id, e) = shared_by[0];
        let cell = &mesh.cells[cell_id];
        let mut edge = Edge {
            kind: cell.kind.edge_kind().unwrap(),
            points: vec![0; cell.kind.edge_nnode()],
        };

        // process points on edge
        for i in 0..edge.points.len() {
            edge.points[i] = cell.points[cell.kind.edge_node_id(e, i)];
            features.points.insert(edge.points[i]);
            for j in 0..mesh.ndim {
                features.min[j] = f64::min(features.min[j], mesh.points[edge.points[i]].coords[j]);
                features.max[j] = f64::max(features.max[j], mesh.points[edge.points[i]].coords[j]);
            }
        }

        // new edge
        features.edges.insert(*edge_key, edge);
    }
    features
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::extract_features_2d;
    use crate::mesh::algorithms::extract_all_2d_edges;
    use crate::mesh::{EdgeKey, Extract, Features, PointId, Samples};
    use crate::util::AsArray2D;
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;

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

    #[test]
    fn extract_features_2d_works() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let edges = extract_all_2d_edges(&mesh);
        let features = extract_features_2d(&mesh, &edges, Extract::Boundary);
        let correct_keys = [(0, 1), (0, 3), (1, 4), (2, 3), (2, 5), (4, 5)];
        let correct_points = [[1, 0], [0, 3], [4, 1], [3, 2], [2, 5], [5, 4]];
        validate_edges(&features, &correct_keys, &correct_points);
        assert_eq!(features.min, &[0.0, 0.0]);
        assert_eq!(features.max, &[2.0, 1.0]);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[0, 1, 2, 3, 4, 5]);
        Ok(())
    }

    #[test]
    fn extract_features_2d_all_works() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let edges = extract_all_2d_edges(&mesh);
        let features = extract_features_2d(&mesh, &edges, Extract::All);
        let correct_keys = [(0, 1), (0, 3), (1, 2), (1, 4), (2, 3), (2, 5), (4, 5)];
        let correct_points = [[1, 0], [0, 3], [2, 1], [4, 1], [3, 2], [2, 5], [5, 4]];
        validate_edges(&features, &correct_keys, &correct_points);
        assert_eq!(features.min, &[0.0, 0.0]);
        assert_eq!(features.max, &[2.0, 1.0]);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[0, 1, 2, 3, 4, 5]);
        Ok(())
    }

    #[test]
    fn extract_features_2d_works_interior() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let edges = extract_all_2d_edges(&mesh);
        let features = extract_features_2d(&mesh, &edges, Extract::Interior);
        let correct_keys = [(1, 2)];
        let correct_points = [[2, 1]];
        validate_edges(&features, &correct_keys, &correct_points);
        assert_eq!(features.min, &[1.0, 0.0]);
        assert_eq!(features.max, &[1.0, 1.0]);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[1, 2]);
        Ok(())
    }

    #[test]
    fn extract_features_2d_mixed_works() -> Result<(), StrError> {
        //           4---------3
        //           |         |
        //           |   [1]   |
        //           |         |
        //  0--------1---------2--------5
        let mesh = Samples::mixed_shapes_2d();
        let edges = extract_all_2d_edges(&mesh);
        let features = extract_features_2d(&mesh, &edges, Extract::Boundary);
        let correct_keys = [(1, 2), (1, 4), (2, 3), (3, 4)];
        let correct_points = [[2, 1], [1, 4], [3, 2], [4, 3]];
        validate_edges(&features, &correct_keys, &correct_points);
        assert_eq!(features.min, &[1.0, 0.0]); // note that Lin is ignored
        assert_eq!(features.max, &[2.0, 1.0]); // note that Lin is ignored
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[1, 2, 3, 4,]); // note that Lin is ignored
        Ok(())
    }

    #[test]
    fn extract_features_2d_block_works() -> Result<(), StrError> {
        // 29---34----31---28---44---42----40
        //  |               |               |
        // 32   39    38   33   48   47    43
        //  |               |               |
        // 35   36    37   30   45   46    41
        //  |               |               |
        //  3---10-----6----2---23---20----17
        //  |               |               |
        //  7   15    14    9   27   26    22
        //  |               |               |
        // 11   12    13    5   24   25    19
        //  |               |               |
        //  0----4-----8----1---18---21----16
        let mesh = Samples::block_2d_four_qua16();
        let edges = extract_all_2d_edges(&mesh);
        let features = extract_features_2d(&mesh, &edges, Extract::Boundary);
        let correct_keys = [(0, 1), (0, 3), (1, 16), (3, 29), (16, 17), (17, 40), (28, 29), (28, 40)];
        let correct_points = [
            [1, 0, 8, 4],
            [0, 3, 11, 7],
            [16, 1, 21, 18],
            [3, 29, 35, 32],
            [17, 16, 22, 19],
            [40, 17, 43, 41],
            [29, 28, 34, 31],
            [28, 40, 44, 42],
        ];
        validate_edges(&features, &correct_keys, &correct_points);
        assert_eq!(features.min, &[0.0, 0.0]);
        assert_eq!(features.max, &[3.0, 3.0]);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(
            points,
            &[0, 1, 3, 4, 7, 8, 11, 16, 17, 18, 19, 21, 22, 28, 29, 31, 32, 34, 35, 40, 41, 42, 43, 44]
        );
        Ok(())
    }

    #[test]
    fn extract_features_2d_all_block_works() -> Result<(), StrError> {
        // 29---34----31---28---44---42----40
        //  |               |               |
        // 32   39    38   33   48   47    43
        //  |               |               |
        // 35   36    37   30   45   46    41
        //  |               |               |
        //  3---10-----6----2---23---20----17
        //  |               |               |
        //  7   15    14    9   27   26    22
        //  |               |               |
        // 11   12    13    5   24   25    19
        //  |               |               |
        //  0----4-----8----1---18---21----16
        let mesh = Samples::block_2d_four_qua16();
        let edges = extract_all_2d_edges(&mesh);
        let features = extract_features_2d(&mesh, &edges, Extract::All);
        let correct_keys = [
            (0, 1),
            (0, 3),
            (1, 2),
            (1, 16),
            (2, 3),
            (2, 17),
            (2, 28),
            (3, 29),
            (16, 17),
            (17, 40),
            (28, 29),
            (28, 40),
        ];
        let correct_points = [
            [1, 0, 8, 4],
            [0, 3, 11, 7],
            [2, 1, 9, 5],
            [16, 1, 21, 18],
            [3, 2, 10, 6],
            [2, 17, 23, 20],
            [28, 2, 33, 30],
            [3, 29, 35, 32],
            [17, 16, 22, 19],
            [40, 17, 43, 41],
            [29, 28, 34, 31],
            [28, 40, 44, 42],
        ];
        validate_edges(&features, &correct_keys, &correct_points);
        assert_eq!(features.min, &[0.0, 0.0]);
        assert_eq!(features.max, &[3.0, 3.0]);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(
            points,
            &[
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 22, 23, 28, 29, 30, 31, 32, 33, 34, 35,
                40, 41, 42, 43, 44
            ]
        );
        Ok(())
    }

    #[test]
    fn extract_features_2d_ring_works() -> Result<(), StrError> {
        // 2.0   14---36--,__11
        //        |          `,-..33
        // 1.75  24   [7]   22     `-,
        //        |         ,  [5]    ,8.
        // 1.5   13--35--10/        20   `*
        //        |       ,`*32    ,'      30
        // 1.25  23 [6] 21     *.7     [3]   *
        //        |     ,  [4]  , *.          5
        // 1.0   12-34-9      19    29     18' *
        //              `31. ,' [2]   *  _,     *
        //                  6.       _.4'        *
        //                   28  _.17   *   [1]  27
        //                     3'  [0]  26        *
        //                     25        *        *
        //        +             0---15---1---16---2
        //
        //                     1.0 1.25  1.5 1.75  2.0
        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        let edges = extract_all_2d_edges(&mesh);
        let features = extract_features_2d(&mesh, &edges, Extract::Boundary);
        let correct_keys = [
            (0, 1),
            (0, 3),
            (1, 2),
            (2, 5),
            (3, 6),
            (5, 8),
            (6, 9),
            (8, 11),
            (9, 12),
            (11, 14),
            (12, 13),
            (13, 14),
        ];
        let correct_points = [
            [1, 0, 15],
            [0, 3, 25],
            [2, 1, 16],
            [5, 2, 27],
            [3, 6, 28],
            [8, 5, 30],
            [6, 9, 31],
            [11, 8, 33],
            [9, 12, 34],
            [14, 11, 36],
            [12, 13, 23],
            [13, 14, 24],
        ];
        validate_edges(&features, &correct_keys, &correct_points);
        assert_vec_approx_eq!(features.min, &[0.0, 0.0], 1e-15);
        assert_vec_approx_eq!(features.max, &[2.0, 2.0], 1e-15);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(
            points,
            &[0, 1, 2, 3, 5, 6, 8, 9, 11, 12, 13, 14, 15, 16, 23, 24, 25, 27, 28, 30, 31, 33, 34, 36,]
        );
        Ok(())
    }

    #[test]
    fn extract_features_2d_ring_interior_works() -> Result<(), StrError> {
        // 2.0   14---36--,__11
        //        |          `,-..33
        // 1.75  24   [7]   22     `-,
        //        |         ,  [5]    ,8.
        // 1.5   13--35--10/        20   `*
        //        |       ,`*32    ,'      30
        // 1.25  23 [6] 21     *.7     [3]   *
        //        |     ,  [4]  , *.          5
        // 1.0   12-34-9      19    29     18' *
        //              `31. ,' [2]   *  _,     *
        //                  6.       _.4'        *
        //                   28  _.17   *   [1]  27
        //                     3'  [0]  26        *
        //                     25        *        *
        //        +             0---15---1---16---2
        //
        //                     1.0 1.25  1.5 1.75  2.0
        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        let edges = extract_all_2d_edges(&mesh);
        let features = extract_features_2d(&mesh, &edges, Extract::Interior);
        let correct_keys = [
            (1, 4),
            (3, 4),
            (4, 5),
            (4, 7),
            (6, 7),
            (7, 8),
            (7, 10),
            (9, 10),
            (10, 11),
            (10, 13),
        ];
        let correct_points = [
            [4, 1, 26],
            [3, 4, 17],
            [4, 5, 18],
            [7, 4, 29],
            [6, 7, 19],
            [7, 8, 20],
            [10, 7, 32],
            [9, 10, 21],
            [10, 11, 22],
            [13, 10, 35],
        ];
        validate_edges(&features, &correct_keys, &correct_points);
        let mut points: Vec<_> = features.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(
            points,
            &[1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 17, 18, 19, 20, 21, 22, 26, 29, 32, 35]
        );
        Ok(())
    }
}
