use super::{all_edges_2d, all_faces_3d, CellId, Edge, EdgeKey, Face, FaceKey, Mesh, PointId};
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
    /// Note: A boundary face is such that it is shared by one 3D cell only (2D cells are ignored)
    pub faces: HashMap<FaceKey, Face>,

    /// The minimum coordinates; len = space_ndim
    pub min: Vec<f64>,

    /// The maximum coordinates; len = space_ndim
    pub max: Vec<f64>,
}

impl Boundary {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, shapes: &Vec<Shape>) -> Result<Self, StrError> {
        let mut boundary = match mesh.space_ndim {
            2 => {
                let edges = all_edges_2d(mesh, shapes)?;
                Boundary::two_dim(mesh, shapes, &edges)?
            }
            3 => {
                let faces = all_faces_3d(mesh, shapes)?;
                Boundary::three_dim(mesh, shapes, &faces)?
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
    use crate::mesh::{alloc_cell_shapes, EdgeKey, FaceKey, PointId, Samples};
    use crate::util::AsArray2D;
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;

    fn validate_edges<'a, T>(
        boundary: &Boundary,
        correct_keys: &[EdgeKey], // sorted
        correct_points: &'a T,
    ) where
        T: AsArray2D<'a, PointId>,
    {
        let mut keys: Vec<_> = boundary.edges.keys().map(|k| *k).collect();
        keys.sort();
        assert_eq!(keys, correct_keys);
        for i in 0..keys.len() {
            assert_eq!(boundary.edges.get(&keys[i]).unwrap().points, correct_points.row(i));
        }
    }

    fn validate_faces<'a, T>(
        boundary: &Boundary,
        correct_keys: &[FaceKey], // sorted
        correct_points: &'a T,
    ) where
        T: AsArray2D<'a, PointId>,
    {
        let mut keys: Vec<_> = boundary.faces.keys().map(|k| *k).collect();
        keys.sort();
        assert_eq!(keys, correct_keys);
        for i in 0..keys.len() {
            assert_eq!(boundary.faces.get(&keys[i]).unwrap().points, correct_points.row(i));
        }
    }

    #[test]
    fn boundary_2d_works() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let shapes = alloc_cell_shapes(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;
        let correct_keys = [(0, 1), (0, 3), (1, 4), (2, 3), (2, 5), (4, 5)];
        let correct_points = [[1, 0], [0, 3], [4, 1], [3, 2], [2, 5], [5, 4]];
        validate_edges(&boundary, &correct_keys, &correct_points);
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
        let shapes = alloc_cell_shapes(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;
        let correct_keys = [(1, 2), (1, 4), (2, 3), (3, 4)];
        let correct_points = [[2, 1], [1, 4], [3, 2], [4, 3]];
        validate_edges(&boundary, &correct_keys, &correct_points);
        assert_eq!(boundary.min, &[0.0, 0.0]);
        assert_eq!(boundary.max, &[2.0, 1.0]);
        let mut points: Vec<_> = boundary.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[0, 1, 2, 3, 4]);
        Ok(())
    }

    #[test]
    fn boundary_2d_block_works() -> Result<(), StrError> {
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
        let shapes = alloc_cell_shapes(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;
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
        validate_edges(&boundary, &correct_keys, &correct_points);
        assert_eq!(boundary.min, &[0.0, 0.0]);
        assert_eq!(boundary.max, &[3.0, 3.0]);
        let mut points: Vec<_> = boundary.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(
            points,
            &[0, 1, 3, 4, 7, 8, 11, 16, 17, 18, 19, 21, 22, 28, 29, 31, 32, 34, 35, 40, 41, 42, 43, 44]
        );
        Ok(())
    }

    #[test]
    fn boundary_2d_ring_works() -> Result<(), StrError> {
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
        let shapes = alloc_cell_shapes(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;
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
        validate_edges(&boundary, &correct_keys, &correct_points);
        assert_vec_approx_eq!(boundary.min, &[0.0, 0.0], 1e-15);
        assert_vec_approx_eq!(boundary.max, &[2.0, 2.0], 1e-15);
        let mut points: Vec<_> = boundary.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(
            points,
            &[0, 1, 2, 3, 5, 6, 8, 9, 11, 12, 13, 14, 15, 16, 23, 24, 25, 27, 28, 30, 31, 33, 34, 36,]
        );
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
        let shapes = alloc_cell_shapes(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;
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
        validate_edges(&boundary, &correct_edge_keys, &correct_edge_points);
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
        validate_faces(&boundary, &correct_face_keys, &correct_face_points);
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
        let shapes = alloc_cell_shapes(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;
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
        validate_edges(&boundary, &correct_edge_keys, &correct_edge_points);
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
        validate_faces(&boundary, &correct_face_keys, &correct_face_points);
        assert_eq!(boundary.min, &[0.0, -1.0, 0.0]);
        assert_eq!(boundary.max, &[1.0, 2.0, 1.0]);
        let mut points: Vec<_> = boundary.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
        Ok(())
    }

    #[test]
    fn boundary_3d_block_works() -> Result<(), StrError> {
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
        let shapes = alloc_cell_shapes(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;
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
        validate_edges(&boundary, &correct_edge_keys, &correct_edge_points);
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
        validate_faces(&boundary, &correct_face_keys, &correct_face_points);
        assert_eq!(boundary.min, &[0.0, 0.0, 0.0]);
        assert_eq!(boundary.max, &[2.0, 2.0, 4.0]);
        let mut points: Vec<_> = boundary.points.iter().map(|id| *id).collect();
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
}
