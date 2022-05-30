use super::{Edge, EdgeKey, EdgesCellsMap2D, Face, FaceKey, FacesCellsMap3D, Mesh, PointId};
use crate::{shapes::Shape, StrError};
use russell_lab::sort2;
use std::collections::{HashMap, HashSet};

/// Holds extracted mesh features such as points, edges and faces on the mesh boundary or interior
///
/// # Examples
///
/// ## Two-dimensional
///
/// ```
/// use gemlab::mesh::{all_edges_2d, allocate_shapes, Cell, Extract, Features, Mesh, Point};
/// use gemlab::StrError;
///
/// fn main() -> Result<(), StrError> {
///     //  3---------2---------5
///     //  |         |         |
///     //  |   [0]   |   [1]   |
///     //  |         |         |
///     //  0---------1---------4
///     let mesh = Mesh {
///         space_ndim: 2,
///         points: vec![
///             Point { id: 0, coords: vec![0.0, 0.0] },
///             Point { id: 1, coords: vec![1.0, 0.0] },
///             Point { id: 2, coords: vec![1.0, 1.0] },
///             Point { id: 3, coords: vec![0.0, 1.0] },
///             Point { id: 4, coords: vec![2.0, 0.0] },
///             Point { id: 5, coords: vec![2.0, 1.0] },
///         ],
///         cells: vec![
///             Cell { id: 0, attribute_id: 1, geo_ndim: 2, points: vec![0, 1, 2, 3] },
///             Cell { id: 1, attribute_id: 2, geo_ndim: 2, points: vec![1, 4, 5, 2] },
///         ],
///     };
///
///     let shapes = allocate_shapes(&mesh)?;
///     let with_internal = false;
///     let edges = all_edges_2d(&mesh, &shapes)?;
///     let boundary = Features::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;
///
///     let mut points: Vec<_> = boundary.points.iter().copied().collect();
///     points.sort();
///     assert_eq!(points, [0, 1, 2, 3, 4, 5]);
///
///     let mut edges: Vec<_> = boundary.edges.keys().copied().collect();
///     edges.sort();
///     assert_eq!(edges, [(0, 1), (0, 3), (1, 4), (2, 3), (2, 5), (4, 5)]);
///     Ok(())
/// }
/// ```
///
/// ## Three-dimensional
///
/// ```
/// use gemlab::mesh::{all_faces_3d, allocate_shapes, Cell, Extract, Features, Mesh, Point};
/// use gemlab::StrError;
///
/// fn main() -> Result<(), StrError> {
///     //          .4--------------7
///     //        ,' |            ,'|
///     //      ,'              ,'  |
///     //    ,'     |        ,'    |
///     //  5'==============6'      |
///     //  |               |       |
///     //  |        |      |       |
///     //  |       ,0- - - | - - - 3
///     //  |     ,'        |     ,'
///     //  |   ,'          |   ,'
///     //  | ,'            | ,'
///     //  1'--------------2'
///     let mesh = Mesh {
///         space_ndim: 3,
///         points: vec![
///             Point { id: 0, coords: vec![0.0, 0.0, 0.0] },
///             Point { id: 1, coords: vec![1.0, 0.0, 0.0] },
///             Point { id: 2, coords: vec![1.0, 1.0, 0.0] },
///             Point { id: 3, coords: vec![0.0, 1.0, 0.0] },
///             Point { id: 4, coords: vec![0.0, 0.0, 1.0] },
///             Point { id: 5, coords: vec![1.0, 0.0, 1.0] },
///             Point { id: 6, coords: vec![1.0, 1.0, 1.0] },
///             Point { id: 7, coords: vec![0.0, 1.0, 1.0] },
///         ],
///         cells: vec![
///             Cell { id: 0, attribute_id: 1, geo_ndim: 3, points: vec![0,1,2,3, 4,5,6,7] },
///         ],
///     };
///
///     let shapes = allocate_shapes(&mesh)?;
///     let faces = all_faces_3d(&mesh, &shapes)?;
///     let boundary = Features::extract(&mesh, &shapes, None, Some(&faces), Extract::Boundary)?;
///
///     let mut points: Vec<_> = boundary.points.iter().copied().collect();
///     points.sort();
///     assert_eq!(points, (0..8).collect::<Vec<_>>());
///
///     let mut edges: Vec<_> = boundary.edges.keys().copied().collect();
///     edges.sort();
///     assert_eq!(
///         edges,
///         [
///             (0, 1),
///             (0, 3),
///             (0, 4),
///             (1, 2),
///             (1, 5),
///             (2, 3),
///             (2, 6),
///             (3, 7),
///             (4, 5),
///             (4, 7),
///             (5, 6),
///             (6, 7)
///         ]
///     );
///
///     let mut faces: Vec<_> = boundary.faces.keys().copied().collect();
///     faces.sort();
///     assert_eq!(
///         faces,
///         [
///             (0, 1, 2, 3),
///             (0, 1, 4, 5),
///             (0, 3, 4, 7),
///             (1, 2, 5, 6),
///             (2, 3, 6, 7),
///             (4, 5, 6, 7),
///         ]
///     );
///     Ok(())
/// }
/// ```
pub struct Features {
    /// Set of points on the boundary edges/faces or on interior edges/faces
    ///
    /// **Notes:**
    ///
    /// 1. A boundary point belongs to a boundary edge or a boundary face
    /// 2. The interior case means that the points on interior edges/faces are stored;
    ///    not all interior points are stored though (e.g. the middle nodes of a Qua9)
    pub points: HashSet<PointId>,

    /// Set of edges on the mesh boundary or interior
    ///
    /// **Notes:**
    ///
    /// 1. In 2D, a boundary edge is such that it is shared by one 2D cell only (1D cells are ignored)
    /// 2. In 3D, a boundary edge belongs to a boundary face
    pub edges: HashMap<EdgeKey, Edge>,

    /// Set of faces on the mesh boundary or interior
    ///
    /// **Notes:**
    ///
    /// 1. A boundary face is such that it is shared by one 3D cell only (2D cells are ignored)
    pub faces: HashMap<FaceKey, Face>,

    /// The minimum coordinates of the features collected here (len = space_ndim)
    pub min: Vec<f64>,

    /// The maximum coordinates of the features collected here (len = space_ndim)
    pub max: Vec<f64>,
}

/// Defines what features are to be extracted
pub enum Extract {
    /// Extracts boundary and interior features
    All,

    /// Extracts boundary features only
    Boundary,

    /// Extracts interior features only
    Interior,
}

impl Features {
    /// Extracts features (points, edges, faces)
    ///
    /// **Note:** The points of rods/shells-in-3D are only extracted with All or Boundary.
    ///
    /// # Input
    ///
    /// * `mesh` -- the Mesh
    /// * `shapes` -- all shapes of cells (len == cells.len())
    /// * `edges` -- 2D only: all edges (boundary and interior)
    /// * `faces` -- 3D only: all edges (boundary and interior)
    /// * `extract` -- option regarding what to extract (All and Boundary will extract points of rods and shells)
    pub fn extract(
        mesh: &Mesh,
        shapes: &Vec<Shape>,
        edges: Option<&EdgesCellsMap2D>,
        faces: Option<&FacesCellsMap3D>,
        extract: Extract,
    ) -> Result<Self, StrError> {
        let mut features = match mesh.space_ndim {
            2 => match edges {
                Some(edges) => extract_features_2d(mesh, shapes, &edges, &extract),
                None => return Err("edges map must be given in 2D"),
            },
            3 => match faces {
                Some(faces) => extract_features_3d(mesh, shapes, &faces, &extract),
                None => return Err("faces map must be given in 3D"),
            },
            _ => return Err("space_ndim must be 2 or 3"),
        };
        match extract {
            Extract::All => extract_points_of_rods_or_shells(mesh, &mut features),
            Extract::Boundary => extract_points_of_rods_or_shells(mesh, &mut features),
            Extract::Interior => (),
        }
        Ok(features)
    }
}

/// Extracts points of rods (2D or 3D) or shells (3D)
#[inline]
fn extract_points_of_rods_or_shells(mesh: &Mesh, features: &mut Features) {
    mesh.cells.iter().for_each(|cell| {
        if cell.geo_ndim == 1 || (cell.geo_ndim == 2 && mesh.space_ndim == 3) {
            cell.points.iter().for_each(|id| {
                features.points.insert(*id);
                for j in 0..mesh.space_ndim {
                    if mesh.points[*id].coords[j] < features.min[j] {
                        features.min[j] = mesh.points[*id].coords[j];
                    }
                    if mesh.points[*id].coords[j] > features.max[j] {
                        features.max[j] = mesh.points[*id].coords[j];
                    }
                }
            });
        }
    });
}

/// Extracts mesh features in 2D
fn extract_features_2d(mesh: &Mesh, shapes: &Vec<Shape>, edges: &EdgesCellsMap2D, extract: &Extract) -> Features {
    assert_eq!(mesh.space_ndim, 2);

    // output
    let mut features = Features {
        points: HashSet::new(),
        edges: HashMap::new(),
        faces: HashMap::new(),
        min: vec![f64::MAX; mesh.space_ndim],
        max: vec![f64::MIN; mesh.space_ndim],
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
        let shape = &shapes[cell_id];
        let mut edge = Edge {
            points: vec![0; shape.edge_nnode],
        };

        // process points on edge
        for i in 0..shape.edge_nnode {
            edge.points[i] = cell.points[shape.edge_node_id(e, i)];
            features.points.insert(edge.points[i]);
            for j in 0..mesh.space_ndim {
                features.min[j] = f64::min(features.min[j], mesh.points[edge.points[i]].coords[j]);
                features.max[j] = f64::max(features.max[j], mesh.points[edge.points[i]].coords[j]);
            }
        }

        // new edge
        features.edges.insert(*edge_key, edge);
    }
    features
}

/// Extracts mesh features in 3D
fn extract_features_3d(mesh: &Mesh, shapes: &Vec<Shape>, faces: &FacesCellsMap3D, extract: &Extract) -> Features {
    assert_eq!(mesh.space_ndim, 3);

    // output
    let mut features = Features {
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
        let shape = &shapes[cell_id];
        let mut face = Face {
            points: vec![0; shape.face_nnode],
        };

        // process points on face
        for i in 0..shape.face_nnode {
            face.points[i] = cell.points[shape.face_node_id(f, i)];
            features.points.insert(face.points[i]);
            for j in 0..mesh.space_ndim {
                features.min[j] = f64::min(features.min[j], mesh.points[face.points[i]].coords[j]);
                features.max[j] = f64::max(features.max[j], mesh.points[face.points[i]].coords[j]);
            }
        }

        // loop over all edges on face
        let face_shape = Shape::new(mesh.space_ndim, 2, shape.face_nnode).unwrap(); // should not fail
        for e in 0..face_shape.nedge {
            // define edge key (sorted point ids)
            let mut edge_key: EdgeKey = (
                face.points[face_shape.edge_node_id(e, 0)],
                face.points[face_shape.edge_node_id(e, 1)],
            );
            sort2(&mut edge_key);

            // skip already handled edge
            if features.edges.contains_key(&edge_key) {
                continue;
            }

            // new edge
            let mut edge = Edge {
                points: vec![0; face_shape.edge_nnode],
            };
            for i in 0..face_shape.edge_nnode {
                edge.points[i] = face.points[face_shape.edge_node_id(e, i)];
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
    use super::{Extract, Features};
    use crate::mesh::{all_edges_2d, all_faces_3d, allocate_shapes, EdgeKey, FaceKey, Mesh, PointId, Samples};
    use crate::util::AsArray2D;
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;

    #[test]
    fn capture_some_wrong_input() {
        let mesh = Mesh {
            space_ndim: 1,
            points: Vec::new(),
            cells: Vec::new(),
        };
        assert_eq!(
            Features::extract(&mesh, &Vec::new(), None, None, Extract::All).err(),
            Some("space_ndim must be 2 or 3")
        );
    }

    fn validate_edges<'a, T>(
        boundary: &Features,
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
        boundary: &Features,
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
        let shapes = allocate_shapes(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let boundary = Features::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;
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
    fn boundary_2d_all_works() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let shapes = allocate_shapes(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let boundary = Features::extract(&mesh, &shapes, Some(&edges), None, Extract::All)?;
        let correct_keys = [(0, 1), (0, 3), (1, 2), (1, 4), (2, 3), (2, 5), (4, 5)];
        let correct_points = [[1, 0], [0, 3], [2, 1], [4, 1], [3, 2], [2, 5], [5, 4]];
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
        //  0--------1---------2--------5
        let mesh = Samples::mixed_shapes_2d();
        let shapes = allocate_shapes(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let boundary = Features::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;
        let correct_keys = [(1, 2), (1, 4), (2, 3), (3, 4)];
        let correct_points = [[2, 1], [1, 4], [3, 2], [4, 3]];
        validate_edges(&boundary, &correct_keys, &correct_points);
        assert_eq!(boundary.min, &[0.0, 0.0]);
        assert_eq!(boundary.max, &[3.0, 1.0]);
        let mut points: Vec<_> = boundary.points.iter().map(|id| *id).collect();
        points.sort();
        assert_eq!(points, &[0, 1, 2, 3, 4, 5]);
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
        let shapes = allocate_shapes(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let boundary = Features::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;
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
    fn boundary_2d_all_block_works() -> Result<(), StrError> {
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
        let shapes = allocate_shapes(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let boundary = Features::extract(&mesh, &shapes, Some(&edges), None, Extract::All)?;
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
        validate_edges(&boundary, &correct_keys, &correct_points);
        assert_eq!(boundary.min, &[0.0, 0.0]);
        assert_eq!(boundary.max, &[3.0, 3.0]);
        let mut points: Vec<_> = boundary.points.iter().map(|id| *id).collect();
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
        let shapes = allocate_shapes(&mesh)?;
        let edges = all_edges_2d(&mesh, &shapes)?;
        let boundary = Features::extract(&mesh, &shapes, Some(&edges), None, Extract::Boundary)?;
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
        let shapes = allocate_shapes(&mesh)?;
        let faces = all_faces_3d(&mesh, &shapes)?;
        let boundary = Features::extract(&mesh, &shapes, None, Some(&faces), Extract::Boundary)?;
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
    fn boundary_3d_all_works() -> Result<(), StrError> {
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
        let faces = all_faces_3d(&mesh, &shapes)?;
        let boundary = Features::extract(&mesh, &shapes, None, Some(&faces), Extract::All)?;
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
        let shapes = allocate_shapes(&mesh)?;
        let faces = all_faces_3d(&mesh, &shapes)?;
        let boundary = Features::extract(&mesh, &shapes, None, Some(&faces), Extract::Boundary)?;
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
        let shapes = allocate_shapes(&mesh)?;
        let faces = all_faces_3d(&mesh, &shapes)?;
        let boundary = Features::extract(&mesh, &shapes, None, Some(&faces), Extract::Boundary)?;
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

    #[test]
    fn boundary_3d_all_block_works() -> Result<(), StrError> {
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
        let shapes = allocate_shapes(&mesh)?;
        let faces = all_faces_3d(&mesh, &shapes)?;
        let boundary = Features::extract(&mesh, &shapes, None, Some(&faces), Extract::All)?;
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
        validate_edges(&boundary, &correct_edge_keys, &correct_edge_points);
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
        validate_faces(&boundary, &correct_face_keys, &correct_face_points);
        assert_eq!(boundary.min, &[0.0, 0.0, 0.0]);
        assert_eq!(boundary.max, &[2.0, 2.0, 4.0]);
        let mut points: Vec<_> = boundary.points.iter().map(|id| *id).collect();
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
