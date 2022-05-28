use super::{Boundary, EdgeKey, FaceKey, Mesh};
use crate::shapes::{Shape, StateOfShape};
use crate::StrError;
use russell_lab::Vector;

/// Holds data to evaluate normal vectors at edge or face
///
/// # Examples
///
/// ## Two-dimensional
///
/// ```
/// use gemlab::mesh::{Boundary, Cell, Mesh, NormalVector, Point, Shapes};
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
///     let shapes = Shapes::new(&mesh)?;
///     let boundary = Boundary::new(&mesh, &shapes)?;
///
///     // the magnitude of the normal vector is equal to
///     // 0.5 = edge_length / 2.0 where 2.0 corresponds to
///     // the edge_length in the reference system
///
///     let mut normal = NormalVector::at_edge(&mesh, &boundary, (0, 1))?;
///     normal.evaluate(&[0.0, 0.0])?;
///     assert_eq!(normal.value.as_data(), &[0.0, -0.5]);
///
///     let mut normal = NormalVector::at_edge(&mesh, &boundary, (4, 5))?;
///     normal.evaluate(&[0.0, 0.0])?;
///     assert_eq!(normal.value.as_data(), &[0.5, 0.0]);
///
///     let mut normal = NormalVector::at_edge(&mesh, &boundary, (2, 5))?;
///     normal.evaluate(&[0.0, 0.0])?;
///     assert_eq!(normal.value.as_data(), &[0.0, 0.5]);
///
///     let mut normal = NormalVector::at_edge(&mesh, &boundary, (0, 3))?;
///     normal.evaluate(&[0.0, 0.0])?;
///     assert_eq!(normal.value.as_data(), &[-0.5, 0.0]);
///     Ok(())
/// }
/// ```
///
/// ## Three-dimensional
///
/// ```
/// use gemlab::mesh::{Boundary, Cell, Mesh, NormalVector, Point, Shapes};
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
///     let shapes = Shapes::new(&mesh)?;
///     let boundary = Boundary::new(&mesh, &shapes)?;
///
///     // the magnitude of the normal vector is equal to
///     // 0.25 = face_area / 4.0 where 4.0 corresponds to
///     // the face_area in the reference system
///
///     let mut normal = NormalVector::at_face(&mesh, &boundary, (0, 3, 4, 7))?;
///     normal.evaluate(&[0.0, 0.0, 0.0])?;
///     assert_eq!(normal.value.as_data(), &[-0.25, 0.0, 0.0]);
///
///     let mut normal = NormalVector::at_face(&mesh, &boundary, (1, 2, 5, 6))?;
///     normal.evaluate(&[0.0, 0.0, 0.0])?;
///     assert_eq!(normal.value.as_data(), &[0.25, 0.0, 0.0]);
///
///     let mut normal = NormalVector::at_face(&mesh, &boundary, (0, 1, 4, 5))?;
///     normal.evaluate(&[0.0, 0.0, 0.0])?;
///     assert_eq!(normal.value.as_data(), &[0.0, -0.25, 0.0]);
///
///     let mut normal = NormalVector::at_face(&mesh, &boundary, (2, 3, 6, 7))?;
///     normal.evaluate(&[0.0, 0.0, 0.0])?;
///     assert_eq!(normal.value.as_data(), &[0.0, 0.25, 0.0]);
///
///     let mut normal = NormalVector::at_face(&mesh, &boundary, (0, 1, 2, 3))?;
///     normal.evaluate(&[0.0, 0.0, 0.0])?;
///     assert_eq!(normal.value.as_data(), &[0.0, 0.0, -0.25]);
///
///     let mut normal = NormalVector::at_face(&mesh, &boundary, (4, 5, 6, 7))?;
///     normal.evaluate(&[0.0, 0.0, 0.0])?;
///     assert_eq!(normal.value.as_data(), &[0.0, 0.0, 0.25]);
///     Ok(())
/// }
/// ```
#[derive(Clone, Debug)]
pub struct NormalVector {
    /// Holds the Shape data for a given boundary edge or face
    pub shape: Shape,

    /// Holds the state for the Shape data corresponding to an edge or face
    pub state: StateOfShape,

    /// Holds the normal vector output from the evaluate function
    pub value: Vector,
}

impl NormalVector {
    /// Allocates data to compute normal vector at edge
    ///
    /// **Important:** You must call [NormalVector::evaluate] to compute the actual values.
    ///
    /// **Note:** This function works in 2D only
    pub fn at_edge(mesh: &Mesh, boundary: &Boundary, edge_key: EdgeKey) -> Result<Self, StrError> {
        if mesh.space_ndim != 2 {
            return Err("normal at_edge works in 2D only");
        }
        const GEO_NDIM: usize = 1;
        let edge = match boundary.edges.get(&edge_key) {
            Some(e) => e,
            None => return Err("edge_key is not present in boundary"),
        };
        let shape = Shape::new(mesh.space_ndim, GEO_NDIM, edge.points.len())?;
        let state = StateOfShape::new(
            GEO_NDIM,
            &edge
                .points
                .iter()
                .map(|id| mesh.points[*id].coords.clone())
                .collect::<Vec<_>>(),
        )?;
        Ok(NormalVector {
            shape,
            state,
            value: Vector::new(mesh.space_ndim),
        })
    }

    /// Allocates data to compute normal vector at face
    ///
    /// **Important:** You must call [NormalVector::evaluate] to compute the actual values.
    ///
    /// **Note:** This function works in 3D only
    pub fn at_face(mesh: &Mesh, boundary: &Boundary, face_key: FaceKey) -> Result<Self, StrError> {
        if mesh.space_ndim != 3 {
            return Err("normal at_face works in 3D only");
        }
        const GEO_NDIM: usize = 2;
        let face = match boundary.faces.get(&face_key) {
            Some(e) => e,
            None => return Err("face_key is not present in boundary"),
        };
        let shape = Shape::new(mesh.space_ndim, GEO_NDIM, face.points.len())?;
        let state = StateOfShape::new(
            GEO_NDIM,
            &face
                .points
                .iter()
                .map(|id| mesh.points[*id].coords.clone())
                .collect::<Vec<_>>(),
        )?;
        Ok(NormalVector {
            shape,
            state,
            value: Vector::new(mesh.space_ndim),
        })
    }

    /// Evaluates boundary normal
    ///
    /// # Input
    ///
    /// * `ksi` -- ξ reference coordinate. The length of ξ must be equal to geo_ndim at least,
    ///            while lengths greater than geo_ndim are allowed (and ignored). In this way,
    ///            we can pass a slice with integration point data such as `[f64; 4]`.
    pub fn evaluate(&mut self, ksi: &[f64]) -> Result<(), StrError> {
        self.shape.calc_boundary_normal(&mut self.value, &mut self.state, ksi)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::NormalVector;
    use crate::mesh::{Boundary, Edge, Face, Samples, Shapes};
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;
    use std::collections::{HashMap, HashSet};

    #[test]
    fn capture_some_wrong_input() {
        let mesh = Samples::two_quads_horizontal();
        let boundary = Boundary {
            points: HashSet::new(),
            edges: HashMap::new(),
            faces: HashMap::new(),
            min: Vec::new(),
            max: Vec::new(),
        };
        assert_eq!(
            NormalVector::at_edge(&mesh, &boundary, (0, 1)).err(),
            Some("edge_key is not present in boundary")
        );
        assert_eq!(
            NormalVector::at_face(&mesh, &boundary, (0, 1, 2, 3)).err(),
            Some("face_key is not present in boundary")
        );

        let boundary = Boundary {
            points: HashSet::new(),
            edges: HashMap::from([((0, 1), Edge { points: Vec::new() })]),
            faces: HashMap::from([((0, 1, 2, 3), Face { points: Vec::new() })]),
            min: Vec::new(),
            max: Vec::new(),
        };
        assert_eq!(
            NormalVector::at_edge(&mesh, &boundary, (0, 1)).err(),
            Some("(geo_ndim,nnode) combination is invalid")
        );
        assert_eq!(
            NormalVector::at_face(&mesh, &boundary, (0, 1, 2, 3)).err(),
            Some("(geo_ndim,nnode) combination is invalid")
        );
    }

    #[test]
    fn at_edge_and_evaluate_work() -> Result<(), StrError> {
        //  3---------2---------5
        //  |         |         |
        //  |   [0]   |   [1]   |
        //  |         |         |
        //  0---------1---------4
        let mesh = Samples::two_quads_horizontal();
        let shapes = Shapes::new(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;

        // the magnitude (l) of the normal vector should be equal to
        // 0.5 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
        let l = 0.5; // magnitude of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 4)], [0.0, -l]),
            // right
            (vec![(4, 5)], [l, 0.0]),
            // top
            (vec![(2, 3), (2, 5)], [0.0, l]),
            // left
            (vec![(0, 3)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let ksi = &[0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let mut n = NormalVector::at_edge(&mesh, &boundary, *edge_key)?;
                n.evaluate(ksi)?;
                assert_vec_approx_eq!(n.value.as_data(), solution, 1e-15);
            }
        }
        Ok(())
    }

    #[test]
    fn at_edge_and_evaluate_work_qua8() -> Result<(), StrError> {
        // 14------16------13------20------18
        //  |               |               |
        //  |               |               |
        // 17              15              19
        //  |               |               |
        //  |               |               |
        //  3-------6-------2------12-------9
        //  |               |               |
        //  |               |               |
        //  7               5              11
        //  |               |               |
        //  |               |               |
        //  0-------4-------1------10-------8
        let mesh = Samples::block_2d_four_qua8();
        let shapes = Shapes::new(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;

        // the magnitude (l) of the normal vector should be equal to
        // 0.5 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
        let l = 0.5; // magnitude of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 8)], [0.0, -l]),
            // right
            (vec![(8, 9), (9, 18)], [l, 0.0]),
            // top
            (vec![(13, 14), (13, 18)], [0.0, l]),
            // left
            (vec![(0, 3), (3, 14)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let ksi = &[0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let mut n = NormalVector::at_edge(&mesh, &boundary, *edge_key)?;
                n.evaluate(ksi)?;
                assert_vec_approx_eq!(n.value.as_data(), solution, 1e-15);
            }
        }
        Ok(())
    }

    #[test]
    fn at_edge_and_evaluate_work_qua9() -> Result<(), StrError> {
        // 16------18------15------23------21
        //  |               |               |
        //  |               |               |
        // 19      20      17      24      22
        //  |               |               |
        //  |               |               |
        //  3-------6-------2------13------10
        //  |               |               |
        //  |               |               |
        //  7       8       5      14      12
        //  |               |               |
        //  |               |               |
        //  0-------4-------1------11-------9
        let mesh = Samples::block_2d_four_qua9();
        let shapes = Shapes::new(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;

        // the magnitude (l) of the normal vector should be equal to
        // 0.5 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
        let l = 0.5; // magnitude of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 9)], [0.0, -l]),
            // right
            (vec![(9, 10), (10, 21)], [l, 0.0]),
            // top
            (vec![(15, 16), (15, 21)], [0.0, l]),
            // left
            (vec![(0, 3), (3, 16)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let ksi = &[0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let mut n = NormalVector::at_edge(&mesh, &boundary, *edge_key)?;
                n.evaluate(ksi)?;
                assert_vec_approx_eq!(n.value.as_data(), solution, 1e-14);
            }
        }
        Ok(())
    }

    #[test]
    fn at_edge_and_evaluate_work_qua12() -> Result<(), StrError> {
        // 21---26---23----20---32---30----28
        //  |               |               |
        // 24              25              31
        //  |               |               |
        // 27              22              29
        //  |               |               |
        //  3---10-----6----2---19---16----13
        //  |               |               |
        //  7               9              18
        //  |               |               |
        // 11               5              15
        //  |               |               |
        //  0----4-----8----1---14---17----12
        let mesh = Samples::block_2d_four_qua12();
        let shapes = Shapes::new(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;

        // the magnitude (l) of the normal vector should be equal to
        // 0.75 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
        let l = 0.75; // magnitude of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 12)], [0.0, -l]),
            // right
            (vec![(12, 13), (13, 28)], [l, 0.0]),
            // top
            (vec![(20, 21), (20, 28)], [0.0, l]),
            // left
            (vec![(0, 3), (3, 21)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let ksi = &[0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let mut n = NormalVector::at_edge(&mesh, &boundary, *edge_key)?;
                n.evaluate(ksi)?;
                assert_vec_approx_eq!(n.value.as_data(), solution, 1e-14);
            }
        }
        Ok(())
    }

    #[test]
    fn at_edge_and_evaluate_work_qua16() -> Result<(), StrError> {
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
        let shapes = Shapes::new(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;

        // the magnitude (l) of the normal vector should be equal to
        // 0.75 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
        let l = 0.75; // magnitude of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 16)], [0.0, -l]),
            // right
            (vec![(16, 17), (17, 40)], [l, 0.0]),
            // top
            (vec![(28, 29), (28, 40)], [0.0, l]),
            // left
            (vec![(0, 3), (3, 29)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let ksi = &[0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let mut n = NormalVector::at_edge(&mesh, &boundary, *edge_key)?;
                n.evaluate(ksi)?;
                assert_vec_approx_eq!(n.value.as_data(), solution, 1e-14);
            }
        }
        Ok(())
    }

    #[test]
    fn at_edge_and_evaluate_work_qua17() -> Result<(), StrError> {
        // 30---38---35---32---29---47---45---43---41
        //  |                   |                   |
        // 33                  37                  46
        //  |                   |                   |
        // 36        40        34        48        44
        //  |                   |                   |
        // 39                  31                  42
        //  |                   |                   |
        //  3---14---10----6----2---27---24---21---18
        //  |                   |                   |
        //  7                  13                  26
        //  |                   |                   |
        // 11        16         9        28        23
        //  |                   |                   |
        // 15                   5                  20
        //  |                   |                   |
        //  0----4----8---12----1---19---22---25---17
        let mesh = Samples::block_2d_four_qua17();
        let shapes = Shapes::new(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;

        // the magnitude (l) of the normal vector should be equal to
        // 1.0 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
        let l = 1.0; // magnitude of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 17)], [0.0, -l]),
            // right
            (vec![(17, 18), (18, 41)], [l, 0.0]),
            // top
            (vec![(29, 30), (29, 41)], [0.0, l]),
            // left
            (vec![(0, 3), (3, 30)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let ksi = &[0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let mut n = NormalVector::at_edge(&mesh, &boundary, *edge_key)?;
                n.evaluate(ksi)?;
                assert_vec_approx_eq!(n.value.as_data(), solution, 1e-14);
            }
        }
        Ok(())
    }

    #[test]
    fn at_face_and_evaluate_work() -> Result<(), StrError> {
        //      8-------------11
        //     /.             /|
        //    / .            / |
        //   /  .           /  |
        //  /   .          /   |
        // 9-------------10    |
        // |    .         |    |
        // |    4---------|----7
        // |   /.         |   /|
        // |  / .         |  / |
        // | /  .         | /  |
        // |/   .         |/   |
        // 5--------------6    |
        // |    .         |    |
        // |    0---------|----3
        // |   /          |   /
        // |  /           |  /
        // | /            | /
        // |/             |/
        // 1--------------2
        let mesh = Samples::two_cubes_vertical();
        let shapes = Shapes::new(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;

        // the magnitude (l) of the normal vector should be equal to
        // 0.25 = face_area / 4.0 where 4.0 corresponds to the face_area in the reference system
        let l = 0.25; // magnitude of normal vector

        // face keys and correct normal vectors (solutions)
        let face_keys_and_solutions = [
            // behind
            (vec![(0, 3, 4, 7), (4, 7, 8, 11)], [-l, 0.0, 0.0]),
            // front
            (vec![(1, 2, 5, 6), (5, 6, 9, 10)], [l, 0.0, 0.0]),
            // left
            (vec![(0, 1, 4, 5), (4, 5, 8, 9)], [0.0, -l, 0.0]),
            // right
            (vec![(2, 3, 6, 7), (6, 7, 10, 11)], [0.0, l, 0.0]),
            // bottom
            (vec![(0, 1, 2, 3)], [0.0, 0.0, -l]),
            // top
            (vec![(8, 9, 10, 11)], [0.0, 0.0, l]),
        ];

        let ksi = &[0.0, 0.0, 0.0];
        for (face_keys, solution) in &face_keys_and_solutions {
            for face_key in face_keys {
                let mut n = NormalVector::at_face(&mesh, &boundary, *face_key)?;
                n.evaluate(ksi)?;
                assert_vec_approx_eq!(n.value.as_data(), solution, 1e-15);
            }
        }
        Ok(())
    }

    #[test]
    fn at_face_and_evaluate_work_hex20() -> Result<(), StrError> {
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
        let shapes = Shapes::new(&mesh)?;
        let boundary = Boundary::new(&mesh, &shapes)?;

        // the magnitude (l) of the normal vector should be equal to
        // face_area / 4.0 where 4.0 corresponds to the face_area in the reference system

        // face keys and correct normal vectors (solutions)
        let face_keys_and_solutions = [
            (vec![(0, 1, 2, 3)], [0.0, 0.0, -1.0 / 4.0]),
            (vec![(52, 53, 63, 64)], [0.0, 0.0, 1.0 / 4.0]),
            (vec![(4, 7, 51, 54)], [-2.0 / 4.0, 0.0, 0.0]),
            (vec![(34, 45, 70, 77)], [0.0, 2.0 / 4.0, 0.0]),
        ];

        let ksi = &[0.0, 0.0, 0.0];
        for (face_keys, solution) in &face_keys_and_solutions {
            for face_key in face_keys {
                let mut n = NormalVector::at_face(&mesh, &boundary, *face_key)?;
                n.evaluate(ksi)?;
                assert_vec_approx_eq!(n.value.as_data(), solution, 1e-15);
            }
        }
        Ok(())
    }

    #[test]
    fn derive_works() {
        let mesh = Samples::two_quads_horizontal();
        let shapes = Shapes::new(&mesh).unwrap();
        let boundary = Boundary::new(&mesh, &shapes).unwrap();
        let n01 = NormalVector::at_edge(&mesh, &boundary, (0, 1)).unwrap();
        let n01_clone = n01.clone();
        assert_eq!(format!("{:?}", n01), "NormalVector { shape: Shape { class: Lin, kind: Lin2, space_ndim: 2, geo_ndim: 1, nnode: 2, nedge: 0, nface: 0, edge_nnode: 0, face_nnode: 0, face_nedge: 0, fn_interp: FnInterp, fn_deriv: FnDeriv }, state: StateOfShape { coords_transp: NumMatrix { nrow: 2, ncol: 2, data: [1.0, 0.0, 0.0, 0.0] }, coords_min: [0.0, 0.0], coords_max: [1.0, 0.0], interp: NumVector { data: [0.0, 0.0] }, deriv: NumMatrix { nrow: 2, ncol: 1, data: [0.0, 0.0] }, jacobian: NumMatrix { nrow: 2, ncol: 1, data: [0.0, 0.0] }, inv_jacobian: NumMatrix { nrow: 0, ncol: 0, data: [] }, gradient: NumMatrix { nrow: 0, ncol: 0, data: [] } }, value: NumVector { data: [0.0, 0.0] } }");
        assert_eq!(n01_clone.value.dim(), n01.value.dim());
    }
}
