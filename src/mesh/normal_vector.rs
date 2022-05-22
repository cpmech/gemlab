use super::{Boundary, EdgeKey, FaceKey, Mesh};
use crate::shapes::{Shape, StateOfShape};
use crate::StrError;
use russell_lab::Vector;

/// Holds data to evaluate normal vectors at edge or face
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
    pub fn at_edge(mesh: &Mesh, boundary: &Boundary, edge_key: EdgeKey) -> Result<Self, StrError> {
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
    pub fn at_face(mesh: &Mesh, boundary: &Boundary, face_key: FaceKey) -> Result<Self, StrError> {
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
