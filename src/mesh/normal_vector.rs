use super::{EdgeKey, FaceKey, Features, Mesh};
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
    pub normal: Vector,
}

impl NormalVector {
    /// Allocates data to compute normal vector at edge
    ///
    /// **Important:** You must call [NormalVector::evaluate] to compute the actual values.
    ///
    /// **Note:** This function works in 2D only
    pub fn at_edge(mesh: &Mesh, features: &Features, edge_key: EdgeKey) -> Result<Self, StrError> {
        if mesh.space_ndim != 2 {
            return Err("normal at_edge works in 2D only");
        }
        const GEO_NDIM: usize = 1;
        let edge = match features.edges.get(&edge_key) {
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
            normal: Vector::new(mesh.space_ndim),
        })
    }

    /// Allocates data to compute normal vector at face
    ///
    /// **Important:** You must call [NormalVector::evaluate] to compute the actual values.
    ///
    /// **Note:** This function works in 3D only
    pub fn at_face(mesh: &Mesh, features: &Features, face_key: FaceKey) -> Result<Self, StrError> {
        if mesh.space_ndim != 3 {
            return Err("normal at_face works in 3D only");
        }
        const GEO_NDIM: usize = 2;
        let face = match features.faces.get(&face_key) {
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
            normal: Vector::new(mesh.space_ndim),
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
        self.shape.calc_boundary_normal(&mut self.normal, &mut self.state, ksi)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {}
