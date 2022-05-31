#![allow(unused_imports)]
#![allow(unused_mut)]
#![allow(unused_variables)]

use super::{allocate_shapes, allocate_state, allocate_states, EdgeKey, Mesh, NormalVector, Region};
use crate::shapes::{Shape, StateOfShape};
use crate::StrError;
use plotpy::{Curve, Plot, Shapes, Text};
use std::collections::HashMap;

pub struct Plotter {}

impl Plotter {
    pub fn with(region: &Region) -> Result<Plot, StrError> {
        /*
        let shapes = allocate_shapes(mesh)?;
        let mut states = allocate_states(mesh)?;
        let mut canvas = Shapes::new();
        for cell_id in 0..mesh.cells.len() {
            let cell = &mesh.cells[cell_id];
            let shape = &shapes[cell_id];
            let state = &states[cell_id];
            if shape.geo_ndim == 2 {
                for e in 0..shape.nedge {}
            }
            let points = canvas.draw_polyline(&[[0.0, 0.0]], true);
        }
        */

        let space_ndim = region.mesh.space_ndim;
        let mut normals: HashMap<EdgeKey, NormalVector> = HashMap::new();
        for (edge_key, edge) in &region.features.edges {
            let normal =
                normals
                    .entry(*edge_key)
                    .or_insert(NormalVector::at_edge(&region.mesh, &region.features, *edge_key)?);
        }

        let mut plot = Plot::new();
        Ok(plot)
    }
}
