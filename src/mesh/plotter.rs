#![allow(unused_imports)]
#![allow(unused_mut)]
#![allow(unused_variables)]

use super::{allocate_shapes, allocate_states, Mesh};
use crate::StrError;
use plotpy::{Curve, Plot, Shapes, Text};

pub struct Plotter {}

impl Plotter {
    pub fn with(mesh: &Mesh) -> Result<Plot, StrError> {
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
        let mut plot = Plot::new();
        Ok(plot)
    }
}
