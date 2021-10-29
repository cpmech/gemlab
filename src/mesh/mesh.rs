use crate::{Cell, Vertex};

pub struct Mesh {
    pub ndim: usize,
    pub vertices: Vec<Vertex>,
    pub cells: Vec<Cell>,
}
