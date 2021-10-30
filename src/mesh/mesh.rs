use crate::{Cell, Point};

pub struct Mesh {
    pub ndim: usize,
    pub vertices: Vec<Point>,
    pub cells: Vec<Cell>,
}
