use crate::{Edge, Face};

#[derive(Clone)]
pub struct Cell {
    pub id: usize,
    pub group: usize,
    pub vertices: Vec<usize>,
    pub edges: Vec<Edge>,
    pub faces: Vec<Face>,
}
