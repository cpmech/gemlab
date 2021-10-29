#[derive(Clone)]
pub struct Vertex {
    pub id: usize,
    pub group: usize,
    pub coords: Vec<f64>,
    pub shared_by_cells: Vec<usize>,
}

/*
impl Clone for Vertex {
    fn clone(&self) -> Self {
        Vertex {
            id: self.id,
            group: self.group,
            coords: self.coords.clone(),
            shared_by_cells: self.shared_by_cells.clone(),
        }
    }
}
*/
