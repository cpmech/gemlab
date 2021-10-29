#[derive(Clone)]
pub struct Edge {
    pub id: usize,
    pub group: usize,
    pub vertices: Vec<usize>,
    pub shared_by_cells: Vec<usize>,
}
