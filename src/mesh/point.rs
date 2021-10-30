#[derive(Clone, Debug)]
pub struct Point {
    pub id: usize,
    pub group: usize,
    pub coords: Vec<f64>,
    pub shared_by_cell_ids: Vec<usize>,
}
