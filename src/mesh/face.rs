#[derive(Clone, Debug)]
pub struct Face {
    pub id: usize,
    pub group: usize,
    pub point_ids: Vec<usize>,
    pub shared_by_cell_ids: Vec<usize>,
}
