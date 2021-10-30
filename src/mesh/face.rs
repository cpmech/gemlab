use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Face {
    pub id: usize,
    pub group: usize,
    pub point_ids: Vec<usize>,
    pub shared_by_cell_ids: Vec<usize>,
}
