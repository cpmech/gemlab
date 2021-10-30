use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Cell {
    pub id: usize,
    pub group: usize,
    pub point_ids: Vec<usize>,
    pub edge_ids: Vec<usize>,
    pub face_ids: Vec<usize>,
}
