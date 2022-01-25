use super::DofIndex;
use crate::mesh::{EdgeKey, PointId};

pub type FnTimeSpace = fn(f64, &[f64]) -> f64;

#[derive(Clone, Copy)]
pub enum Bc {
    Essential,
    Natural,
}

#[derive(Clone, Copy)]
pub struct PointBc {
    pub point_id: PointId,
    pub dof_index: DofIndex,
    pub condition: Bc,
    pub f: FnTimeSpace,
}

#[derive(Clone, Copy)]
pub struct EdgeBc {
    pub edge_key: EdgeKey,
    pub dof_index: DofIndex,
    pub condition: Bc,
    pub f: FnTimeSpace,
}
