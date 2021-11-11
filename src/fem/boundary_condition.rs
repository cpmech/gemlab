use super::degree_of_freedom::Dof;
use crate::mesh::{EdgeKey, Index};

pub type FnTimeSpace = fn(f64, &[f64]) -> f64;

#[derive(Clone, Copy)]
pub enum Bc {
    Essential,
    Natural,
}

#[derive(Clone, Copy)]
pub struct PointBc {
    pub bc: Bc,
    pub dof: Dof,
    pub f: FnTimeSpace,
    pub point_id: Index,
}

#[derive(Clone, Copy)]
pub struct EdgeBc {
    pub bc: Bc,
    pub dof: Dof,
    pub f: FnTimeSpace,
    pub edge_key: EdgeKey,
}
