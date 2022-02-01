#![allow(dead_code, unused_mut, unused_variables)]

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

pub struct BoundaryCondition {
    pub point_id: Option<PointId>,
    pub edge_key: Option<EdgeKey>,
    pub dof_index: DofIndex,
    pub condition: Bc,
    pub f: FnTimeSpace,
}

pub struct PointBoundaryConditions {
    // TODO
}

pub struct EdgeBoundaryConditions {
    // TODO
}

pub struct BoundaryConditions {
    conditions: Vec<BoundaryCondition>,
}

/*
pub fn add_point_bc(&mut self, group: &str, condition: Bc, dof_index: DofIndex, f: FnTimeSpace) -> &mut Self {
    let ids = self.mesh.get_boundary_point_ids_sorted(group);
    for point_id in ids {
        self.point_bcs.push(PointBc {
            point_id,
            dof_index,
            condition,
            f,
        });
    }
    self
}

pub fn add_edge_bc(&mut self, group: &str, condition: Bc, dof_index: DofIndex, f: FnTimeSpace) -> &mut Self {
    let keys = self.mesh.get_boundary_edge_keys_sorted(group);
    for edge_key in keys {
        self.edge_bcs.push(EdgeBc {
            edge_key,
            dof_index,
            condition,
            f,
        });
    }
    self
}
*/
