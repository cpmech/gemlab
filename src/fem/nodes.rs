use super::{EquationNumber, PointDofKey};
use std::collections::HashMap;

pub struct Nodes {
    dof_map: HashMap<PointDofKey, EquationNumber>,
}
