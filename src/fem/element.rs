use super::{ElementKind, ElementSolid, PointDofs};
use crate::mesh::Mesh;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::SparseTriplet;

pub trait Element {
    fn get_dofs(&self) -> Result<Vec<PointDofs>, StrError>;
    fn compute_ke(&mut self) -> Result<(), StrError>;
    fn add_ke_to_kk(&self, kk: &mut SparseTriplet) -> Result<(), StrError>;
    fn add_fe_to_ff(&self, ff: &mut Vector) -> Result<(), StrError>;
}

pub fn new_element(kind: ElementKind, mesh: &Mesh, cell_id: usize) -> Result<Box<dyn Element>, StrError> {
    match kind {
        ElementKind::Solid => Ok(Box::new(ElementSolid::new(mesh, cell_id)?)),
        _ => panic!("Element kind {:?} is not available yet", kind),
    }
}
