use super::{ElementSolid, EquationNumbers};
use crate::mesh::Mesh;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::SparseTriplet;

#[derive(Clone, Copy, Debug)]
pub enum ElementKind {
    Diffusion,
    Truss,
    Beam,
    Solid,
    SeepageLiq,
    SeepageLiqGas,
    PorousLiq,
    PorousLiqGas,
    PorousLiqGasTemp,
    PorousLiqVel,
    PorousLiqGasVel,
}

pub trait Element {
    fn activate_equation_numbers(&self, equation_numbers: &mut EquationNumbers);

    fn get_nnz(&self) -> usize;

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
