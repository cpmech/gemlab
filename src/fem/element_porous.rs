use super::{Element, EquationNumbers};
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::SparseTriplet;

pub struct ElementPorous {
    // todo
}

impl Element for ElementPorous {
    fn activate_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize {
        0
    }

    fn compute_local_rhs_vector(&mut self) -> Result<(), StrError> {
        Err("TODO")
    }

    fn compute_local_k_matrix(&mut self, first_iteration: bool) -> Result<(), StrError> {
        Err("TODO")
    }

    fn assemble_rhs_vector(&self, rhs: &mut Vector) -> Result<(), StrError> {
        Err("TODO")
    }

    fn assemble_k_matrix(&self, kk: &mut SparseTriplet) -> Result<(), StrError> {
        Err("TODO")
    }
}
