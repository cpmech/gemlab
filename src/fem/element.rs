use super::EquationNumbers;
use crate::mesh::{CellAttributeId, Mesh};
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::SparseTriplet;
use std::collections::HashMap;

/// Defines a trait for (finite) elements
pub trait Element {
    /// Activates an equation number, if not set yet
    fn activate_equation_numbers(&self, equation_numbers: &mut EquationNumbers) -> usize;

    /// Computes the element RHS-vector
    fn compute_local_rhs_vector(&mut self) -> Result<(), StrError>;

    /// Computes the element K-matrix
    fn compute_local_k_matrix(&mut self, first_iteration: bool) -> Result<(), StrError>;

    /// Assembles local right-hand side (RHS) vector into global RHS-vector
    fn assemble_rhs_vector(&self, rhs: &mut Vector) -> Result<(), StrError>;

    /// Assembles the local K-matrix into the global K-matrix
    fn assemble_k_matrix(&self, kk: &mut SparseTriplet) -> Result<(), StrError>;
}

/// Defines a trait for element configuration
pub trait ElementConfig {
    /// Tells whether the element is active or not
    fn is_active(&self) -> bool;

    /// Allocates a new element
    fn allocate(&self, mesh: &Mesh, cell_id: usize) -> Result<Box<dyn Element + '_>, StrError>;
}

/// Maps element attributes to configuration
pub struct ElementAttributes {
    attributes: HashMap<CellAttributeId, Box<dyn ElementConfig>>,
}

impl ElementAttributes {
    /// Returns a new ElementAttributes instance
    pub fn new() -> Self {
        ElementAttributes {
            attributes: HashMap::new(),
        }
    }

    /// Sets a new element attribute
    pub fn set(&mut self, attribute_id: CellAttributeId, attribute: Box<dyn ElementConfig>) -> Result<(), StrError> {
        if self.attributes.contains_key(&attribute_id) {
            return Err("attribute already set");
        }
        self.attributes.insert(attribute_id, attribute);
        Ok(())
    }

    /// Returns an existent element attribute
    pub fn get(&self, attribute_id: CellAttributeId) -> Result<&Box<dyn ElementConfig>, StrError> {
        match self.attributes.get(&attribute_id) {
            Some(attribute) => Ok(attribute),
            None => Err("cannot find an element attribute"),
        }
    }
}
