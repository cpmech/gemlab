use crate::StrError;
use russell_lab::Vector;
use russell_tensor::{Tensor2, Tensor4};

pub struct ModelSolidState {
    pub stress: Tensor2, // stress
    pub strain: Tensor2, // strain
    pub ivs: Vector,     // internal values
}

pub trait ModelSolid {
    /// Calculates internal values and initializes state
    fn initialize_state(&self, stress_ini: &Tensor2) -> Result<ModelSolidState, StrError>;

    /// Updates state for given delta_strain
    fn update_state(
        &self,
        state_new: &mut ModelSolidState,
        state: &ModelSolidState,
        delta_strain: &Tensor2,
    ) -> Result<(), StrError>;

    /// Calculates the consistent tangent modulus at new state
    fn consistent_modulus(
        &self,
        dd_new: &mut Tensor4,
        state_new: &ModelSolidState,
        first_iteration: bool,
    ) -> Result<(), StrError>;
}
