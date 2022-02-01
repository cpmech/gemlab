#![allow(dead_code, unused_mut, unused_variables)]

use super::{ModelSolid, ModelSolidState};
use crate::StrError;
use russell_lab::{add_vectors, copy_matrix, copy_vector, Vector};
use russell_tensor::{t4_ddot_t2, LinElasticity, Tensor2, Tensor4};

pub struct ModelSolidLinearElastic {
    two_dim: bool,
    lin_elasticity: LinElasticity,
}

impl ModelSolidLinearElastic {
    pub fn new(young: f64, poisson: f64, two_dim: bool, plane_stress: bool) -> Self {
        ModelSolidLinearElastic {
            two_dim,
            lin_elasticity: LinElasticity::new(young, poisson, two_dim, plane_stress),
        }
    }
}

impl ModelSolid for ModelSolidLinearElastic {
    /// Calculates internal values and initializes state
    fn initialize_state(&self, stress_ini: &Tensor2) -> Result<ModelSolidState, StrError> {
        let mut stress = Tensor2::new(true, self.two_dim);
        copy_vector(&mut stress.vec, &stress_ini.vec)?;
        Ok(ModelSolidState {
            stress,
            strain: Tensor2::new(true, self.two_dim),
            ivs: Vector::new(0),
        })
    }

    /// Updates state for given delta_strain
    fn update_state(
        &self,
        state_new: &mut ModelSolidState,
        state: &ModelSolidState,
        delta_strain: &Tensor2,
    ) -> Result<(), StrError> {
        // update strain
        let eps_new = &mut state_new.strain.vec;
        let eps = &state.strain.vec;
        let deps = &delta_strain.vec;
        add_vectors(eps_new, 1.0, eps, 1.0, deps)?;

        // update stress
        let dd = self.lin_elasticity.get_modulus();
        t4_ddot_t2(&mut state_new.stress, 1.0, &dd, &state_new.strain)
    }

    /// Calculates the consistent tangent modulus at new state
    fn consistent_modulus(
        &self,
        dd_new: &mut Tensor4,
        state_new: &ModelSolidState,
        first_iteration: bool,
    ) -> Result<(), StrError> {
        let dd_elastic = self.lin_elasticity.get_modulus();
        copy_matrix(&mut dd_new.mat, &dd_elastic.mat)
    }
}
