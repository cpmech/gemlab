use crate::{mesh::PointId, StrError};
use russell_lab::GenericMatrix;
use std::fmt;

/// Alias for DOF index
pub type DofIndex = usize;

/// DOF index: Displacement along the first dimension
pub const DOF_UX: DofIndex = 0;

/// DOF index: Displacement along the second dimension
pub const DOF_UY: DofIndex = 1;

/// DOF index: Displacement along the third dimension
pub const DOF_UZ: DofIndex = 2;

/// DOF index: Rotation around the first axis
pub const DOF_RX: DofIndex = 3;

/// DOF index: Rotation around the second axis
pub const DOF_RY: DofIndex = 4;

/// DOF index: Rotation around the third axis
pub const DOF_RZ: DofIndex = 5;

/// DOF index: Temperature
pub const DOF_T: DofIndex = 6;

/// DOF index: Liquid pressure
pub const DOF_PL: DofIndex = 7;

/// DOF index: Gas pressure
pub const DOF_PG: DofIndex = 8;

/// DOF index: Free-surface-output (fso) enrichment
pub const DOF_FSO: DofIndex = 9;

/// Total number of available DOFs
pub const DOF_TOTAL: usize = 10;

/// Holds equation numbers (DOF numbers)
pub struct EquationNumbers {
    /// Total number of equations
    count: i32,

    /// Equation numbers matrix [point][dof]
    numbers: GenericMatrix<i32>,
}

impl EquationNumbers {
    /// Creates a new Equation Numbers object
    pub fn new(npoint: usize) -> Self {
        EquationNumbers {
            count: 0,
            numbers: GenericMatrix::filled(npoint, DOF_TOTAL, -1),
        }
    }

    /// Activates equation corresponding to a point-DOF pair
    ///
    /// Note: Also increments the number of equations count
    ///       if the equation does not exist yet
    pub fn activate_equation(&mut self, point_id: PointId, dof_index: DofIndex) {
        if self.numbers[point_id][dof_index] < 0 {
            self.numbers[point_id][dof_index] = self.count;
            self.count += 1;
        }
    }

    /// Returns the current total number of equations (DOFs)
    pub fn get_number_of_equations(&self) -> usize {
        self.count as usize
    }

    /// Returns the equation number corresponding to a point-DOF pair
    pub fn get_equation_number(&self, point_id: PointId, dof_index: DofIndex) -> Result<usize, StrError> {
        if self.numbers[point_id][dof_index] < 0 {
            return Err("equation number has not been set");
        }
        Ok(self.numbers[point_id][dof_index] as usize)
    }
}

impl fmt::Display for EquationNumbers {
    /// Generates a string representation of the EquationNumbers
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.numbers)
    }
}
