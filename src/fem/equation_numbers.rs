use crate::{mesh::PointId, StrError};
use std::cmp;
use std::fmt::{self, Write};

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
    numbers: Vec<Vec<i32>>,
}

impl EquationNumbers {
    /// Creates a new Equation Numbers object
    pub fn new(npoint: usize) -> Self {
        EquationNumbers {
            count: 0,
            numbers: vec![vec![-1; DOF_TOTAL]; npoint],
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
        // handle empty matrix
        if self.numbers.len() < 1 {
            return Ok(());
        }
        let nrow = self.numbers.len();
        let ncol = self.numbers[0].len();
        // find largest width
        let mut width = 0;
        let mut buf = String::new();
        for i in 0..nrow {
            for j in 0..ncol {
                let val = self.numbers[i][j];
                match f.precision() {
                    Some(v) => write!(&mut buf, "{:.1$}", val, v)?,
                    None => write!(&mut buf, "{}", val)?,
                }
                width = cmp::max(buf.chars().count(), width);
                buf.clear();
            }
        }
        // draw matrix
        width += 1;
        write!(f, "┌{:1$}┐\n", " ", width * ncol + 1)?;
        for i in 0..nrow {
            if i > 0 {
                write!(f, " │\n")?;
            }
            for j in 0..ncol {
                if j == 0 {
                    write!(f, "│")?;
                }
                let val = self.numbers[i][j];
                match f.precision() {
                    Some(v) => write!(f, "{:>1$.2$}", val, width, v)?,
                    None => write!(f, "{:>1$}", val, width)?,
                }
            }
        }
        write!(f, " │\n")?;
        write!(f, "└{:1$}┘", " ", width * ncol + 1)?;
        Ok(())
    }
}
