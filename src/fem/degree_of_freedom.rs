use std::collections::HashMap;

/// Defines the type of degree-of-freedom
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Dof {
    T,     // temperature or any other scalar
    Pl,    // liquid pressure
    Pg,    // gas pressure
    Ux,    // displacement along x
    Uy,    // displacement along y
    Uz,    // displacement along z
    Rx,    // rotation around x
    Ry,    // rotation around y
    Rz,    // rotation around z
    Wlx,   // relative liquid velocity along x
    Wly,   // relative liquid velocity along y
    Wlz,   // relative liquid velocity along z
    Extra, // extra degree of freedom (e.g., enriched)
}

/// Defines a key using the point index (ID) and the DOF type
type PointIdDofKey = (usize, Dof);

/// Corresponds to the position of a DOF in the global system of equations
type EquationNumber = usize;

pub struct SystemDofs {
    /// Maps (point_id, dof) to equation number
    map: HashMap<PointIdDofKey, EquationNumber>,

    /// Maps equation number to (point_id, dof)
    ///
    /// The length of this map equals the total number of equations.
    reverse_map: Vec<PointIdDofKey>,
}

impl SystemDofs {
    pub fn new() -> Self {
        SystemDofs {
            map: HashMap::new(),
            reverse_map: Vec::new(),
        }
    }

    pub fn update(&mut self, point_id: usize, dof: Dof) -> &mut Self {
        let key: PointIdDofKey = (point_id, dof);
        if !self.map.contains_key(&key) {
            let nequation = self.reverse_map.len();
            self.map.insert(key, nequation);
            self.reverse_map.push(key);
        }
        self
    }
}
