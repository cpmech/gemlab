use crate::mesh::{EdgeKey, Index};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug)]
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
pub type PointDofKey = (usize, Dof);

/// Holds the number of a DOF in the global system of equations
pub type EquationNumber = usize;

pub type PointDofs = (Index, Vec<Dof>);

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

pub type FnTimeSpace = fn(f64, &[f64]) -> f64;

#[derive(Clone, Copy)]
pub enum Bc {
    Essential,
    Natural,
}

#[derive(Clone, Copy)]
pub(super) struct PointBc {
    pub(super) bc: Bc,
    pub(super) dof: Dof,
    pub(super) f: FnTimeSpace,
    pub(super) id: Index,
}

#[derive(Clone, Copy)]
pub(super) struct EdgeBc {
    pub(super) bc: Bc,
    pub(super) dof: Dof,
    pub(super) f: FnTimeSpace,
    pub(super) key: EdgeKey,
}

#[derive(Clone, Debug)]
pub struct Attribute {
    pub(super) kind: ElementKind,
    pub(super) inactive: bool,
    pub(super) parameters: HashMap<String, f64>,
    pub(super) flags: HashMap<String, bool>,
}

impl Attribute {
    pub fn new(kind: ElementKind) -> Self {
        Attribute {
            kind,
            inactive: false,
            parameters: HashMap::new(),
            flags: HashMap::new(),
        }
    }

    pub fn set_parameter(&mut self, name: &str, value: f64) -> &mut Self {
        if let Some(parameter) = self.parameters.get_mut(name) {
            *parameter = value;
        }
        self
    }

    pub fn set_flag(&mut self, name: &str, value: bool) -> &mut Self {
        if let Some(flag) = self.flags.get_mut(name) {
            *flag = value;
        }
        self
    }
}

pub(super) fn possible_dofs(ndim: usize, kind: &ElementKind) -> Vec<Dof> {
    match kind {
        ElementKind::Diffusion => {
            vec![Dof::T]
        }
        ElementKind::Truss => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            }
        }
        ElementKind::Beam => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy, Dof::Rz]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz]
            }
        }
        ElementKind::Solid => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            }
        }
        ElementKind::SeepageLiq => {
            vec![Dof::Pl, Dof::Extra]
        }
        ElementKind::SeepageLiqGas => {
            vec![Dof::Pl, Dof::Pg, Dof::Extra]
        }
        ElementKind::PorousLiq => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy, Dof::Pl, Dof::Extra]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Pl, Dof::Extra]
            }
        }
        ElementKind::PorousLiqGas => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy, Dof::Pl, Dof::Pg, Dof::Extra]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Pl, Dof::Pg, Dof::Extra]
            }
        }
        ElementKind::PorousLiqGasTemp => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy, Dof::Pl, Dof::Pg, Dof::Extra, Dof::T]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Pl, Dof::Pg, Dof::Extra, Dof::T]
            }
        }
        ElementKind::PorousLiqVel => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy, Dof::Pl, Dof::Extra, Dof::Wlx, Dof::Wly]
            } else {
                vec![
                    Dof::Ux,
                    Dof::Uy,
                    Dof::Uz,
                    Dof::Pl,
                    Dof::Extra,
                    Dof::Wlx,
                    Dof::Wly,
                    Dof::Wlz,
                ]
            }
        }
        ElementKind::PorousLiqGasVel => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy, Dof::Pl, Dof::Pg, Dof::Extra, Dof::Wlx, Dof::Wly]
            } else {
                vec![
                    Dof::Ux,
                    Dof::Uy,
                    Dof::Uz,
                    Dof::Pl,
                    Dof::Pg,
                    Dof::Extra,
                    Dof::Wlx,
                    Dof::Wly,
                    Dof::Wlz,
                ]
            }
        }
    }
}
