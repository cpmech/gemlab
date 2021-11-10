use crate::mesh::{EdgeKey, Index};
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::SparseTriplet;
use std::collections::HashMap;

#[derive(Clone, Copy)]
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

pub enum ElementType {
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

pub struct Node {
    // dofs: Vec<(DOF, usize)>, // (dof, equation_number)
}

pub struct Element {
    //
}

impl Element {
    pub fn compute_ke(&mut self) -> Result<(), StrError> {
        Ok(())
    }
    pub fn add_ke_to_kk(&self, kk: &mut SparseTriplet) -> Result<(), StrError> {
        Ok(())
    }
    pub fn add_fe_to_ff(&self, ff: &mut Vector) -> Result<(), StrError> {
        Ok(())
    }
}

pub struct Attribute {
    pub inactive: bool,
    pub element_type: ElementType,
    pub properties: HashMap<String, f64>,
}

pub(super) fn possible_dofs(ndim: usize, element_type: &ElementType) -> Vec<Dof> {
    match element_type {
        ElementType::Diffusion => {
            vec![Dof::T]
        }
        ElementType::Truss => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            }
        }
        ElementType::Beam => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy, Dof::Rz]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz]
            }
        }
        ElementType::Solid => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz]
            }
        }
        ElementType::SeepageLiq => {
            vec![Dof::Pl, Dof::Extra]
        }
        ElementType::SeepageLiqGas => {
            vec![Dof::Pl, Dof::Pg, Dof::Extra]
        }
        ElementType::PorousLiq => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy, Dof::Pl, Dof::Extra]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Pl, Dof::Extra]
            }
        }
        ElementType::PorousLiqGas => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy, Dof::Pl, Dof::Pg, Dof::Extra]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Pl, Dof::Pg, Dof::Extra]
            }
        }
        ElementType::PorousLiqGasTemp => {
            if ndim == 2 {
                vec![Dof::Ux, Dof::Uy, Dof::Pl, Dof::Pg, Dof::Extra, Dof::T]
            } else {
                vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Pl, Dof::Pg, Dof::Extra, Dof::T]
            }
        }
        ElementType::PorousLiqVel => {
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
        ElementType::PorousLiqGasVel => {
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
