#![allow(dead_code, unused_mut, unused_variables)]

use super::{Element, ElementConfig, ModelSolid, ModelSolidState};
use crate::fem::{DOF_UX, DOF_UY, DOF_UZ};
use crate::mesh::Mesh;
use crate::shapes::Shape;
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_sparse::SparseTriplet;
use russell_tensor::{Tensor2, Tensor4};

/// Implements the element for solid mechanics simulations
pub struct ElementSolid<'a> {
    /// Shape object to perform integrations
    shape: Shape,

    /// Map nodes ids to global point ids
    point_ids: Vec<usize>,

    /// Local DOF indices
    dof_indices: Vec<usize>,

    /// Local RHS-vector (neq)
    rhs: Vector,

    /// Local K-matrix (neq,neq)
    kk: Matrix,

    /// Auxiliary stress-strain modulus
    dd_aux: Tensor4,

    /// Auxiliary stress tensor
    sig_aux: Tensor2,

    /// Material model
    model: &'a Box<dyn ModelSolid>,

    /// Thickness
    thickness: f64,
    // Stress-strain-ivs state
    // state: Vec<ModelSolidState>,
}

/// Holds configuration data for a group of ElementSolid
pub struct ElementSolidConfig {
    /// Material model
    model: Box<dyn ModelSolid>,

    /// Thickness
    thickness: f64,
    // Use default integration points
    // int_points_default: bool,

    // /// Number of integration points
    // int_points_number: usize,

    // /// Flag for "edge" integration points
    // int_points_edge: bool,
}

impl ElementSolidConfig {
    /// Returns a new ElementSolidGroup instance
    pub fn new(model: Box<dyn ModelSolid>, thickness: f64) -> Self {
        ElementSolidConfig { model, thickness }
    }

    // Set integration points
    // pub fn set_int_points(&mut self, nip: usize, edge: bool, ws: bool) {
    // todo
    // }
}

impl ElementConfig for ElementSolidConfig {
    /// Tells whether the element group is active or not
    fn is_active(&self) -> bool {
        true
    }

    /// Allocates a new element belonging to this group
    fn allocate(&self, mesh: &Mesh, cell_id: usize) -> Result<Box<dyn Element + '_>, StrError> {
        // geometry
        let cell = &mesh.cells[cell_id];
        let geo_ndim = cell.geo_ndim;
        let space_ndim = mesh.space_ndim;
        let npoint = cell.points.len();
        let two_dim = space_ndim == 2;

        // check
        if geo_ndim < 2 {
            return Err("geometry ndim must be 2 or 3");
        }

        // new shape and set node coordinates
        let mut shape = Shape::new(space_ndim, geo_ndim, npoint)?;
        mesh.extract_coords(&mut shape, cell_id)?;

        // DOF indices and number of local equations
        let dof_indices = if space_ndim == 2 {
            vec![DOF_UX, DOF_UY]
        } else {
            vec![DOF_UX, DOF_UY, DOF_UZ]
        };
        let neq = shape.nnode * dof_indices.len();

        // new element
        let element = ElementSolid {
            shape,
            point_ids: cell.points.clone(),
            dof_indices,
            rhs: Vector::new(neq),
            kk: Matrix::new(neq, neq),
            dd_aux: Tensor4::new(true, two_dim),
            sig_aux: Tensor2::new(true, two_dim),
            model: &self.model,
            thickness: self.thickness,
            // state:vec![State]
        };
        Ok(Box::new(element))
    }
}

impl Element for ElementSolid<'_> {
    fn activate_equation_numbers(&self, equation_numbers: &mut super::EquationNumbers) -> usize {
        for point_id in &self.point_ids {
            for dof_index in &self.dof_indices {
                equation_numbers.activate_equation(*point_id, *dof_index);
            }
        }
        let neq = self.shape.nnode * self.dof_indices.len();
        neq * neq
    }

    fn compute_local_rhs_vector(&mut self) -> Result<(), StrError> {
        let fn_sig = |sig: &mut Tensor2, _: usize| {
            // todo
        };
        self.shape
            .integ_vec_d_tg(&mut self.rhs, fn_sig, &mut self.sig_aux, self.thickness)
    }

    fn compute_local_k_matrix(&mut self, first_iteration: bool) -> Result<(), StrError> {
        let two_dim = self.shape.space_ndim == 2;
        let state_new = ModelSolidState {
            stress: Tensor2::new(true, two_dim),
            strain: Tensor2::new(true, two_dim),
            ivs: Vector::new(0),
        };
        let fn_dd = |dd: &mut Tensor4, _: usize| {
            self.model.consistent_modulus(dd, &state_new, first_iteration).unwrap();
        };
        self.shape
            .integ_mat_10_gdg(&mut self.kk, fn_dd, &mut self.dd_aux, self.thickness)
    }

    fn assemble_rhs_vector(&self, rhs: &mut Vector) -> Result<(), StrError> {
        Ok(())
    }

    fn assemble_k_matrix(&self, kk: &mut SparseTriplet) -> Result<(), StrError> {
        Ok(())
    }
}
