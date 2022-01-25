use super::Element;
use crate::fem::{DOF_UX, DOF_UY, DOF_UZ};
use crate::mesh::Mesh;
use crate::shapes::Shape;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::SparseTriplet;

#[allow(dead_code)]

pub struct ElementSolid {
    shape: Shape,
    point_ids: Vec<usize>,
    dof_indices: Vec<usize>,
}

impl ElementSolid {
    pub fn new(mesh: &Mesh, cell_id: usize) -> Result<Self, StrError> {
        // shape
        let cell = &mesh.cells[cell_id];
        let space_ndim = mesh.space_ndim;
        let geo_ndim = cell.geo_ndim;
        let npoint = cell.points.len();
        let mut shape = Shape::new(space_ndim, geo_ndim, npoint)?;
        mesh.extract_coords(&mut shape, cell_id)?;
        // element
        Ok(ElementSolid {
            shape,
            point_ids: cell.points.clone(),
            dof_indices: if space_ndim == 2 {
                vec![DOF_UX, DOF_UY]
            } else {
                vec![DOF_UX, DOF_UY, DOF_UZ]
            },
        })
    }
}

impl Element for ElementSolid {
    fn activate_equation_numbers(&self, equation_numbers: &mut super::EquationNumbers) {
        for point_id in &self.point_ids {
            for dof_index in &self.dof_indices {
                equation_numbers.activate_equation(*point_id, *dof_index);
            }
        }
    }

    fn get_nnz(&self) -> usize {
        0
    }

    fn compute_ke(&mut self) -> Result<(), StrError> {
        Ok(())
    }

    fn add_ke_to_kk(&self, _: &mut SparseTriplet) -> Result<(), StrError> {
        Ok(())
    }

    fn add_fe_to_ff(&self, _: &mut Vector) -> Result<(), StrError> {
        Ok(())
    }
}
