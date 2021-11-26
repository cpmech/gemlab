use super::{Dof, Element, SystemDofs};
use crate::mesh::Mesh;
use crate::shapes::Shape;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::SparseTriplet;

pub struct ElementSolid {
    shape: Shape,
    point_ids: Vec<usize>,
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
        })
    }
}

impl Element for ElementSolid {
    fn assign_dofs(&self, dofs: &mut SystemDofs) {
        let dof_per_point = if self.shape.space_ndim == 2 {
            vec![Dof::Ux, Dof::Uy]
        } else {
            vec![Dof::Ux, Dof::Uy, Dof::Uz]
        };
        let npoint = self.shape.npoint;
        for m in 0..npoint {
            for dof in &dof_per_point {
                dofs.update(self.point_ids[m], *dof);
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
