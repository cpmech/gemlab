use super::{Dof, Element, SystemDofs};
use crate::mesh::Mesh;
use crate::shapes::Shape;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::SparseTriplet;

pub struct ElementSolid {
    space_ndim: usize,
    shape: Shape,
    point_ids: Vec<usize>,
}

impl ElementSolid {
    pub fn new(mesh: &Mesh, cell_id: usize) -> Result<Self, StrError> {
        // shape
        let cell = &mesh.cells[cell_id];
        let space_ndim = mesh.space_ndim;
        let shape_ndim = cell.shape_ndim;
        let npoint = cell.points.len();
        let mut shape = Shape::new(space_ndim, shape_ndim, npoint)?;
        for m in 0..npoint {
            for i in 0..space_ndim {
                shape.set_coords(m, i, mesh.points[cell.points[m]].coords[i]);
            }
        }

        // done
        Ok(ElementSolid {
            space_ndim,
            shape,
            point_ids: cell.points.clone(),
        })
    }
}

impl Element for ElementSolid {
    fn assign_dofs(&self, dofs: &mut SystemDofs) {
        let dof_per_point = if self.space_ndim == 2 {
            vec![Dof::Ux, Dof::Uy]
        } else {
            vec![Dof::Ux, Dof::Uy, Dof::Uz]
        };
        let npoint = self.shape.get_npoint();
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

    fn add_ke_to_kk(&self, kk: &mut SparseTriplet) -> Result<(), StrError> {
        Ok(())
    }

    fn add_fe_to_ff(&self, ff: &mut Vector) -> Result<(), StrError> {
        Ok(())
    }
}
