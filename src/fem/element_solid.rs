use super::{Dof, Element, PointDofs};
use crate::mesh::Mesh;
use crate::shapes::{new_shape, Shape};
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::SparseTriplet;

pub struct Node {
    id: usize,
    dofs: Vec<Dof>,
    coords: Vec<f64>,
}

pub struct ElementSolid {
    ndim: usize,
    shape: Box<dyn Shape>,
    nodes: Vec<Node>,
}

impl ElementSolid {
    pub fn new(mesh: &Mesh, cell_id: usize) -> Result<Self, StrError> {
        // shape
        let cell = &mesh.cells[cell_id];
        let ndim = mesh.ndim;
        let npoint = cell.points.len();
        let shape = new_shape(ndim, npoint)?;

        // node DOFs
        let node_dofs = if ndim == 2 {
            vec![Dof::Ux, Dof::Uy]
        } else {
            vec![Dof::Ux, Dof::Uy, Dof::Uz]
        };

        // nodes
        let nodes: Vec<_> = cell
            .points
            .iter()
            .map(|id| Node {
                id: *id,
                dofs: node_dofs.clone(),
                coords: mesh.points[*id].coords.clone(),
            })
            .collect();

        // done
        Ok(ElementSolid { ndim, shape, nodes })
    }
}

impl Element for ElementSolid {
    fn get_dofs(&self) -> Result<Vec<PointDofs>, StrError> {
        Ok(Vec::new())
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
