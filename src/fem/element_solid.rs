use super::{Dof, Element, Nodes, PointDofs};
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
    space_ndim: usize,
    shape: Box<dyn Shape>,
    // nodes: Vec<Node>,
}

impl ElementSolid {
    pub fn new(mesh: &Mesh, cell_id: usize) -> Result<Self, StrError> {
        // shape
        let cell = &mesh.cells[cell_id];
        let space_ndim = mesh.space_ndim;
        let shape_ndim = cell.shape_ndim;
        let npoint = cell.points.len();
        let mut shape = new_shape(space_ndim, shape_ndim, npoint)?;
        for m in 0..npoint {
            for i in 0..space_ndim {
                shape.set_coords(m, i, mesh.points[cell.points[m]].coords[i]);
            }
        }

        // node DOFs
        let node_dofs = if space_ndim == 2 {
            vec![Dof::Ux, Dof::Uy]
        } else {
            vec![Dof::Ux, Dof::Uy, Dof::Uz]
        };

        /*
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
            */

        // let coords:Vec<_>=cell.points.iter().map(|id| mesh.points[*id].coords)

        // done
        Ok(ElementSolid { space_ndim, shape })
    }

    fn assign_dofs(&self, nodes: &mut Nodes) {}
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
