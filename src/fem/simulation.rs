use super::definitions::{Attribute, Bc, Dof, EdgeBc, Element, FnTimeSpace, PointBc};
use super::possible_dofs;
use crate::mesh::Mesh;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::{SparseTriplet, Symmetry};
use std::collections::HashMap;

pub struct Simulation {
    mesh: Mesh,
    point_bcs: Vec<PointBc>,
    edge_bcs: Vec<EdgeBc>,
    elements: Vec<Element>,
    attributes: HashMap<usize, Attribute>,
}

impl Simulation {
    pub fn new(mesh: Mesh) -> Result<Self, StrError> {
        Ok(Simulation {
            mesh,
            point_bcs: Vec::new(),
            edge_bcs: Vec::new(),
            elements: Vec::new(),
            attributes: HashMap::new(),
        })
    }

    pub fn add_attribute(&mut self, attribute_id: usize, a: Attribute) -> &mut Self {
        self
    }

    pub fn add_point_bc(&mut self, group: &str, bc: Bc, dof: Dof, f: FnTimeSpace) -> &mut Self {
        let ids = self.mesh.get_boundary_point_ids_sorted(group);
        for id in ids {
            self.point_bcs.push(PointBc { bc, dof, f, id });
        }
        self
    }

    pub fn add_edge_bc(&mut self, group: &str, bc: Bc, dof: Dof, f: FnTimeSpace) -> &mut Self {
        let keys = self.mesh.get_boundary_edge_keys_sorted(group);
        for key in keys {
            self.edge_bcs.push(EdgeBc { bc, dof, f, key });
        }
        self
    }

    pub fn initialize(&mut self) -> Result<(), StrError> {
        // allocate all elements and assign numbers to the DOFs
        let ndof = 0_usize;
        for cell in &self.mesh.cells {
            let props = match self.attributes.get(&cell.attribute_id) {
                Some(a) => {
                    if a.inactive {
                        continue;
                    }
                    let dofs = possible_dofs(self.mesh.ndim, &a.element_type);
                    &a.properties
                }
                None => return Err("cannot find cell with a specific attribute id"),
            };
            // how to get tne DOFs for each element/problem?
        }
        Ok(())
    }

    pub fn run(&mut self) -> Result<(), StrError> {
        let time = 0.0; // time

        let n = 2;
        let nnz = 2;
        let mut uu = Vector::new(n);
        let mut ff = Vector::new(n);
        let mut kk = SparseTriplet::new(n, n, nnz, Symmetry::No)?;

        // assemble element Ke matrices (can be done in parallel)
        for element in &mut self.elements {
            element.compute_ke()?;
        }

        // add element Ke matrix to global K matrix (must be serial)
        for element in &self.elements {
            element.add_ke_to_kk(&mut kk)?;
        }

        // add element Fe vector to global F vector (must be serial)
        for element in &self.elements {
            element.add_fe_to_ff(&mut ff)?;
        }

        let i = 0;

        // handle point boundary conditions
        for item in &self.point_bcs {
            let x = &self.mesh.points[item.id].coords;
            let val = (item.f)(time, x);
            match item.bc {
                Bc::Essential => {
                    kk.put(i, i, 1.0);
                    uu[i] = val;
                    ff[i] = val;
                }
                Bc::Natural => {
                    ff[i] = val;
                }
            }
        }

        // handle edge boundary conditions
        for item in &self.edge_bcs {
            let edge = match self.mesh.boundary_edges.get(&item.key) {
                Some(e) => e,
                None => return Err("cannot find boundary edge"),
            };
            match item.bc {
                Bc::Essential => {
                    for id in &edge.points {
                        let x = &self.mesh.points[*id].coords;
                        let val = (item.f)(time, x);
                        kk.put(i, i, 1.0);
                        uu[i] = val;
                        ff[i] = val;
                    }
                }
                Bc::Natural => {
                    // perform integration
                }
            }
        }

        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::At;

    #[test]
    fn run_works() -> Result<(), StrError> {
        // create mesh
        let mut mesh = Mesh::from_text(
            r"
            # header
            # ndim npoint ncell
            2 6 2

            # points
            # id x y
            0  0.0 0.0
            1  1.0 0.0
            2  1.0 1.0
            3  0.0 1.0
            4  2.0 0.0
            5  2.0 1.0

            # cells
            # id attribute ndim npoint point_ids...
            0 1  2 4  0 1 2 3
            1 0  2 4  1 4 5 2
            ",
        )?;

        // set groups of points and edges at boundaries
        mesh.set_points("origin", At::XY(0.0, 0.0))?
            .set_edges("left", At::X(0.0))?
            .set_edges("right", At::X(1.0))?
            .set_edges("bottom", At::Y(0.0))?
            .set_edges("top", At::Y(1.0))?;

        // create simulation
        let mut sim = Simulation::new(mesh)?;

        // set boundary conditions
        sim.add_point_bc("origin", Bc::Essential, Dof::Ux, |_, _| 0.0)
            .add_point_bc("origin", Bc::Essential, Dof::Uy, |_, _| 0.0)
            .add_edge_bc("left", Bc::Essential, Dof::Ux, |_, _| 0.0)
            .add_edge_bc("right", Bc::Essential, Dof::Ux, |_, _| 0.0)
            .add_edge_bc("bottom", Bc::Essential, Dof::Uy, |_, _| 0.0);

        // run simulation
        sim.run();

        Ok(())
    }
}
