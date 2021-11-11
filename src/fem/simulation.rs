use super::{new_element, Attribute, Bc, Dof, EdgeBc, Element, FnTimeSpace, PointBc, SystemDofs};
use crate::mesh::Mesh;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry};
use std::collections::HashMap;
use std::fmt;

pub struct Simulation {
    mesh: Mesh,
    point_bcs: Vec<PointBc>,
    edge_bcs: Vec<EdgeBc>,
    attributes: HashMap<usize, Attribute>,
    elements: Vec<Box<dyn Element>>,
    system_dofs: SystemDofs,
    system_kk: SparseTriplet,
    system_uu: Vector,
    system_ff: Vector,
}

impl Simulation {
    pub fn new(mesh: Mesh) -> Result<Self, StrError> {
        Ok(Simulation {
            mesh,
            point_bcs: Vec::new(),
            edge_bcs: Vec::new(),
            attributes: HashMap::new(),
            elements: Vec::new(),
            system_dofs: SystemDofs::new(),
            system_kk: SparseTriplet::new(0, 0, 0, Symmetry::No)?,
            system_uu: Vector::new(0),
            system_ff: Vector::new(0),
        })
    }

    pub fn add_point_bc(&mut self, group: &str, bc: Bc, dof: Dof, f: FnTimeSpace) -> &mut Self {
        let ids = self.mesh.get_boundary_point_ids_sorted(group);
        for id in ids {
            self.point_bcs.push(PointBc {
                bc,
                dof,
                f,
                point_id: id,
            });
        }
        self
    }

    pub fn add_edge_bc(&mut self, group: &str, bc: Bc, dof: Dof, f: FnTimeSpace) -> &mut Self {
        let keys = self.mesh.get_boundary_edge_keys_sorted(group);
        for key in keys {
            self.edge_bcs.push(EdgeBc {
                bc,
                dof,
                f,
                edge_key: key,
            });
        }
        self
    }

    pub fn set_attribute(&mut self, id: usize, a: Attribute) -> &mut Self {
        self.attributes.insert(id, a);
        self
    }

    pub fn initialize(&mut self) -> Result<(), StrError> {
        // allocate all elements and DOFs and count the number of non-zeros in the kk matrix
        let mut nnz = 0;
        for cell in &self.mesh.cells {
            match self.attributes.get(&cell.attribute_id) {
                Some(attribute) => {
                    if attribute.inactive {
                        continue;
                    }
                    let element = new_element(attribute.kind, &self.mesh, cell.id)?;
                    element.assign_dofs(&mut self.system_dofs);
                    nnz += element.get_nnz();
                    self.elements.push(element);
                }
                None => return Err("cannot find cell with a specific attribute id"),
            };
        }

        // allocate system vector and matrix
        let nequation = self.system_dofs.get_number_of_equations();
        self.system_ff = Vector::new(nequation);
        self.system_kk = SparseTriplet::new(nequation, nequation, nnz, Symmetry::No)?;
        Ok(())
    }

    pub fn run(&mut self) -> Result<(), StrError> {
        let time = 0.0; // time

        // assemble element Ke matrices (can be done in parallel)
        for element in &mut self.elements {
            element.compute_ke()?;
        }

        // add element Ke matrix to global K matrix (must be serial)
        for element in &self.elements {
            element.add_ke_to_kk(&mut self.system_kk)?;
        }

        // add element Fe vector to global F vector (must be serial)
        for element in &self.elements {
            element.add_fe_to_ff(&mut self.system_ff)?;
        }

        // handle edge boundary conditions
        for item in &self.edge_bcs {
            let edge = match self.mesh.boundary_edges.get(&item.edge_key) {
                Some(e) => e,
                None => return Err("cannot find boundary edge"),
            };
            match item.bc {
                Bc::Essential => {
                    for point_id in &edge.points {
                        let x = &self.mesh.points[*point_id].coords;
                        let val = (item.f)(time, x);
                        let eq = self.system_dofs.get_equation(*point_id, item.dof)?;
                        self.system_kk.put(eq, eq, 1.0);
                        self.system_uu[eq] = val;
                        self.system_ff[eq] = val;
                    }
                }
                Bc::Natural => {
                    // perform integration
                }
            }
        }

        // handle point boundary conditions
        for item in &self.point_bcs {
            let x = &self.mesh.points[item.point_id].coords;
            let val = (item.f)(time, x);
            let eq = self.system_dofs.get_equation(item.point_id, item.dof)?;
            match item.bc {
                Bc::Essential => {
                    self.system_kk.put(eq, eq, 1.0);
                    self.system_uu[eq] = val;
                    self.system_ff[eq] = val;
                }
                Bc::Natural => {
                    self.system_ff[eq] = val;
                }
            }
        }

        // solve system
        let config = ConfigSolver::new();
        let mut solver = Solver::new(config)?;
        solver.initialize(&self.system_kk)?;
        solver.factorize()?;
        solver.solve(&mut self.system_uu, &self.system_ff)?;

        Ok(())
    }
}

impl fmt::Display for Simulation {
    /// Prints the simulation data (may be large)
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "system_dofs:\n{}", self.system_dofs).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fem::ElementKind;
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
            0 0  2 4  0 1 2 3
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

        // add boundary conditions
        sim.add_point_bc("origin", Bc::Essential, Dof::Ux, |_, _| 0.0)
            .add_point_bc("origin", Bc::Essential, Dof::Uy, |_, _| 0.0)
            .add_edge_bc("left", Bc::Essential, Dof::Ux, |_, _| 0.0)
            .add_edge_bc("right", Bc::Essential, Dof::Ux, |_, _| 0.0)
            .add_edge_bc("bottom", Bc::Essential, Dof::Uy, |_, _| 0.0);

        // define attributes
        let mut att0 = Attribute::new(ElementKind::Solid);
        att0.set_parameter("Young modulus", 1000.0)
            .set_parameter("Poisson's coefficient", 0.25)
            .set_flag("Plane-stress", true);

        // set attributes
        sim.set_attribute(0, att0);

        // initialize simulation
        sim.initialize()?;

        // run simulation
        sim.run()?;

        println!("{}", sim);

        Ok(())
    }
}
