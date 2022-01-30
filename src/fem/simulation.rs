use super::{Bc, DofIndex, EdgeBc, Element, ElementAttributes, EquationNumbers, FnTimeSpace, PointBc};
use crate::mesh::Mesh;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::{ConfigSolver, Solver, SparseTriplet, Symmetry};
use std::fmt;

pub struct Simulation {
    mesh: Mesh,
    elements: Vec<Box<dyn Element>>,
    equation_numbers: EquationNumbers,
    system_kk: SparseTriplet,
    system_xx: Vector,
    system_rhs: Vector,
    point_bcs: Vec<PointBc>,
    edge_bcs: Vec<EdgeBc>,
}

impl Simulation {
    pub fn new(mesh: Mesh, attributes: &ElementAttributes) -> Result<Self, StrError> {
        // elements and equation numbers (DOFs)
        let npoint = mesh.points.len();
        let mut elements = Vec::<Box<dyn Element>>::new();
        let mut equation_numbers = EquationNumbers::new(npoint);

        // allocate all elements and DOFs and estimate the max number of non-zeros in the K-matrix
        let mut nnz_max = 0;
        for cell in &mesh.cells {
            let attribute = attributes.get(cell.attribute_id)?;
            if attribute.is_active() {
                // let element = new_element(attribute.kind, &mesh, cell.id)?;
                // nnz_max += element.activate_equation_numbers(&mut equation_numbers);
                // elements.push(element);
            }
        }

        // return simulation data
        let neq = equation_numbers.get_number_of_equations();
        Ok(Simulation {
            mesh,
            elements,
            equation_numbers,
            system_kk: SparseTriplet::new(neq, neq, nnz_max, Symmetry::No)?,
            system_xx: Vector::new(neq),
            system_rhs: Vector::new(neq),
            point_bcs: Vec::new(),
            edge_bcs: Vec::new(),
        })
    }

    pub fn run(&mut self) -> Result<(), StrError> {
        let time = 0.0; // time

        // assemble element Ke matrices (can be done in parallel)
        let first_iteration = true;
        for element in &mut self.elements {
            element.compute_local_k_matrix(first_iteration)?;
        }

        // assemble RHS-vector
        for element in &self.elements {
            element.assemble_rhs_vector(&mut self.system_rhs)?;
        }

        // assemble K-matrix
        for element in &self.elements {
            element.assemble_k_matrix(&mut self.system_kk)?;
        }

        /*
        // handle edge boundary conditions
        for item in &self.edge_bcs {
            let edge = match self.mesh.boundary_edges.get(&item.edge_key) {
                Some(e) => e,
                None => return Err("cannot find boundary edge"),
            };
            match item.condition {
                Bc::Essential => {
                    for point_id in &edge.points {
                        let x = &self.mesh.points[*point_id].coords;
                        let val = (item.f)(time, x);
                        let eq = self.equation_numbers.get_equation_number(*point_id, item.dof_index)?;
                        // self.system_kk.put(eq, eq, 1.0);
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
            let eq = self
                .equation_numbers
                .get_equation_number(item.point_id, item.dof_index)?;
            match item.condition {
                Bc::Essential => {
                    // self.system_kk.put(eq, eq, 1.0);
                    self.system_uu[eq] = val;
                    self.system_ff[eq] = val;
                }
                Bc::Natural => {
                    self.system_ff[eq] = val;
                }
            }
        }
        */

        // solve system
        let config = ConfigSolver::new();
        let mut solver = Solver::new(config)?;
        solver.initialize(&self.system_kk)?;
        solver.factorize()?;
        solver.solve(&mut self.system_xx, &self.system_rhs)?;
        Ok(())
    }
}

impl fmt::Display for Simulation {
    /// Prints the simulation data (may be large)
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "equation numbers:\n{}", self.equation_numbers).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fem::{
        ElementAttributes, ElementSolidGroup, ModelSeepageStandard, ModelSolidLinearElastic, DOF_UX, DOF_UY,
    };
    use crate::mesh::At;

    #[test]
    fn run_works() -> Result<(), StrError> {
        // create mesh
        let mut mesh = Mesh::from_text(
            r"
            # sizes
            # ndim npoint ncell
                 2     13     2

            # points
            # id   x   y
               0 0.0 0.0
               1 0.5 0.0
               2 1.0 0.0
               3 1.5 0.0
               4 2.0 0.0
               5 0.0 0.5
               6 1.0 0.5
               7 2.0 0.5
               8 0.0 1.0
               9 0.5 1.0
              10 1.0 1.0
              11 1.5 1.0
              12 2.0 1.0

            # cells
            # id attribute ndim npoint point_ids...
               0         1    2      8 0 2 10 8 1 6 9 5
               1         2    2      8 2 4 12 10 3 7 11 6
            ",
        )?;

        // set groups of points and edges at boundaries
        mesh.set_points("origin", At::XY(0.0, 0.0))?
            .set_edges("left", At::X(0.0))?
            .set_edges("right", At::X(1.0))?
            .set_edges("bottom", At::Y(0.0))?
            .set_edges("top", At::Y(1.0))?;

        // define material models
        let two_dim = mesh.space_ndim == 2;
        let plane_stress = true;
        let model_1 = ModelSolidLinearElastic::new(5_000.0, 0.2, two_dim, plane_stress);
        let model_2 = ModelSolidLinearElastic::new(10_000.0, 0.2, two_dim, plane_stress);
        let model_3 = ModelSeepageStandard::new(1e-4, 1e-4, 1e-4);

        // define attributes
        let att_1 = ElementSolidGroup::new(Box::new(model_1), 0.25);
        let att_2 = ElementSolidGroup::new(Box::new(model_2), 0.25);
        let mut attributes = ElementAttributes::new();
        attributes.set(1, Box::new(att_1))?;
        attributes.set(2, Box::new(att_2))?;

        // let mut attributes = ElementAttributes::new();
        // attributes.set(-1)
        // let mut att0 = ElementAttribute::new(ElementKind::Solid);
        // att0.set_parameter("Young modulus", 1000.0)
        //     .set_parameter("Poisson's coefficient", 0.25)
        //     .set_flag("Plane-stress", true);

        // create simulation
        // let mut sim = Simulation::new(mesh, attributes)?;

        /*
        // add boundary conditions
        sim.add_point_bc("origin", Bc::Essential, DOF_UX, |_, _| 0.0)
            .add_point_bc("origin", Bc::Essential, DOF_UY, |_, _| 0.0)
            .add_edge_bc("left", Bc::Essential, DOF_UX, |_, _| 0.0)
            .add_edge_bc("right", Bc::Essential, DOF_UX, |_, _| 0.0)
            .add_edge_bc("bottom", Bc::Essential, DOF_UY, |_, _| 0.0);

        // set attributes
        sim.set_attribute(0, att0);
        */

        // run simulation
        // sim.run()?;

        // println!("{}", sim);

        Ok(())
    }
}
