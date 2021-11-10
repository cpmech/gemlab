//! Implements a simple linear finite element solver

// #![allow(dead_code)]
#![allow(unused_variables)]

macro_rules! print_bc {
    ($b:expr, $c:expr, $k:expr, $g:expr, $t:expr, $x:expr, $f:expr) => {{
        println!("{} {} k:{:?} g:{} p:{} t:{} x:{:?} f:{}", $b, $c, $k, $g, 0, $t, $x, $f);
    }};
}

use crate::mesh::{Group, Mesh};
use crate::StrError;

pub type FnTimeSpace = fn(f64, &[f64]) -> f64;

pub enum Condition {
    Essential(FnTimeSpace),
    Natural(FnTimeSpace),
}

pub enum Boundary {
    Point(Group, Condition), // String is Group
    Edge(Group, Condition),  // String is Group
    Face(Group, Condition),  // String is Group
}

pub struct Simulation {
    mesh: Mesh,
    boundary_conditions: Vec<Boundary>,
}

impl Simulation {
    pub fn new(mesh: Mesh) -> Result<Self, StrError> {
        Ok(Simulation {
            mesh,
            boundary_conditions: Vec::new(),
        })
    }

    pub fn add_bc(&mut self, bc: Boundary) -> &mut Self {
        self.boundary_conditions.push(bc);
        self
    }

    pub fn run(&mut self) {
        let time = 0.0; // time
        for boundary_condition in &self.boundary_conditions {
            match boundary_condition {
                Boundary::Point(group, condition) => match condition {
                    Condition::Essential(f) => {
                        let x = &[0.0, 0.0];
                        println!("point {} essential {}", group, f(0.0, x));
                    }
                    Condition::Natural(f) => {
                        let x = &[0.0, 0.0];
                        println!("point {} natural {}", group, f(0.0, x));
                    }
                },
                Boundary::Edge(group, condition) => match condition {
                    Condition::Essential(f) => {
                        if let Some(edges) = self.mesh.groups_boundary_edges.get(&group) {
                            for key in edges {
                                if let Some(edge) = self.mesh.boundary_edges.get(key) {
                                    for id in &edge.points {
                                        let point = &self.mesh.points[*id];
                                        let x = &point.coords;
                                        let f_val = f(time, x);
                                        print_bc!("edge", "essential", key, group, time, x, f_val);
                                    }
                                }
                            }
                        }
                    }
                    Condition::Natural(f) => {
                        let x = &[0.0, 0.0];
                        println!("edge {} natural {}", group, f(0.0, x));
                    }
                },
                Boundary::Face(group, condition) => match condition {
                    Condition::Essential(f) => {
                        let x = &[0.0, 0.0];
                        println!("face {} essential {}", group, f(0.0, x));
                    }
                    Condition::Natural(f) => {
                        let x = &[0.0, 0.0];
                        println!("face {} natural {}", group, f(0.0, x));
                    }
                },
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_works() -> Result<(), StrError> {
        let mesh = Mesh::from_text(
            r"
            # header
            # ndim npoint ncell
            2 6 2

            # points
            # id group x y
            0 0  0.0 0.0
            1 0  1.0 0.0
            2 1  1.0 1.0
            3 0  0.0 1.0
            4 0  2.0 0.0
            5 0  2.0 1.0

            # cells
            # id group ndim npoint point_ids...
            0 1  2 4  0 1 2 3
            1 0  2 4  1 4 5 2
            ",
        )?;

        let mut sim = Simulation::new(mesh)?;

        sim.add_bc(Boundary::Point(1, Condition::Essential(|t, x| 1.0)));
        sim.add_bc(Boundary::Edge(2, Condition::Essential(|t, x| 2.0)));
        sim.add_bc(Boundary::Edge(3, Condition::Natural(|t, x| 3.0)));
        sim.add_bc(Boundary::Face(4, Condition::Essential(|t, x| 4.0)));

        sim.run();

        Ok(())
    }
}
