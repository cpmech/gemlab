//! Implements a simple linear finite element solver

#![allow(dead_code)]
#![allow(unused_variables)]

macro_rules! print_bc {
    ($b:expr, $c:expr, $g:expr, $t:expr, $x:expr, $f:expr) => {{
        println!("{} {} g:{} p:{} t:{} x:{:?} f:{}", $b, $c, $g, 0, $t, $x, $f);
    }};
}

use crate::mesh::{EdgeKey, Group, Index, Mesh};
use crate::StrError;

type FnTimeSpace = fn(f64, &[f64]) -> f64;

enum Condition {
    Essential(FnTimeSpace),
    Natural(FnTimeSpace),
}

enum Boundary {
    Point(Group, Condition), // String is Group
    Edge(Group, Condition),  // String is Group
    Face(Group, Condition),  // String is Group
}

pub struct Simulation {
    pub mesh: Mesh,
}

impl Simulation {
    pub fn new(ndim: usize) -> Result<Self, StrError> {
        Ok(Simulation { mesh: Mesh::new(ndim)? })
    }

    pub fn run(&mut self) {
        let boundary_conditions = vec![
            Boundary::Point(1, Condition::Essential(|t, x| 1.0)),
            Boundary::Edge(2, Condition::Essential(|t, x| 2.0)),
            Boundary::Edge(3, Condition::Natural(|t, x| 3.0)),
            Boundary::Face(4, Condition::Essential(|t, x| 4.0)),
        ];
        let time = 0.0; // time
        for boundary_condition in boundary_conditions {
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
                        if let Some(edges) = self.mesh.edge_groups.get(&group) {
                            for key in edges {
                                print_bc!("edge", "essential", group, time, 0, 0);
                            }
                        }
                        let x = &[0.0, 0.0];
                        println!("edge {} essential {}", group, f(0.0, x));
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

    fn print_edge_bc(&self, group: Group, condition: Condition, key: &EdgeKey, time: f64) {
        let point_ids: [Index; 2] = [key.0, key.1];
        for index in point_ids {
            let point = &self.mesh.points[index];
            let x = &point.coords;
            // println!(
            //     "edge {} g:{} p:{} t:{} x:{:?} f:{}",
            //     condition,
            //     group,
            //     index,
            //     time,
            //     x,
            //     f(time, x)
            // );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::parse_mesh;
    // use crate::Mesh;

    #[test]
    fn run_works() -> Result<(), StrError> {
        let mut simulation = Simulation::new(2)?;

        simulation.mesh = parse_mesh(
            r"
            # header
            # ndim npoint ncell
            2 6 2

            # points
            # id group x y
            0 1  0.0 0.0
            1 1  1.0 0.0
            2 11 1.0 1.0
            3 1  0.0 1.0
            4 1  2.0 0.0
            5 1  2.0 1.0

            # cells
            # id group ndim npoint point_ids...
            0 1  2 4  0 1 2 3
            1 8  2 4  1 4 5 2
            ",
        )?;

        simulation.run();

        Ok(())
    }
}
