//! Implements a simple linear finite element solver

// #![allow(dead_code)]
#![allow(unused_variables)]

macro_rules! print_bc {
    ($b:expr, $c:expr, $k:expr, $g:expr, $t:expr, $x:expr, $f:expr) => {{
        println!("{} {} k:{:?} g:{} p:{} t:{} x:{:?} f:{}", $b, $c, $k, $g, 0, $t, $x, $f);
    }};
}

use crate::mesh::{At, EdgeKey, Group, Index, Mesh};
use crate::StrError;

pub type FnTimeSpace = fn(f64, &[f64]) -> f64;

#[derive(Clone, Copy)]
pub enum BC {
    Essential(FnTimeSpace),
    Natural(FnTimeSpace),
}

#[derive(Clone, Copy)]
pub struct PointBC {
    id: Index,
    condition: BC,
}

#[derive(Clone, Copy)]
pub struct EdgeBC {
    key: EdgeKey,
    condition: BC,
}

pub struct Simulation {
    mesh: Mesh,
    point_bcs: Vec<PointBC>,
    edge_bcs: Vec<EdgeBC>,
}

impl Simulation {
    pub fn new(mesh: Mesh) -> Result<Self, StrError> {
        Ok(Simulation {
            mesh,
            point_bcs: Vec::new(),
            edge_bcs: Vec::new(),
        })
    }

    pub fn add_point_bc(&mut self, group: &str, condition: BC) -> &mut Self {
        let ids = self.mesh.get_boundary_point_ids_sorted(group);
        for id in ids {
            self.point_bcs.push(PointBC { id, condition });
        }
        self
    }

    pub fn add_edge_bc(&mut self, group: &str, condition: BC) -> &mut Self {
        let keys = self.mesh.get_boundary_edge_keys_sorted(group);
        for key in keys {
            self.edge_bcs.push(EdgeBC { key, condition });
        }
        self
    }

    pub fn run(&mut self) {
        let time = 0.0; // time
                        /*
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
                        */
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_works() -> Result<(), StrError> {
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

        mesh.set_group_points("origin", At::XY(0.0, 0.0))?
            .set_group_edges("left", At::X(0.0))?
            .set_group_edges("right", At::X(1.0))?
            .set_group_edges("bottom", At::Y(0.0))?
            .set_group_edges("top", At::Y(1.0))?;

        let mut sim = Simulation::new(mesh)?;

        sim.add_point_bc("origin", BC::Essential(|t, x| 1.0))
            .add_edge_bc("left", BC::Essential(|t, x| 2.0))
            .add_edge_bc("right", BC::Essential(|t, x| 3.0))
            .add_edge_bc("bottom", BC::Essential(|t, x| 3.0))
            .add_edge_bc("top", BC::Natural(|t, x| 3.0));

        sim.run();

        Ok(())
    }
}
