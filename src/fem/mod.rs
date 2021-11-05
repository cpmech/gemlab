#![allow(dead_code)]
#![allow(unused_variables)]

use crate::{Group, Mesh};

type FnTimeSpace = fn(f64, &[f64]) -> f64;

enum Condition {
    Essential(FnTimeSpace),
    Natural(FnTimeSpace),
}

enum Boundary {
    Point(Group, Condition),
    Edge(Group, Condition),
    Face(Group, Condition),
}

pub fn run(mesh: &Mesh) {
    let boundary_conditions: Vec<Boundary> = vec![
        Boundary::Point(1, Condition::Essential(|t, x| 1.0)),
        Boundary::Edge(1, Condition::Essential(|t, x| 2.0)),
        Boundary::Edge(2, Condition::Natural(|t, x| 3.0)),
        Boundary::Face(3, Condition::Essential(|t, x| 4.0)),
    ];
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
                    if let Some(edges) = mesh.edge_groups.get(&group) {
                        for key in edges {
                            println!("{}, {}", key.0, key.1);
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Mesh;

    #[test]
    fn run_works() -> Result<(), &'static str> {
        let mesh = Mesh::new_zeroed(2, 4, 1)?;
        run(&mesh);
        Ok(())
    }
}
