use crate::{Cell, Edge, Face, Point};
use std::collections::HashMap;

#[derive(Hash, Eq, PartialEq, Debug)]
pub struct KeyPoint {
    pub x: u64,
    pub y: u64,
    pub z: u64,
}

#[derive(Hash, Eq, PartialEq, Debug)]
pub struct KeyEdge {
    pub a: usize,
    pub b: usize,
}

#[derive(Hash, Eq, PartialEq, Debug)]
pub struct KeyFace {
    pub a: usize,
    pub b: usize,
    pub c: usize,
}

pub struct Mesh {
    pub ndim: usize,
    pub points: HashMap<KeyPoint, Point>,
    pub edges: HashMap<KeyEdge, Edge>,
    pub faces: HashMap<KeyFace, Face>,
    pub cells: Vec<Cell>,
}

impl Mesh {
    pub fn new(ndim: usize) -> Self {
        Mesh {
            ndim,
            points: HashMap::<KeyPoint, Point>::new(),
            edges: HashMap::<KeyEdge, Edge>::new(),
            faces: HashMap::<KeyFace, Face>::new(),
            cells: Vec::<Cell>::new(),
        }
    }

    pub fn print(&self) {
        for point in self.points.values() {
            println!("{:?}", point);
        }
        println!();
        for edge in self.edges.values() {
            println!("{:?}", edge);
        }
        println!();
        for face in self.faces.values() {
            println!("{:?}", face);
        }
        println!();
        for cell in &self.cells {
            println!("{:?}", cell);
        }
    }
}
