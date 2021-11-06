#![allow(unused)]

use crate::{GridSearch, Mesh, StrError};

/// Defines the location of points
pub enum At {
    /// At (x,y)
    XY(f64, f64),

    /// At (x,y,z)
    XYZ(f64, f64, f64),

    /// Horizontal line passing at (x,y)
    Horizontal(f64, f64),

    /// Vertical line passing at (x,y)
    Vertical(f64, f64),

    /// Circumference of circle at (x,y,radius)
    Circle(f64, f64, f64),

    /// Circumference of circle at (x,y,z,radius)
    Circle3D(f64, f64, f64, f64),

    /// Surface of cylinder parallel to x (x,y,z,radius)
    CylinderX(f64, f64, f64, f64),

    /// Surface of cylinder parallel to y (x,y,z,radius)
    CylinderY(f64, f64, f64, f64),

    /// Surface of cylinder parallel to z (x,y,z,radius)
    CylinderZ(f64, f64, f64, f64),

    /// At plane with normal along x and passing at (x,y,z)
    PlaneNormalX(f64, f64, f64),

    /// At plane with normal along y and passing at (x,y,z)
    PlaneNormalY(f64, f64, f64),

    /// At plane with normal along z and passing at (x,y,z)
    PlaneNormalZ(f64, f64, f64),
}

/// Defines geometric features by locating individual points or points along edges or faces
pub enum Geo {
    Point(At),
    Edge(At),
    Face(At),
}

pub struct Domain {
    mesh: Mesh,
    grid: GridSearch,
}

impl Domain {
    pub fn new(mesh: Mesh) -> Self {
        const GRID_NDIV: usize = 20;
        const GRID_TOL: f64 = 1e-4;
        let ndim = mesh.ndim;
        let min = &mesh.min.clone();
        let max = &mesh.min.clone();
        let ndiv = &vec![GRID_NDIV; ndim];
        let tol = &vec![GRID_TOL; ndim];
        Domain {
            mesh,
            grid: GridSearch::new(ndiv, min, max, tol).unwrap(),
        }
    }

    pub fn set_group(&mut self, name: &str, geo: Geo) -> Result<(), StrError> {
        let mut points: Vec<usize> = Vec::new();
        match geo {
            Geo::Point(at) => {
                for index in 0..self.mesh.points.len() {
                    if self.point_belongs_to(index, &at) {
                        points.push(index);
                    }
                }
            }
            Geo::Edge(at) => {}
            Geo::Face(at) => {}
        };
        Ok(())
    }

    fn point_belongs_to(&self, index: usize, at: &At) -> bool {
        let coords = &self.mesh.points[index].coords;
        match at {
            At::XY(x, y) => false,
            At::XYZ(x, y, z) => false,
            At::Horizontal(x, y) => false,
            At::Vertical(x, y) => false,
            At::Circle(x, y, z) => false,
            At::Circle3D(x, y, z, r) => false,
            At::CylinderX(x, y, z, r) => false,
            At::CylinderY(x, y, z, r) => false,
            At::CylinderZ(x, y, z, r) => false,
            At::PlaneNormalX(x, y, z) => false,
            At::PlaneNormalY(x, y, z) => false,
            At::PlaneNormalZ(x, y, z) => false,
        }
    }

    fn near_points(&self, a: &[f64], b: &[f64]) -> bool {
        let tol = 1e-15;
        for i in 0..a.len() {
            if f64::abs(a[i] - b[i]) > tol {
                return false;
            }
        }
        true
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{parse_mesh, At, Domain, Geo, StrError};

    #[test]
    fn test() -> Result<(), StrError> {
        return Ok(());

        let mut mesh = parse_mesh(
            r"
            2 4 1
            0 1 0.0 0.0
            1 1 1.0 0.0
            2 1 1.0 1.0
            3 1 0.0 1.0
            0 1  2 4  0 1 2 3
        ",
        )?;

        let mut domain = Domain::new(mesh);

        domain.set_group("origin", Geo::Point(At::XY(0.0, 0.0)))?;
        domain.set_group("origin", Geo::Point(At::Horizontal(0.0, 0.0)))?;
        domain.set_group("origin", Geo::Point(At::Vertical(0.0, 0.0)))?;
        domain.set_group("origin", Geo::Point(At::XYZ(0.0, 0.0, 0.0)))?;
        domain.set_group("bedrock", Geo::Edge(At::Horizontal(0.0, 0.0)))?;
        domain.set_group("left", Geo::Edge(At::Vertical(0.0, 0.0)))?;
        domain.set_group("right", Geo::Edge(At::Vertical(0.0, 0.0)))?;
        domain.set_group("tunnel", Geo::Face(At::Circle(0.0, 0.0, 0.5)))?;
        domain.set_group("tunnel", Geo::Face(At::Circle3D(0.0, 0.0, 0.0, 0.5)))?;
        domain.set_group("tunnel", Geo::Face(At::CylinderX(0.0, 0.0, 0.0, 0.5)))?;
        domain.set_group("bedrock", Geo::Face(At::PlaneNormalZ(0.0, 0.0, 0.0)))?;

        println!("{}", domain.mesh);

        Ok(())
    }
}
