use super::{ArgsRing, Block, Mesh};
use crate::shapes::GeoKind;
use crate::util::PI;
use crate::StrError;

/// Groups generators of structured meshes (Qua and Hex -- sometimes Tri)
pub struct Structured {}

/// TODO Groups generators of unstructured meshes (Tri and Tet only)
pub struct Unstructured {}

impl Structured {
    /// Generates a mesh representing a quarter of a ring in 2D
    ///
    /// ```text
    ///   ^
    ///   |
    ///   ***=---__
    ///   |        '*._
    ///   |            *._
    ///   |               *.
    ///   ***=-__           *.
    ///   |      '-.          *
    ///   |         *.         *
    ///   |           *         *
    ///   |            *         *
    ///   |             *         *
    ///   |             #         #
    ///   o ----------- # ------- # --> r
    ///               rmin       rmax
    /// ```
    ///
    /// # Input
    ///
    /// * `rmin` -- inner radius
    /// * `rmax` -- outer radius
    /// * `nr` -- number of divisions along the radius (must be > 0)
    /// * `na` -- number of divisions along alpha (must be > 0)
    /// * `target` -- Qua shapes only
    pub fn quarter_ring_2d(rmin: f64, rmax: f64, nr: usize, na: usize, target: GeoKind) -> Result<Mesh, StrError> {
        let mut block = Block::new_square(1.0);
        block.set_ndiv(&[nr, na])?;
        block.set_transform_into_ring(Some(ArgsRing {
            amin: 0.0,
            amax: PI / 2.0,
            rmin,
            rmax,
            zmin: 0.0,
            zmax: 1.0,
        }))?;
        block.subdivide(target)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Structured;
    use crate::mesh::draw_mesh;
    use crate::shapes::GeoKind;
    use crate::StrError;

    #[test]
    fn quarter_ring_2d_works() -> Result<(), StrError> {
        let mesh = Structured::quarter_ring_2d(3.0, 6.0, 1, 1, GeoKind::Qua16)?;
        assert_eq!(mesh.points.len(), 16);
        assert_eq!(mesh.cells.len(), 1);
        if false {
            draw_mesh(mesh, true, "/tmp/gemlab/test_quarter_ring_2d.svg")?;
        }
        Ok(())
    }
}
