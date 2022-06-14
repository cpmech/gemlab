use super::{join_meshes, ArgsRing, Block, Constraint2D, Mesh};
use crate::shapes::GeoKind;
use crate::util::{COS_PI_BY_8, ONE_BY_SQRT_2, PI, SIN_PI_BY_8, SQRT_2};
use crate::StrError;

/// Groups generators of structured meshes (Qua and Hex -- sometimes Tri)
pub struct Structured {}

/// TODO Groups generators of unstructured meshes (Tri and Tet only)
pub struct Unstructured {}

impl Structured {
    /// Generates a mesh representing a quarter of a ring in 2D
    ///
    /// ```text
    /// y ^
    ///   |
    ///   ***=---__
    ///   |        '*._
    ///   |            *._
    ///   |               *.
    ///   ***=-__           *.
    ///   |      '-.          *
    ///             *.         *
    ///   |           *         *
    ///                *         *
    ///   |             *         *
    ///                 #         #
    ///   o -   -   -   # ------- # --> x
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

    /// Generates a mesh representing a quarter of a ring in 3D (extruded along the z-direction)
    ///
    /// ```text
    /// y ^
    ///   |
    ///   ***=---__
    ///   |        '*._
    ///   |            *._
    ///   |               *.
    ///   ***=-__           *.
    ///   |      '-.          *
    ///             *.         *
    ///   |           *         *
    ///                *         *
    ///   |             *         *
    ///                 #         #
    /// z o -   -   -   # ------- # --> x
    ///               rmin       rmax
    /// ```
    ///
    /// # Input
    ///
    /// * `rmin` -- inner radius
    /// * `rmax` -- outer radius
    /// * `zmax` -- max z = thickness, since zmin = 0.0
    /// * `nr` -- number of divisions along the radius (must be > 0)
    /// * `na` -- number of divisions along alpha (must be > 0)
    /// * `nz` -- number of divisions along z (must be > 0)
    /// * `target` -- Qua shapes only
    pub fn quarter_ring_3d(
        rmin: f64,
        rmax: f64,
        zmax: f64,
        nr: usize,
        na: usize,
        nz: usize,
        target: GeoKind,
    ) -> Result<Mesh, StrError> {
        let mut block = Block::new_cube(1.0);
        block.set_ndiv(&[nr, na, nz])?;
        block.set_transform_into_ring(Some(ArgsRing {
            amin: 0.0,
            amax: PI / 2.0,
            rmin,
            rmax,
            zmin: 0.0,
            zmax,
        }))?;
        block.subdivide(target)
    }

    /// Generates a mesh representing a quarter of a disk in 2D (A-version)
    ///
    /// ```text
    /// y ^
    ///   .
    ///   #**=---__
    ///   |        '*._
    ///   |            *._
    ///   |               *.
    ///   |                 *.
    ///   |                   *
    ///   |                    *
    ///   |                     *
    ///   |                      *
    ///   |                       *
    ///   |                       #
    ///   o-----------------------# --> x
    ///                        radius
    /// ```
    ///
    /// # Input
    ///
    /// * `radius` -- the radius
    /// * `ndiv` -- number of divisions (must be > 0)
    /// * `target` -- Qua shapes only
    pub fn quarter_disk_2d_a(radius: f64, ndiv: usize, target: GeoKind) -> Result<Mesh, StrError> {
        let m = radius / 2.0;
        let n = radius / SQRT_2;
        let p = 1.15 * m / SQRT_2;
        let mut block_1 = Block::new(&[[0.0, 0.0], [m, 0.0], [p, p], [0.0, m]])?;
        let mut block_2 = Block::new(&[[m, 0.0], [radius, 0.0], [n, n], [p, p]])?;
        let mut block_3 = Block::new(&[[0.0, m], [p, p], [n, n], [0.0, radius]])?;
        block_1.set_ndiv(&[ndiv, ndiv])?;
        block_2.set_ndiv(&[ndiv, ndiv])?;
        block_3.set_ndiv(&[ndiv, ndiv])?;
        let ct = Constraint2D::Circle(0.0, 0.0, radius);
        block_2.set_edge_constraint(1, Some(ct.clone()))?;
        block_3.set_edge_constraint(2, Some(ct))?;
        let mesh_1 = block_1.subdivide(target)?;
        let mesh_2 = block_2.subdivide(target)?;
        let mesh_3 = block_3.subdivide(target)?;
        let mesh_1_2 = join_meshes(&mesh_1, &mesh_2)?;
        join_meshes(&mesh_1_2, &mesh_3)
    }

    /// Generates a mesh representing a quarter of a disk in 2D (B-version)
    ///
    /// ```text
    /// y ^
    ///   .
    ///   #**=---__
    ///   |        '*._
    ///   |            *._
    ///   |               *.
    ///   |                 *.
    ///   |                   *
    ///   |                    *
    ///   |                     *
    ///   |                      *
    ///   |                       *
    ///   |                       #
    ///   o-----------------------# --> x
    ///                        radius
    /// ```
    ///
    /// # Input
    ///
    /// * `radius` -- the radius
    /// * `ndiv` -- number of divisions (must be > 0)
    /// * `target` -- Qua shapes only
    pub fn quarter_disk_2d_b(a: f64, r: f64, ndiv_a: usize, ndiv_b: usize, target: GeoKind) -> Result<Mesh, StrError> {
        if a >= r {
            return Err("'a' must be smaller than 'r'");
        }
        const COS_PI_BY_4: f64 = ONE_BY_SQRT_2;
        let b = r - a;
        let c1 = a * SIN_PI_BY_8;
        let c2 = a * COS_PI_BY_4;
        let c3 = a * COS_PI_BY_8;
        let d1 = r * SIN_PI_BY_8;
        let d2 = r * COS_PI_BY_4;
        let d3 = r * COS_PI_BY_8;
        let m = a / 2.0;
        let n = a + b / 2.0;
        let e = n * COS_PI_BY_4;
        let mut block_1 = Block::new(&[
            [0.0, 0.0],
            [a, 0.0],
            [c2, c2],
            [0.0, a],
            [m, 0.0],
            [c3, c1],
            [c1, c3],
            [0.0, m],
        ])?;
        let mut block_2 = Block::new(&[
            [a, 0.0],
            [r, 0.0],
            [d2, d2],
            [c2, c2],
            [n, 0.0],
            [d3, d1],
            [e, e],
            [c3, c1],
        ])?;
        let mut block_3 = Block::new(&[
            [0.0, a],
            [c2, c2],
            [d2, d2],
            [0.0, r],
            [c1, c3],
            [e, e],
            [d1, d3],
            [0.0, n],
        ])?;
        block_1.set_ndiv(&[ndiv_a, ndiv_a])?;
        block_2.set_ndiv(&[ndiv_b, ndiv_a])?;
        block_3.set_ndiv(&[ndiv_a, ndiv_b])?;
        let ct = Constraint2D::Circle(0.0, 0.0, r);
        block_2.set_edge_constraint(1, Some(ct.clone()))?;
        block_3.set_edge_constraint(2, Some(ct))?;
        let mesh_1 = block_1.subdivide(target)?;
        let mesh_2 = block_2.subdivide(target)?;
        let mesh_3 = block_3.subdivide(target)?;
        let mesh_1_2 = join_meshes(&mesh_1, &mesh_2)?;
        join_meshes(&mesh_1_2, &mesh_3)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Structured;
    use crate::geometry::point_point_distance;
    use crate::mesh::{check_overlapping_points, draw_mesh};
    use crate::shapes::GeoKind;
    use crate::StrError;
    use russell_chk::assert_approx_eq;

    #[test]
    fn quarter_ring_2d_works() -> Result<(), StrError> {
        let mesh = Structured::quarter_ring_2d(3.0, 6.0, 1, 1, GeoKind::Qua16)?;
        check_overlapping_points(&mesh, 0.02)?;
        assert_eq!(mesh.points.len(), 16);
        assert_eq!(mesh.cells.len(), 1);
        for p in [0, 11, 7, 3] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0])?;
            assert_approx_eq!(d, 3.0, 1e-15);
        }
        for p in [4, 12, 15, 10] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0])?;
            assert_approx_eq!(d, 4.0, 1e-15);
        }
        for p in [8, 13, 14, 6] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0])?;
            assert_approx_eq!(d, 5.0, 1e-15);
        }
        for p in [1, 5, 9, 2] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0])?;
            assert_approx_eq!(d, 6.0, 1e-15);
        }
        if false {
            draw_mesh(mesh, true, "/tmp/gemlab/test_quarter_ring_2d.svg")?;
        }
        Ok(())
    }

    #[test]
    fn quarter_ring_3d_works() -> Result<(), StrError> {
        let mesh = Structured::quarter_ring_3d(3.0, 6.0, 2.0, 1, 2, 1, GeoKind::Hex32)?;
        check_overlapping_points(&mesh, 0.02)?;
        assert_eq!(mesh.points.len(), 52);
        assert_eq!(mesh.cells.len(), 2);
        for p in [0, 15, 14, 3, 41, 40, 33] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0])?;
            assert_approx_eq!(d, 3.0, 1e-15);
        }
        for p in [1, 10, 11, 2, 36, 37, 32] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0])?;
            assert_approx_eq!(d, 6.0, 1e-15);
        }
        if false {
            draw_mesh(mesh, true, "/tmp/gemlab/test_quarter_ring_3d.svg")?;
        }
        Ok(())
    }

    #[test]
    fn quarter_disk_2d_a_works_qua8() -> Result<(), StrError> {
        let mesh = Structured::quarter_disk_2d_a(6.0, 1, GeoKind::Qua8)?;
        check_overlapping_points(&mesh, 0.02)?;
        assert_eq!(mesh.points.len(), 16);
        assert_eq!(mesh.cells.len(), 3);
        for p in [8, 11, 9, 14, 13] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0])?;
            assert_approx_eq!(d, 6.0, 1e-15);
        }
        if false {
            draw_mesh(mesh, true, "/tmp/gemlab/test_quarter_disk_2d_a_qua8.svg")?;
        }
        Ok(())
    }

    #[test]
    fn quarter_disk_2d_a_works_qua8_finer() -> Result<(), StrError> {
        let mesh = Structured::quarter_disk_2d_a(6.0, 3, GeoKind::Qua8)?;
        check_overlapping_points(&mesh, 0.02)?;
        assert_eq!(mesh.points.len(), 100);
        assert_eq!(mesh.cells.len(), 27);
        for p in [50, 53, 51, 62, 61, 71, 70, 99, 96, 98, 91, 94, 92] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0])?;
            assert_approx_eq!(d, 6.0, 1e-15);
        }
        if false {
            draw_mesh(mesh, true, "/tmp/gemlab/test_quarter_disk_2d_a_qua8_finer.svg")?;
        }
        Ok(())
    }

    #[test]
    fn quarter_disk_2d_a_works_qua16() -> Result<(), StrError> {
        let mesh = Structured::quarter_disk_2d_a(6.0, 1, GeoKind::Qua16)?;
        check_overlapping_points(&mesh, 0.02)?;
        assert_eq!(mesh.points.len(), 37);
        assert_eq!(mesh.cells.len(), 3);
        for p in [16, 19, 22, 17, 29, 31, 28] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0])?;
            assert_approx_eq!(d, 6.0, 1e-15);
        }
        if false {
            draw_mesh(mesh, true, "/tmp/gemlab/test_quarter_disk_2d_a_qua16.svg")?;
        }
        Ok(())
    }

    #[test]
    fn quarter_disk_2d_b_works_qua8() -> Result<(), StrError> {
        let mesh = Structured::quarter_disk_2d_b(3.0, 6.0, 1, 1, GeoKind::Qua8)?;
        check_overlapping_points(&mesh, 0.02)?;
        assert_eq!(mesh.points.len(), 16);
        assert_eq!(mesh.cells.len(), 3);
        for p in [8, 11, 9, 14, 13] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0])?;
            assert_approx_eq!(d, 6.0, 1e-15);
        }
        if false {
            draw_mesh(mesh, true, "/tmp/gemlab/test_quarter_disk_2d_b_qua8.svg")?;
        }
        Ok(())
    }

    #[test]
    fn quarter_disk_2d_b_works_qua8_finer() -> Result<(), StrError> {
        let mesh = Structured::quarter_disk_2d_b(3.0, 6.0, 3, 3, GeoKind::Qua8)?;
        check_overlapping_points(&mesh, 0.02)?;
        assert_eq!(mesh.points.len(), 100);
        assert_eq!(mesh.cells.len(), 27);
        for p in [50, 53, 51, 62, 61, 71, 70, 99, 96, 98, 91, 94, 92] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0])?;
            assert_approx_eq!(d, 6.0, 1e-15);
        }
        if false {
            draw_mesh(mesh, true, "/tmp/gemlab/test_quarter_disk_2d_b_qua8_finer.svg")?;
        }
        Ok(())
    }

    #[test]
    fn quarter_disk_2d_b_works_qua16() -> Result<(), StrError> {
        let mesh = Structured::quarter_disk_2d_b(3.0, 6.0, 1, 1, GeoKind::Qua16)?;
        check_overlapping_points(&mesh, 0.02)?;
        assert_eq!(mesh.points.len(), 37);
        assert_eq!(mesh.cells.len(), 3);
        for p in [16, 19, 22, 17, 29, 31, 28] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0])?;
            assert_approx_eq!(d, 6.0, 1e-15);
        }
        if false {
            draw_mesh(mesh, true, "/tmp/gemlab/test_quarter_disk_2d_b_qua16.svg")?;
        }
        Ok(())
    }
}
