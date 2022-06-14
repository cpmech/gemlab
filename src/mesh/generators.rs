use super::{join_meshes, ArgsRing, Block, Constraint2D, Constraint3D, Mesh};
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
    ///   .      '-.          *
    ///             *.         *
    ///   .           *         *
    ///                *         *
    ///   .             *         *
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
    /// * `target` -- [crate::shapes::GeoClass::Qua] shapes only
    pub fn quarter_ring_2d(rmin: f64, rmax: f64, nr: usize, na: usize, target: GeoKind) -> Result<Mesh, StrError> {
        let mut block = Block::new_square(1.0);
        block.set_ndiv(&[nr, na])?;
        block.set_transform_into_ring(Some(ArgsRing {
            amin: 0.0,
            amax: PI / 2.0,
            rmin,
            rmax,
            zmin: 0.0, // ignored
            zmax: 1.0, // ignored
        }))?;
        block.subdivide(target)
    }

    /// Generates a mesh representing a quarter of a ring in 3D (extrusion along the z-direction)
    ///
    /// ```text
    /// y ^
    ///   |
    ///   ***=---__
    ///   |        '*._
    ///   |            *._
    ///   |               *.
    ///   ***=-__           *.
    ///   .      '-.          *
    ///             *.         *
    ///   .           *         *
    ///                *         *
    ///   .             *         *
    ///                 #         #
    /// z o -   -   -   # ------- # --> x
    ///               rmin       rmax
    /// ```
    ///
    /// # Input
    ///
    /// * `rmin` -- inner radius
    /// * `rmax` -- outer radius
    /// * `z` -- thickness (zmin = 0.0)
    /// * `nr` -- number of divisions along the radius (must be > 0)
    /// * `na` -- number of divisions along alpha (must be > 0)
    /// * `nz` -- number of divisions along z (thickness) (must be > 0)
    /// * `target` -- [crate::shapes::GeoClass::Qua] shapes only
    pub fn quarter_ring_3d(
        rmin: f64,
        rmax: f64,
        z: f64,
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
            zmax: z,
        }))?;
        block.subdivide(target)
    }

    /// Generates a mesh representing a quarter of a disk in 2D (A-version)
    ///
    /// ```text
    /// y ^
    ///   .
    /// r #**=---__
    ///   |        '*._
    ///   |            *._
    ///   |              .*.
    ///   |            .'   *.
    ///   |          .'       *
    ///   |---------'          *
    ///   |         |           *
    ///   |         |            *
    ///   |         |             *
    ///   |         |             #
    ///   o-----------------------# --> x
    ///                           r
    ///   |<---a--->|<-----b----->|
    ///   |<----------r---------->|
    /// ```
    ///
    /// # Input
    ///
    /// * `r` -- the radius
    /// * `na` -- number of divisions along 'a' (must be > 0)
    /// * `nb` -- number of divisions along 'b' (must be > 0)
    /// * `target` -- [crate::shapes::GeoClass::Qua] shapes only
    pub fn quarter_disk_2d_a(r: f64, na: usize, nb: usize, target: GeoKind) -> Result<Mesh, StrError> {
        let m = r / 2.0;
        let n = r / SQRT_2;
        let p = 1.15 * m / SQRT_2;
        let mut block_1 = Block::new(&[[0.0, 0.0], [m, 0.0], [p, p], [0.0, m]])?;
        let mut block_2 = Block::new(&[[m, 0.0], [r, 0.0], [n, n], [p, p]])?;
        let mut block_3 = Block::new(&[[0.0, m], [p, p], [n, n], [0.0, r]])?;
        block_1.set_ndiv(&[na, na])?;
        block_2.set_ndiv(&[nb, na])?;
        block_3.set_ndiv(&[na, nb])?;
        let ct = Constraint2D::Circle(0.0, 0.0, r);
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
    /// r #**=---__
    ///   |        '*._
    ///   |            *._
    ///   |               ..
    ///   |             .'  *.
    ///   |**=-__     .'      *
    ///   |      '*..'         *
    ///   |         \           *
    ///   |          '           *
    ///   |           |           *
    ///   |           |           #
    ///   o-----------------------# --> x
    ///                           r
    ///   |<----a---->|<----b---->|
    ///   |<----------r---------->|
    /// ```
    ///
    /// # Input
    ///
    /// * `a` -- distance such that the outer cells belong to a ring (must be < r)
    /// * `r` -- the radius
    /// * `na` -- number of divisions along 'a' (must be > 0)
    /// * `nb` -- number of divisions along 'b' (must be > 0)
    /// * `target` -- [crate::shapes::GeoClass::Qua] shapes only
    pub fn quarter_disk_2d_b(a: f64, r: f64, na: usize, nb: usize, target: GeoKind) -> Result<Mesh, StrError> {
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
        #[rustfmt::skip]
        let mut block_1 = Block::new(&[
            [0.0, 0.0],
            [  a, 0.0],
            [ c2,  c2],
            [0.0,   a],
            [  m, 0.0],
            [ c3,  c1],
            [ c1,  c3],
            [0.0,   m],
        ])?;
        #[rustfmt::skip]
        let mut block_2 = Block::new(&[
            [ a, 0.0],
            [ r, 0.0],
            [d2,  d2],
            [c2,  c2],
            [ n, 0.0],
            [d3,  d1],
            [ e,   e],
            [c3,  c1],
        ])?;
        #[rustfmt::skip]
        let mut block_3 = Block::new(&[
            [0.0,  a],
            [ c2, c2],
            [ d2, d2],
            [0.0,  r],
            [ c1, c3],
            [ e,   e],
            [ d1, d3],
            [0.0,  n],
        ])?;
        block_1.set_ndiv(&[na, na])?;
        block_2.set_ndiv(&[nb, na])?;
        block_3.set_ndiv(&[na, nb])?;
        let ct = Constraint2D::Circle(0.0, 0.0, r);
        block_2.set_edge_constraint(1, Some(ct.clone()))?;
        block_3.set_edge_constraint(2, Some(ct))?;
        let mesh_1 = block_1.subdivide(target)?;
        let mesh_2 = block_2.subdivide(target)?;
        let mesh_3 = block_3.subdivide(target)?;
        let mesh_1_2 = join_meshes(&mesh_1, &mesh_2)?;
        join_meshes(&mesh_1_2, &mesh_3)
    }

    /// Generates a mesh representing a quarter of a disk in 3D (extrusion along the z-direction) (A-version)
    ///
    /// ```text
    /// y ^
    ///   .
    /// r #**=---__
    ///   |        '*._
    ///   |            *._
    ///   |              .*.
    ///   |            .'   *.
    ///   |          .'       *
    ///   |---------'          *
    ///   |         |           *
    ///   |         |            *
    ///   |         |             *
    ///   |         |             #
    /// z o-----------------------# --> x
    ///                           r
    ///   |<---a--->|<-----b----->|
    ///   |<----------r---------->|
    /// ```
    ///
    /// # Input
    ///
    /// * `r` -- radius
    /// * `z` -- thickness (zmin = 0.0)
    /// * `na` -- number of divisions along 'a' on the x-y plane (must be > 0)
    /// * `nb` -- number of divisions along 'b' on the x-y plane (must be > 0)
    /// * `nz` -- number of divisions along 'z' (thickness) (must be > 0)
    /// * `target` -- [crate::shapes::GeoClass::Hex] shapes only
    pub fn quarter_disk_3d_a(
        r: f64,
        z: f64,
        na: usize,
        nb: usize,
        nz: usize,
        target: GeoKind,
    ) -> Result<Mesh, StrError> {
        let m = r / 2.0;
        let n = r / SQRT_2;
        let p = 1.15 * m / SQRT_2;
        #[rustfmt::skip]
        let mut block_1 = Block::new(&[
            [0.0, 0.0, 0.0],
            [  m, 0.0, 0.0],
            [  p,   p, 0.0],
            [0.0,   m, 0.0],
            [0.0, 0.0,   z],
            [  m, 0.0,   z],
            [  p,   p,   z],
            [0.0,   m,   z],
        ])?;
        #[rustfmt::skip]
        let mut block_2 = Block::new(&[
            [m, 0.0, 0.0],
            [r, 0.0, 0.0],
            [n,   n, 0.0],
            [p,   p, 0.0],
            [m, 0.0,   z],
            [r, 0.0,   z],
            [n,   n,   z],
            [p,   p,   z],
        ])?;
        #[rustfmt::skip]
        let mut block_3 = Block::new(&[
            [0.0, m, 0.0],
            [  p, p, 0.0],
            [  n, n, 0.0],
            [0.0, r, 0.0],
            [0.0, m,   z],
            [  p, p,   z],
            [  n, n,   z],
            [0.0, r,   z],
        ])?;
        block_1.set_ndiv(&[na, na, nz])?;
        block_2.set_ndiv(&[nb, na, nz])?;
        block_3.set_ndiv(&[na, nb, nz])?;
        let ct = Constraint3D::CylinderZ(0.0, 0.0, r);
        block_2.set_face_constraint(1, Some(ct.clone()))?;
        block_3.set_face_constraint(3, Some(ct))?;
        let mesh_1 = block_1.subdivide(target)?;
        let mesh_2 = block_2.subdivide(target)?;
        let mesh_3 = block_3.subdivide(target)?;
        let mesh_1_2 = join_meshes(&mesh_1, &mesh_2)?;
        join_meshes(&mesh_1_2, &mesh_3)
    }

    /// Generates a mesh representing a quarter of a disk in 3D (extrusion along the z-direction) (B-version)
    ///
    /// ```text
    /// y ^
    ///   .
    /// r #**=---__
    ///   |        '*._
    ///   |            *._
    ///   |               ..
    ///   |             .'  *.
    ///   |**=-__     .'      *
    ///   |      '*..'         *
    ///   |         \           *
    ///   |          '           *
    ///   |           |           *
    ///   |           |           #
    /// z o-----------------------# --> x
    ///                           r
    ///   |<----a---->|<----b---->|
    ///   |<----------r---------->|
    /// ```
    ///
    /// # Input
    ///
    /// * `a` -- distance such that the outer cells belong to a ring (must be < r)
    /// * `r` -- the radius
    /// * `z` -- thickness
    /// * `na` -- number of divisions along 'a' on the x-y plane (must be > 0)
    /// * `nb` -- number of divisions along 'b' on the x-y plane (must be > 0)
    /// * `nz` -- number of divisions along 'z' (thickness) (must be > 0)
    /// * `target` -- [crate::shapes::GeoClass::Hex] shapes only
    pub fn quarter_disk_3d_b(
        a: f64,
        r: f64,
        z: f64,
        na: usize,
        nb: usize,
        nz: usize,
        target: GeoKind,
    ) -> Result<Mesh, StrError> {
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
        let f = z / 2.0;
        #[rustfmt::skip]
        let mut block_1 = Block::new(&[
            [0.0, 0.0, 0.0],
            [  a, 0.0, 0.0],
            [ c2,  c2, 0.0],
            [0.0,   a, 0.0],
            [0.0, 0.0,   z],
            [  a, 0.0,   z],
            [ c2,  c2,   z],
            [0.0,   a,   z],
            [  m, 0.0, 0.0],
            [ c3,  c1, 0.0],
            [ c1,  c3, 0.0],
            [0.0,   m, 0.0],
            [  m, 0.0,   z],
            [ c3,  c1,   z],
            [ c1,  c3,   z],
            [0.0,   m,   z],
            [0.0, 0.0,   f],
            [  a, 0.0,   f],
            [ c2,  c2,   f],
            [0.0,   a,   f],
        ])?;
        #[rustfmt::skip]
        let mut block_2 = Block::new(&[
            [ a, 0.0, 0.0],
            [ r, 0.0, 0.0],
            [d2,  d2, 0.0],
            [c2,  c2, 0.0],
            [ a, 0.0,   z],
            [ r, 0.0,   z],
            [d2,  d2,   z],
            [c2,  c2,   z],
            [ n, 0.0, 0.0],
            [d3,  d1, 0.0],
            [ e,   e, 0.0],
            [c3,  c1, 0.0],
            [ n, 0.0,   z],
            [d3,  d1,   z],
            [ e,   e,   z],
            [c3,  c1,   z],
            [ a, 0.0,   f],
            [ r, 0.0,   f],
            [d2,  d2,   f],
            [c2,  c2,   f],
        ])?;
        #[rustfmt::skip]
        let mut block_3 = Block::new(&[
            [0.0,  a, 0.0],
            [ c2, c2, 0.0],
            [ d2, d2, 0.0],
            [0.0,  r, 0.0],
            [0.0,  a,   z],
            [ c2, c2,   z],
            [ d2, d2,   z],
            [0.0,  r,   z],
            [ c1, c3, 0.0],
            [ e,   e, 0.0],
            [ d1, d3, 0.0],
            [0.0,  n, 0.0],
            [ c1, c3,   z],
            [ e,   e,   z],
            [ d1, d3,   z],
            [0.0,  n,   z],
            [0.0,  a,   f],
            [ c2, c2,   f],
            [ d2, d2,   f],
            [0.0,  r,   f],
        ])?;
        block_1.set_ndiv(&[na, na, nz])?;
        block_2.set_ndiv(&[nb, na, nz])?;
        block_3.set_ndiv(&[na, nb, nz])?;
        let ct = Constraint3D::CylinderZ(0.0, 0.0, r);
        block_2.set_face_constraint(1, Some(ct.clone()))?;
        block_3.set_face_constraint(3, Some(ct))?;
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
        let mesh = Structured::quarter_disk_2d_a(6.0, 1, 1, GeoKind::Qua8)?;
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
        let mesh = Structured::quarter_disk_2d_a(6.0, 3, 3, GeoKind::Qua8)?;
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
        let mesh = Structured::quarter_disk_2d_a(6.0, 1, 1, GeoKind::Qua16)?;
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

    #[test]
    fn quarter_disk_3d_a_works_hex32() -> Result<(), StrError> {
        let mesh = Structured::quarter_disk_3d_a(6.0, 1.5, 1, 1, 1, GeoKind::Hex32)?;
        check_overlapping_points(&mesh, 0.02)?;
        assert_eq!(mesh.cells.len(), 3);
        for p in [
            32, 38, 39, 33, 54, 55, 52, // z-min
            48, 50, 62, // z-1/3
            49, 51, 63, // z-2/3
            34, 44, 45, 35, 58, 59, 53, // z-max
        ] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0])?;
            assert_approx_eq!(d, 6.0, 1e-15);
        }
        if false {
            draw_mesh(mesh, true, "/tmp/gemlab/test_quarter_disk_3d_a_hex32.svg")?;
        }
        Ok(())
    }

    #[test]
    fn quarter_disk_3d_b_works_hex32() -> Result<(), StrError> {
        let mesh = Structured::quarter_disk_3d_b(3.0, 6.0, 1.5, 1, 1, 1, GeoKind::Hex32)?;
        check_overlapping_points(&mesh, 0.02)?;
        assert_eq!(mesh.cells.len(), 3);
        for p in [
            32, 38, 39, 33, 54, 55, 52, // z-min
            48, 50, 62, // z-1/3
            49, 51, 63, // z-2/3
            34, 44, 45, 35, 58, 59, 53, // z-max
        ] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0])?;
            assert_approx_eq!(d, 6.0, 1e-15);
        }
        if false {
            draw_mesh(mesh, true, "/tmp/gemlab/test_quarter_disk_3d_b_hex32.svg")?;
        }
        Ok(())
    }
}
