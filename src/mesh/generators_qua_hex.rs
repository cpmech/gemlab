use super::{join_meshes, ArgsRing, Block, Constraint2D, Constraint3D, Graph, Mesh};
use crate::shapes::GeoKind;
use crate::StrError;
use russell_lab::math::{COS_PI_BY_8, ONE_BY_SQRT_2, PI, SIN_PI_BY_8, SQRT_2};

/// Groups generators of structured meshes (Qua and Hex -- sometimes Tri)
pub struct Structured {}

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
    /// * `renumber` -- renumbers the points to minimize the bandwidth of the associated graph's matrix
    pub fn quarter_ring_2d(
        rmin: f64,
        rmax: f64,
        nr: usize,
        na: usize,
        target: GeoKind,
        renumber: bool,
    ) -> Result<Mesh, StrError> {
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
        let mut mesh = block.subdivide(target)?;
        if renumber {
            Graph::renumber_mesh(&mut mesh, false)?;
        }
        Ok(mesh)
    }

    /// Generates a mesh representing a quarter of a ring in 3D (extrusion along z)
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
    /// * `renumber` -- renumbers the points to minimize the bandwidth of the associated graph's matrix
    pub fn quarter_ring_3d(
        rmin: f64,
        rmax: f64,
        z: f64,
        nr: usize,
        na: usize,
        nz: usize,
        target: GeoKind,
        renumber: bool,
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
        let mut mesh = block.subdivide(target)?;
        if renumber {
            Graph::renumber_mesh(&mut mesh, false)?;
        }
        Ok(mesh)
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
    /// * `renumber` -- renumbers the points to minimize the bandwidth of the associated graph's matrix
    pub fn quarter_disk_2d_a(r: f64, na: usize, nb: usize, target: GeoKind, renumber: bool) -> Result<Mesh, StrError> {
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
        let mut mesh = join_meshes(&[&mesh_1, &mesh_2, &mesh_3])?;
        if renumber {
            Graph::renumber_mesh(&mut mesh, false)?;
        }
        Ok(mesh)
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
    /// * `renumber` -- renumbers the points to minimize the bandwidth of the associated graph's matrix
    pub fn quarter_disk_2d_b(
        a: f64,
        r: f64,
        na: usize,
        nb: usize,
        target: GeoKind,
        renumber: bool,
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
        let mut mesh = join_meshes(&[&mesh_1, &mesh_2, &mesh_3])?;
        if renumber {
            Graph::renumber_mesh(&mut mesh, false)?;
        }
        Ok(mesh)
    }

    /// Generates a mesh representing a quarter of a disk in 3D (extrusion along z) (A-version)
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
    /// * `renumber` -- renumbers the points to minimize the bandwidth of the associated graph's matrix
    pub fn quarter_disk_3d_a(
        r: f64,
        z: f64,
        na: usize,
        nb: usize,
        nz: usize,
        target: GeoKind,
        renumber: bool,
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
        let mut mesh = join_meshes(&[&mesh_1, &mesh_2, &mesh_3])?;
        if renumber {
            Graph::renumber_mesh(&mut mesh, false)?;
        }
        Ok(mesh)
    }

    /// Generates a mesh representing a quarter of a disk in 3D (extrusion along z) (B-version)
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
    /// * `renumber` -- renumbers the points to minimize the bandwidth of the associated graph's matrix
    pub fn quarter_disk_3d_b(
        a: f64,
        r: f64,
        z: f64,
        na: usize,
        nb: usize,
        nz: usize,
        target: GeoKind,
        renumber: bool,
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
        let mut mesh = join_meshes(&[&mesh_1, &mesh_2, &mesh_3])?;
        if renumber {
            Graph::renumber_mesh(&mut mesh, false)?;
        }
        Ok(mesh)
    }

    /// Generates a mesh representing a quarter of a plate with a hole in 2D
    ///
    /// ```text
    /// y ^
    ///   ---------------------------------
    ///   |                             .'|
    ///   |                           .'  |
    ///   |                         .'    |
    ///   |                       .'      |
    ///   #**=---__             .'        |
    ///   |        '*._       .'          |
    ///   |            *._  .'            |
    ///   |               .'              |
    ///   |             .'  *.            |
    ///   #**=-__     .'      *           |
    ///          '*..'         *          |
    ///   .         \           *         |
    ///              '           *        |
    ///   .           |          *        |
    ///               |          |        |
    ///   o   -   -   ---------------------  --> x
    ///   |<----r---->|<----a--->|<--b--->|
    ///   |<----------------l------------>|
    /// ```
    ///
    /// # Input
    ///
    /// * `r` -- the radius of the hole
    /// * `a` -- distance to make a ring around the hole
    /// * `b` -- the difference `l-(r+a)` with 'l' being the length of the square plate
    /// * `na` -- number of divisions along 'a' (must be > 0)
    /// * `nb` -- number of divisions along 'b' (must be > 0)
    /// * `n45` -- number of divisions along the 45 degrees angle
    /// * `target` -- [crate::shapes::GeoClass::Qua] shapes only
    /// * `renumber` -- renumbers the points to minimize the bandwidth of the associated graph's matrix
    pub fn quarter_plate_hole_2d(
        r: f64,
        a: f64,
        b: f64,
        na: usize,
        nb: usize,
        n45: usize,
        target: GeoKind,
        renumber: bool,
    ) -> Result<Mesh, StrError> {
        const COS_PI_BY_4: f64 = ONE_BY_SQRT_2;
        let ra = r + a;
        let c1 = r * SIN_PI_BY_8;
        let c2 = r * COS_PI_BY_4;
        let c3 = r * COS_PI_BY_8;
        let d1 = ra * SIN_PI_BY_8;
        let d2 = ra * COS_PI_BY_4;
        let d3 = ra * COS_PI_BY_8;
        let m = r + a / 2.0;
        let n = ra + b / 2.0;
        let l = ra + b;
        let hl = l / 2.0;
        let e = (c2 + d2) / 2.0;
        let f = (d2 + l) / 2.0;
        #[rustfmt::skip]
        let mut block_1 = Block::new(&[
            [ r, 0.0],
            [ra, 0.0],
            [d2,  d2],
            [c2,  c2],
            [ m, 0.0],
            [d3,  d1],
            [ e,   e],
            [c3,  c1],
        ])?;
        #[rustfmt::skip]
        let mut block_2 = Block::new(&[
            [ c2, c2],
            [ d2, d2],
            [0.0, ra],
            [0.0,  r],
            [  e,  e],
            [ d1, d3],
            [0.0,  m],
            [ c1, c3],
        ])?;
        #[rustfmt::skip]
        let mut block_3 = Block::new(&[
            [ ra, 0.0],
            [  l, 0.0],
            [  l,   l],
            [ d2,  d2],
            [  n, 0.0],
            [  l,  hl],
            [  f,   f],
            [ d3,  d1],
        ])?;
        #[rustfmt::skip]
        let mut block_4 = Block::new(&[
            [ d2, d2],
            [  l,  l],
            [0.0,  l],
            [0.0, ra],
            [  f,  f],
            [ hl,  l],
            [0.0,  n],
            [ d1, d3],
        ])?;
        let ct_r = Constraint2D::Circle(0.0, 0.0, r);
        let ct_ra = Constraint2D::Circle(0.0, 0.0, ra);
        block_1.set_edge_constraint(3, Some(ct_r.clone()))?;
        block_1.set_edge_constraint(1, Some(ct_ra.clone()))?;
        block_1.set_ndiv(&[na, n45])?;
        block_2.set_edge_constraint(3, Some(ct_r))?;
        block_2.set_edge_constraint(1, Some(ct_ra.clone()))?;
        block_2.set_ndiv(&[na, n45])?;
        block_3.set_edge_constraint(3, Some(ct_ra.clone()))?;
        block_3.set_ndiv(&[nb, n45])?;
        block_4.set_edge_constraint(3, Some(ct_ra.clone()))?;
        block_4.set_ndiv(&[nb, n45])?;
        let mesh_1 = block_1.subdivide(target)?;
        let mesh_2 = block_2.subdivide(target)?;
        let mesh_3 = block_3.subdivide(target)?;
        let mesh_4 = block_4.subdivide(target)?;
        let mut mesh = join_meshes(&[&mesh_1, &mesh_2, &mesh_3, &mesh_4])?;
        if renumber {
            Graph::renumber_mesh(&mut mesh, false)?;
        }
        Ok(mesh)
    }

    /// Generates a mesh representing a quarter of a plate with a hole in 3D (extrusion along z)
    ///
    /// ```text
    /// y ^
    ///   ---------------------------------
    ///   |                             .'|
    ///   |                           .'  |
    ///   |                         .'    |
    ///   |                       .'      |
    ///   #**=---__             .'        |
    ///   |        '*._       .'          |
    ///   |            *._  .'            |
    ///   |               .'              |
    ///   |             .'  *.            |
    ///   #**=-__     .'      *           |
    ///          '*..'         *          |
    ///   .         \           *         |
    ///              '           *        |
    ///   .           |          *        |
    ///               |          |        |
    /// z o   -   -   ---------------------  --> x
    ///   |<----r---->|<----a--->|<--b--->|
    ///   |<----------------l------------>|
    /// ```
    ///
    /// # Input
    ///
    /// * `r` -- the radius of the hole
    /// * `a` -- distance to make a ring around the hole
    /// * `b` -- the difference `l-(r+a)` with 'l' being the length of the square plate
    /// * `z` -- thickness (zmin = 0)
    /// * `na` -- number of divisions along 'a' (must be > 0)
    /// * `nb` -- number of divisions along 'b' (must be > 0)
    /// * `n45` -- number of divisions along the 45 degrees angle
    /// * `nz` -- number of divisions along 'z' (must be > 0)
    /// * `target` -- [crate::shapes::GeoClass::Hex] shapes only
    /// * `renumber` -- renumbers the points to minimize the bandwidth of the associated graph's matrix
    pub fn quarter_plate_hole_3d(
        r: f64,
        a: f64,
        b: f64,
        z: f64,
        na: usize,
        nb: usize,
        n45: usize,
        nz: usize,
        target: GeoKind,
        renumber: bool,
    ) -> Result<Mesh, StrError> {
        const COS_PI_BY_4: f64 = ONE_BY_SQRT_2;
        let ra = r + a;
        let c1 = r * SIN_PI_BY_8;
        let c2 = r * COS_PI_BY_4;
        let c3 = r * COS_PI_BY_8;
        let d1 = ra * SIN_PI_BY_8;
        let d2 = ra * COS_PI_BY_4;
        let d3 = ra * COS_PI_BY_8;
        let m = r + a / 2.0;
        let n = ra + b / 2.0;
        let l = ra + b;
        let hl = l / 2.0;
        let e = (c2 + d2) / 2.0;
        let f = (d2 + l) / 2.0;
        let hz = z / 2.0;
        #[rustfmt::skip]
        let mut block_1 = Block::new(&[
            [ r, 0.0, 0.0],
            [ra, 0.0, 0.0],
            [d2,  d2, 0.0],
            [c2,  c2, 0.0],
            [ r, 0.0,   z],
            [ra, 0.0,   z],
            [d2,  d2,   z],
            [c2,  c2,   z],
            [ m, 0.0, 0.0],
            [d3,  d1, 0.0],
            [ e,   e, 0.0],
            [c3,  c1, 0.0],
            [ m, 0.0,   z],
            [d3,  d1,   z],
            [ e,   e,   z],
            [c3,  c1,   z],
            [ r, 0.0,  hz],
            [ra, 0.0,  hz],
            [d2,  d2,  hz],
            [c2,  c2,  hz],
        ])?;
        #[rustfmt::skip]
        let mut block_2 = Block::new(&[
            [ c2, c2, 0.0],
            [ d2, d2, 0.0],
            [0.0, ra, 0.0],
            [0.0,  r, 0.0],
            [ c2, c2,   z],
            [ d2, d2,   z],
            [0.0, ra,   z],
            [0.0,  r,   z],
            [  e,  e, 0.0],
            [ d1, d3, 0.0],
            [0.0,  m, 0.0],
            [ c1, c3, 0.0],
            [  e,  e,   z],
            [ d1, d3,   z],
            [0.0,  m,   z],
            [ c1, c3,   z],
            [ c2, c2,  hz],
            [ d2, d2,  hz],
            [0.0, ra,  hz],
            [0.0,  r,  hz],
        ])?;
        #[rustfmt::skip]
        let mut block_3 = Block::new(&[
            [ ra, 0.0, 0.0],
            [  l, 0.0, 0.0],
            [  l,   l, 0.0],
            [ d2,  d2, 0.0],
            [ ra, 0.0,   z],
            [  l, 0.0,   z],
            [  l,   l,   z],
            [ d2,  d2,   z],
            [  n, 0.0, 0.0],
            [  l,  hl, 0.0],
            [  f,   f, 0.0],
            [ d3,  d1, 0.0],
            [  n, 0.0,   z],
            [  l,  hl,   z],
            [  f,   f,   z],
            [ d3,  d1,   z],
            [ ra, 0.0,  hz],
            [  l, 0.0,  hz],
            [  l,   l,  hz],
            [ d2,  d2,  hz],
        ])?;
        #[rustfmt::skip]
        let mut block_4 = Block::new(&[
            [ d2, d2, 0.0],
            [  l,  l, 0.0],
            [0.0,  l, 0.0],
            [0.0, ra, 0.0],
            [ d2, d2,   z],
            [  l,  l,   z],
            [0.0,  l,   z],
            [0.0, ra,   z],
            [  f,  f, 0.0],
            [ hl,  l, 0.0],
            [0.0,  n, 0.0],
            [ d1, d3, 0.0],
            [  f,  f,   z],
            [ hl,  l,   z],
            [0.0,  n,   z],
            [ d1, d3,   z],
            [ d2, d2,  hz],
            [  l,  l,  hz],
            [0.0,  l,  hz],
            [0.0, ra,  hz],
        ])?;
        let ct_r = Constraint3D::CylinderZ(0.0, 0.0, r);
        let ct_ra = Constraint3D::CylinderZ(0.0, 0.0, ra);
        block_1.set_face_constraint(0, Some(ct_r.clone()))?;
        block_1.set_face_constraint(1, Some(ct_ra.clone()))?;
        block_1.set_ndiv(&[na, n45, nz])?;
        block_2.set_face_constraint(0, Some(ct_r))?;
        block_2.set_face_constraint(1, Some(ct_ra.clone()))?;
        block_2.set_ndiv(&[na, n45, nz])?;
        block_3.set_face_constraint(0, Some(ct_ra.clone()))?;
        block_3.set_ndiv(&[nb, n45, nz])?;
        block_4.set_face_constraint(0, Some(ct_ra.clone()))?;
        block_4.set_ndiv(&[nb, n45, nz])?;
        let mesh_1 = block_1.subdivide(target)?;
        let mesh_2 = block_2.subdivide(target)?;
        let mesh_3 = block_3.subdivide(target)?;
        let mesh_4 = block_4.subdivide(target)?;
        let mut mesh = join_meshes(&[&mesh_1, &mesh_2, &mesh_3, &mesh_4])?;
        if renumber {
            Graph::renumber_mesh(&mut mesh, false)?;
        }
        Ok(mesh)
    }

    /// Generates a rectangle with horizontal layers
    ///
    /// ```text
    /// n is the number of horizontal levels
    /// n_layer = n - 1
    ///
    ///         xa        xb             xc
    /// y[n-1]  o---------o---------------o
    ///         |         |               |
    ///         |         |               | << ny[n-2]
    ///         |         |               |
    /// y[n-2]  o---------o---------------o
    ///         |         |               |
    ///                . . . . . . .        << ny[1]
    ///         |         |               |
    /// y[1]    o---------o---------------o
    ///         |         |               | << ny[0]
    /// y[0]    o---------o---------------o
    ///         xa        xb             xc
    ///                   na           nb
    /// ```
    ///
    /// # Input
    ///
    /// * `xa` -- minimum x-coordinate
    /// * `xb` -- (optional) a transitional x-coordinate (e.g., to accommodate a "footing" width)
    /// * `xc` -- maximum x-coordinate
    /// * `na` -- If `xb` is Some value, `na` is the number of spacings between `xa` and `xb`.
    ///   Otherwise, if `xb` is None, `na` is the number of spacings between `xa` and `xc`.
    ///   `na` must be greater than 0.
    /// * `nb` -- If `xb` is Some value, `nb` is the number of spacings between `xb` and `xc`.
    ///   Otherwise, if `xb` is None, `nb` is ignored.
    ///   `nb` must be greater than 0.
    /// * `y` -- is a list of increasing y-coordinates, including the first and last y-coordinate
    ///   The length of `y` must be greater than or equal to 2.
    /// * `ny` -- is the number of spacings between `y[i]` and `y[i+1]`. The length of `ny`, which
    ///   is equal to the number of layers, must be equal to the length of `y` minus one, i.e.,
    ///   `n_layer = ny.len() = y.len() - 1`
    /// * `attributes` -- is a list of attributes for each layer; it's length is equal to `n_layer`
    /// * `target` -- is the resulting GeoKind and must be a [crate::shapes::GeoClass::Qua]
    /// * `renumber` -- renumbers the points to minimize the bandwidth of the associated graph's matrix
    pub fn rectangle(
        xa: f64,
        xb: Option<f64>,
        xc: f64,
        na: usize,
        nb: usize,
        y: &[f64],
        ny: &[usize],
        attributes: &[usize],
        target: GeoKind,
        renumber: bool,
    ) -> Result<Mesh, StrError> {
        if xc <= xa {
            return Err("xc must be > xa");
        }
        if let Some(xxb) = xb {
            if xxb <= xa || xxb >= xc {
                return Err("xb must satisfy: xa < xb < xc");
            }
        }
        if na < 1 || nb < 1 {
            return Err("na and nb must be > 0");
        }
        if y.len() < 2 {
            return Err("y.len() must be ≥ 2");
        }
        let n_layer = y.len() - 1;
        if ny.len() != n_layer {
            return Err("ny.len() must be equal to n_layer = y.len() - 1");
        }
        if attributes.len() != n_layer {
            return Err("attributes.len() must be equal to n_layer = y.len() - 1");
        }
        let mut meshes = Vec::new();
        let mut ya = y[0];
        if let Some(xxb) = xb {
            for l in 0..n_layer {
                let yb = y[l + 1];
                if yb <= ya {
                    return Err("y values must be sorted ascending");
                }
                let mut ba = Block::new(&[[xa, ya], [xxb, ya], [xxb, yb], [xa, yb]]).unwrap();
                let mut bb = Block::new(&[[xxb, ya], [xc, ya], [xc, yb], [xxb, yb]]).unwrap();
                ba.set_ndiv(&[na, ny[l]]).unwrap();
                bb.set_ndiv(&[nb, ny[l]]).unwrap();
                ba.set_attribute(attributes[l]);
                bb.set_attribute(attributes[l]);
                meshes.push(ba.subdivide(target)?);
                meshes.push(bb.subdivide(target).unwrap());
                ya = yb;
            }
        } else {
            for l in 0..n_layer {
                let yb = y[l + 1];
                if yb <= ya {
                    return Err("y values must be sorted ascending");
                }
                let mut block = Block::new(&[[xa, ya], [xc, ya], [xc, yb], [xa, yb]]).unwrap();
                block.set_ndiv(&[na, ny[l]]).unwrap();
                block.set_attribute(attributes[l]);
                let mesh = block.subdivide(target)?;
                if n_layer == 1 {
                    return Ok(mesh);
                }
                meshes.push(mesh);
                ya = yb;
            }
        }
        let mut mesh = join_meshes(&meshes.iter().collect::<Vec<_>>())?;
        if renumber {
            Graph::renumber_mesh(&mut mesh, false)?;
        }
        Ok(mesh)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Structured;
    use crate::geometry::point_point_distance;
    use crate::mesh::{Figure, Graph, Mesh};
    use crate::shapes::GeoKind;
    use russell_lab::{approx_eq, array_approx_eq};

    const SAVE_FIGURE: bool = false;
    const MAX_NPOINT_PRINT: usize = 200;

    fn print_bandwidth(mesh: &mut Mesh) {
        let graph = Graph::new(&mesh, false).unwrap();
        if mesh.points.len() < MAX_NPOINT_PRINT {
            graph.print_non_zero_pattern();
        }
        Graph::renumber_mesh(mesh, false).unwrap();
        let graph_after = Graph::new(&mesh, false).unwrap();
        if mesh.points.len() < MAX_NPOINT_PRINT {
            graph_after.print_non_zero_pattern();
        }
        println!("bandwidth (before) = {}", graph.calc_bandwidth());
        println!("bandwidth (after)  = {}", graph_after.calc_bandwidth());
    }

    fn draw(mesh: &Mesh, larger: bool, filename: &str) {
        let mut fig = Figure::new();
        fig.cell_ids = true;
        fig.point_ids = true;
        if larger {
            fig.figure_size = Some((600.0, 600.0));
        }
        mesh.draw(Some(fig), filename, |_, _| {}).unwrap();
    }

    #[test]
    fn quarter_ring_2d_captures_errors() {
        assert_eq!(
            Structured::quarter_ring_2d(3.0, 6.0, 0, 1, GeoKind::Qua16, false).err(),
            Some("ndiv must be ≥ 1")
        );
    }

    #[test]
    fn quarter_ring_2d_works() {
        let mesh = Structured::quarter_ring_2d(3.0, 6.0, 1, 1, GeoKind::Qua16, false).unwrap();
        mesh.check_overlapping_points(0.02).unwrap();
        assert_eq!(mesh.points.len(), 16);
        assert_eq!(mesh.cells.len(), 1);
        for p in [0, 11, 7, 3] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 3.0, 1e-15);
        }
        for p in [4, 12, 15, 10] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 4.0, 1e-15);
        }
        for p in [8, 13, 14, 6] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 5.0, 1e-15);
        }
        for p in [1, 5, 9, 2] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_quarter_ring_2d.svg");
        }
    }

    #[test]
    fn quarter_ring_3d_captures_errors() {
        assert_eq!(
            Structured::quarter_ring_3d(3.0, 6.0, 2.0, 0, 1, 1, GeoKind::Hex8, false).err(),
            Some("ndiv must be ≥ 1")
        );
    }

    #[test]
    fn quarter_ring_3d_works() {
        let mut mesh = Structured::quarter_ring_3d(3.0, 6.0, 2.0, 1, 2, 1, GeoKind::Hex32, false).unwrap();
        mesh.check_overlapping_points(0.02).unwrap();
        assert_eq!(mesh.points.len(), 52);
        assert_eq!(mesh.cells.len(), 2);
        for p in [0, 15, 14, 3, 41, 40, 33] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, 3.0, 1e-15);
        }
        for p in [1, 10, 11, 2, 36, 37, 32] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
        if SAVE_FIGURE {
            draw(&mesh, true, "/tmp/gemlab/test_quarter_ring_3d.svg");
        }
        if false {
            print_bandwidth(&mut mesh);
        }
    }

    #[test]
    fn quarter_disk_2d_a_works_qua8() {
        let mesh = Structured::quarter_disk_2d_a(6.0, 1, 1, GeoKind::Qua8, false).unwrap();
        mesh.check_overlapping_points(0.02).unwrap();
        assert_eq!(mesh.points.len(), 16);
        assert_eq!(mesh.cells.len(), 3);
        for p in [8, 11, 9, 14, 13] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_quarter_disk_2d_a_qua8.svg");
        }
    }

    #[test]
    fn quarter_disk_2d_a_works_qua8_finer() {
        let mesh = Structured::quarter_disk_2d_a(6.0, 3, 3, GeoKind::Qua8, false).unwrap();
        mesh.check_overlapping_points(0.02).unwrap();
        assert_eq!(mesh.points.len(), 100);
        assert_eq!(mesh.cells.len(), 27);
        for p in [50, 53, 51, 62, 61, 71, 70, 99, 96, 98, 91, 94, 92] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_quarter_disk_2d_a_qua8_finer.svg");
        }
    }

    #[test]
    fn quarter_disk_2d_a_works_qua16() {
        let mesh = Structured::quarter_disk_2d_a(6.0, 1, 1, GeoKind::Qua16, false).unwrap();
        mesh.check_overlapping_points(0.02).unwrap();
        assert_eq!(mesh.points.len(), 37);
        assert_eq!(mesh.cells.len(), 3);
        for p in [16, 19, 22, 17, 29, 31, 28] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_quarter_disk_2d_a_qua16.svg");
        }
    }

    #[test]
    fn quarter_disk_2d_b_captures_errors() {
        assert_eq!(
            Structured::quarter_disk_2d_b(6.1, 6.0, 1, 1, GeoKind::Qua8, false).err(),
            Some("'a' must be smaller than 'r'")
        );
    }

    #[test]
    fn quarter_disk_2d_b_works_qua8() {
        let mesh = Structured::quarter_disk_2d_b(3.0, 6.0, 1, 1, GeoKind::Qua8, false).unwrap();
        mesh.check_overlapping_points(0.02).unwrap();
        assert_eq!(mesh.points.len(), 16);
        assert_eq!(mesh.cells.len(), 3);
        for p in [8, 11, 9, 14, 13] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_quarter_disk_2d_b_qua8.svg");
        }
    }

    #[test]
    fn quarter_disk_2d_b_works_qua8_finer() {
        let mesh = Structured::quarter_disk_2d_b(3.0, 6.0, 3, 3, GeoKind::Qua8, false).unwrap();
        mesh.check_overlapping_points(0.02).unwrap();
        assert_eq!(mesh.points.len(), 100);
        assert_eq!(mesh.cells.len(), 27);
        for p in [50, 53, 51, 62, 61, 71, 70, 99, 96, 98, 91, 94, 92] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
        if SAVE_FIGURE {
            draw(&mesh, true, "/tmp/gemlab/test_quarter_disk_2d_b_qua8_finer.svg");
        }
    }

    #[test]
    fn quarter_disk_2d_b_works_qua16() {
        let mesh = Structured::quarter_disk_2d_b(3.0, 6.0, 1, 1, GeoKind::Qua16, false).unwrap();
        mesh.check_overlapping_points(0.02).unwrap();
        assert_eq!(mesh.points.len(), 37);
        assert_eq!(mesh.cells.len(), 3);
        for p in [16, 19, 22, 17, 29, 31, 28] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_quarter_disk_2d_b_qua16.svg");
        }
    }

    #[test]
    fn quarter_disk_3d_a_works_hex32() {
        let mesh = Structured::quarter_disk_3d_a(6.0, 1.5, 1, 1, 1, GeoKind::Hex32, false).unwrap();
        mesh.check_overlapping_points(0.02).unwrap();
        assert_eq!(mesh.cells.len(), 3);
        for p in [
            32, 38, 39, 33, 54, 55, 52, // z-min
            48, 50, 62, // z-1/3
            49, 51, 63, // z-2/3
            34, 44, 45, 35, 58, 59, 53, // z-max
        ] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
        if SAVE_FIGURE {
            draw(&mesh, true, "/tmp/gemlab/test_quarter_disk_3d_a_hex32.svg");
        }
    }

    #[test]
    fn quarter_disk_3d_b_captures_errors() {
        assert_eq!(
            Structured::quarter_disk_3d_b(6.1, 6.0, 1.5, 1, 1, 1, GeoKind::Hex32, false).err(),
            Some("'a' must be smaller than 'r'")
        );
    }

    #[test]
    fn quarter_disk_3d_b_works_hex32() {
        let mesh = Structured::quarter_disk_3d_b(3.0, 6.0, 1.5, 1, 1, 1, GeoKind::Hex32, false).unwrap();
        mesh.check_overlapping_points(0.02).unwrap();
        assert_eq!(mesh.cells.len(), 3);
        for p in [
            32, 38, 39, 33, 54, 55, 52, // z-min
            48, 50, 62, // z-1/3
            49, 51, 63, // z-2/3
            34, 44, 45, 35, 58, 59, 53, // z-max
        ] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
        if SAVE_FIGURE {
            draw(&mesh, true, "/tmp/gemlab/test_quarter_disk_3d_b_hex32.svg");
        }
    }

    #[test]
    fn quarter_plate_hole_2d_works() {
        let mesh = Structured::quarter_plate_hole_2d(1.0, 1.0, 1.0, 1, 1, 1, GeoKind::Qua12, false).unwrap();
        assert_eq!(mesh.points.len(), 33);
        assert_eq!(mesh.cells.len(), 4);
        mesh.check_overlapping_points(0.01).unwrap();
        for p in [0, 11, 7, 3, 19, 16, 13] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 1.0, 1e-15);
        }
        for p in [1, 5, 9, 2, 14, 17, 12] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 2.0, 1e-15);
        }
        assert_eq!(mesh.points[21].coords, &[3.0, 3.0]);
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_quarter_plate_hole_2d.svg");
        }
    }

    #[test]
    fn quarter_plate_hole_3d_works() {
        let mesh = Structured::quarter_plate_hole_3d(1.0, 1.0, 1.0, 1.5, 1, 1, 1, 1, GeoKind::Hex32, false).unwrap();
        assert_eq!(mesh.points.len(), 66 + 18);
        assert_eq!(mesh.cells.len(), 4);
        mesh.check_overlapping_points(0.01).unwrap();
        // z-min
        for p in [0, 1, 2, 32, 33, 52, 53, 72] {
            assert_eq!(mesh.points[p].coords[2], 0.0);
        }
        // z-max
        for p in [4, 5, 6, 7, 34, 35, 54, 55, 73] {
            assert_eq!(mesh.points[p].coords[2], 1.5);
        }
        // hole/inner
        for p in [0, 15, 14, 3, 41, 40, 33] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, 1.0, 1e-15);
        }
        // hole/outer (ring)
        for p in [1, 10, 11, 2, 36, 37, 32] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, 2.0, 1e-15);
        }
        array_approx_eq(&mesh.points[26].coords, &[2.0, 0.0, 0.5], 1e-15);
        array_approx_eq(&mesh.points[49].coords, &[0.0, 2.0, 1.0], 1e-15);
        array_approx_eq(&mesh.points[70].coords, &[3.0, 3.0, 0.5], 1e-15);
        if SAVE_FIGURE {
            draw(&mesh, true, "/tmp/gemlab/test_quarter_plate_hole_3d.svg");
        }
    }

    #[test]
    fn rectangle_handles_errors() {
        assert_eq!(
            Structured::rectangle(0.0, None, -1.0, 1, 1, &[1.0, 2.0], &[1], &[10], GeoKind::Qua4, false).err(),
            Some("xc must be > xa")
        );
        assert_eq!(
            Structured::rectangle(
                0.0,
                Some(0.0),
                1.0,
                1,
                1,
                &[1.0, 2.0],
                &[1],
                &[10],
                GeoKind::Qua4,
                false
            )
            .err(),
            Some("xb must satisfy: xa < xb < xc")
        );
        assert_eq!(
            Structured::rectangle(0.0, None, 1.0, 0, 1, &[1.0, 2.0], &[1], &[10], GeoKind::Qua4, false).err(),
            Some("na and nb must be > 0")
        );
        assert_eq!(
            Structured::rectangle(0.0, None, 1.0, 1, 0, &[1.0, 2.0], &[1], &[10], GeoKind::Qua4, false).err(),
            Some("na and nb must be > 0")
        );
        assert_eq!(
            Structured::rectangle(0.0, None, 1.0, 1, 1, &[1.0], &[1], &[10], GeoKind::Qua4, false).err(),
            Some("y.len() must be ≥ 2")
        );
        assert_eq!(
            Structured::rectangle(0.0, None, 1.0, 1, 1, &[1.0, 2.0], &[], &[10], GeoKind::Qua4, false).err(),
            Some("ny.len() must be equal to n_layer = y.len() - 1")
        );
        assert_eq!(
            Structured::rectangle(0.0, None, 1.0, 1, 1, &[1.0, 2.0], &[1], &[], GeoKind::Qua4, false).err(),
            Some("attributes.len() must be equal to n_layer = y.len() - 1")
        );
        assert_eq!(
            Structured::rectangle(0.0, None, 1.0, 1, 1, &[2.0, 2.0], &[1], &[10], GeoKind::Qua4, false).err(),
            Some("y values must be sorted ascending")
        );
        assert_eq!(
            Structured::rectangle(0.0, None, 1.0, 1, 1, &[1.0, 2.0], &[1], &[10], GeoKind::Tri3, false).err(),
            Some("in 2D, the GeoClass of target must be Qua")
        );
        assert_eq!(
            Structured::rectangle(
                0.0,
                Some(0.5),
                1.0,
                1,
                1,
                &[2.0, 2.0],
                &[1],
                &[10],
                GeoKind::Qua4,
                false
            )
            .err(),
            Some("y values must be sorted ascending")
        );
        assert_eq!(
            Structured::rectangle(
                0.0,
                Some(0.5),
                1.0,
                1,
                1,
                &[1.0, 2.0],
                &[1],
                &[10],
                GeoKind::Tri3,
                false
            )
            .err(),
            Some("in 2D, the GeoClass of target must be Qua")
        );
    }

    #[test]
    fn rectangle_works() {
        let (xa, xc, na, nb) = (1.0, 3.0, 1, 1);
        let target = GeoKind::Qua4;

        // one column / one layer = single cell -------------------------------

        let mesh = Structured::rectangle(xa, None, xc, na, nb, &[2.0, 5.0], &[1], &[10], target, false).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_layered_rectangle_1.svg");
        }
        mesh.check_overlapping_points(1e-2).unwrap();
        mesh.check_all().unwrap();
        assert_eq!(
            format!("{}", mesh),
            "# header\n\
             # ndim npoint ncell\n\
             2 4 1\n\
             \n\
             # points\n\
             # id marker x y {z}\n\
             0 0 1.0 2.0\n\
             1 0 3.0 2.0\n\
             2 0 3.0 5.0\n\
             3 0 1.0 5.0\n\
             \n\
             # cells\n\
             # id attribute kind points\n\
             0 10 qua4 0 1 2 3\n"
        );

        // two columns / one layer = two cells -------------------------------

        let mesh = Structured::rectangle(xa, Some(1.5), xc, na, nb, &[2.0, 5.0], &[1], &[20], target, false).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_layered_rectangle_2.svg");
        }
        mesh.check_overlapping_points(1e-2).unwrap();
        mesh.check_all().unwrap();
        assert_eq!(
            format!("{}", mesh),
            "# header\n\
             # ndim npoint ncell\n\
             2 6 2\n\
             \n\
             # points\n\
             # id marker x y {z}\n\
             0 0 1.0 2.0\n\
             1 0 1.5 2.0\n\
             2 0 1.5 5.0\n\
             3 0 1.0 5.0\n\
             4 0 3.0 2.0\n\
             5 0 3.0 5.0\n\
             \n\
             # cells\n\
             # id attribute kind points\n\
             0 20 qua4 0 1 2 3\n\
             1 20 qua4 1 4 5 2\n"
        );

        // one column / two layers = two cells -------------------------------

        let mesh = Structured::rectangle(
            xa,
            None,
            xc,
            na,
            nb,
            &[2.0, 3.0, 5.0],
            &[1, 1],
            &[10, 20],
            target,
            false,
        )
        .unwrap();
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_layered_rectangle_3.svg");
        }
        mesh.check_overlapping_points(1e-2).unwrap();
        mesh.check_all().unwrap();
        assert_eq!(
            format!("{}", mesh),
            "# header\n\
             # ndim npoint ncell\n\
             2 6 2\n\
             \n\
             # points\n\
             # id marker x y {z}\n\
             0 0 1.0 2.0\n\
             1 0 3.0 2.0\n\
             2 0 3.0 3.0\n\
             3 0 1.0 3.0\n\
             4 0 3.0 5.0\n\
             5 0 1.0 5.0\n\
             \n\
             # cells\n\
             # id attribute kind points\n\
             0 10 qua4 0 1 2 3\n\
             1 20 qua4 3 2 4 5\n"
        );

        // two columns / two layers = four cells -------------------------------

        let mut mesh = Structured::rectangle(
            xa,
            Some(1.5),
            xc,
            na,
            nb,
            &[2.0, 3.0, 5.0],
            &[1, 1],
            &[10, 20],
            target,
            false,
        )
        .unwrap();
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_layered_rectangle_4.svg");
        }
        mesh.check_overlapping_points(1e-2).unwrap();
        mesh.check_all().unwrap();
        assert_eq!(
            format!("{}", mesh),
            "# header\n\
             # ndim npoint ncell\n\
             2 9 4\n\
             \n\
             # points\n\
             # id marker x y {z}\n\
             0 0 1.0 2.0\n\
             1 0 1.5 2.0\n\
             2 0 1.5 3.0\n\
             3 0 1.0 3.0\n\
             4 0 3.0 2.0\n\
             5 0 3.0 3.0\n\
             6 0 1.5 5.0\n\
             7 0 1.0 5.0\n\
             8 0 3.0 5.0\n\
             \n\
             # cells\n\
             # id attribute kind points\n\
             0 10 qua4 0 1 2 3\n\
             1 10 qua4 1 4 5 2\n\
             2 20 qua4 3 2 6 7\n\
             3 20 qua4 2 5 8 6\n"
        );
        if false {
            print_bandwidth(&mut mesh);
        }
    }
}
