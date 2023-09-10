use super::{Cell, Mesh, Point};
use crate::shapes::GeoKind;
use crate::StrError;
use russell_lab::math::PI;
use tritet::Trigen;

/// Groups generators of unstructured meshes (Tri and Tet only)
pub struct Unstructured {}

impl Unstructured {
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
    /// * `o2` -- generate middle nodes
    pub fn quarter_ring_2d(
        rmin: f64,
        rmax: f64,
        nr: usize,
        na: usize,
        o2: bool,
        global_max_area: Option<f64>,
    ) -> Result<Mesh, StrError> {
        // check
        if nr < 1 {
            return Err("number of divisions along the radius must be > 0");
        }
        if na < 1 {
            return Err("number of divisions along alpha must be > 0");
        }

        // allocate data
        let npoint = 2 * (nr + 1) + 2 * (na - 1);
        let nsegment = npoint;
        let nregion = 1;
        let mut trigen = Trigen::new(npoint, Some(nsegment), Some(nregion), None)?;

        // constants
        const AMIN: f64 = 0.0;
        const AMAX: f64 = PI / 2.0;
        let dr = (rmax - rmin) / (nr as f64);
        let da = (AMAX - AMIN) / (na as f64);

        // counters
        let mut index = 0;

        // horizontal line
        trigen.set_point(index, rmin, 0.0)?;
        index += 1;
        for i in 1..(nr + 1) {
            let r = rmin + (i as f64) * dr;
            trigen.set_point(index, r, 0.0)?;
            trigen.set_segment(index - 1, index - 1, index)?;
            index += 1;
        }

        // outer circle
        for i in 1..(na + 1) {
            let a = AMIN + (i as f64) * da;
            let x = rmax * f64::cos(a);
            let y = rmax * f64::sin(a);
            trigen.set_point(index, x, y)?;
            trigen.set_segment(index - 1, index - 1, index)?;
            index += 1;
        }

        // vertical line
        for i in 1..(nr + 1) {
            let r = rmin + ((nr - i) as f64) * dr;
            trigen.set_point(index, 0.0, r)?;
            trigen.set_segment(index - 1, index - 1, index)?;
            index += 1;
        }

        // inner circle
        for i in 1..na {
            let a = AMIN + ((na - i) as f64) * da;
            let x = rmin * f64::cos(a);
            let y = rmin * f64::sin(a);
            trigen.set_point(index, x, y)?;
            trigen.set_segment(index - 1, index - 1, index)?;
            index += 1;
        }
        trigen.set_segment(index - 1, index - 1, 0)?;

        // region
        trigen.set_region(0, rmin + 1e-4, 1e-4, 1, None)?;

        // generate mesh
        trigen.generate_mesh(false, o2, false, global_max_area, None)?;

        // allocate data
        const NDIM: usize = 2;
        let npoint = trigen.npoint();
        let ncell = trigen.ntriangle();
        let nnode = trigen.nnode();
        let kind = if o2 { GeoKind::Tri6 } else { GeoKind::Tri3 };
        let zero_point = Point {
            id: 0,
            coords: vec![0.0; NDIM],
        };
        let zero_cell = Cell {
            id: 0,
            attribute_id: 1,
            kind,
            points: vec![0; nnode],
        };
        let mut mesh = Mesh {
            ndim: NDIM,
            points: vec![zero_point; npoint],
            cells: vec![zero_cell; ncell],
        };

        // set mesh data
        for i in 0..npoint {
            mesh.points[i].id = i;
            mesh.points[i].coords[0] = trigen.point(i, 0);
            mesh.points[i].coords[1] = trigen.point(i, 1);
        }
        for i in 0..ncell {
            mesh.cells[i].id = i;
            mesh.cells[i].attribute_id = trigen.triangle_attribute(i);
            for m in 0..nnode {
                mesh.cells[i].points[m] = trigen.triangle_node(i, m);
            }
        }

        // results
        Ok(mesh)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Unstructured;
    use crate::geometry::point_point_distance;
    use crate::mesh::{check_all, check_overlapping_points};
    use russell_chk::approx_eq;

    #[allow(unused_imports)]
    use crate::mesh::draw_mesh;

    #[test]
    fn tri_quarter_ring_2d_captures_errors() {
        assert_eq!(
            Unstructured::quarter_ring_2d(3.0, 6.0, 0, 4, false, None).err(),
            Some("number of divisions along the radius must be > 0")
        );
        assert_eq!(
            Unstructured::quarter_ring_2d(3.0, 6.0, 2, 0, false, None).err(),
            Some("number of divisions along alpha must be > 0")
        );
    }

    #[test]
    fn tri_quarter_ring_2d_works() {
        let mesh = Unstructured::quarter_ring_2d(3.0, 6.0, 2, 4, false, None).unwrap();
        // draw_mesh(&mesh, false, true, false, "/tmp/gemlab/test_tri_quarter_ring_2d.svg").unwrap();
        assert_eq!(mesh.points.len(), 14);
        assert_eq!(mesh.cells.len(), 14);
        check_all(&mesh).unwrap();
        check_overlapping_points(&mesh, 0.18).unwrap();
        for p in [0, 11, 10, 9, 8] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, 3.0, 1e-15);
        }
        for p in [2, 3, 4, 5, 6] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
    }

    #[test]
    fn tri_quarter_ring_2d_global_max_area_works() {
        let global_max_area = 0.4;
        let mesh = Unstructured::quarter_ring_2d(3.0, 6.0, 2, 4, false, Some(global_max_area)).unwrap();
        // let svg = "/tmp/gemlab/test_tri_quarter_ring_2d_global_max_area.svg";
        // draw_mesh(&mesh, false, true, false, svg).unwrap();
        assert_eq!(mesh.points.len(), 40);
        assert_eq!(mesh.cells.len(), 66);
        check_all(&mesh).unwrap();
        check_overlapping_points(&mesh, 0.18).unwrap();
        for p in [0, 11, 10, 9, 8] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, 3.0, 1e-15);
        }
        for p in [2, 3, 4, 5, 6] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
    }
}
