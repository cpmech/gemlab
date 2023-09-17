use super::{Cell, Mesh, Point};
use crate::shapes::{GeoClass, GeoKind};
use crate::StrError;
use russell_lab::math::PI;
use tritet::Trigen;

/// Groups generators of unstructured meshes (Tri and Tet only)
pub struct Unstructured {}

fn apply_constraints(mesh: &mut Mesh, rmin: f64, rmax: f64) {
    //  -3 =---__
    //   |        '*._
    //   | -40        *._
    //   |               *.  -20
    //  -4 ==-__           *.
    //          '-.          *
    //             *.         *
    //          -10  *         *
    //                *         *
    //                 *         *
    //                 #   -30   #
    //                -1 ------- -2

    // apply constraints
    const TOL: f64 = 1e-15;
    for i in 0..mesh.points.len() {
        // points on inner circle
        if mesh.points[i].marker == -10 {
            let dx = mesh.points[i].coords[0] - 0.0;
            let dy = mesh.points[i].coords[1] - 0.0;
            let d = f64::sqrt(dx * dx + dy * dy);
            let gap = rmin - d;
            if f64::abs(gap) > TOL {
                let move_x = gap * dx / d;
                let move_y = gap * dy / d;
                mesh.points[i].coords[0] += move_x;
                mesh.points[i].coords[1] += move_y;
            }
        }
        // points on outer circle
        if mesh.points[i].marker == -20 {
            let dx = mesh.points[i].coords[0] - 0.0;
            let dy = mesh.points[i].coords[1] - 0.0;
            let d = f64::sqrt(dx * dx + dy * dy);
            let gap = rmax - d;
            if f64::abs(gap) > TOL {
                let move_x = gap * dx / d;
                let move_y = gap * dy / d;
                mesh.points[i].coords[0] += move_x;
                mesh.points[i].coords[1] += move_y;
            }
        }
    }
}

impl Unstructured {
    /// Generates a mesh representing a quarter of a ring in 2D
    ///
    /// ```text
    /// Geometry:
    ///
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
    /// ```text
    /// Point markers:
    ///
    ///  -3 ==---__
    ///   |        '*._
    ///   | -40        *._
    ///   |               *.  -20
    ///  -4 ==-__           *.
    ///          '-.          *
    ///             *.         *
    ///          -10  *         *
    ///                *         *
    ///                 *         *
    ///                 #   -30   #
    ///                -1 ------- -2
    /// ```
    ///
    /// # Input
    ///
    /// * `rmin` -- inner radius
    /// * `rmax` -- outer radius
    /// * `nr` -- number of divisions along the radius (must be > 0)
    /// * `na` -- number of divisions along alpha (must be > 0)
    /// * `target` -- a triangle GeoKind
    /// * `global_max_area` -- max area allowed for all triangles
    /// * `trigen_vtu_filename` -- an optional filename for the VTU generated by trigen (may be "")
    pub fn quarter_ring_2d(
        rmin: f64,
        rmax: f64,
        nr: usize,
        na: usize,
        target: GeoKind,
        global_max_area: Option<f64>,
    ) -> Result<Mesh, StrError> {
        // check
        if nr < 1 {
            return Err("number of divisions along the radius must be > 0");
        }
        if na < 1 {
            return Err("number of divisions along alpha must be > 0");
        }
        if target.class() != GeoClass::Tri {
            return Err("the GeoClass of target must be Tri");
        }

        // generate o2 triangles (i.e., Tri6) (need to use o2 for others too, e.g.,
        // Tri10, Tri15, because the the middle-edge markers will be replicated by trigen
        let o2 = if target.nnode() > 3 { true } else { false };

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
        trigen.set_point(index, -1, rmin, 0.0)?;
        index += 1;
        for i in 1..(nr + 1) {
            let marker = if i == nr { -2 } else { 0 };
            let r = rmin + (i as f64) * dr;
            trigen.set_point(index, marker, r, 0.0)?;
            trigen.set_segment(index - 1, -30, index - 1, index)?;
            index += 1;
        }

        // outer circle
        for i in 1..(na + 1) {
            let marker = if i == na { -3 } else { 0 };
            let a = AMIN + (i as f64) * da;
            let x = rmax * f64::cos(a);
            let y = rmax * f64::sin(a);
            trigen.set_point(index, marker, x, y)?;
            trigen.set_segment(index - 1, -20, index - 1, index)?;
            index += 1;
        }

        // vertical line
        for i in 1..(nr + 1) {
            let marker = if i == nr { -4 } else { 0 };
            let r = rmin + ((nr - i) as f64) * dr;
            trigen.set_point(index, marker, 0.0, r)?;
            trigen.set_segment(index - 1, -40, index - 1, index)?;
            index += 1;
        }

        // inner circle
        for i in 1..na {
            let a = AMIN + ((na - i) as f64) * da;
            let x = rmin * f64::cos(a);
            let y = rmin * f64::sin(a);
            trigen.set_point(index, 0, x, y)?;
            trigen.set_segment(index - 1, -10, index - 1, index)?;
            index += 1;
        }
        trigen.set_segment(index - 1, -10, index - 1, 0)?;

        // region
        trigen.set_region(0, 1, rmin + 1e-4, 1e-4, None)?;

        // generate mesh
        trigen.generate_mesh(false, o2, true, global_max_area, None)?;

        // allocate data
        const NDIM: usize = 2;
        let npoint = trigen.out_npoint();
        let ncell = trigen.out_ncell();
        let nnode = trigen.out_cell_npoint();
        let kind = if o2 { GeoKind::Tri6 } else { GeoKind::Tri3 };
        let zero_point = Point {
            id: 0,
            marker: 0,
            coords: vec![0.0; NDIM],
        };
        let zero_cell = Cell {
            id: 0,
            attribute: 1,
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
            let marker = trigen.out_point_marker(i);
            mesh.points[i].id = i;
            mesh.points[i].marker = marker;
            mesh.points[i].coords[0] = trigen.out_point(i, 0);
            mesh.points[i].coords[1] = trigen.out_point(i, 1);
        }
        for i in 0..ncell {
            mesh.cells[i].id = i;
            mesh.cells[i].attribute = trigen.out_cell_attribute(i);
            for m in 0..nnode {
                mesh.cells[i].points[m] = trigen.out_cell_point(i, m);
            }
        }

        // apply constraints (need to be done before the upgrade because
        // Steiner points may be added even for Tri3)
        apply_constraints(&mut mesh, rmin, rmax);

        // upgrade mesh (need to apply constraints again because
        // new mid-edge points may be created)
        if target.nnode() > 6 {
            let mut new_mesh = mesh.convert_2d(target)?;
            apply_constraints(&mut new_mesh, rmin, rmax);
            return Ok(new_mesh);
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
    use crate::mesh::{At, Figure, Find, Mesh};
    use crate::shapes::GeoKind;
    use crate::util::any_x;
    use russell_chk::approx_eq;

    const RMIN: f64 = 3.0;
    const RMAX: f64 = 6.0;
    const SAVE_FIGURE: bool = false;

    fn draw(mesh: &Mesh, larger: bool, filename: &str) {
        let mut fig = Figure::new();
        fig.param_cell_ids = true;
        fig.param_point_ids = true;
        if larger {
            fig.param_figure_size = Some((800.0, 800.0));
        } else {
            fig.param_figure_size = Some((600.0, 600.0));
        }
        mesh.draw(Some(fig), filename).unwrap();
    }

    fn check_constraints(mesh: &Mesh) {
        let inner: Vec<_> = mesh.points.iter().filter(|p| p.marker == -10).collect();
        let outer: Vec<_> = mesh.points.iter().filter(|p| p.marker == -20).collect();
        for point in inner {
            let d = point_point_distance(&point.coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for point in outer {
            let d = point_point_distance(&point.coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
    }

    fn check_corner_markers(mesh: &Mesh) {
        let find = Find::new(mesh, None);
        let res = find.point_ids(At::XY(RMIN, 0.0), any_x).unwrap();
        assert_eq!(res.len(), 1);
        assert_eq!(mesh.points[res[0]].marker, -1);

        let res = find.point_ids(At::XY(RMAX, 0.0), any_x).unwrap();
        assert_eq!(res.len(), 1);
        assert_eq!(mesh.points[res[0]].marker, -2);

        let res = find.point_ids(At::XY(0.0, RMAX), any_x).unwrap();
        assert_eq!(res.len(), 1);
        assert_eq!(mesh.points[res[0]].marker, -3);

        let res = find.point_ids(At::XY(0.0, RMIN), any_x).unwrap();
        assert_eq!(res.len(), 1);
        assert_eq!(mesh.points[res[0]].marker, -4);
    }

    fn check_edge_point_markers(mesh: &Mesh, inner: &[usize], outer: &[usize]) {
        for p in inner {
            assert_eq!(mesh.points[*p].marker, -10);
        }
        for p in outer {
            assert_eq!(mesh.points[*p].marker, -20);
        }
    }

    #[test]
    fn tri_quarter_ring_2d_captures_errors() {
        assert_eq!(
            Unstructured::quarter_ring_2d(RMIN, RMAX, 0, 4, GeoKind::Tri3, None).err(),
            Some("number of divisions along the radius must be > 0")
        );
        assert_eq!(
            Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 0, GeoKind::Tri3, None).err(),
            Some("number of divisions along alpha must be > 0")
        );
        assert_eq!(
            Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 2, GeoKind::Qua4, None).err(),
            Some("the GeoClass of target must be Tri")
        );
    }

    #[test]
    fn tri_quarter_ring_2d_works() {
        let mesh = Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 4, GeoKind::Tri3, None).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_tri_quarter_ring_2d.svg");
        }
        assert_eq!(mesh.points.len(), 14);
        assert_eq!(mesh.cells.len(), 14);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.18).unwrap();
        check_constraints(&mesh);
        check_corner_markers(&mesh);
        check_edge_point_markers(&mesh, &[11, 10, 9], &[3, 4, 5]);
        for p in [0, 11, 10, 9, 8] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for p in [2, 3, 4, 5, 6] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
    }

    #[test]
    fn tri_quarter_ring_2d_o2_works() {
        let mesh = Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 4, GeoKind::Tri6, None).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_tri_quarter_ring_2d_o2.svg");
        }
        assert_eq!(mesh.points.len(), 41);
        assert_eq!(mesh.cells.len(), 14);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.18).unwrap();
        check_constraints(&mesh);
        check_corner_markers(&mesh);
        check_edge_point_markers(&mesh, &[32, 11, 36, 10, 21, 9, 22], &[34, 3, 40, 4, 17, 5, 26]);
        for p in [0, 32, 11, 36, 10, 21, 9, 22, 8] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for p in [2, 34, 3, 40, 4, 17, 5, 26, 6] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
    }

    #[test]
    fn tri_quarter_ring_2d_global_max_area_works() {
        let global_max_area = Some(0.4);
        let mesh = Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 4, GeoKind::Tri3, global_max_area).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, true, "/tmp/gemlab/test_tri_quarter_ring_2d_global_max_area.svg");
        }
        assert_eq!(mesh.points.len(), 50);
        assert_eq!(mesh.cells.len(), 78);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.18).unwrap();
        check_constraints(&mesh);
        check_corner_markers(&mesh);
        check_edge_point_markers(&mesh, &[11, 10, 17, 9], &[40, 3, 25, 4, 24, 5, 33]);
        for p in [0, 11, 10, 17, 9, 8] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for p in [2, 40, 3, 25, 4, 24, 5, 33, 6] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
    }

    #[test]
    fn tri_quarter_ring_2d_o2_global_max_area_works() {
        let global_max_area = Some(0.4);
        let mesh = Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 4, GeoKind::Tri6, global_max_area).unwrap();
        if SAVE_FIGURE {
            draw(
                &mesh,
                true,
                "/tmp/gemlab/test_tri_quarter_ring_2d_o2_global_max_area.svg",
            );
        }
        assert_eq!(mesh.points.len(), 177);
        assert_eq!(mesh.cells.len(), 78);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.1).unwrap();
        check_constraints(&mesh);
        check_corner_markers(&mesh);
        check_edge_point_markers(
            &mesh,
            &[150, 11, 81, 10, 58, 17, 99, 9, 146],
            &[154, 40, 176, 3, 174, 25, 130, 4, 118, 24, 136, 5, 124, 33, 134],
        );
        for p in [0, 150, 11, 81, 10, 58, 17, 99, 9, 146, 8] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for p in [2, 154, 40, 176, 3, 174, 25, 130, 4, 118, 24, 136, 5, 124, 33, 134, 6] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
    }

    #[test]
    fn tri_quarter_ring_2d_tri10_works() {
        let mesh = Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 4, GeoKind::Tri10, None).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_tri_quarter_ring_2d_tri10.svg");
        }
        assert_eq!(mesh.points.len(), 82);
        assert_eq!(mesh.cells.len(), 14);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.1).unwrap();
        check_constraints(&mesh);
        check_corner_markers(&mesh);
        check_edge_point_markers(
            &mesh,
            &[57, 54, 51, 70, 68, 2, 23, 21, 19, 30, 27],
            &[62, 65, 60, 79, 80, 10, 12, 15, 11, 36, 39],
        );
        for p in [49, 57, 54, 51, 70, 68, 2, 23, 21, 19, 30, 27, 25] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for p in [59, 62, 65, 60, 79, 80, 10, 12, 15, 11, 36, 39, 34] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
    }

    #[test]
    fn tri_quarter_ring_2d_tri15_works() {
        let mesh = Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 4, GeoKind::Tri15, None).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, true, "/tmp/gemlab/test_tri_quarter_ring_2d_tri15.svg");
        }
        assert_eq!(mesh.points.len(), 137);
        assert_eq!(mesh.cells.len(), 14);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.1).unwrap();
        check_constraints(&mesh);
        check_corner_markers(&mesh);
        check_edge_point_markers(
            &mesh,
            &[92, 86, 91, 83, 113, 110, 112, 2, 35, 31, 34, 29, 45, 41, 44],
            &[103, 99, 104, 97, 132, 131, 133, 15, 20, 17, 21, 16, 59, 55, 60],
        );
        for p in [81, 92, 86, 91, 83, 113, 110, 112, 2, 35, 31, 34, 29, 45, 41, 44, 39] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for p in [96, 103, 99, 104, 97, 132, 131, 133, 15, 20, 17, 21, 16, 59, 55, 60, 53] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
    }
}
