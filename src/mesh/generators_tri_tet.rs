#![allow(unused)]

use super::{Cell, Mesh, Point};
use crate::shapes::{GeoClass, GeoKind};
use crate::StrError;
use russell_lab::math::PI;
use tritet::{Tetgen, Trigen};

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
    /// * `nr` -- number of divisions along the radius (must be ≥ 1)
    /// * `na` -- number of divisions along alpha (must be ≥ 1)
    /// * `target` -- [crate::shapes::GeoClass::Tri] shapes only
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
            return Err("number of divisions along the radius must be ≥ 1");
        }
        if na < 1 {
            return Err("number of divisions along alpha must be ≥ 1");
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

    /// Generates a mesh representing a quarter of a ring in 3D
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
    /// * `thickness` -- thickness (zmin = 0.0)
    /// * `nr` -- number of divisions along the radius (must be ≥ 1)
    /// * `na` -- number of divisions along alpha (must be ≥ 1)
    /// * `nz` -- number of divisions along z (thickness) (must be ≥ 1)
    /// * `target` -- [crate::shapes::GeoClass::Tet] shapes only
    /// * `global_max_volume` -- max volume allowed for all tetrahedra
    /// * `trigen_vtu_filename` -- an optional filename for the VTU generated by trigen (may be "")
    pub fn quarter_ring_3d(
        rmin: f64,
        rmax: f64,
        thickness: f64,
        nr: usize,
        na: usize,
        target: GeoKind,
        global_max_volume: Option<f64>,
    ) -> Result<Mesh, StrError> {
        // check
        if nr < 1 {
            return Err("number of divisions along the radius must be ≥ 1");
        }
        if na < 1 {
            return Err("number of divisions along alpha must be ≥ 1");
        }
        if thickness <= 0.0 {
            return Err("thickness must be > 0.0");
        }
        if target.class() != GeoClass::Tet {
            return Err("the GeoClass of target must be Tet");
        }

        // generate o2 triangles (i.e., Tet10) (need to use o2 for others too, e.g.,
        // Tet20, because the the middle-edge markers will be replicated by tetgen
        let o2 = if target.nnode() > 4 { true } else { false };

        // allocate data
        let nz = 1;
        let npoint_cap_facets = 2 * (nr + 1) + 2 * (na - 1);
        let npoint_sym_facets = 2 * (nr + 1); // + 2 * (nz - 1)
        let npoint_cyl_facets = 2 * (na + 1); // + 2 * (nz - 1)
        let npoint = 2 * npoint_cap_facets;
        let nfacet = 2 + 2 + 2 * na;
        let mut facets_npoint = vec![0; nfacet];
        facets_npoint[0] = npoint_cap_facets; // cap facet @ z = 0
        facets_npoint[1] = npoint_cap_facets; // cap facet @ z = thickness
        facets_npoint[2] = npoint_sym_facets; // sym facet parallel to x
        facets_npoint[3] = npoint_sym_facets; // sym facet parallel to y
        for i in 0..na {
            facets_npoint[4 + i] = 4; // inner cyl facets
            facets_npoint[4 + na + i] = 4; // outer cyl facets
        }
        let nregion = 1;
        // let mut tetgen = Tetgen::new(npoint, None, None, None)?;
        let mut tetgen = Tetgen::new(npoint, Some(facets_npoint), Some(nregion), None)?;

        // constants
        const AMIN: f64 = 0.0;
        const AMAX: f64 = PI / 2.0;
        let dr = (rmax - rmin) / (nr as f64);
        let da = (AMAX - AMIN) / (na as f64);
        let dz = thickness / (nz as f64);

        // global index of point
        let mut p = 0;

        for k in 0..(nz + 1) {
            // calc z
            let z = (k as f64) * dz;

            // horizontal line
            tetgen.set_point(p, -1, rmin, 0.0, z)?;
            println!("p={} x={:?}, y={:?}, z={:?}", p, rmin, 0.0, z);
            p += 1;
            for i in 1..(nr + 1) {
                let marker = if i == nr { -2 } else { 0 };
                let r = rmin + (i as f64) * dr;
                println!("p={} x={:?}, y={:?}, z={:?}", p, r, 0.0, z);
                tetgen.set_point(p, marker, r, 0.0, z)?;
                p += 1;
            }

            // outer circle
            for i in 1..(na + 1) {
                let marker = if i == na { -3 } else { 0 };
                let a = AMIN + (i as f64) * da;
                let (x, y) = if i == na {
                    (0.0, rmax)
                } else {
                    (rmax * f64::cos(a), rmax * f64::sin(a))
                };
                println!("p={} x={:?}, y={:?}, z={:?}", p, x, y, z);
                tetgen.set_point(p, marker, x, y, z)?;
                p += 1;
            }

            // vertical line
            for i in 1..(nr + 1) {
                let marker = if i == nr { -4 } else { 0 };
                let r = rmin + ((nr - i) as f64) * dr;
                println!("p={} x={:?}, y={:?}, z={:?}", p, 0.0, r, z);
                tetgen.set_point(p, marker, 0.0, r, z)?;
                p += 1;
            }

            // inner circle
            for i in 1..na {
                let a = AMIN + ((na - i) as f64) * da;
                let x = rmin * f64::cos(a);
                let y = rmin * f64::sin(a);
                println!("p={} x={:?}, y={:?}, z={:?}", p, x, y, z);
                tetgen.set_point(p, 0, x, y, z)?;
                p += 1;
            }
        }

        // auxiliary indices
        let mut m; // local index of point on facet
        let p_pivot_a = 0; // pivot point @ (rmin, 0, 0)
        let p_pivot_b = p_pivot_a + nr + na + nr; // pivot point @ (0, rmin, 0)
        let p_pivot_c = p_pivot_a + nr; // pivot point @ (rmax, 0, 0)
        let p_pivot_d = p_pivot_b - nr; // pivot point @ (0, rmax, 0)
        let p_pivot_aa = p_pivot_a + npoint_cap_facets; // pivot point @ (rmin, 0, thickness)
        let p_pivot_bb = p_pivot_b + npoint_cap_facets; // pivot point @ (0, rmin, thickness)
        let p_pivot_cc = p_pivot_aa + nr; // pivot point @ (rmax, 0, thickness)
        let p_pivot_dd = p_pivot_bb - nr; // pivot point @ (rmax, 0, thickness)

        // On page 52 of TetGen's manual:
        // * Each polygon of a facet is described by an ordered list of vertices.
        // * The order of the vertices can be in either clockwise or counterclockwise order.

        // cap facet @ z = 0
        println!("cap facet @ z = 0");
        for i in 0..npoint_cap_facets {
            m = i;
            p = p_pivot_a + i;
            println!("m={} p={}", m, p);
            tetgen.set_facet_point(0, m, p)?;
        }

        // cap facet @ z = thickness
        println!("cap facet @ z = thickness");
        for i in 0..npoint_cap_facets {
            m = i;
            p = p_pivot_aa + i;
            println!("m={} p={}", m, p);
            tetgen.set_facet_point(1, m, p)?;
        }

        // sym facet parallel to x
        println!("sym facet parallel to x");
        for i in 0..(nr + 1) {
            m = i;
            p = p_pivot_a + i;
            print!("m={} p={}", m, p);
            tetgen.set_facet_point(2, m, p)?;
            m = nr + 1 + nr - i;
            p = p_pivot_aa + i;
            println!("   m={} p={}", m, p);
            tetgen.set_facet_point(2, m, p)?;
        }

        // sym facet parallel to y
        println!("sym facet parallel to y");
        for i in 0..(nr + 1) {
            m = i;
            p = p_pivot_b - i;
            print!("m={} p={}", m, p);
            tetgen.set_facet_point(3, m, p)?;
            m = nr + 1 + nr - i;
            p = p_pivot_bb - i;
            println!("   m={} p={}", m, p);
            tetgen.set_facet_point(3, m, p)?;
        }

        // inner cyl facets
        println!("inner cyl facet");
        let pp = [p_pivot_b, p_pivot_b + 1, p_pivot_bb + 1, p_pivot_bb];
        for i in 0..(na - 1) {
            for m in 0..4 {
                p = pp[m] + i;
                println!("m={} p={}", m, p);
                tetgen.set_facet_point(4 + i, m, p)?;
            }
        }
        tetgen.set_facet_point(4 + na - 1, 0, p_pivot_b + na - 1)?;
        tetgen.set_facet_point(4 + na - 1, 1, p_pivot_a)?;
        tetgen.set_facet_point(4 + na - 1, 2, p_pivot_aa)?;
        tetgen.set_facet_point(4 + na - 1, 3, p_pivot_bb + na - 1)?;

        // outer cyl facets
        println!("outer cyl facet");
        let pp = [p_pivot_c, p_pivot_c + 1, p_pivot_cc + 1, p_pivot_cc];
        for i in 0..(na - 1) {
            for m in 0..4 {
                p = pp[m] + i;
                println!("m={} p={}", m, p);
                tetgen.set_facet_point(4 + na + i, m, p)?;
            }
        }
        tetgen.set_facet_point(4 + 2 * na - 1, 0, p_pivot_c + na - 1)?;
        tetgen.set_facet_point(4 + 2 * na - 1, 1, p_pivot_d)?;
        tetgen.set_facet_point(4 + 2 * na - 1, 2, p_pivot_dd)?;
        tetgen.set_facet_point(4 + 2 * na - 1, 3, p_pivot_cc + na - 1)?;

        // region
        tetgen.set_region(0, 1, rmin + 1e-4, 1e-4, 1e-4, None)?;

        // generate mesh
        tetgen.generate_mesh(false, o2, global_max_volume, None)?;
        // tetgen.generate_delaunay(true)?;

        // allocate data
        const NDIM: usize = 3;
        let npoint = tetgen.out_npoint();
        let ncell = tetgen.out_ncell();
        let nnode = tetgen.out_cell_npoint();
        let kind = if o2 { GeoKind::Tet10 } else { GeoKind::Tet4 };
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
            let marker = tetgen.out_point_marker(i);
            mesh.points[i].id = i;
            mesh.points[i].marker = marker;
            mesh.points[i].coords[0] = tetgen.out_point(i, 0);
            mesh.points[i].coords[1] = tetgen.out_point(i, 1);
            mesh.points[i].coords[2] = tetgen.out_point(i, 2);
        }
        for i in 0..ncell {
            mesh.cells[i].id = i;
            mesh.cells[i].attribute = tetgen.out_cell_attribute(i);
            for m in 0..nnode {
                mesh.cells[i].points[m] = tetgen.out_cell_point(i, m);
            }
        }

        // apply constraints (need to be done before the upgrade because
        // Steiner points may be added even for Tri3)
        // apply_constraints(&mut mesh, rmin, rmax);

        // upgrade mesh (need to apply constraints again because
        // new mid-edge points may be created)
        // if target.nnode() > 10 {
        //     let mut new_mesh = mesh.convert_2d(target)?;
        //     apply_constraints(&mut new_mesh, rmin, rmax);
        //     return Ok(new_mesh);
        // }

        // results
        Ok(mesh)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Unstructured;
    use crate::geometry::point_point_distance;
    use crate::mesh::{At, Features, Figure, Mesh};
    use crate::shapes::GeoKind;
    use crate::util::any_x;
    use plotpy::Surface;
    use russell_lab::approx_eq;

    const RMIN: f64 = 3.0;
    const RMAX: f64 = 6.0;
    const SAVE_FIGURE: bool = false;

    fn draw(mesh: &Mesh, larger: bool, filename: &str) {
        let mut fig = Figure::new();
        fig.cell_ids = true;
        fig.point_ids = true;
        if larger {
            fig.figure_size = Some((800.0, 800.0));
        } else {
            fig.figure_size = Some((600.0, 600.0));
        }
        mesh.draw(Some(fig), filename, |_, _| {}).unwrap();
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
        let feat = Features::new(mesh, false);
        let res = feat.search_point_ids(At::XY(RMIN, 0.0), any_x).unwrap();
        assert_eq!(res.len(), 1);
        assert_eq!(mesh.points[res[0]].marker, -1);

        let res = feat.search_point_ids(At::XY(RMAX, 0.0), any_x).unwrap();
        assert_eq!(res.len(), 1);
        assert_eq!(mesh.points[res[0]].marker, -2);

        let res = feat.search_point_ids(At::XY(0.0, RMAX), any_x).unwrap();
        assert_eq!(res.len(), 1);
        assert_eq!(mesh.points[res[0]].marker, -3);

        let res = feat.search_point_ids(At::XY(0.0, RMIN), any_x).unwrap();
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

    #[test]
    fn tri_quarter_ring_3d_works() {
        let mesh = Unstructured::quarter_ring_3d(RMIN, RMAX, 1.0, 2, 4, GeoKind::Tet4, None).unwrap();
        if true {
            let mut cylin_in = Surface::new();
            let mut cylin_out = Surface::new();
            cylin_in.set_solid_color("#ff000020");
            cylin_out.set_solid_color("#ff000020");
            cylin_in
                .draw_cylinder(&[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], RMIN, 5, 81)
                .unwrap();
            cylin_out
                .draw_cylinder(&[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], RMAX, 5, 81)
                .unwrap();

            let mut fig = Figure::new();
            fig.figure_size = Some((1000.0, 1000.0));
            fig.point_ids = true;
            mesh.draw(Some(fig), "/tmp/gemlab/test_tri_quarter_ring_3d.svg", |plot, before| {
                if before {
                    plot.add(&cylin_in).add(&cylin_out);
                }
            });
        }
        // assert_eq!(mesh.points.len(), 14);
        // assert_eq!(mesh.cells.len(), 14);
        // mesh.check_all().unwrap();
        // mesh.check_overlapping_points(0.18).unwrap();
        // check_constraints(&mesh);
        // check_corner_markers(&mesh);
        // check_edge_point_markers(&mesh, &[11, 10, 9], &[3, 4, 5]);
        // for p in [0, 11, 10, 9, 8] {
        //     let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
        //     approx_eq(d, RMIN, 1e-15);
        // }
        // for p in [2, 3, 4, 5, 6] {
        //     let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
        //     approx_eq(d, RMAX, 1e-15);
        // }
    }
}
