use super::{CellMarker, Mesh};
use crate::graph::GraphUnd;
use crate::mesh::Features;
use crate::shapes::{GeoClass, GeoKind};
use crate::StrError;
use russell_lab::math::PI;
use std::collections::HashMap;
use tritet::{Tetgen, Trigen};

/// Groups generators of unstructured meshes (Tri and Tet only)
pub struct Unstructured {}

const MARKER_A: i32 = -1;
const MARKER_B: i32 = -2;
const MARKER_C: i32 = -3;
const MARKER_D: i32 = -4;
const MARKER_INNER: i32 = -5;
const MARKER_OUTER: i32 = -6;
const MARKER_SYM_X: i32 = -7;
const MARKER_SYM_Y: i32 = -8;
const MARKER_ZMIN: i32 = -9;
const MARKER_ZMAX: i32 = -10;

fn apply_constraints(mesh: &mut Mesh, rmin: f64, rmax: f64) {
    // Point markers:
    //                              higher marker values
    //  D ==---__                   have PRIORITY over
    //  |        '*._               lower marker values
    //  | SYM_Y      *._
    //  |               *.  OUTER   A: -1    INNER: -5
    //  B ==-__           *.        B: -2    OUTER: -6
    //         '-.          *       C: -3    SYM_X: -7
    //            *.         *      D: -4    SYM_Y: -8
    //       INNER  *         *
    //               *         *    ZMIN:  -9
    //                *         *   ZMAX: -10
    //                #  SYM_X  #
    //                A ------- C

    // apply constraints
    const TOL: f64 = 1e-15;
    for i in 0..mesh.points.len() {
        // points on inner circle
        if mesh.points[i].marker == MARKER_INNER {
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
        if mesh.points[i].marker == MARKER_OUTER {
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
    /// Allocates a mesh from the data stored in a Trigen instance
    pub fn from_trigen(trigen: &Trigen) -> Mesh {
        // allocate data
        const NDIM: usize = 2;
        let npoint = trigen.out_npoint();
        let ncell = trigen.out_ncell();
        let nnode = trigen.out_cell_npoint();
        let kind = if nnode == 6 { GeoKind::Tri6 } else { GeoKind::Tri3 };
        let mut mesh = Mesh::new_zero_homogeneous(NDIM, npoint, ncell, kind).unwrap();

        // set essential mesh data
        for i in 0..npoint {
            let marker = trigen.out_point_marker(i);
            mesh.points[i].id = i;
            mesh.points[i].marker = marker;
            mesh.points[i].coords[0] = trigen.out_point(i, 0);
            mesh.points[i].coords[1] = trigen.out_point(i, 1);
        }
        for i in 0..ncell {
            mesh.cells[i].id = i;
            mesh.cells[i].marker = trigen.out_cell_marker(i);
            for m in 0..nnode {
                mesh.cells[i].points[m] = trigen.out_cell_point(i, m);
            }
        }

        // set marked edges
        for i in 0..trigen.out_nsegment() {
            let marker = trigen.out_segment_marker(i);
            if marker != 0 {
                let a = trigen.out_segment_point(i, 0);
                let b = trigen.out_segment_point(i, 1);
                mesh.marked_edges.push((marker, a, b));
            }
        }
        mesh
    }

    /// Allocates a mesh from the data stored in a Tetgen instance
    ///
    /// # Notes
    ///
    /// 1. Zero face markers are ignored.
    /// 2. TetGen automatically assigns the marker 1 for points on the boundary
    /// thus, we cannot use the marker 1 to identify corner points
    pub fn from_tetgen(tetgen: &Tetgen) -> Mesh {
        // allocate data
        const NDIM: usize = 3;
        let npoint = tetgen.out_npoint();
        let ncell = tetgen.out_ncell();
        let nnode = tetgen.out_cell_npoint();
        let kind = if nnode == 10 { GeoKind::Tet10 } else { GeoKind::Tet4 };
        let mut mesh = Mesh::new_zero_homogeneous(NDIM, npoint, ncell, kind).unwrap();

        // set mesh data
        for i in 0..npoint {
            // note: TetGen automatically assigns the marker 1 for points on the boundary
            // thus, we cannot use the marker 1 to identify corner points
            let tet_mark = tetgen.out_point_marker(i);
            let marker = if tet_mark == 1 { 0 } else { tet_mark };
            mesh.points[i].id = i;
            mesh.points[i].marker = marker;
            mesh.points[i].coords[0] = tetgen.out_point(i, 0);
            mesh.points[i].coords[1] = tetgen.out_point(i, 1);
            mesh.points[i].coords[2] = tetgen.out_point(i, 2);
        }
        for i in 0..ncell {
            mesh.cells[i].id = i;
            mesh.cells[i].marker = tetgen.out_cell_marker(i);
            for m in 0..nnode {
                mesh.cells[i].points[m] = tetgen.out_cell_point(i, m);
            }
        }

        // set marked faces
        let mut pp = [0i32; 6];
        for i in 0..tetgen.out_n_marked_face() {
            let (marker, _) = tetgen.out_marked_face(i, &mut pp);
            if marker != 0 {
                mesh.marked_faces
                    .push((marker, pp[0] as usize, pp[1] as usize, pp[2] as usize, usize::MAX));
            }
        }
        mesh
    }

    /// Generates a triangular mesh from a Planar Straight Line Graph (PSLG) defined by a Mesh
    pub fn call_trigen(
        pslg: &Mesh,
        holes: &Vec<(f64, f64)>,
        o2: bool,
        max_areas: Option<HashMap<CellMarker, f64>>,
        global_max_area: Option<f64>,
        global_min_angle: Option<f64>,
        renumber: bool,
    ) -> Result<Mesh, StrError> {
        // check
        if pslg.ndim != 2 {
            return Err("the PSLG mesh must be 2D");
        }

        // constants
        let npoint = pslg.points.len();
        let nregion = pslg.cells.len();
        let nhole = holes.len();
        let extract_all = true; // we need interior edges as well
        let features = Features::new(&pslg, extract_all);
        let nsegment = features.edges.len();

        // allocate trigen structure
        let mut trigen = Trigen::new(npoint, Some(nsegment), Some(nregion), Some(nhole))?;

        // set points
        for i in 0..npoint {
            trigen.set_point(
                i,
                pslg.points[i].marker,
                pslg.points[i].coords[0],
                pslg.points[i].coords[1],
            )?;
        }

        // set segments
        let mut sorted_edge_keys: Vec<_> = features.edges.keys().cloned().collect();
        sorted_edge_keys.sort();
        for i in 0..sorted_edge_keys.len() {
            let edge = &features.edges.get(&sorted_edge_keys[i]).unwrap();
            if edge.points.len() != 2 {
                return Err("edges must have 2 points");
            }
            trigen.set_segment(i, edge.marker, edge.points[0], edge.points[1])?;
        }

        // set regions
        for cell in &pslg.cells {
            let mut cen_x = 0.0;
            let mut cen_y = 0.0;
            let nnode = cell.points.len();
            for p in &cell.points {
                cen_x += pslg.points[*p].coords[0];
                cen_y += pslg.points[*p].coords[1];
            }
            cen_x /= nnode as f64;
            cen_y /= nnode as f64;
            let max_area = max_areas.as_ref().and_then(|ma| ma.get(&cell.marker)).cloned();
            trigen.set_region(cell.id, cell.marker, cen_x, cen_y, max_area)?;
        }

        // set holes
        for (i, hole) in holes.iter().enumerate() {
            trigen.set_hole(i, hole.0, hole.1)?;
        }

        // generate mesh
        trigen.generate_mesh(false, o2, true, global_max_area, global_min_angle)?;
        let mut mesh = Unstructured::from_trigen(&trigen);

        // results
        if renumber {
            GraphUnd::renumber_mesh(&mut mesh, false)?;
        }
        Ok(mesh)
    }

    /// Generates a tetrahedral mesh from a Piecewise Linear Complex (PLC) defined by a Mesh
    ///
    /// The PLC must contain only faces (shells) and each face must be either a Tri3 or a Qua4.
    pub fn call_tetgen(
        plc: &Mesh,
        regions: &Vec<(i32, f64, f64, f64)>,
        holes: &Vec<(f64, f64, f64)>,
        o2: bool,
        max_volumes: Option<HashMap<CellMarker, f64>>,
        global_max_volume: Option<f64>,
        global_min_angle: Option<f64>,
        renumber: bool,
    ) -> Result<Mesh, StrError> {
        // check
        if plc.ndim != 3 {
            return Err("the PLC mesh must be 3D");
        }

        // constants
        let npoint = plc.points.len();
        let nregion = regions.len();
        let nhole = holes.len();

        // extract features
        let extract_all = true; // we need all facets
        let features = Features::new(&plc, extract_all);
        let nfacet = features.shells.len() + features.faces.len();

        // shells: counter the number of points on each facet
        let mut facet_npoint = Vec::with_capacity(nfacet);
        for cell_id in &features.shells {
            facet_npoint.push(plc.cells[*cell_id].points.len());
        }

        // faces: counter the number of points on each facet
        let mut sorted_face_keys: Vec<_> = features.faces.keys().cloned().collect();
        sorted_face_keys.sort();
        for i in 0..sorted_face_keys.len() {
            let face = &features.faces.get(&sorted_face_keys[i]).unwrap();
            facet_npoint.push(face.points.len());
        }

        // allocate trigen structure
        let mut tetgen = Tetgen::new(npoint, Some(facet_npoint), Some(nregion), Some(nhole))?;

        // set points
        for i in 0..npoint {
            tetgen.set_point(
                i,
                plc.points[i].marker,
                plc.points[i].coords[0],
                plc.points[i].coords[1],
                plc.points[i].coords[2],
            )?;
        }

        // shells: set facets
        for (i, cell_id) in features.shells.iter().enumerate() {
            let facet = &plc.cells[*cell_id];
            tetgen.set_facet_marker(i, facet.marker)?;
            for (m, p) in facet.points.iter().enumerate() {
                tetgen.set_facet_point(i, m, *p)?;
            }
        }

        // faces: set facets
        let start = features.shells.len();
        for (i, face_key) in sorted_face_keys.iter().enumerate() {
            let facet = &features.faces.get(face_key).unwrap();
            tetgen.set_facet_marker(start + i, facet.marker)?;
            for (m, p) in facet.points.iter().enumerate() {
                tetgen.set_facet_point(start + i, m, *p)?;
            }
        }

        // set regions
        for (i, region) in regions.iter().enumerate() {
            let max_volume = max_volumes.as_ref().and_then(|ma| ma.get(&region.0)).cloned();
            tetgen.set_region(i, region.0, region.1, region.2, region.3, max_volume)?;
        }

        // set holes
        for (i, hole) in holes.iter().enumerate() {
            tetgen.set_hole(i, hole.0, hole.1, hole.2)?;
        }

        // generate mesh
        tetgen.generate_mesh(false, o2, global_max_volume, global_min_angle)?;
        let mut mesh = Unstructured::from_tetgen(&tetgen);

        // results
        if renumber {
            GraphUnd::renumber_mesh(&mut mesh, false)?;
        }
        Ok(mesh)
    }

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
    /// (-4)                          higher marker values
    ///   D ==---__                   have PRIORITY over
    ///   |        '*._               lower marker values
    ///   | SYM_Y(-8)   *._   (-6)
    ///   |               *.  OUTER   A: -1    INNER: -5
    ///   B ==-__           *.        B: -2    OUTER: -6
    /// (-2)     '-.          *       C: -3    SYM_X: -7
    ///             *.         *      D: -4    SYM_Y: -8
    ///        INNER  *         *
    ///         (-5)    *        *    ZMIN:  -9
    ///                 *  (-7)   *   ZMAX: -10
    ///                 #  SYM_X  #
    ///                 A ------- C
    ///               (-1)      (-3)
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
        renumber: bool,
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

        // global index of point
        let mut p = 0;

        // y=0 line (horizontal line)
        trigen.set_point(p, MARKER_A, rmin, 0.0)?;
        p += 1;
        for i in 1..(nr + 1) {
            let marker = if i == nr { MARKER_C } else { 0 };
            let r = rmin + (i as f64) * dr;
            trigen.set_point(p, marker, r, 0.0)?;
            trigen.set_segment(p - 1, MARKER_SYM_X, p - 1, p)?;
            p += 1;
        }

        // outer circle
        for i in 1..(na + 1) {
            let marker = if i == na { MARKER_D } else { 0 };
            let a = AMIN + (i as f64) * da;
            let (x, y) = if i == na {
                (0.0, rmax)
            } else {
                (rmax * f64::cos(a), rmax * f64::sin(a))
            };
            trigen.set_point(p, marker, x, y)?;
            trigen.set_segment(p - 1, MARKER_OUTER, p - 1, p)?;
            p += 1;
        }

        // x=0 line (vertical line)
        for i in 1..(nr + 1) {
            let marker = if i == nr { MARKER_B } else { 0 };
            let r = rmin + ((nr - i) as f64) * dr;
            trigen.set_point(p, marker, 0.0, r)?;
            trigen.set_segment(p - 1, MARKER_SYM_Y, p - 1, p)?;
            p += 1;
        }

        // inner circle
        for i in 1..na {
            let a = AMIN + ((na - i) as f64) * da;
            let x = rmin * f64::cos(a);
            let y = rmin * f64::sin(a);
            trigen.set_point(p, 0, x, y)?;
            trigen.set_segment(p - 1, MARKER_INNER, p - 1, p)?;
            p += 1;
        }
        trigen.set_segment(p - 1, MARKER_INNER, p - 1, 0)?;

        // region
        trigen.set_region(0, 1, rmin + 1e-4, 1e-4, None)?;

        // generate mesh
        trigen.generate_mesh(false, o2, true, global_max_area, None)?;

        // convert trigen to Mesh
        let mut mesh = Unstructured::from_trigen(&trigen);

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
        if renumber {
            GraphUnd::renumber_mesh(&mut mesh, false)?;
        }
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
    /// (-4)                          higher marker values
    ///   D ==---__                   have PRIORITY over
    ///   |        '*._               lower marker values
    ///   | SYM_Y(-8)   *._   (-6)
    ///   |               *.  OUTER   A: -1    INNER: -5
    ///   B ==-__           *.        B: -2    OUTER: -6
    /// (-2)     '-.          *       C: -3    SYM_X: -7
    ///             *.         *      D: -4    SYM_Y: -8
    ///        INNER  *         *
    ///         (-5)    *        *    ZMIN:  -9
    ///                 *  (-7)   *   ZMAX: -10
    ///                 #  SYM_X  #
    ///                 A ------- C
    ///               (-1)      (-3)
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
        renumber: bool,
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
        if target != GeoKind::Tet4 && target != GeoKind::Tet10 {
            return Err("only Tet4 and Tet10 are available currently");
        }

        // generate o2 triangles (i.e., Tet10) (need to use o2 for others too, e.g.,
        // Tet20, because the the middle-edge markers will be replicated by tetgen
        let o2 = if target.nnode() > 4 { true } else { false };

        // allocate data
        let nz = 1;
        let npoint_cap_facets = 2 * (nr + 1) + 2 * (na - 1);
        let npoint_sym_facets = 2 * (nr + 1); // + 2 * (nz - 1)
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

            // y=0 line
            tetgen.set_point(p, MARKER_A, rmin, 0.0, z)?;
            p += 1;
            for i in 1..(nr + 1) {
                let marker = if i == nr { MARKER_C } else { 0 };
                let r = rmin + (i as f64) * dr;
                tetgen.set_point(p, marker, r, 0.0, z)?;
                p += 1;
            }

            // outer circle
            for i in 1..(na + 1) {
                let marker = if i == na { MARKER_D } else { 0 };
                let a = AMIN + (i as f64) * da;
                let (x, y) = if i == na {
                    (0.0, rmax)
                } else {
                    (rmax * f64::cos(a), rmax * f64::sin(a))
                };
                tetgen.set_point(p, marker, x, y, z)?;
                p += 1;
            }

            // x=0 line
            for i in 1..(nr + 1) {
                let marker = if i == nr { MARKER_B } else { 0 };
                let r = rmin + ((nr - i) as f64) * dr;
                tetgen.set_point(p, marker, 0.0, r, z)?;
                p += 1;
            }

            // inner circle
            for i in 1..na {
                let a = AMIN + ((na - i) as f64) * da;
                let x = rmin * f64::cos(a);
                let y = rmin * f64::sin(a);
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
        let p_pivot_cc = p_pivot_c + npoint_cap_facets; // pivot point @ (rmax, 0, thickness)
        let p_pivot_dd = p_pivot_d + npoint_cap_facets; // pivot point @ (0, rmax, thickness)

        // On page 52 of TetGen's manual:
        // * Each polygon of a facet is described by an ordered list of vertices.
        // * The order of the vertices can be in either clockwise or counterclockwise order.

        // cap facet @ z = 0
        tetgen.set_facet_marker(0, MARKER_ZMIN)?;
        for i in 0..npoint_cap_facets {
            m = i;
            p = p_pivot_a + i;
            tetgen.set_facet_point(0, m, p)?;
        }

        // cap facet @ z = thickness
        tetgen.set_facet_marker(1, MARKER_ZMAX)?;
        for i in 0..npoint_cap_facets {
            m = i;
            p = p_pivot_aa + i;
            tetgen.set_facet_point(1, m, p)?;
        }

        // sym facet parallel to x (y=0)
        tetgen.set_facet_marker(2, MARKER_SYM_X)?;
        for i in 0..(nr + 1) {
            m = i;
            p = p_pivot_a + i;
            tetgen.set_facet_point(2, m, p)?;
            m = nr + 1 + nr - i;
            p = p_pivot_aa + i;
            tetgen.set_facet_point(2, m, p)?;
        }

        // sym facet parallel to y
        tetgen.set_facet_marker(3, MARKER_SYM_Y)?;
        for i in 0..(nr + 1) {
            m = i;
            p = p_pivot_b - i;
            tetgen.set_facet_point(3, m, p)?;
            m = nr + 1 + nr - i;
            p = p_pivot_bb - i;
            tetgen.set_facet_point(3, m, p)?;
        }

        // inner cyl facets
        let pp = [p_pivot_b, p_pivot_b + 1, p_pivot_bb + 1, p_pivot_bb];
        for i in 0..(na - 1) {
            tetgen.set_facet_marker(4 + i, MARKER_INNER)?;
            for m in 0..4 {
                p = pp[m] + i;
                tetgen.set_facet_point(4 + i, m, p)?;
            }
        }
        let last = 4 + na - 1;
        tetgen.set_facet_marker(last, MARKER_INNER)?;
        tetgen.set_facet_point(last, 0, p_pivot_b + na - 1)?;
        tetgen.set_facet_point(last, 1, p_pivot_a)?;
        tetgen.set_facet_point(last, 2, p_pivot_aa)?;
        tetgen.set_facet_point(last, 3, p_pivot_bb + na - 1)?;

        // outer cyl facets
        let pp = [p_pivot_c, p_pivot_c + 1, p_pivot_cc + 1, p_pivot_cc];
        for i in 0..(na - 1) {
            tetgen.set_facet_marker(4 + na + i, MARKER_OUTER)?;
            for m in 0..4 {
                p = pp[m] + i;
                tetgen.set_facet_point(4 + na + i, m, p)?;
            }
        }
        let last = 4 + na + na - 1;
        tetgen.set_facet_marker(last, MARKER_OUTER)?;
        tetgen.set_facet_point(last, 0, p_pivot_c + na - 1)?;
        tetgen.set_facet_point(last, 1, p_pivot_d)?;
        tetgen.set_facet_point(last, 2, p_pivot_dd)?;
        tetgen.set_facet_point(last, 3, p_pivot_cc + na - 1)?;

        // region
        tetgen.set_region(0, 1, rmin + 1e-4, 1e-4, 1e-4, None)?;

        // generate mesh
        tetgen.generate_mesh(false, o2, global_max_volume, None)?;

        // convert trigen to Mesh
        let mut mesh = Unstructured::from_tetgen(&tetgen);

        // set markers of points on marked faces
        let n_face_point = if o2 { 6 } else { 3 };
        let mut face_points = [0; 6];
        for i in 0..tetgen.out_n_marked_face() {
            let (marker, _) = tetgen.out_marked_face(i, &mut face_points);
            for m in 0..n_face_point {
                let p = face_points[m] as usize;
                if mesh.points[p].marker == 0 || mesh.points[p].marker < marker {
                    // set only if it has a higher priority (higher marker == higher priority)
                    mesh.points[p].marker = marker;
                }
            }
        }

        // apply constraints
        apply_constraints(&mut mesh, rmin, rmax);

        // results
        if renumber {
            GraphUnd::renumber_mesh(&mut mesh, false)?;
        }
        Ok(mesh)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::point_point_distance;
    use crate::mesh::{At, Draw, Features, Samples};
    use crate::util::any_x;
    use plotpy::Surface;
    use russell_lab::approx_eq;

    const RMIN: f64 = 3.0;
    const RMAX: f64 = 6.0;
    const SAVE_FIGURE: bool = false;
    const MAX_NPOINT_PRINT: usize = 200;

    fn print_bandwidth(mesh: &mut Mesh) {
        let graph = GraphUnd::from_mesh(&mesh, true, false).unwrap();
        if mesh.points.len() < MAX_NPOINT_PRINT {
            graph.print_non_zero_pattern();
        }
        GraphUnd::renumber_mesh(mesh, false).unwrap();
        let graph_after = GraphUnd::from_mesh(&mesh, true, false).unwrap();
        if mesh.points.len() < MAX_NPOINT_PRINT {
            graph_after.print_non_zero_pattern();
        }
        println!("bandwidth (before) = {}", graph.calc_bandwidth());
        println!("bandwidth (after)  = {}", graph_after.calc_bandwidth());
    }

    fn draw(mesh: &Mesh, larger: bool, filename: &str) {
        let mut draw = Draw::new();
        draw.show_cell_ids(true)
            .show_point_ids(true)
            .show_cell_marker(true)
            .show_point_marker(true)
            .show_edge_markers(true)
            .show_face_markers(true)
            .set_view_flag(false);
        if larger {
            draw.set_size(800.0, 800.0);
        } else {
            draw.set_size(600.0, 600.0);
        }
        draw.all(&mesh, filename).unwrap();
    }

    fn check_corner_markers(mesh: &Mesh, o2_3d: bool) {
        let (m, n) = if mesh.ndim == 3 {
            if o2_3d {
                (3, 2)
            } else {
                (2, 2)
            }
        } else {
            (1, 1)
        };

        let feat = Features::new(mesh, false);

        let res = feat.search_point_ids(At::XY(RMIN, 0.0), any_x).unwrap();
        assert_eq!(res.len(), m);
        for i in 0..n {
            assert_eq!(mesh.points[res[i]].marker, MARKER_A);
        }

        let res = feat.search_point_ids(At::XY(RMAX, 0.0), any_x).unwrap();
        assert_eq!(res.len(), m);
        for i in 0..n {
            assert_eq!(mesh.points[res[i]].marker, MARKER_C);
        }

        let res = feat.search_point_ids(At::XY(0.0, RMAX), any_x).unwrap();
        assert_eq!(res.len(), m);
        for i in 0..n {
            assert_eq!(mesh.points[res[i]].marker, MARKER_D);
        }

        let res = feat.search_point_ids(At::XY(0.0, RMIN), any_x).unwrap();
        assert_eq!(res.len(), m);
        for i in 0..n {
            assert_eq!(mesh.points[res[i]].marker, MARKER_B);
        }
    }

    fn check_point_markers(mesh: &Mesh, inner: &[usize], outer: &[usize], sym_x: &[usize], sym_y: &[usize]) {
        for p in inner {
            assert_eq!(mesh.points[*p].marker, MARKER_INNER);
        }
        for p in outer {
            assert_eq!(mesh.points[*p].marker, MARKER_OUTER);
        }
        for p in sym_x {
            assert_eq!(mesh.points[*p].marker, MARKER_SYM_X);
        }
        for p in sym_y {
            assert_eq!(mesh.points[*p].marker, MARKER_SYM_Y);
        }
    }

    fn check_constraints(mesh: &Mesh) {
        let inner: Vec<_> = mesh.points.iter().filter(|p| p.marker == MARKER_INNER).collect();
        let outer: Vec<_> = mesh.points.iter().filter(|p| p.marker == MARKER_OUTER).collect();
        for point in inner {
            let d = point_point_distance(&point.coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for point in outer {
            let d = point_point_distance(&point.coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
    }

    #[test]
    fn from_trigen_works() {
        let mut trigen = Trigen::new(3, Some(3), Some(1), None).unwrap();
        trigen.set_point(0, -100, 0.0, 0.0).unwrap();
        trigen.set_point(1, -200, 1.0, 0.0).unwrap();
        trigen.set_point(2, -300, 0.0, 1.0).unwrap();
        trigen.set_segment(0, -10, 0, 1).unwrap();
        trigen.set_segment(1, -20, 1, 2).unwrap();
        trigen.set_segment(2, -30, 2, 0).unwrap();
        trigen.set_region(0, 8, 0.1, 0.1, None).unwrap();
        trigen.generate_mesh(false, false, false, None, None).unwrap();
        let mesh = Unstructured::from_trigen(&trigen);
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_from_trigen.svg");
        }
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 3);
        assert_eq!(mesh.cells.len(), 1);
        assert_eq!(mesh.points[0].id, 0);
        assert_eq!(mesh.points[0].marker, -100);
        assert_eq!(mesh.points[0].coords, vec![0.0, 0.0]);
        assert_eq!(mesh.points[1].id, 1);
        assert_eq!(mesh.points[1].marker, -200);
        assert_eq!(mesh.points[1].coords, vec![1.0, 0.0]);
        assert_eq!(mesh.points[2].id, 2);
        assert_eq!(mesh.points[2].marker, -300);
        assert_eq!(mesh.points[2].coords, vec![0.0, 1.0]);
        assert_eq!(mesh.cells[0].id, 0);
        assert_eq!(mesh.cells[0].marker, 8);
        assert_eq!(mesh.cells[0].kind, GeoKind::Tri3);
        assert_eq!(mesh.cells[0].points, vec![0, 1, 2]);
    }

    #[test]
    fn from_tetgen_works() {
        // allocate data for 4 points
        let mut tetgen = Tetgen::new(4, Some(vec![3, 3, 3, 3]), Some(1), None).unwrap();

        // set points
        tetgen.set_point(0, 0, 0.0, 1.0, 0.0).unwrap();
        tetgen.set_point(1, 0, 0.0, 0.0, 0.0).unwrap();
        tetgen.set_point(2, 0, 1.0, 1.0, 0.0).unwrap();
        tetgen.set_point(3, 0, 0.0, 1.0, 1.0).unwrap();

        // set facets
        // 0
        tetgen.set_facet_point(0, 0, 0).unwrap();
        tetgen.set_facet_point(0, 1, 2).unwrap();
        tetgen.set_facet_point(0, 2, 1).unwrap();
        // 1
        tetgen.set_facet_point(1, 0, 0).unwrap();
        tetgen.set_facet_point(1, 1, 1).unwrap();
        tetgen.set_facet_point(1, 2, 3).unwrap();
        // 2
        tetgen.set_facet_point(2, 0, 0).unwrap();
        tetgen.set_facet_point(2, 1, 3).unwrap();
        tetgen.set_facet_point(2, 2, 2).unwrap();
        // 3
        tetgen.set_facet_point(3, 0, 1).unwrap();
        tetgen.set_facet_point(3, 1, 2).unwrap();
        tetgen.set_facet_point(3, 2, 3).unwrap();

        // set region
        tetgen.set_region(0, 1, 0.1, 0.9, 0.1, None).unwrap();

        // generate mesh
        tetgen.generate_mesh(false, false, None, None).unwrap();

        let mesh = Unstructured::from_tetgen(&tetgen);
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_from_tetgen.svg");
        }
        assert_eq!(mesh.ndim, 3);
        assert_eq!(mesh.points.len(), 5);
        assert_eq!(mesh.cells.len(), 2);
    }

    #[test]
    fn from_tetgen_works_o2() {
        // allocate data for 4 points
        let mut tetgen = Tetgen::new(4, Some(vec![3, 3, 3, 3]), Some(1), None).unwrap();

        // set points
        tetgen.set_point(0, 0, 0.0, 1.0, 0.0).unwrap();
        tetgen.set_point(1, 0, 0.0, 0.0, 0.0).unwrap();
        tetgen.set_point(2, 0, 1.0, 1.0, 0.0).unwrap();
        tetgen.set_point(3, 0, 0.0, 1.0, 1.0).unwrap();

        // set facets
        // 0
        tetgen.set_facet_point(0, 0, 0).unwrap();
        tetgen.set_facet_point(0, 1, 2).unwrap();
        tetgen.set_facet_point(0, 2, 1).unwrap();
        // 1
        tetgen.set_facet_point(1, 0, 0).unwrap();
        tetgen.set_facet_point(1, 1, 1).unwrap();
        tetgen.set_facet_point(1, 2, 3).unwrap();
        // 2
        tetgen.set_facet_point(2, 0, 0).unwrap();
        tetgen.set_facet_point(2, 1, 3).unwrap();
        tetgen.set_facet_point(2, 2, 2).unwrap();
        // 3
        tetgen.set_facet_point(3, 0, 1).unwrap();
        tetgen.set_facet_point(3, 1, 2).unwrap();
        tetgen.set_facet_point(3, 2, 3).unwrap();

        // set region
        tetgen.set_region(0, 1, 0.1, 0.9, 0.1, None).unwrap();

        // generate mesh
        tetgen.generate_mesh(false, true, None, None).unwrap();

        let mesh = Unstructured::from_tetgen(&tetgen);
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_from_tetgen_o2.svg");
        }
        assert_eq!(mesh.ndim, 3);
        assert_eq!(mesh.points.len(), 14);
        assert_eq!(mesh.cells.len(), 2);
    }

    #[test]
    fn from_trigen_works_o2() {
        let mut trigen = Trigen::new(3, Some(3), Some(1), None).unwrap();
        trigen.set_point(0, -100, 0.0, 0.0).unwrap();
        trigen.set_point(1, -200, 1.0, 0.0).unwrap();
        trigen.set_point(2, -300, 0.0, 1.0).unwrap();
        trigen.set_segment(0, -10, 0, 1).unwrap();
        trigen.set_segment(1, -20, 1, 2).unwrap();
        trigen.set_segment(2, -30, 2, 0).unwrap();
        trigen.set_region(0, 8, 0.1, 0.1, None).unwrap();
        trigen.generate_mesh(false, true, false, None, None).unwrap();
        let mesh = Unstructured::from_trigen(&trigen);
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_from_trigen_o2.svg");
        }
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 6);
        assert_eq!(mesh.cells.len(), 1);
        assert_eq!(mesh.points[0].id, 0);
        assert_eq!(mesh.points[0].marker, -100);
        assert_eq!(mesh.points[0].coords, vec![0.0, 0.0]);
        assert_eq!(mesh.points[1].id, 1);
        assert_eq!(mesh.points[1].marker, -200);
        assert_eq!(mesh.points[1].coords, vec![1.0, 0.0]);
        assert_eq!(mesh.points[2].id, 2);
        assert_eq!(mesh.points[2].marker, -300);
        assert_eq!(mesh.points[2].coords, vec![0.0, 1.0]);
        assert_eq!(mesh.points[3].id, 3);
        assert_eq!(mesh.points[3].marker, -10);
        assert_eq!(mesh.points[3].coords, vec![0.5, 0.0]);
        assert_eq!(mesh.points[4].id, 4);
        assert_eq!(mesh.points[4].marker, -20);
        assert_eq!(mesh.points[4].coords, vec![0.5, 0.5]);
        assert_eq!(mesh.points[5].id, 5);
        assert_eq!(mesh.points[5].marker, -30);
        assert_eq!(mesh.points[5].coords, vec![0.0, 0.5]);
        assert_eq!(mesh.cells[0].id, 0);
        assert_eq!(mesh.cells[0].marker, 8);
        assert_eq!(mesh.cells[0].kind, GeoKind::Tri6);
        assert_eq!(mesh.cells[0].points, vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn call_trigen_works_1() {
        let pslg = Samples::two_qua4();
        let holes = Vec::new();

        // convert two_qua4 to triangles
        let mesh_1 = Unstructured::call_trigen(&pslg, &holes, false, None, None, None, false).unwrap();
        mesh_1.check_all().unwrap();
        if SAVE_FIGURE {
            draw(&mesh_1, false, "/tmp/gemlab/test_call_trigen_works_1.svg");
        }
        let correct_mesh_1 = "# header\n\
             # ndim npoint ncell nmarked_edge nmarked_face\n\
             2 6 4 6 0\n\
             \n\
             # points\n\
             # id marker x y {z}\n\
             0 -1 0.0 0.0\n\
             1 -2 1.0 0.0\n\
             2 -3 1.0 1.0\n\
             3 -4 0.0 1.0\n\
             4 -5 2.0 0.0\n\
             5 -6 2.0 1.0\n\
             \n\
             # cells\n\
             # id marker kind points\n\
             0 1 tri3 0 1 3\n\
             1 1 tri3 3 1 2\n\
             2 2 tri3 2 4 5\n\
             3 2 tri3 4 2 1\n\
             \n\
             # marked edges\n\
             # marker p1 p2\n\
             -400 1 0\n\
             -300 0 3\n\
             -400 4 1\n\
             -100 3 2\n\
             -200 2 5\n\
             -300 5 4\n";
        assert_eq!(format!("{}", mesh_1), correct_mesh_1);

        // convert again (will be the same)
        let mesh_2 = Unstructured::call_trigen(&mesh_1, &holes, false, None, None, None, false).unwrap();
        if SAVE_FIGURE {
            draw(&mesh_2, false, "/tmp/gemlab/test_call_trigen_works_2.svg");
        }
        assert_eq!(format!("{}", mesh_2), correct_mesh_1);

        // now convert but use o2 elements
        let mesh_3 = Unstructured::call_trigen(&mesh_2, &holes, true, None, None, None, false).unwrap();
        if SAVE_FIGURE {
            draw(&mesh_3, false, "/tmp/gemlab/test_call_trigen_works_3.svg");
        }
        let correct_mesh_3 = "# header\n\
            # ndim npoint ncell nmarked_edge nmarked_face\n\
            2 15 4 6 0\n\
            \n\
            # points\n\
            # id marker x y {z}\n\
            0 -1 0.0 0.0\n\
            1 -2 1.0 0.0\n\
            2 -3 1.0 1.0\n\
            3 -4 0.0 1.0\n\
            4 -5 2.0 0.0\n\
            5 -6 2.0 1.0\n\
            6 -400 0.5 0.0\n\
            7 0 0.5 0.5\n\
            8 -300 0.0 0.5\n\
            9 0 1.0 0.5\n\
            10 -100 0.5 1.0\n\
            11 0 1.5 0.5\n\
            12 -300 2.0 0.5\n\
            13 -200 1.5 1.0\n\
            14 -400 1.5 0.0\n\
            \n\
            # cells\n\
            # id marker kind points\n\
            0 1 tri6 0 1 3 6 7 8\n\
            1 1 tri6 3 1 2 7 9 10\n\
            2 2 tri6 2 4 5 11 12 13\n\
            3 2 tri6 4 2 1 11 9 14\n\
            \n\
            # marked edges\n\
            # marker p1 p2\n\
            -400 1 0\n\
            -300 0 3\n\
            -400 4 1\n\
            -100 3 2\n\
            -200 2 5\n\
            -300 5 4\n";
        assert_eq!(format!("{}", mesh_3), correct_mesh_3);
    }

    #[test]
    fn call_trigen_works_4() {
        let pslg = Samples::two_qua4();
        let holes = Vec::new();

        // generate again with more elements
        let max_areas = HashMap::from([(1, 1.0), (2, 0.1)]);
        let mesh = Unstructured::call_trigen(&pslg, &holes, false, Some(max_areas), None, None, false).unwrap();
        assert_eq!(mesh.points.len(), 15);
        assert_eq!(mesh.cells.len(), 19);
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_call_trigen_works_4.svg");
        }

        // generate again with more elements
        let mesh = Unstructured::call_trigen(&pslg, &holes, false, None, Some(0.1), None, false).unwrap();
        assert_eq!(mesh.points.len(), 22);
        assert_eq!(mesh.cells.len(), 30);
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_call_trigen_works_5.svg");
        }
    }

    #[test]
    fn call_tetgen_works_1() {
        let plc = Samples::two_hex8();
        let regions = vec![(1, 0.5, 0.5, 0.5), (2, 0.5, 0.5, 1.5)];
        let holes = Vec::new();
        let mesh = Unstructured::call_tetgen(&plc, &regions, &holes, false, None, None, None, false).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, true, "/tmp/gemlab/test_call_tetgen_works_1.svg");
        }
        assert_eq!(mesh.ndim, 3);
        assert_eq!(mesh.points.len(), 12);
        assert_eq!(mesh.cells.len(), 12);
        println!("{}", mesh);
        let correct = "# header\n\
            # ndim npoint ncell nmarked_edge nmarked_face\n\
            3 12 12 0 8\n\
            \n\
            # points\n\
            # id marker x y {z}\n\
            0 0 0.0 0.0 0.0\n\
            1 -1 1.0 0.0 0.0\n\
            2 -1 1.0 1.0 0.0\n\
            3 0 0.0 1.0 0.0\n\
            4 0 0.0 0.0 1.0\n\
            5 0 1.0 0.0 1.0\n\
            6 0 1.0 1.0 1.0\n\
            7 -1 0.0 1.0 1.0\n\
            8 0 0.0 0.0 2.0\n\
            9 0 1.0 0.0 2.0\n\
            10 0 1.0 1.0 2.0\n\
            11 0 0.0 1.0 2.0\n\
            \n\
            # cells\n\
            # id marker kind points\n\
            0 1 tet4 0 3 7 2\n\
            1 1 tet4 0 2 6 1\n\
            2 1 tet4 0 6 5 1\n\
            3 2 tet4 8 9 4 10\n\
            4 2 tet4 8 4 11 10\n\
            5 1 tet4 0 7 4 6\n\
            6 1 tet4 7 0 2 6\n\
            7 1 tet4 4 0 6 5\n\
            8 2 tet4 4 7 11 6\n\
            9 2 tet4 4 11 10 6\n\
            10 2 tet4 9 4 10 5\n\
            11 2 tet4 10 4 6 5\n\
            \n\
            # marked faces\n\
            # marker p1 p2 p3 {p4}\n\
            -8 2 7 3\n\
            -10 1 6 2\n\
            -10 1 5 6\n\
            -11 1 0 5\n\
            -9 9 8 10\n\
            -9 10 8 11\n\
            -8 6 7 2\n\
            -11 0 4 5\n";
        assert_eq!(format!("{}", mesh), correct);
    }

    #[test]
    fn call_tetgen_works_2() {
        let plc = Samples::one_hex8();
        let regions = vec![(1, 0.5, 0.5, 0.5)];
        let holes = Vec::new();
        let max_vol = HashMap::from([(1, 0.3)]);
        let mesh = Unstructured::call_tetgen(&plc, &regions, &holes, false, Some(max_vol), None, None, false).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, true, "/tmp/gemlab/test_call_tetgen_works_2.svg");
        }
        assert_eq!(mesh.ndim, 3);
        assert_eq!(mesh.points.len(), 14);
        assert_eq!(mesh.cells.len(), 24);
    }

    #[test]
    fn tri_quarter_ring_2d_captures_errors() {
        assert_eq!(
            Unstructured::quarter_ring_2d(RMIN, RMAX, 0, 4, GeoKind::Tri3, None, false).err(),
            Some("number of divisions along the radius must be ≥ 1")
        );
        assert_eq!(
            Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 0, GeoKind::Tri3, None, false).err(),
            Some("number of divisions along alpha must be ≥ 1")
        );
        assert_eq!(
            Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 2, GeoKind::Qua4, None, false).err(),
            Some("the GeoClass of target must be Tri")
        );
    }

    #[test]
    fn tri_quarter_ring_2d_works() {
        let mut mesh = Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 4, GeoKind::Tri3, None, false).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_tri_quarter_ring_2d.svg");
        }
        assert_eq!(mesh.points.len(), 14);
        assert_eq!(mesh.cells.len(), 14);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.01).unwrap();
        check_corner_markers(&mesh, false);
        check_point_markers(&mesh, &[11, 10, 9], &[3, 4, 5], &[1], &[7]);
        check_constraints(&mesh);
        for p in [0, 11, 10, 9, 8] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for p in [2, 3, 4, 5, 6] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
        if false {
            print_bandwidth(&mut mesh);
        }
    }

    #[test]
    fn tri_quarter_ring_2d_o2_works() {
        let mut mesh = Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 4, GeoKind::Tri6, None, false).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_tri_quarter_ring_2d_o2.svg");
        }
        assert_eq!(mesh.points.len(), 41);
        assert_eq!(mesh.cells.len(), 14);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.01).unwrap();
        check_corner_markers(&mesh, false);
        check_point_markers(
            &mesh,
            &[32, 11, 36, 10, 21, 9, 22],
            &[34, 3, 40, 4, 17, 5, 26],
            &[30, 1, 33],
            &[24, 7, 27],
        );
        check_constraints(&mesh);
        for p in [0, 32, 11, 36, 10, 21, 9, 22, 8] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for p in [2, 34, 3, 40, 4, 17, 5, 26, 6] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
        if false {
            print_bandwidth(&mut mesh);
        }
    }

    #[test]
    fn tri_quarter_ring_2d_global_max_area_works() {
        let global_max_area = Some(0.4);
        let mut mesh = Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 4, GeoKind::Tri3, global_max_area, false).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, true, "/tmp/gemlab/test_tri_quarter_ring_2d_global_max_area.svg");
        }
        assert_eq!(mesh.points.len(), 50);
        assert_eq!(mesh.cells.len(), 78);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.01).unwrap();
        check_corner_markers(&mesh, false);
        check_point_markers(
            &mesh,
            &[11, 10, 17, 9],
            &[40, 3, 25, 4, 24, 5, 33],
            &[14, 1, 42],
            &[13, 7],
        );
        check_constraints(&mesh);
        for p in [0, 11, 10, 17, 9, 8] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for p in [2, 40, 3, 25, 4, 24, 5, 33, 6] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
        if false {
            print_bandwidth(&mut mesh);
        }
    }

    #[test]
    fn tri_quarter_ring_2d_o2_global_max_area_works() {
        let global_max_area = Some(0.4);
        let mesh = Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 4, GeoKind::Tri6, global_max_area, false).unwrap();
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
        mesh.check_overlapping_points(0.01).unwrap();
        check_corner_markers(&mesh, false);
        check_point_markers(
            &mesh,
            &[150, 11, 81, 10, 58, 17, 99, 9, 146],
            &[154, 40, 176, 3, 174, 25, 130, 4, 118, 24, 136, 5, 124, 33, 134],
            &[89, 14, 78, 1, 167, 42, 163],
            &[67, 13, 140, 7, 145],
        );
        check_constraints(&mesh);
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
        let mut mesh = Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 4, GeoKind::Tri10, None, false).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, false, "/tmp/gemlab/test_tri_quarter_ring_2d_tri10.svg");
        }
        assert_eq!(mesh.points.len(), 82);
        assert_eq!(mesh.cells.len(), 14);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.01).unwrap();
        check_corner_markers(&mesh, false);
        check_point_markers(
            &mesh,
            &[57, 54, 51, 70, 68, 2, 23, 21, 19, 30, 27],
            &[62, 65, 60, 79, 80, 10, 12, 15, 11, 36, 39],
            &[52, 55, 50, 61, 64],
            &[32, 29, 26, 40, 37],
        );
        check_constraints(&mesh);
        for p in [49, 57, 54, 51, 70, 68, 2, 23, 21, 19, 30, 27, 25] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for p in [59, 62, 65, 60, 79, 80, 10, 12, 15, 11, 36, 39, 34] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
        if false {
            print_bandwidth(&mut mesh);
        }
    }

    #[test]
    fn tri_quarter_ring_2d_tri15_works() {
        let mesh = Unstructured::quarter_ring_2d(RMIN, RMAX, 2, 4, GeoKind::Tri15, None, false).unwrap();
        if SAVE_FIGURE {
            draw(&mesh, true, "/tmp/gemlab/test_tri_quarter_ring_2d_tri15.svg");
        }
        assert_eq!(mesh.points.len(), 137);
        assert_eq!(mesh.cells.len(), 14);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.01).unwrap();
        check_corner_markers(&mesh, false);
        check_point_markers(
            &mesh,
            &[92, 86, 91, 83, 113, 110, 112, 2, 35, 31, 34, 29, 45, 41, 44],
            &[103, 99, 104, 97, 132, 131, 133, 15, 20, 17, 21, 16, 59, 55, 60],
            &[87, 84, 88, 82, 101, 98, 102],
            &[49, 43, 48, 40, 62, 56, 61],
        );
        check_constraints(&mesh);
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
    fn tri_quarter_ring_3d_handles_errors() {
        assert_eq!(
            Unstructured::quarter_ring_3d(RMIN, RMAX, 1.0, 0, 4, GeoKind::Tet4, None, false).err(),
            Some("number of divisions along the radius must be ≥ 1")
        );
        assert_eq!(
            Unstructured::quarter_ring_3d(RMIN, RMAX, 1.0, 2, 0, GeoKind::Tet4, None, false).err(),
            Some("number of divisions along alpha must be ≥ 1")
        );
        assert_eq!(
            Unstructured::quarter_ring_3d(RMIN, RMAX, 1.0, 2, 2, GeoKind::Qua4, None, false).err(),
            Some("the GeoClass of target must be Tet")
        );
        assert_eq!(
            Unstructured::quarter_ring_3d(RMIN, RMAX, 1.0, 2, 4, GeoKind::Tet20, None, false).err(),
            Some("only Tet4 and Tet10 are available currently")
        )
    }

    fn draw_ring_3d_with_cylin(mesh: &Mesh, filename: &str) {
        let mut cylin_in = Surface::new();
        let mut cylin_out = Surface::new();
        cylin_in.set_surf_color("#ff000020");
        cylin_out.set_surf_color("#ff000020");
        cylin_in
            .draw_cylinder(&[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], RMIN, 5, 81)
            .unwrap();
        cylin_out
            .draw_cylinder(&[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], RMAX, 5, 81)
            .unwrap();

        let mut draw = Draw::new();
        draw.set_size(800.0, 800.0).show_point_ids(true).show_point_dots(true);
        draw.extra(|plot, before| {
            if before {
                plot.add(&cylin_in).add(&cylin_out);
            }
        })
        .all(&mesh, filename)
        .unwrap();
    }

    #[test]
    fn tri_quarter_ring_3d_works() {
        let mesh = Unstructured::quarter_ring_3d(RMIN, RMAX, 1.0, 2, 4, GeoKind::Tet4, None, false).unwrap();
        if SAVE_FIGURE {
            draw_ring_3d_with_cylin(&mesh, "/tmp/gemlab/test_tri_quarter_ring_3d.svg");
        }
        assert_eq!(mesh.points.len(), 24);
        assert_eq!(mesh.cells.len(), 30);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.01).unwrap();
        check_corner_markers(&mesh, false);
        check_point_markers(
            &mesh,
            &[9, 10, 11, 21, 22, 23],
            &[3, 4, 5, 15, 16, 17],
            &[1, 13],
            &[7, 19],
        );
        check_constraints(&mesh);
        for p in [0, 11, 10, 9, 8, 20, 21, 22, 23, 12] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMIN, 1e-15);
        }
        for p in [2, 3, 4, 5, 6, 14, 15, 16, 17, 18] {
            let d = point_point_distance(&mesh.points[p].coords[0..2], &[0.0, 0.0]).unwrap();
            approx_eq(d, RMAX, 1e-15);
        }
    }

    #[test]
    fn tri_quarter_ring_3d_o2_works() {
        let mesh = Unstructured::quarter_ring_3d(RMIN, RMAX, 1.0, 2, 4, GeoKind::Tet10, None, false).unwrap();
        if SAVE_FIGURE {
            draw_ring_3d_with_cylin(&mesh, "/tmp/gemlab/test_tri_quarter_ring_3d_o2.svg");
        }
        assert_eq!(mesh.points.len(), 99);
        assert_eq!(mesh.cells.len(), 30);
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.01).unwrap();
        check_corner_markers(&mesh, true);
        check_point_markers(
            &mesh,
            &[
                9, 10, 11, 21, 22, 23, 64, 45, 30, 62, 58, 57, 67, 60, 98, 97, 96, 72, 74, 28, 76,
            ],
            &[
                3, 4, 5, 15, 16, 17, 85, 83, 90, 51, 69, 50, 55, 71, 79, 80, 81, 49, 36, 43, 47,
            ],
            &[1, 13, 24, 25, 77, 29, 84, 87, 88],
            &[7, 19, 44, 34, 31, 32, 37, 39, 42],
        );
        for p in &[73, 82, 95, 93, 78, 59, 61, 40, 63] {
            assert_eq!(mesh.points[*p].marker, MARKER_ZMIN);
        }
        for p in &[26, 89, 92, 52, 53, 68, 65, 48, 35] {
            assert_eq!(mesh.points[*p].marker, MARKER_ZMAX);
        }
        check_constraints(&mesh);
    }

    #[test]
    fn tri_quarter_ring_3d_o2_max_vol_works() {
        let global_max_volume = Some(0.5);
        let mut mesh =
            Unstructured::quarter_ring_3d(RMIN, RMAX, 1.0, 2, 4, GeoKind::Tet10, global_max_volume, false).unwrap();
        if SAVE_FIGURE {
            draw_ring_3d_with_cylin(&mesh, "/tmp/gemlab/test_tri_quarter_ring_3d_o2_max_vol.svg");
        }
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.01).unwrap();
        check_constraints(&mesh);
        if false {
            print_bandwidth(&mut mesh);
        }
    }
}
