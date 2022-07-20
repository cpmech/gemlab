use super::{Features, Mesh, PointId, Region};
use crate::shapes::{op, GeoKind, Scratchpad};
use crate::util::ONE_BY_3;
use crate::StrError;
use plotpy::{Canvas, Curve, Plot, PolyCode, Text};
use russell_lab::Vector;
use std::collections::HashMap;

/// Adds (Polyline or Bezier) curve to the canvas
fn add_curve(
    canvas: &mut Canvas,
    mesh: &Mesh,
    lin_kind: GeoKind,
    lin_points: &[PointId],
    begin: bool,
    end: bool,
    pads: &mut HashMap<GeoKind, Scratchpad>,
) -> Result<(), StrError> {
    let pp = &mesh.points;
    let ii = lin_points;
    let mut x = Vector::new(mesh.ndim);
    match lin_kind {
        GeoKind::Lin2 => {
            if mesh.ndim == 2 {
                if begin {
                    canvas.polycurve_begin();
                    canvas.polycurve_add(pp[ii[0]].coords[0], pp[ii[0]].coords[1], PolyCode::MoveTo);
                }
                canvas.polycurve_add(pp[ii[1]].coords[0], pp[ii[1]].coords[1], PolyCode::LineTo);
                if end {
                    canvas.polycurve_end(false);
                }
            } else {
                if begin {
                    canvas.polyline_3d_begin();
                    canvas.polyline_3d_add(pp[ii[0]].coords[0], pp[ii[0]].coords[1], pp[ii[0]].coords[2]);
                }
                canvas.polyline_3d_add(pp[ii[1]].coords[0], pp[ii[1]].coords[1], pp[ii[1]].coords[2]);
                if end {
                    canvas.polyline_3d_end();
                }
            }
        }
        GeoKind::Lin3 => {
            if mesh.ndim == 2 {
                // add poly-curve points or control points (q..) for quadratic Bezier
                // (see file bezier-curves-math.pdf under data/derivations)
                let (xa, ya) = (pp[ii[0]].coords[0], pp[ii[0]].coords[1]);
                let (xb, yb) = (pp[ii[1]].coords[0], pp[ii[1]].coords[1]);
                let (xc, yc) = (pp[ii[2]].coords[0], pp[ii[2]].coords[1]);
                let qx = (-xa - xb + 4.0 * xc) / 2.0;
                let qy = (-ya - yb + 4.0 * yc) / 2.0;
                if begin {
                    canvas.polycurve_begin();
                    canvas.polycurve_add(xa, ya, PolyCode::MoveTo);
                }
                canvas
                    .polycurve_add(qx, qy, PolyCode::Curve3)
                    .polycurve_add(xb, yb, PolyCode::Curve3);
                if end {
                    canvas.polycurve_end(false);
                }
            } else {
                if begin {
                    canvas.polyline_3d_begin();
                    canvas.polyline_3d_add(pp[ii[0]].coords[0], pp[ii[0]].coords[1], pp[ii[0]].coords[2]);
                }
                canvas
                    .polyline_3d_add(pp[ii[2]].coords[0], pp[ii[2]].coords[1], pp[ii[2]].coords[2])
                    .polyline_3d_add(pp[ii[1]].coords[0], pp[ii[1]].coords[1], pp[ii[1]].coords[2]);
                if end {
                    canvas.polyline_3d_end();
                }
            }
        }
        GeoKind::Lin4 => {
            if mesh.ndim == 2 {
                // add poly-curve points or control points (q..) for cubic Bezier
                // (see file bezier-curves-math.pdf under data/derivations)
                let (xa, ya) = (pp[ii[0]].coords[0], pp[ii[0]].coords[1]);
                let (xb, yb) = (pp[ii[1]].coords[0], pp[ii[1]].coords[1]);
                let (xc, yc) = (pp[ii[2]].coords[0], pp[ii[2]].coords[1]);
                let (xd, yd) = (pp[ii[3]].coords[0], pp[ii[3]].coords[1]);
                let qcx = (-5.0 * xa + 2.0 * xb + 18.0 * xc - 9.0 * xd) / 6.0;
                let qcy = (-5.0 * ya + 2.0 * yb + 18.0 * yc - 9.0 * yd) / 6.0;
                let qdx = (2.0 * xa - 5.0 * xb - 9.0 * xc + 18.0 * xd) / 6.0;
                let qdy = (2.0 * ya - 5.0 * yb - 9.0 * yc + 18.0 * yd) / 6.0;
                if begin {
                    canvas.polycurve_begin();
                    canvas.polycurve_add(xa, ya, PolyCode::MoveTo);
                }
                canvas
                    .polycurve_add(qcx, qcy, PolyCode::Curve4)
                    .polycurve_add(qdx, qdy, PolyCode::Curve4)
                    .polycurve_add(xb, yb, PolyCode::Curve4);
                if end {
                    canvas.polycurve_end(false);
                }
            } else {
                if begin {
                    canvas.polyline_3d_begin();
                    canvas.polyline_3d_add(pp[ii[0]].coords[0], pp[ii[0]].coords[1], pp[ii[0]].coords[2]);
                }
                canvas
                    .polyline_3d_add(pp[ii[2]].coords[0], pp[ii[2]].coords[1], pp[ii[2]].coords[2])
                    .polyline_3d_add(pp[ii[3]].coords[0], pp[ii[3]].coords[1], pp[ii[3]].coords[2])
                    .polyline_3d_add(pp[ii[1]].coords[0], pp[ii[1]].coords[1], pp[ii[1]].coords[2]);
                if end {
                    canvas.polyline_3d_end();
                }
            }
        }
        GeoKind::Lin5 => {
            if mesh.ndim == 2 {
                // add poly-curve points or control points (q..) for cubic Bezier
                // (see file bezier-curves-math.pdf under data/derivations)
                // in this case, we have to interpolate the Lin5 to make it a Lin4
                let mut pad = pads.entry(lin_kind).or_insert(Scratchpad::new(mesh.ndim, lin_kind)?);
                for m in 0..lin_points.len() {
                    for j in 0..mesh.ndim {
                        pad.set_xx(m, j, pp[ii[m]].coords[j]);
                    }
                }
                let (xa, ya) = (pp[ii[0]].coords[0], pp[ii[0]].coords[1]);
                let (xb, yb) = (pp[ii[1]].coords[0], pp[ii[1]].coords[1]);
                op::calc_coords(&mut x, &mut pad, &[-ONE_BY_3, 0.0]).unwrap();
                let (xc, yc) = (x[0], x[1]);
                op::calc_coords(&mut x, &mut pad, &[ONE_BY_3, 0.0]).unwrap();
                let (xd, yd) = (x[0], x[1]);
                let qcx = (-5.0 * xa + 2.0 * xb + 18.0 * xc - 9.0 * xd) / 6.0;
                let qcy = (-5.0 * ya + 2.0 * yb + 18.0 * yc - 9.0 * yd) / 6.0;
                let qdx = (2.0 * xa - 5.0 * xb - 9.0 * xc + 18.0 * xd) / 6.0;
                let qdy = (2.0 * ya - 5.0 * yb - 9.0 * yc + 18.0 * yd) / 6.0;
                if begin {
                    canvas.polycurve_begin();
                    canvas.polycurve_add(xa, ya, PolyCode::MoveTo);
                }
                canvas
                    .polycurve_add(qcx, qcy, PolyCode::Curve4)
                    .polycurve_add(qdx, qdy, PolyCode::Curve4)
                    .polycurve_add(xb, yb, PolyCode::Curve4);
                if end {
                    canvas.polycurve_end(false);
                }
            } else {
                if begin {
                    canvas.polyline_3d_begin();
                    canvas.polyline_3d_add(pp[ii[0]].coords[0], pp[ii[0]].coords[1], pp[ii[0]].coords[2]);
                }
                canvas
                    .polyline_3d_add(pp[ii[3]].coords[0], pp[ii[3]].coords[1], pp[ii[3]].coords[2])
                    .polyline_3d_add(pp[ii[2]].coords[0], pp[ii[2]].coords[1], pp[ii[2]].coords[2])
                    .polyline_3d_add(pp[ii[4]].coords[0], pp[ii[4]].coords[1], pp[ii[4]].coords[2])
                    .polyline_3d_add(pp[ii[1]].coords[0], pp[ii[1]].coords[1], pp[ii[1]].coords[2]);
                if end {
                    canvas.polyline_3d_end();
                }
            }
        }
        _ => return Err("lin_kind is not Lin"),
    }
    Ok(())
}

/// Draws cell
#[rustfmt::skip]
fn draw_cell(
    canvas: &mut Canvas,
    mesh: &Mesh,
    cell_kind: GeoKind,
    cell_points: &Vec<PointId>,
    pads: &mut HashMap<GeoKind, Scratchpad>,
) -> Result<(), StrError> {
    let ii = cell_points;
    match cell_kind {
        // Lin
        GeoKind::Lin2 => {
            add_curve(canvas, mesh, cell_kind, cell_points, true, true, pads)?;
        }
        GeoKind::Lin3 => {
            add_curve(canvas, mesh, cell_kind, cell_points, true, true, pads)?;
        }
        GeoKind::Lin4 => {
            add_curve(canvas, mesh, cell_kind, cell_points, true, true, pads)?;
        }
        GeoKind::Lin5 => {
            add_curve(canvas, mesh, cell_kind, cell_points, true, true, pads)?;
        }
        // Tri
        GeoKind::Tri3 => {
            add_curve(canvas, mesh, GeoKind::Lin2, &[ii[0], ii[1]], true, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin2, &[ii[1], ii[2]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin2, &[ii[2], ii[0]], false, true, pads)?;
        }
        GeoKind::Tri6 => {
            add_curve(canvas, mesh, GeoKind::Lin3, &[ii[0], ii[1], ii[3]], true, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin3, &[ii[1], ii[2], ii[4]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin3, &[ii[2], ii[0], ii[5]], false, true, pads)?;
        }
        GeoKind::Tri10 => {
            add_curve(canvas, mesh, GeoKind::Lin4, &[ii[0], ii[1], ii[3], ii[6]], true, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin4, &[ii[1], ii[2], ii[4], ii[7]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin4, &[ii[2], ii[0], ii[5], ii[8]], false, true, pads)?;
        }
        GeoKind::Tri15 => {
            add_curve(canvas, mesh, GeoKind::Lin5, &[ii[0], ii[1], ii[3], ii[6], ii[7]], true, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin5, &[ii[1], ii[2], ii[4], ii[8], ii[9]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin5, &[ii[2], ii[0], ii[5], ii[10], ii[11]], false, true, pads)?;
        }
        // Qua
        GeoKind::Qua4 => {
            add_curve(canvas, mesh, GeoKind::Lin2, &[ii[0], ii[1]], true, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin2, &[ii[1], ii[2]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin2, &[ii[2], ii[3]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin2, &[ii[3], ii[0]], false, true, pads)?;
        }
        GeoKind::Qua8 => {
            add_curve(canvas, mesh, GeoKind::Lin3, &[ii[0], ii[1], ii[4]], true, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin3, &[ii[1], ii[2], ii[5]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin3, &[ii[2], ii[3], ii[6]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin3, &[ii[3], ii[0], ii[7]], false, true, pads)?;
        }
        GeoKind::Qua9 => {
            add_curve(canvas, mesh, GeoKind::Lin3, &[ii[0], ii[1], ii[4]], true, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin3, &[ii[1], ii[2], ii[5]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin3, &[ii[2], ii[3], ii[6]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin3, &[ii[3], ii[0], ii[7]], false, true, pads)?;
        }
        GeoKind::Qua12 => {
            add_curve(canvas, mesh, GeoKind::Lin4, &[ii[0], ii[1], ii[4], ii[8]], true, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin4, &[ii[1], ii[2], ii[5], ii[9]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin4, &[ii[2], ii[3], ii[6], ii[10]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin4, &[ii[3], ii[0], ii[7], ii[11]], false, true, pads)?;
        }
        GeoKind::Qua16 => {
            add_curve(canvas, mesh, GeoKind::Lin4, &[ii[0], ii[1], ii[4], ii[8]], true, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin4, &[ii[1], ii[2], ii[5], ii[9]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin4, &[ii[2], ii[3], ii[6], ii[10]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin4, &[ii[3], ii[0], ii[7], ii[11]], false, true, pads)?;
        }
        GeoKind::Qua17 => {
            add_curve(canvas, mesh, GeoKind::Lin5, &[ii[0], ii[1], ii[4], ii[9], ii[10]], true, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin5, &[ii[1], ii[2], ii[5], ii[11], ii[12]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin5, &[ii[2], ii[3], ii[6], ii[13], ii[14]], false, false, pads)?;
            add_curve(canvas, mesh, GeoKind::Lin5, &[ii[3], ii[0], ii[7], ii[15], ii[16]], false, true, pads)?;
        }
        // Tet
        GeoKind::Tet4 => {
            for e in 0..cell_kind.nedge() {
                let a = cell_kind.edge_node_id(e, 0);
                let b = cell_kind.edge_node_id(e, 1);
                add_curve(canvas, mesh, GeoKind::Lin2, &[ii[a], ii[b]], true, true, pads)?;
            }
        }
        GeoKind::Tet10 => {
            for e in 0..cell_kind.nedge() {
                let a = cell_kind.edge_node_id(e, 0);
                let b = cell_kind.edge_node_id(e, 1);
                let c = cell_kind.edge_node_id(e, 2);
                add_curve(canvas, mesh, GeoKind::Lin2, &[ii[a], ii[b], ii[c]], true, true, pads)?;
            }
        }
        GeoKind::Tet20 => {
            for e in 0..cell_kind.nedge() {
                let a = cell_kind.edge_node_id(e, 0);
                let b = cell_kind.edge_node_id(e, 1);
                let c = cell_kind.edge_node_id(e, 2);
                let d = cell_kind.edge_node_id(e, 3);
                add_curve(canvas, mesh, GeoKind::Lin3, &[ii[a], ii[b], ii[c], ii[d]], true, true, pads)?;
            }
        }
        // Hex
        GeoKind::Hex8 => {
            for e in 0..cell_kind.nedge() {
                let a = cell_kind.edge_node_id(e, 0);
                let b = cell_kind.edge_node_id(e, 1);
                add_curve(canvas, mesh, GeoKind::Lin2, &[ii[a], ii[b]], true, true, pads)?;
            }
        }
        GeoKind::Hex20 => {
            for e in 0..cell_kind.nedge() {
                let a = cell_kind.edge_node_id(e, 0);
                let b = cell_kind.edge_node_id(e, 1);
                let c = cell_kind.edge_node_id(e, 2);
                add_curve(canvas, mesh, GeoKind::Lin2, &[ii[a], ii[b], ii[c]], true, true, pads)?;
            }
        }
        GeoKind::Hex32 => {
            for e in 0..cell_kind.nedge() {
                let a = cell_kind.edge_node_id(e, 0);
                let b = cell_kind.edge_node_id(e, 1);
                let c = cell_kind.edge_node_id(e, 2);
                let d = cell_kind.edge_node_id(e, 3);
                add_curve(canvas, mesh, GeoKind::Lin3, &[ii[a], ii[b], ii[c], ii[d]], true, true, pads)?;
            }
        }
    }
    Ok(())
}

/// Implements functions to draw cells, edges and faces
pub struct Draw {
    /// Canvas to draw edges
    pub canvas_edges: Canvas,

    /// Canvas to draw points (markers)
    pub canvas_points: Curve,

    /// Canvas to draw point ids
    pub canvas_point_ids: Text,

    /// Canvas to draw cell ids
    pub canvas_cell_ids: Text,

    /// Canvas to draw cells
    pub canvas_cells: Canvas,

    /// Canvas to draw lin cells
    pub canvas_lin_cells: Canvas,
}

impl Draw {
    /// Allocates a new instance
    pub fn new() -> Self {
        let mut canvas_edges = Canvas::new();
        let mut canvas_points = Curve::new();
        let mut canvas_point_ids = Text::new();
        let mut canvas_cell_ids = Text::new();
        let mut canvas_cells = Canvas::new();
        let mut canvas_lin_cells = Canvas::new();
        canvas_edges
            .set_stop_clip(true)
            .set_face_color("None")
            .set_line_width(2.0)
            .set_edge_color("#2440cd");
        canvas_points
            .set_stop_clip(true)
            .set_marker_color("black")
            .set_marker_line_color("white")
            .set_marker_style("o")
            .set_line_style("None");
        canvas_point_ids
            .set_color("red")
            .set_fontsize(8.0)
            .set_align_horizontal("center")
            .set_align_vertical("center")
            .set_bbox(true)
            .set_bbox_facecolor("white")
            .set_bbox_edgecolor("None")
            .set_bbox_style("circle,pad=0.15");
        canvas_cell_ids
            .set_color("#22971f")
            .set_fontsize(9.0)
            .set_align_horizontal("center")
            .set_align_vertical("center")
            .set_bbox(true)
            .set_bbox_facecolor("white")
            .set_bbox_edgecolor("#b7b7b7")
            .set_bbox_style("square,pad=0.15");
        canvas_cells
            .set_stop_clip(true)
            .set_face_color("#e3f3ff")
            .set_edge_color("#0055d4")
            .set_line_width(1.5);
        canvas_lin_cells
            .set_stop_clip(true)
            .set_face_color("None")
            .set_edge_color("#cd0000")
            .set_line_width(3.0);
        Draw {
            canvas_edges,
            canvas_points,
            canvas_point_ids,
            canvas_cell_ids,
            canvas_cells,
            canvas_lin_cells,
        }
    }

    /// Draws all points (markers)
    ///
    /// # Input
    ///
    /// * `plot` -- the plot instance (to be updated)
    /// * `mesh` -- the mesh
    pub fn points(&mut self, plot: &mut Plot, mesh: &Mesh) {
        if mesh.ndim == 2 {
            self.canvas_points.points_begin();
            mesh.points.iter().for_each(|point| {
                self.canvas_points.points_add(point.coords[0], point.coords[1]);
            });
            self.canvas_points.points_end();
        } else {
            self.canvas_points.points_3d_begin();
            mesh.points.iter().for_each(|point| {
                self.canvas_points
                    .points_3d_add(point.coords[0], point.coords[1], point.coords[2]);
            });
            self.canvas_points.points_3d_end();
        }
        plot.add(&self.canvas_points);
    }

    /// Draws all point ids (labels)
    ///
    /// # Input
    ///
    /// * `plot` -- the plot instance (to be updated)
    /// * `mesh` -- the mesh
    pub fn point_ids(&mut self, plot: &mut Plot, mesh: &Mesh) {
        if mesh.ndim == 2 {
            mesh.points.iter().for_each(|point| {
                self.canvas_point_ids
                    .draw(point.coords[0], point.coords[1], format!("{}", point.id).as_str());
            });
        } else {
            mesh.points.iter().for_each(|point| {
                self.canvas_point_ids.draw_3d(
                    point.coords[0],
                    point.coords[1],
                    point.coords[2],
                    format!("{}", point.id).as_str(),
                );
            });
        }
        plot.add(&self.canvas_point_ids);
    }

    /// Draw cells
    pub fn cells(&mut self, plot: &mut Plot, mesh: &Mesh, set_range: bool) -> Result<(), StrError> {
        // limits
        let mut xmin = vec![f64::MAX; mesh.ndim];
        let mut xmax = vec![f64::MIN; mesh.ndim];

        // memoization of scratchpads
        let mut pads = HashMap::new();

        // loop over cells
        for cell in &mesh.cells {
            let canvas = if cell.kind.is_lin() {
                &mut self.canvas_lin_cells
            } else {
                &mut self.canvas_cells
            };
            draw_cell(canvas, &mesh, cell.kind, &cell.points, &mut pads)?;
            for point_id in &cell.points {
                for j in 0..mesh.ndim {
                    xmin[j] = f64::min(xmin[j], mesh.points[*point_id].coords[j]);
                    xmax[j] = f64::max(xmax[j], mesh.points[*point_id].coords[j]);
                }
            }
        }

        // add to plot
        plot.add(&self.canvas_cells);
        plot.add(&self.canvas_lin_cells);
        if set_range {
            plot.set_range(xmin[0], xmax[0], xmin[1], xmax[1]);
        }
        Ok(())
    }

    /// Draws ids and attributes of cells
    ///
    /// # Input
    ///
    /// * `plot` -- the plot instance (to be updated)
    /// * `mesh` -- the mesh
    pub fn cell_ids(&mut self, plot: &mut Plot, mesh: &Mesh) -> Result<(), StrError> {
        // auxiliary
        let mut x = Vector::new(mesh.ndim);

        // loop over all cells
        for cell_id in 0..mesh.cells.len() {
            // compute coordinates at the center of cell
            let cell = &mesh.cells[cell_id];
            x.fill(0.0);
            for m in 0..cell.points.len() {
                for i in 0..mesh.ndim {
                    x[i] += mesh.points[cell.points[m]].coords[i];
                }
            }
            for i in 0..mesh.ndim {
                x[i] /= cell.points.len() as f64;
            }

            // add label
            let msg = format!("{} a {}", cell.id, cell.attribute_id);
            if mesh.ndim == 2 {
                self.canvas_cell_ids.draw(x[0], x[1], msg.as_str());
            } else {
                self.canvas_cell_ids.draw_3d(x[0], x[1], x[2], msg.as_str());
            }
        }

        // add to plot
        plot.add(&self.canvas_cell_ids);
        Ok(())
    }

    /// Draws edges
    ///
    /// # Input
    ///
    /// * `plot` -- the plot instance (to be updated)
    /// * `region` -- the region containing the mesh and features
    /// * `set_range` -- sets the range of `plot` to the limits of region; otherwise do not modifies the range/limits.
    pub fn edges(&mut self, plot: &mut Plot, region: &Region, set_range: bool) -> Result<(), StrError> {
        if region.mesh.ndim == 2 {
            self.edges_2d(plot, &region.mesh, &region.features, set_range)
        } else {
            self.edges_3d(plot, &region.mesh, &region.features, set_range)
        }
    }

    /// Draws 2D edges
    fn edges_2d(&mut self, plot: &mut Plot, mesh: &Mesh, features: &Features, set_range: bool) -> Result<(), StrError> {
        // space dimension
        let space_ndim = mesh.ndim;
        assert_eq!(space_ndim, 2);

        // middle points (x..) on Bezier curve and control points (q..)
        let mut xc = Vector::new(space_ndim); // first middle point on Bezier curve with t=⅓(ξ=-⅓) or t=½(ξ=0)
        let mut xd = Vector::new(space_ndim); // second middle point on Bezier curve with t=⅔(ξ=+⅓)
        let mut qc = Vector::new(space_ndim); // control point corresponding to pc
        let mut qd = Vector::new(space_ndim); // control point corresponding to pd

        // reference coordinate
        let mut ksi = vec![0.0; space_ndim];

        // memoization of scratchpads
        let mut pads_memo: HashMap<GeoKind, Scratchpad> = HashMap::new();

        // begin poly-curve
        self.canvas_edges.polycurve_begin();

        // loop over edges (GeoKind::Lin)
        for (_, edge) in &features.edges {
            // retrieve pad or allocate new
            let mut pad = pads_memo
                .entry(edge.kind)
                .or_insert(Scratchpad::new(space_ndim, edge.kind)?);

            // set coordinates of GeoKind::Lin
            let nnode = edge.points.len();
            for m in 0..nnode {
                for j in 0..space_ndim {
                    pad.set_xx(m, j, mesh.points[edge.points[m]].coords[j]);
                }
            }
            let (xa, ya) = (pad.xxt[0][0], pad.xxt[1][0]);
            let (xb, yb) = (pad.xxt[0][1], pad.xxt[1][1]);

            // add poly-curve points (or control points for quadratic/cubic Bezier)
            // (see file bezier-curves-math.pdf under data/derivations)
            match nnode {
                2 => {
                    self.canvas_edges
                        .polycurve_add(xa, ya, PolyCode::MoveTo)
                        .polycurve_add(xb, yb, PolyCode::LineTo);
                }
                3 => {
                    ksi[0] = 0.0; // middle
                    op::calc_coords(&mut xc, &mut pad, &ksi)?;
                    qc[0] = (-xa - xb + 4.0 * xc[0]) / 2.0;
                    qc[1] = (-ya - yb + 4.0 * xc[1]) / 2.0;
                    self.canvas_edges
                        .polycurve_add(xa, ya, PolyCode::MoveTo)
                        .polycurve_add(qc[0], qc[1], PolyCode::Curve3)
                        .polycurve_add(xb, yb, PolyCode::Curve3);
                }
                _ => {
                    ksi[0] = -1.0 / 3.0; // => t=⅓(ξ=-⅓)
                    op::calc_coords(&mut xc, &mut pad, &ksi)?;
                    ksi[0] = 1.0 / 3.0; // => t=⅔(ξ=+⅓)
                    op::calc_coords(&mut xd, &mut pad, &ksi)?;
                    qc[0] = (-5.0 * xa + 2.0 * xb + 18.0 * xc[0] - 9.0 * xd[0]) / 6.0;
                    qc[1] = (-5.0 * ya + 2.0 * yb + 18.0 * xc[1] - 9.0 * xd[1]) / 6.0;
                    qd[0] = (2.0 * xa - 5.0 * xb - 9.0 * xc[0] + 18.0 * xd[0]) / 6.0;
                    qd[1] = (2.0 * ya - 5.0 * yb - 9.0 * xc[1] + 18.0 * xd[1]) / 6.0;
                    self.canvas_edges
                        .polycurve_add(xa, ya, PolyCode::MoveTo)
                        .polycurve_add(qc[0], qc[1], PolyCode::Curve4)
                        .polycurve_add(qd[0], qd[1], PolyCode::Curve4)
                        .polycurve_add(xb, yb, PolyCode::Curve4);
                }
            }
        }

        // end polycurve and add to plot
        self.canvas_edges.polycurve_end(false);
        plot.add(&self.canvas_edges);
        if set_range {
            plot.set_range(features.min[0], features.max[0], features.min[1], features.max[1]);
        }
        Ok(())
    }

    /// Draws 3D edges
    fn edges_3d(&mut self, plot: &mut Plot, mesh: &Mesh, features: &Features, set_range: bool) -> Result<(), StrError> {
        // space dimension
        let space_ndim = mesh.ndim;
        assert_eq!(space_ndim, 3);

        // middle points (x..) on edge
        let mut xc = Vector::new(space_ndim);
        let mut xd = Vector::new(space_ndim);

        // reference coordinate
        let mut ksi = vec![0.0; space_ndim];

        // memoization of scratchpads
        let mut pads_memo: HashMap<GeoKind, Scratchpad> = HashMap::new();

        // loop over edges (GeoKind::Lin)
        for (_, edge) in &features.edges {
            // retrieve pad or allocate new
            let mut pad = pads_memo
                .entry(edge.kind)
                .or_insert(Scratchpad::new(space_ndim, edge.kind)?);

            // set coordinates of GeoKind::Lin
            let nnode = edge.points.len();
            for m in 0..nnode {
                for j in 0..space_ndim {
                    pad.set_xx(m, j, mesh.points[edge.points[m]].coords[j]);
                }
            }
            let (xa, ya, za) = (pad.xxt[0][0], pad.xxt[1][0], pad.xxt[2][0]);
            let (xb, yb, zb) = (pad.xxt[0][1], pad.xxt[1][1], pad.xxt[2][1]);

            // add poly-line points
            self.canvas_edges.polyline_3d_begin();
            match nnode {
                2 => {
                    self.canvas_edges
                        .polyline_3d_add(xa, ya, za)
                        .polyline_3d_add(xb, yb, zb);
                }
                3 => {
                    ksi[0] = 0.0; // middle
                    op::calc_coords(&mut xc, &mut pad, &ksi)?;
                    self.canvas_edges
                        .polyline_3d_add(xa, ya, za)
                        .polyline_3d_add(xc[0], xc[1], xc[2])
                        .polyline_3d_add(xb, yb, zb);
                }
                4 => {
                    ksi[0] = -1.0 / 3.0; // middle-left
                    op::calc_coords(&mut xc, &mut pad, &ksi)?;
                    ksi[0] = 1.0 / 3.0; // middle-right
                    op::calc_coords(&mut xd, &mut pad, &ksi)?;
                    self.canvas_edges
                        .polyline_3d_add(xa, ya, za)
                        .polyline_3d_add(xc[0], xc[1], xc[2])
                        .polyline_3d_add(xd[0], xd[1], xd[2])
                        .polyline_3d_add(xb, yb, zb);
                }
                _ => return Err("drawing of 3D edge with more than 4 nodes is not available"),
            }
            self.canvas_edges.polyline_3d_end();
        }

        // add to plot
        plot.add(&self.canvas_edges);
        if set_range {
            plot.set_range_3d(
                features.min[0],
                features.max[0],
                features.min[1],
                features.max[1],
                features.min[2],
                features.max[2],
            );
        }
        Ok(())
    }
}

/// Draws mesh
///
/// **Note:** This is a high-level function calling Draw and others
pub fn draw_mesh(mesh: &Mesh, with_ids: bool, filename: &str) -> Result<(), StrError> {
    let mut plot = Plot::new();
    let mut draw = Draw::new();
    draw.cells(&mut plot, &mesh, true)?;
    if with_ids {
        draw.cell_ids(&mut plot, mesh)?;
        draw.point_ids(&mut plot, mesh);
    } else {
        draw.points(&mut plot, mesh);
    }
    if mesh.ndim == 2 {
        plot.grid_and_labels("x", "y");
    }
    plot.set_equal_axes(true)
        .set_figure_size_points(600.0, 600.0)
        .save(filename)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Draw;
    use crate::mesh::{Extract, Mesh, Region, Samples};
    use crate::StrError;
    use plotpy::{Canvas, GraphMaker, Plot, Text};

    fn draw_cell_local_ids(plot: &mut Plot, labels: &mut Text, mesh: &Mesh, dx: f64, dy: f64, dz: f64) {
        for cell in &mesh.cells {
            for m in 0..cell.points.len() {
                let (x, y) = (
                    mesh.points[cell.points[m]].coords[0],
                    mesh.points[cell.points[m]].coords[1],
                );
                if mesh.ndim == 2 {
                    labels.draw(x + dx, y + dy, format!("{}", m).as_str());
                } else {
                    let z = mesh.points[cell.points[m]].coords[2];
                    labels.draw_3d(x + dx, y + dy, z + dz, format!("{}", m).as_str());
                }
            }
        }
        plot.add(labels);
    }

    fn draw_reference_circles() -> Plot {
        let mut plot = Plot::new();
        let mut circle_in = Canvas::new();
        let mut circle_mi = Canvas::new();
        let mut circle_ou = Canvas::new();
        circle_in
            .set_face_color("None")
            .set_edge_color("#bfbfbf")
            .set_line_width(7.0)
            .draw_circle(0.0, 0.0, 1.0);
        circle_mi
            .set_face_color("None")
            .set_edge_color("#bfbfbf")
            .set_line_width(7.0)
            .draw_circle(0.0, 0.0, 1.5);
        circle_ou
            .set_face_color("None")
            .set_edge_color("#bfbfbf")
            .set_line_width(7.0)
            .draw_circle(0.0, 0.0, 2.0);
        plot.add(&circle_in);
        plot.add(&circle_mi);
        plot.add(&circle_ou);
        plot
    }

    #[test]
    fn draw_cells_and_points_work() -> Result<(), StrError> {
        // labels for cell local ids
        let mut labels = Text::new();
        labels
            .set_color("#17af14")
            .set_align_vertical("bottom")
            .set_align_horizontal("center");

        // captions such as Lin3, Lin4, Tri3, ...
        let mut caption = Text::new();
        caption
            .set_fontsize(12.0)
            .set_align_horizontal("center")
            .set_align_vertical("top")
            .set_bbox(true)
            .set_bbox_style("round4,pad=0.3,rounding_size=0.15")
            .set_bbox_edgecolor("None")
            .set_bbox_facecolor("#57e986");

        // lin cells
        let mesh = Samples::lin_cells();
        let mut plot = Plot::new();
        let mut draw = Draw::new();
        draw.cells(&mut plot, &mesh, true)?;
        draw.points(&mut plot, &mesh);
        if false {
            draw_cell_local_ids(&mut plot, &mut labels, &mesh, 0.0, 0.05, 0.0);
            caption.draw(0.6, -0.05, "Lin2");
            caption.draw(1.4 + 0.6, -0.05, "Lin3");
            caption.draw(2.8 + 0.6, -0.05, "Lin4");
            caption.draw(4.2 + 0.6, -0.05, "Lin5");
            plot.add(&caption);
            plot.set_figure_size_points(600.0, 600.0)
                .set_frame_borders(false)
                .set_hide_axes(true)
                .set_equal_axes(true)
                .set_range(-0.1, 5.6, -0.3, 1.5)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_1_lin.svg")?;
        }

        // tri cells
        let mesh = Samples::tri_cells();
        let mut plot = Plot::new();
        let mut draw = Draw::new();
        draw.cells(&mut plot, &mesh, true)?;
        draw.points(&mut plot, &mesh);
        if false {
            labels.clear_buffer();
            caption.clear_buffer();
            draw_cell_local_ids(&mut plot, &mut labels, &mesh, 0.0, 0.02, 0.0);
            caption.draw(0.5, -0.1, "Tri3");
            caption.draw(1.7, -0.1, "Tri6");
            caption.draw(0.5, 1.1, "Tri10");
            caption.draw(1.7, 1.1, "Tri15");
            plot.add(&caption);
            plot.set_figure_size_points(600.0, 600.0)
                .set_frame_borders(false)
                .set_hide_axes(true)
                .set_equal_axes(true)
                .set_range(-0.1, 2.3, -0.2, 2.2)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_2_tri.svg")?;
        }

        let mesh = Samples::qua_cells();
        let mut plot = Plot::new();
        let mut draw = Draw::new();
        draw.cells(&mut plot, &mesh, true)?;
        draw.points(&mut plot, &mesh);
        if false {
            labels.clear_buffer();
            caption.clear_buffer();
            draw_cell_local_ids(&mut plot, &mut labels, &mesh, -0.02, 0.02, 0.0);
            caption.draw(0.4, -0.06, "Qua4");
            caption.draw(1.5, -0.06, "Qua8");
            caption.draw(2.6, -0.06, "Qua9");
            caption.draw(0.4, 1.1, "Qua12");
            caption.draw(1.5, 1.1, "Qua16");
            caption.draw(2.6, 1.1, "Qua17");
            plot.add(&caption);
            plot.set_figure_size_points(600.0, 600.0)
                .set_frame_borders(false)
                .set_hide_axes(true)
                .set_equal_axes(true)
                .set_range(-0.15, 3.1, -0.25, 2.15)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_3_qua.svg")?;
        }

        let mesh = Samples::tet_cells();
        let mut plot = Plot::new();
        let mut draw = Draw::new();
        draw.cells(&mut plot, &mesh, true)?;
        draw.points(&mut plot, &mesh);
        if false {
            labels.clear_buffer();
            caption.clear_buffer();
            draw_cell_local_ids(&mut plot, &mut labels, &mesh, 0.03, 0.0, 0.03);
            caption.draw_3d(0.5, -0.15, 0.0, "Tet4");
            caption.draw_3d(1.7, -0.15, 0.0, "Tet10");
            caption.draw_3d(1.1, 1.05, 0.8, "Tet20");
            plot.add(&caption);
            plot.set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range_3d(0.0, 2.2, 0.0, 2.2, 0.0, 2.2)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_4_tet.svg")?;
        }

        let mesh = Samples::hex_cells();
        let mut plot = Plot::new();
        let mut draw = Draw::new();
        draw.cells(&mut plot, &mesh, true)?;
        draw.points(&mut plot, &mesh);
        if false {
            labels.clear_buffer();
            caption.clear_buffer();
            draw_cell_local_ids(&mut plot, &mut labels, &mesh, 0.03, 0.0, 0.03);
            caption.draw_3d(0.5, -0.15, 0.0, "Hex8");
            caption.draw_3d(2.3, -0.15, 0.0, "Hex20");
            caption.draw_3d(1.1, 1.1, 1.3, "Hex32");
            plot.add(&caption);
            plot.set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range_3d(0.0, 3.0, 0.0, 3.0, 0.0, 2.5)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_5_hex.svg")?;
        }

        let mut plot = draw_reference_circles();
        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        let mut draw = Draw::new();
        draw.cells(&mut plot, &mesh, true)?;
        draw.points(&mut plot, &mesh);
        if false {
            plot.set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range(-0.1, 2.1, -0.1, 2.1)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_6_ring.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_edges_points_and_ids_work_2d_ring() -> Result<(), StrError> {
        // draw reference circles
        let mut plot = draw_reference_circles();

        // draw edges
        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        let region = Region::with(&mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.edges(&mut plot, &region, false)?;

        // draw points and point ids
        draw.canvas_points.set_marker_color("red");
        draw.canvas_point_ids
            .set_align_horizontal("left")
            .set_align_vertical("bottom")
            .set_color("gold")
            .set_fontsize(9.0)
            .set_bbox_facecolor("black")
            .set_bbox_alpha(0.6);
        draw.points(&mut plot, &region.mesh);
        draw.point_ids(&mut plot, &region.mesh);

        // draw cell ids
        draw.cell_ids(&mut plot, &region.mesh)?;

        // save figure
        if false {
            plot.set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range(-0.1, 2.1, -0.1, 2.1)
                .save("/tmp/gemlab/test_draw_edges_points_and_ids_work_2d_ring.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_edges_and_ids_work_2d_qua12() -> Result<(), StrError> {
        let mut plot = Plot::new();
        let mesh = Samples::block_2d_four_qua12();
        let region = Region::with(&mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.edges(&mut plot, &region, true)?;
        draw.point_ids(&mut plot, &region.mesh);
        draw.cell_ids(&mut plot, &region.mesh)?;
        if false {
            plot.set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_and_ids_work_2d_qua12.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_edges_and_ids_work_2d_qua16() -> Result<(), StrError> {
        let mut plot = Plot::new();
        let mesh = Samples::block_2d_four_qua16();
        let region = Region::with(&mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.edges(&mut plot, &region, true)?;
        draw.point_ids(&mut plot, &region.mesh);
        draw.cell_ids(&mut plot, &region.mesh)?;
        if false {
            plot.set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_and_ids_work_2d_qua16.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_edges_and_ids_work_2d_qua17() -> Result<(), StrError> {
        let mut plot = Plot::new();
        let mesh = Samples::block_2d_four_qua17();
        let region = Region::with(&mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.edges(&mut plot, &region, true)?;
        draw.point_ids(&mut plot, &region.mesh);
        if false {
            plot.set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_and_ids_work_2d_qua17.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_edges_and_ids_work_2d_mixed() -> Result<(), StrError> {
        let mut plot = Plot::new();
        let mesh = Samples::mixed_shapes_2d();
        let region = Region::with(&mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.edges(&mut plot, &region, true)?;
        draw.point_ids(&mut plot, &region.mesh);
        draw.cell_ids(&mut plot, &region.mesh)?;
        if false {
            plot.set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_and_ids_work_2d_mixed.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_edges_points_and_ids_work_3d() -> Result<(), StrError> {
        let mut plot = Plot::new();
        let mesh = Samples::two_hex8();
        let region = Region::with(&mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.edges(&mut plot, &region, true)?;
        draw.canvas_point_ids
            .set_align_horizontal("left")
            .set_align_vertical("bottom")
            .set_color("black")
            .set_fontsize(10.0)
            .set_bbox_facecolor("gold")
            .set_bbox_alpha(0.5);
        draw.points(&mut plot, &region.mesh);
        draw.point_ids(&mut plot, &region.mesh);
        draw.cell_ids(&mut plot, &region.mesh)?;
        if false {
            plot.set_figure_size_points(500.0, 500.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_points_and_ids_work_3d.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_edges_and_ids_work_3d_1() -> Result<(), StrError> {
        let mut plot = Plot::new();
        let mesh = Samples::block_3d_eight_hex20();
        let region = Region::with(&mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.edges(&mut plot, &region, true)?;
        draw.point_ids(&mut plot, &region.mesh);
        draw.cell_ids(&mut plot, &region.mesh)?;
        if false {
            plot.set_figure_size_points(800.0, 800.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_and_ids_work_3d_1.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_edges_and_ids_work_3d_2() -> Result<(), StrError> {
        let mut plot = Plot::new();
        let mesh = Samples::mixed_shapes_3d();
        let region = Region::with(&mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.edges(&mut plot, &region, true)?;
        draw.point_ids(&mut plot, &region.mesh);
        draw.cell_ids(&mut plot, &region.mesh)?;
        if false {
            plot.set_figure_size_points(800.0, 800.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_and_ids_work_3d_2.svg")?;
        }
        Ok(())
    }
}
