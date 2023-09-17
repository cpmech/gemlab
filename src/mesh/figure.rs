use super::{Features, Mesh, PointId};
use crate::shapes::{GeoKind, Scratchpad};
use crate::StrError;
use plotpy::{Canvas, Curve, Plot, PolyCode, Text};
use russell_lab::math::ONE_BY_3;
use russell_lab::Vector;
use std::collections::HashMap;
use std::ffi::OsStr;

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
                let pad = pads.entry(lin_kind).or_insert(Scratchpad::new(mesh.ndim, lin_kind)?);
                for m in 0..lin_points.len() {
                    for j in 0..mesh.ndim {
                        pad.set_xx(m, j, pp[ii[m]].coords[j]);
                    }
                }
                let (xa, ya) = (pp[ii[0]].coords[0], pp[ii[0]].coords[1]);
                let (xb, yb) = (pp[ii[1]].coords[0], pp[ii[1]].coords[1]);
                pad.calc_coords(&mut x, &[-ONE_BY_3, 0.0]).unwrap();
                let (xc, yc) = (x[0], x[1]);
                pad.calc_coords(&mut x, &[ONE_BY_3, 0.0]).unwrap();
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
pub fn draw_cell(
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
pub struct Figure {
    /// The plotpy structure to draw figures (plots)
    pub plot: Plot,

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

    /// Parameter: draw cell ids
    pub param_cell_ids: bool,

    /// Parameter: draw point ids
    pub param_point_ids: bool,

    /// Parameter: draw point dots
    pub param_point_dots: bool,

    /// Parameter: plot without equal axes
    pub param_not_equal_exes: bool,

    /// Parameter: specifies the plot range (xmin, xmax, ymin, ymax)
    pub param_range_2d: Option<(f64, f64, f64, f64)>,

    /// Parameter: specifies the plot range (xmin, xmax, ymin, ymax, zmin, zmax)
    pub param_range_3d: Option<(f64, f64, f64, f64, f64, f64)>,

    /// Parameter: figure size in points
    pub param_figure_size: Option<(f64, f64)>,
}

impl Figure {
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
            .set_bbox_style("round,pad=0.15");
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
        Figure {
            plot: Plot::new(),
            canvas_edges,
            canvas_points,
            canvas_point_ids,
            canvas_cell_ids,
            canvas_cells,
            canvas_lin_cells,
            param_cell_ids: false,
            param_point_ids: false,
            param_point_dots: false,
            param_not_equal_exes: false,
            param_range_2d: None,
            param_range_3d: None,
            param_figure_size: None,
        }
    }
}

impl Mesh {
    /// Draws all points (dots)
    pub fn draw_point_dots(&self, fig: &mut Figure) {
        if self.ndim == 2 {
            fig.canvas_points.points_begin();
            self.points.iter().for_each(|point| {
                fig.canvas_points.points_add(point.coords[0], point.coords[1]);
            });
            fig.canvas_points.points_end();
        } else {
            fig.canvas_points.points_3d_begin();
            self.points.iter().for_each(|point| {
                fig.canvas_points
                    .points_3d_add(point.coords[0], point.coords[1], point.coords[2]);
            });
            fig.canvas_points.points_3d_end();
        }
        fig.plot.add(&fig.canvas_points);
    }

    /// Draws all point ids (labels)
    ///
    /// **Note:** Negative point markers are shown within parentheses.
    pub fn draw_point_ids(&self, fig: &mut Figure) {
        if self.ndim == 2 {
            self.points.iter().for_each(|point| {
                let msg = if point.marker < 0 {
                    format!("{}({})", point.id, point.marker)
                } else {
                    format!("{}", point.id)
                };
                fig.canvas_point_ids
                    .draw(point.coords[0], point.coords[1], msg.as_str());
            });
        } else {
            self.points.iter().for_each(|point| {
                fig.canvas_point_ids.draw_3d(
                    point.coords[0],
                    point.coords[1],
                    point.coords[2],
                    format!("{}", point.id).as_str(),
                );
            });
        }
        fig.plot.add(&fig.canvas_point_ids);
    }

    /// Draws cells
    pub fn draw_cells(&self, fig: &mut Figure, set_range: bool) -> Result<(), StrError> {
        // limits
        let mut xmin = vec![f64::MAX; self.ndim];
        let mut xmax = vec![f64::MIN; self.ndim];

        // memoization of scratchpads
        let mut pads = HashMap::new();

        // loop over cells
        for cell in &self.cells {
            let canvas = if cell.kind.is_lin() {
                &mut fig.canvas_lin_cells
            } else {
                &mut fig.canvas_cells
            };
            draw_cell(canvas, self, cell.kind, &cell.points, &mut pads)?;
            for point_id in &cell.points {
                for j in 0..self.ndim {
                    xmin[j] = f64::min(xmin[j], self.points[*point_id].coords[j]);
                    xmax[j] = f64::max(xmax[j], self.points[*point_id].coords[j]);
                }
            }
        }

        // add to plot
        fig.plot.add(&fig.canvas_cells);
        fig.plot.add(&fig.canvas_lin_cells);
        if set_range {
            let dx = xmax[0] - xmin[0];
            let dy = xmax[1] - xmin[1];
            if dx > 0.0 && dy > 0.0 {
                let gx = 0.05 * dx;
                let gy = 0.05 * dy;
                fig.plot
                    .set_range(xmin[0] - gx, xmax[0] + gx, xmin[1] - gy, xmax[1] + gy);
            }
        }
        Ok(())
    }

    /// Draws ids and attributes of cells
    ///
    /// **Note:** Cell attributes are shown within parentheses.
    pub fn draw_cell_ids(&self, fig: &mut Figure) -> Result<(), StrError> {
        // auxiliary
        let mut x = Vector::new(self.ndim);

        // loop over all cells
        for cell_id in 0..self.cells.len() {
            // compute coordinates of the label
            let cell = &self.cells[cell_id];
            x.fill(0.0);
            if cell.kind == GeoKind::Qua9 || cell.kind == GeoKind::Qua17 {
                for i in 0..self.ndim {
                    x[i] = (3.0 * self.points[cell.points[0]].coords[i] + self.points[cell.points[2]].coords[i]) / 4.0;
                }
            } else if cell.kind == GeoKind::Tri10 {
                for i in 0..self.ndim {
                    x[i] = (self.points[cell.points[0]].coords[i] + self.points[cell.points[9]].coords[i]) / 2.0;
                }
            } else {
                for m in 0..cell.points.len() {
                    for i in 0..self.ndim {
                        x[i] += self.points[cell.points[m]].coords[i];
                    }
                }
                for i in 0..self.ndim {
                    x[i] /= cell.points.len() as f64;
                }
            }

            // add label
            let msg = format!("{}({})", cell.id, cell.attribute);
            if self.ndim == 2 {
                fig.canvas_cell_ids.draw(x[0], x[1], msg.as_str());
            } else {
                fig.canvas_cell_ids.draw_3d(x[0], x[1], x[2], msg.as_str());
            }
        }

        // add to plot
        fig.plot.add(&fig.canvas_cell_ids);
        Ok(())
    }

    /// Draws edges
    ///
    /// # Input
    ///
    /// * `plot` -- the plot instance (to be updated)
    /// * `mesh` -- the mesh
    /// * `features` -- structure with edges and faces
    /// * `set_range` -- sets the range of `plot` to the limits of region; otherwise do not modifies the range/limits.
    pub fn draw_edges(&self, fig: &mut Figure, features: &Features, set_range: bool) -> Result<(), StrError> {
        let mut pads = HashMap::new();
        for (_, edge) in &features.edges {
            draw_cell(&mut &mut fig.canvas_edges, self, edge.kind, &edge.points, &mut pads)?;
        }
        fig.plot.add(&fig.canvas_edges);
        if set_range {
            if self.ndim == 2 {
                fig.plot
                    .set_range(features.min[0], features.max[0], features.min[1], features.max[1]);
            } else {
                fig.plot.set_range_3d(
                    features.min[0],
                    features.max[0],
                    features.min[1],
                    features.max[1],
                    features.min[2],
                    features.max[2],
                );
            }
        }
        Ok(())
    }

    /// Draws the cells, points, attributes, and markers
    ///
    /// # Input
    ///
    /// * `fig` -- the Figure struct (optional => use default configuration)
    /// * `filepath` -- may be a String, &str, or Path
    pub fn draw<P>(&self, fig: Option<Figure>, filepath: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let mut figure = if let Some(f) = fig { f } else { Figure::new() };
        self.draw_cells(&mut figure, true)?;
        if figure.param_cell_ids {
            self.draw_cell_ids(&mut figure)?;
        }
        if figure.param_point_dots {
            self.draw_point_dots(&mut figure);
        }
        if figure.param_point_ids {
            self.draw_point_ids(&mut figure);
        }
        if self.ndim == 2 {
            figure.plot.grid_and_labels("x", "y");
        }
        if !figure.param_not_equal_exes {
            figure.plot.set_equal_axes(true);
        }
        if self.ndim == 2 {
            if let Some((xmin, xmax, ymin, ymax)) = figure.param_range_2d {
                figure.plot.set_range(xmin, xmax, ymin, ymax);
            }
        } else {
            if let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = figure.param_range_3d {
                figure.plot.set_range_3d(xmin, xmax, ymin, ymax, zmin, zmax);
            }
        }
        if let Some((width, height)) = figure.param_figure_size {
            figure.plot.set_figure_size_points(width, height);
        }
        figure.plot.save(filepath)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{add_curve, Figure};
    use crate::mesh::{Extract, Features, Mesh, Samples};
    use crate::shapes::GeoKind;
    use plotpy::{Canvas, Plot, Text};
    use std::collections::HashMap;

    const SAVE_FIGURE: bool = false;

    fn labels_and_caption() -> (Text, Text) {
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
        (labels, caption)
    }

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

    fn draw_reference_circles(plot: &mut Plot) {
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
    }

    #[test]
    fn add_curve_catches_errors() {
        let mesh = Samples::one_lin2();
        let mut pads = HashMap::new();
        let mut canvas = Canvas::new();
        assert_eq!(
            add_curve(&mut canvas, &mesh, GeoKind::Tri3, &[], true, true, &mut pads).err(),
            Some("lin_kind is not Lin")
        );
    }

    #[test]
    fn draw_cells_and_points_work() {
        // lin cells
        let mesh = Samples::lin_cells();
        let mut fig = Figure::new();
        mesh.draw_cells(&mut fig, true).unwrap();
        mesh.draw_point_dots(&mut fig);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut fig.plot, &mut labels, &mesh, 0.0, 0.05, 0.0);
            caption.draw(0.6, -0.05, "Lin2");
            caption.draw(1.4 + 0.6, -0.05, "Lin3");
            caption.draw(2.8 + 0.6, -0.05, "Lin4");
            caption.draw(4.2 + 0.6, -0.05, "Lin5");
            fig.plot.add(&caption);
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_frame_borders(false)
                .set_hide_axes(true)
                .set_equal_axes(true)
                .set_range(-0.1, 5.6, -0.3, 1.5)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_1_lin.svg")
                .unwrap();
        }

        // lin cells in 3d
        let mesh = Samples::lin_cells_3d();
        let mut fig = Figure::new();
        mesh.draw_cells(&mut fig, true).unwrap();
        mesh.draw_point_dots(&mut fig);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut fig.plot, &mut labels, &mesh, 0.0, 0.05, 0.0);
            caption.draw_3d(1.2, 1.33, 1.33, "Lin2");
            caption.draw_3d(1.4 + 1.2, 1.33, 1.33, "Lin3");
            caption.draw_3d(2.8 + 1.2, 1.33, 1.33, "Lin4");
            caption.draw_3d(4.2 + 1.2, 1.33, 1.33, "Lin5");
            fig.plot.add(&caption);
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range_3d(0.0, 5.3, 0.0, 1.2, 0.0, 1.2)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_1_lin_3d.svg")
                .unwrap();
        }

        // tri cells
        let mesh = Samples::tri_cells();
        let mut fig = Figure::new();
        mesh.draw_cells(&mut fig, true).unwrap();
        mesh.draw_point_dots(&mut fig);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut fig.plot, &mut labels, &mesh, 0.0, 0.02, 0.0);
            caption.draw(0.5, -0.1, "Tri3");
            caption.draw(1.7, -0.1, "Tri6");
            caption.draw(0.5, 1.1, "Tri10");
            caption.draw(1.7, 1.1, "Tri15");
            fig.plot.add(&caption);
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_frame_borders(false)
                .set_hide_axes(true)
                .set_equal_axes(true)
                .set_range(-0.1, 2.3, -0.2, 2.2)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_2_tri.svg")
                .unwrap();
        }

        let mesh = Samples::qua_cells();
        let mut fig = Figure::new();
        mesh.draw_cells(&mut fig, true).unwrap();
        mesh.draw_point_dots(&mut fig);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut fig.plot, &mut labels, &mesh, -0.02, 0.02, 0.0);
            caption.draw(0.4, -0.06, "Qua4");
            caption.draw(1.5, -0.06, "Qua8");
            caption.draw(2.6, -0.06, "Qua9");
            caption.draw(0.4, 1.1, "Qua12");
            caption.draw(1.5, 1.1, "Qua16");
            caption.draw(2.6, 1.1, "Qua17");
            fig.plot.add(&caption);
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_frame_borders(false)
                .set_hide_axes(true)
                .set_equal_axes(true)
                .set_range(-0.15, 3.1, -0.25, 2.15)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_3_qua.svg")
                .unwrap();
        }

        let mesh = Samples::tet_cells();
        let mut fig = Figure::new();
        mesh.draw_cells(&mut fig, true).unwrap();
        mesh.draw_point_dots(&mut fig);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut fig.plot, &mut labels, &mesh, 0.03, 0.0, 0.03);
            caption.draw_3d(0.5, -0.15, 0.0, "Tet4");
            caption.draw_3d(1.7, -0.15, 0.0, "Tet10");
            caption.draw_3d(1.1, 1.05, 0.8, "Tet20");
            fig.plot.add(&caption);
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range_3d(0.0, 2.2, 0.0, 2.2, 0.0, 2.2)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_4_tet.svg")
                .unwrap();
        }

        let mesh = Samples::hex_cells();
        let mut fig = Figure::new();
        mesh.draw_cells(&mut fig, true).unwrap();
        mesh.draw_point_dots(&mut fig);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut fig.plot, &mut labels, &mesh, 0.03, 0.0, 0.03);
            caption.draw_3d(0.5, -0.15, 0.0, "Hex8");
            caption.draw_3d(2.3, -0.15, 0.0, "Hex20");
            caption.draw_3d(1.1, 1.1, 1.3, "Hex32");
            fig.plot.add(&caption);
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range_3d(0.0, 3.0, 0.0, 3.0, 0.0, 2.5)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_5_hex.svg")
                .unwrap();
        }

        draw_reference_circles(&mut fig.plot);
        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        let mut fig = Figure::new();
        mesh.draw_cells(&mut fig, true).unwrap();
        mesh.draw_point_dots(&mut fig);

        if SAVE_FIGURE {
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range(-0.1, 2.1, -0.1, 2.1)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_6_ring.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_edges_points_and_ids_work_2d_ring() {
        // draw reference circles and edges
        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        let features = Features::new(&mesh, Extract::All);
        let mut fig = Figure::new();
        draw_reference_circles(&mut fig.plot);
        mesh.draw_edges(&mut fig, &features, false).unwrap();

        // draw points and point ids
        fig.canvas_points.set_marker_color("red");
        fig.canvas_point_ids
            .set_align_horizontal("left")
            .set_align_vertical("bottom")
            .set_color("gold")
            .set_fontsize(9.0)
            .set_bbox_facecolor("black")
            .set_bbox_alpha(0.6);
        mesh.draw_point_dots(&mut fig);
        mesh.draw_point_ids(&mut fig);

        // draw cell ids
        mesh.draw_cell_ids(&mut fig).unwrap();

        if SAVE_FIGURE {
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range(-0.1, 2.1, -0.1, 2.1)
                .save("/tmp/gemlab/test_draw_edges_points_and_ids_work_2d_ring.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_edges_and_ids_work_2d_qua12() {
        let mesh = Samples::block_2d_four_qua12();
        let features = Features::new(&mesh, Extract::All);
        let mut fig = Figure::new();
        mesh.draw_edges(&mut fig, &features, true).unwrap();
        mesh.draw_point_ids(&mut fig);
        mesh.draw_cell_ids(&mut fig).unwrap();

        if SAVE_FIGURE {
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_and_ids_work_2d_qua12.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_edges_and_ids_work_2d_qua16() {
        let mesh = Samples::block_2d_four_qua16();
        let features = Features::new(&mesh, Extract::All);
        let mut fig = Figure::new();
        mesh.draw_edges(&mut fig, &features, true).unwrap();
        mesh.draw_point_ids(&mut fig);
        mesh.draw_cell_ids(&mut fig).unwrap();

        if SAVE_FIGURE {
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_and_ids_work_2d_qua16.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_edges_and_ids_work_2d_qua17() {
        let mesh = Samples::block_2d_four_qua17();
        let features = Features::new(&mesh, Extract::All);
        let mut fig = Figure::new();
        mesh.draw_edges(&mut fig, &features, true).unwrap();
        mesh.draw_point_ids(&mut fig);

        if SAVE_FIGURE {
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_and_ids_work_2d_qua17.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_edges_and_ids_work_2d_mixed() {
        let mesh = Samples::mixed_shapes_2d();
        let features = Features::new(&mesh, Extract::All);
        let mut fig = Figure::new();
        mesh.draw_edges(&mut fig, &features, true).unwrap();
        mesh.draw_point_ids(&mut fig);
        mesh.draw_cell_ids(&mut fig).unwrap();

        if SAVE_FIGURE {
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_and_ids_work_2d_mixed.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_edges_points_and_ids_work_3d() {
        let mesh = Samples::two_hex8();
        let features = Features::new(&mesh, Extract::All);
        let mut fig = Figure::new();
        mesh.draw_edges(&mut fig, &features, true).unwrap();
        fig.canvas_point_ids
            .set_align_horizontal("left")
            .set_align_vertical("bottom")
            .set_color("black")
            .set_fontsize(10.0)
            .set_bbox_facecolor("gold")
            .set_bbox_alpha(0.5);
        mesh.draw_point_dots(&mut fig);
        mesh.draw_point_ids(&mut fig);
        mesh.draw_cell_ids(&mut fig).unwrap();

        if SAVE_FIGURE {
            fig.plot
                .set_figure_size_points(500.0, 500.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_points_and_ids_work_3d.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_edges_and_ids_work_3d_1() {
        let mesh = Samples::block_3d_eight_hex20();
        let features = Features::new(&mesh, Extract::All);
        let mut fig = Figure::new();
        mesh.draw_edges(&mut fig, &features, true).unwrap();
        mesh.draw_point_ids(&mut fig);
        mesh.draw_cell_ids(&mut fig).unwrap();

        if SAVE_FIGURE {
            fig.plot
                .set_figure_size_points(800.0, 800.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_and_ids_work_3d_1.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_edges_and_ids_work_3d_2() {
        let mesh = Samples::mixed_shapes_3d();
        let features = Features::new(&mesh, Extract::All);
        let mut fig = Figure::new();
        mesh.draw_edges(&mut fig, &features, true).unwrap();
        mesh.draw_point_ids(&mut fig);
        mesh.draw_cell_ids(&mut fig).unwrap();

        if SAVE_FIGURE {
            fig.plot
                .set_figure_size_points(800.0, 800.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edges_and_ids_work_3d_2.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_point_markers_work_2d() {
        let mesh = Samples::qua8_tri6_lin2();
        let features = Features::new(&mesh, Extract::All);
        let mut fig = Figure::new();
        mesh.draw_edges(&mut fig, &features, true).unwrap();
        mesh.draw_point_ids(&mut fig);

        if SAVE_FIGURE {
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_range(-0.1, 2.0, -0.1, 1.1)
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_point_markers_work_2d.svg")
                .unwrap();
        }
    }
}
