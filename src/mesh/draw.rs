use super::{Extract, Features, Mesh, Region};
use crate::shapes::{op, GeoKind, Scratchpad};
use crate::StrError;
use plotpy::{Canvas, Curve, Plot, PolyCode, Text};
use russell_lab::Vector;
use std::collections::HashMap;

/// Implements functions to draw edges and faces
pub struct Draw {
    /// Canvas to draw edges
    pub canvas_edges: Canvas,

    /// Canvas to draw points (markers)
    pub canvas_points: Curve,

    /// Canvas to draw point ids
    pub canvas_point_ids: Text,

    /// Canvas to draw cell ids
    pub canvas_cell_ids: Text,
}

impl Draw {
    /// Allocates a new instance
    pub fn new() -> Self {
        let mut canvas_edges = Canvas::new();
        let mut canvas_points = Curve::new();
        let mut canvas_point_ids = Text::new();
        let mut canvas_cell_ids = Text::new();
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
        Draw {
            canvas_edges,
            canvas_points,
            canvas_point_ids,
            canvas_cell_ids,
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
            if mesh.ndim == 2 {
                self.canvas_cell_ids
                    .draw(x[0], x[1], format!("{}({})", cell.id, cell.attribute_id).as_str());
            } else {
                self.canvas_cell_ids
                    .draw_3d(x[0], x[1], x[2], format!("{}({})", cell.id, cell.attribute_id).as_str());
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
#[no_coverage]
pub fn draw_mesh(mesh: &Mesh, with_ids: bool, filename: &str) -> Result<(), StrError> {
    let region = Region::with(mesh, Extract::All)?;
    let mut plot = Plot::new();
    let mut draw = Draw::new();
    draw.edges(&mut plot, &region, true)?;
    if with_ids {
        draw.cell_ids(&mut plot, &region.mesh)?;
        draw.point_ids(&mut plot, &region.mesh);
    } else {
        draw.points(&mut plot, &region.mesh);
    }
    plot.set_equal_axes(true)
        .grid_and_labels("x", "y")
        .set_figure_size_points(600.0, 600.0)
        .save(filename)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Draw;
    use crate::mesh::{Extract, Region, Samples};
    use crate::StrError;
    use plotpy::{Canvas, Plot};

    #[test]
    fn draw_works_2d_ring() -> Result<(), StrError> {
        // draw reference circles
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
            plot.set_figure_size_points(400.0, 400.0)
                .set_equal_axes(true)
                .set_range(-0.1, 2.1, -0.1, 2.1)
                .save("/tmp/gemlab/draw_works_2d_ring.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_works_2d_qua12() -> Result<(), StrError> {
        let mut plot = Plot::new();
        let mesh = Samples::block_2d_four_qua12();
        let region = Region::with(&mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.edges(&mut plot, &region, true)?;
        draw.point_ids(&mut plot, &region.mesh);
        draw.cell_ids(&mut plot, &region.mesh)?;
        if false {
            plot.set_figure_size_points(400.0, 400.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/draw_works_2d_qua12.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_works_2d_qua16() -> Result<(), StrError> {
        let mut plot = Plot::new();
        let mesh = Samples::block_2d_four_qua16();
        let region = Region::with(&mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.edges(&mut plot, &region, true)?;
        draw.point_ids(&mut plot, &region.mesh);
        draw.cell_ids(&mut plot, &region.mesh)?;
        if false {
            plot.set_figure_size_points(400.0, 400.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/draw_works_2d_qua16.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_works_2d_qua17() -> Result<(), StrError> {
        let mut plot = Plot::new();
        let mesh = Samples::block_2d_four_qua17();
        let region = Region::with(&mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.edges(&mut plot, &region, true)?;
        draw.point_ids(&mut plot, &region.mesh);
        if false {
            plot.set_figure_size_points(400.0, 400.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/draw_works_2d_qua17.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_works_2d_mixed() -> Result<(), StrError> {
        let mut plot = Plot::new();
        let mesh = Samples::mixed_shapes_2d();
        let region = Region::with(&mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.edges(&mut plot, &region, true)?;
        draw.point_ids(&mut plot, &region.mesh);
        draw.cell_ids(&mut plot, &region.mesh)?;
        if false {
            plot.set_figure_size_points(400.0, 400.0)
                .set_equal_axes(true)
                .save("/tmp/gemlab/draw_works_2d_mixed.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_works_3d() -> Result<(), StrError> {
        let mut plot = Plot::new();
        let mesh = Samples::two_cubes_vertical();
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
                .save("/tmp/gemlab/draw_works_3d.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_works_3d_eight_hex20() -> Result<(), StrError> {
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
                .save("/tmp/gemlab/draw_works_3d_eight_hex20.svg")?;
        }
        Ok(())
    }

    #[test]
    fn draw_works_3d_mixed() -> Result<(), StrError> {
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
                .save("/tmp/gemlab/draw_works_3d_mixed.svg")?;
        }
        Ok(())
    }
}
