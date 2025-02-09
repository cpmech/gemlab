use super::Mesh;
use crate::shapes::GeoKind;
use crate::StrError;
use plotpy::{Canvas, Curve, InsetAxes, Plot, Text};
use russell_lab::Vector;
use std::ffi::OsStr;

/// Implements functions to draw cells, edges and faces
pub struct Figure {
    /// The plotpy structure to draw figures (plots)
    plot: Plot,

    /// Canvas to draw edges
    canvas_edges: Canvas,

    /// Canvas to draw points (markers)
    canvas_points: Curve,

    /// Canvas to draw point ids
    canvas_point_ids: Text,

    /// Canvas to draw cell ids
    canvas_cell_ids: Text,

    /// Canvas to draw cells
    canvas_cells: Canvas,

    /// Canvas to draw lin cells
    canvas_lin_cells: Canvas,

    /// Shows cell ids
    show_cell_ids: bool,

    /// Shows cell attribute within parenthesis
    show_cell_att: bool,

    /// Shows point ids
    show_point_ids: bool,

    /// Shows point marker within parenthesis (if not zero)
    show_point_marker: bool,

    /// Shows point dots
    show_point_dots: bool,

    /// Generates the plot without equal axes
    unequal_exes: bool,

    /// Specifies the plot range (xmin, xmax, ymin, ymax)
    range_2d: Option<(f64, f64, f64, f64)>,

    /// Specifies the plot range (xmin, xmax, ymin, ymax, zmin, zmax)
    range_3d: Option<(f64, f64, f64, f64, f64, f64)>,

    /// Specifies the figure size in points
    figure_size: Option<(f64, f64)>,

    /// Enables zooming a region of the plot (2D meshes only)
    ///
    /// Holds `((xmin, xmax, ymin, ymax), (u0, u0, w, h))` for the zoomed region; where:
    ///
    /// * `(xmin, xmax, ymin, ymax)` -- the region to be zoomed
    /// * `(u0, v0, w, h)` -- the position (normalized coordinates) and size (width,height) of the zoomed region
    zoom_2d: Option<((f64, f64, f64, f64), (f64, f64, f64, f64))>,

    /// Configures some characteristics of the zoom indicator
    ///
    /// Holds `(color, alpha, linewidth)` for the zoom indicator.
    zoom_indicator_config: (String, f64, f64),
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
            .set_face_color("None")
            .set_line_width(1.0)
            .set_edge_color("#2440cd");
        canvas_points
            .set_marker_color("black")
            .set_marker_line_color("white")
            .set_marker_style("o")
            .set_line_style("None");
        canvas_point_ids
            .set_extra("clip_on=True")
            .set_color("red")
            .set_fontsize(7.0)
            .set_align_horizontal("center")
            .set_align_vertical("center")
            .set_bbox(true)
            .set_bbox_facecolor("white")
            .set_bbox_edgecolor("None")
            .set_bbox_style("round,pad=0.15");
        canvas_cell_ids
            .set_extra("clip_on=True")
            .set_color("#22971f")
            .set_fontsize(8.0)
            .set_align_horizontal("center")
            .set_align_vertical("center")
            .set_bbox(true)
            .set_bbox_facecolor("white")
            .set_bbox_edgecolor("None")
            .set_bbox_style("square,pad=0.15");
        canvas_cells
            .set_face_color("#e3f3ff")
            .set_edge_color("#0055d4")
            .set_line_width(1.0);
        canvas_lin_cells
            .set_face_color("None")
            .set_edge_color("#cd0000")
            .set_line_width(2.0);
        Figure {
            plot: Plot::new(),
            canvas_edges,
            canvas_points,
            canvas_point_ids,
            canvas_cell_ids,
            canvas_cells,
            canvas_lin_cells,
            show_cell_ids: false,
            show_cell_att: true,
            show_point_ids: false,
            show_point_marker: false,
            show_point_dots: false,
            unequal_exes: false,
            range_2d: None,
            range_3d: None,
            figure_size: None,
            zoom_2d: None,
            zoom_indicator_config: ("#f9d835".to_string(), 1.0, 2.0),
        }
    }

    /// Get a mutable reference to the canvas edges
    pub fn canvas_edges(&mut self) -> &mut Canvas {
        &mut self.canvas_edges
    }

    /// Get a mutable reference to the canvas points
    pub fn canvas_points(&mut self) -> &mut Curve {
        &mut self.canvas_points
    }

    /// Get a mutable reference to the canvas point ids
    pub fn canvas_point_ids(&mut self) -> &mut Text {
        &mut self.canvas_point_ids
    }

    /// Get a mutable reference to the canvas cell ids
    pub fn canvas_cell_ids(&mut self) -> &mut Text {
        &mut self.canvas_cell_ids
    }

    /// Get a mutable reference to the canvas cells
    pub fn canvas_cells(&mut self) -> &mut Canvas {
        &mut self.canvas_cells
    }

    /// Get a mutable reference to the canvas lin cells
    pub fn canvas_lin_cells(&mut self) -> &mut Canvas {
        &mut self.canvas_lin_cells
    }

    /// Shows cell ids
    pub fn show_cell_ids(&mut self, value: bool) -> &mut Self {
        self.show_cell_ids = value;
        self
    }

    /// Shows cell attribute within parenthesis
    pub fn show_cell_att(&mut self, value: bool) -> &mut Self {
        self.show_cell_att = value;
        self
    }

    /// Shows point ids
    pub fn show_point_ids(&mut self, value: bool) -> &mut Self {
        self.show_point_ids = value;
        self
    }

    /// Shows point marker within parenthesis (if not zero)
    pub fn show_point_marker(&mut self, value: bool) -> &mut Self {
        self.show_point_marker = value;
        self
    }

    /// Shows point dots
    pub fn show_point_dots(&mut self, value: bool) -> &mut Self {
        self.show_point_dots = value;
        self
    }

    /// Generates the plot without equal axes
    pub fn unequal_exes(&mut self, value: bool) -> &mut Self {
        self.unequal_exes = value;
        self
    }

    /// Specifies the plot range in 2D
    pub fn range_2d(&mut self, xmin: f64, xmax: f64, ymin: f64, ymax: f64) -> &mut Self {
        self.range_2d = Some((xmin, xmax, ymin, ymax));
        self
    }

    /// Specifies the plot range in 3D
    pub fn range_3d(&mut self, xmin: f64, xmax: f64, ymin: f64, ymax: f64, zmin: f64, zmax: f64) -> &mut Self {
        self.range_3d = Some((xmin, xmax, ymin, ymax, zmin, zmax));
        self
    }

    /// Specifies the figure size in points
    pub fn size(&mut self, width: f64, height: f64) -> &mut Self {
        self.figure_size = Some((width, height));
        self
    }

    /// Enables zooming a region of the plot (2D meshes only)
    ///
    /// # Input
    ///
    /// * `xmin, xmax, ymin, ymax` - the region to be zoomed
    /// * `u0, v0, w, h` - the position (normalized coordinates) and size (width,height) of the zoomed region
    pub fn zoom_2d(
        &mut self,
        xmin: f64,
        xmax: f64,
        ymin: f64,
        ymax: f64,
        u0: f64,
        v0: f64,
        w: f64,
        h: f64,
    ) -> &mut Self {
        self.zoom_2d = Some(((xmin, xmax, ymin, ymax), (u0, v0, w, h)));
        self
    }

    /// Disables the mesh zooming
    pub fn disable_zoom_2d(&mut self) -> &mut Self {
        self.zoom_2d = None;
        self
    }

    /// Configures some characteristics of the zoom indicator
    ///
    /// # Input
    ///
    /// * `color` - color of the zoom indicator
    /// * `alpha` - transparency of the zoom indicator
    /// * `linewidth` - width of the zoom indicator line
    pub fn zoom_indicator(&mut self, color: &str, alpha: f64, linewidth: f64) -> &mut Self {
        self.zoom_indicator_config = (color.to_string(), alpha, linewidth);
        self
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
    /// **Note:** Non-zero point markers are shown within parentheses.
    pub fn draw_point_ids(&self, fig: &mut Figure) {
        if self.ndim == 2 {
            self.points.iter().for_each(|point| {
                let msg = if point.marker != 0 && fig.show_point_marker {
                    format!("{}({})", point.id, point.marker)
                } else {
                    format!("{}", point.id)
                };
                fig.canvas_point_ids.draw(point.coords[0], point.coords[1], &msg);
            });
        } else {
            self.points.iter().for_each(|point| {
                let msg = if point.marker != 0 && fig.show_point_marker {
                    format!("{}({})", point.id, point.marker)
                } else {
                    format!("{}", point.id)
                };
                fig.canvas_point_ids
                    .draw_3d(point.coords[0], point.coords[1], point.coords[2], &msg);
            });
        }
        fig.plot.add(&fig.canvas_point_ids);
    }

    /// Draws cells
    pub fn draw_cells(&self, fig: &mut Figure, set_range: bool) -> Result<(), StrError> {
        // limits
        let mut xmin = vec![f64::MAX; self.ndim];
        let mut xmax = vec![f64::MIN; self.ndim];

        // loop over cells
        for cell_id in 0..self.cells.len() {
            let canvas = if self.cells[cell_id].kind.is_lin() {
                &mut fig.canvas_lin_cells
            } else {
                &mut fig.canvas_cells
            };
            self.draw_cell(canvas, cell_id)?;
            if set_range {
                for point_id in &self.cells[cell_id].points {
                    for j in 0..self.ndim {
                        xmin[j] = f64::min(xmin[j], self.points[*point_id].coords[j]);
                        xmax[j] = f64::max(xmax[j], self.points[*point_id].coords[j]);
                    }
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
            let msg = if fig.show_cell_att {
                format!("{}({})", cell.id, cell.attribute)
            } else {
                format!("{}", cell.id)
            };
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

    /// Draws the cells, points, attributes, and markers
    ///
    /// # Input
    ///
    /// * `fig` -- the Figure struct (optional => use default configuration)
    /// * `filepath` -- may be a String, &str, or Path
    /// * `extra` -- is a function `|plot, before| {}` to perform some {pre,post}-drawing on the plot area.
    ///   The two arguments of this function are:
    ///     * `plot: &mut Plot` -- the `plot` reference that can be used perform some extra drawings.
    ///     * `before: bool` -- **true** indicates that the function is being called before all other
    ///       drawing functions. Otherwise, **false* indicates that the function is being called after
    ///       all other drawing functions, and just before the `plot.save` call.
    ///   For example, use `|_, _| {}` to do nothing.
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::mesh::Samples;
    /// use gemlab::StrError;
    /// use plotpy::Canvas;
    ///
    /// const SAVE_FIGURE: bool = false;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     // use sample mesh
    ///     let mesh = Samples::ring_eight_qua8_rad1_thick1();
    ///
    ///     // configure circles
    ///     let mut circle_in = Canvas::new();
    ///     let mut circle_out = Canvas::new();
    ///     circle_in
    ///         .set_face_color("None")
    ///         .set_edge_color("#ff000080")
    ///         .set_line_width(7.0)
    ///         .draw_circle(0.0, 0.0, 1.0);
    ///     circle_out
    ///         .set_face_color("None")
    ///         .set_edge_color("#0000ff80")
    ///         .set_line_width(7.0)
    ///         .draw_circle(0.0, 0.0, 2.0);
    ///
    ///     // draw mesh elements and circles
    ///     if SAVE_FIGURE {
    ///         let filename = "/tmp/gemlab/doc_example_mesh_draw.svg";
    ///         mesh.draw(None, filename, |plot, before| {
    ///             if !before {
    ///                 plot.add(&circle_in);
    ///                 plot.add(&circle_out);
    ///             }
    ///         })?;
    ///     }
    ///     Ok(())
    /// }
    /// ```
    ///
    /// ![doc_example_mesh_draw](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/doc_example_mesh_draw.svg)
    pub fn draw<P, F>(&self, fig: Option<Figure>, filepath: &P, mut extra: F) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
        F: FnMut(&mut Plot, bool),
    {
        let mut figure = if let Some(f) = fig { f } else { Figure::new() };
        extra(&mut figure.plot, true);
        self.draw_cells(&mut figure, true)?;
        if figure.show_cell_ids {
            self.draw_cell_ids(&mut figure)?;
        }
        if figure.show_point_dots {
            self.draw_point_dots(&mut figure);
        }
        if figure.show_point_ids {
            self.draw_point_ids(&mut figure);
        }
        if self.ndim == 2 {
            figure.plot.grid_and_labels("x", "y");
        }
        if !figure.unequal_exes {
            figure.plot.set_equal_axes(true);
        }
        if self.ndim == 2 {
            if let Some((xmin, xmax, ymin, ymax)) = figure.range_2d {
                figure.plot.set_range(xmin, xmax, ymin, ymax);
            }
        } else {
            if let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = figure.range_3d {
                figure.plot.set_range_3d(xmin, xmax, ymin, ymax, zmin, zmax);
            }
        }
        if let Some((width, height)) = figure.figure_size {
            figure.plot.set_figure_size_points(width, height);
        }
        extra(&mut figure.plot, false);
        if let Some(((xmin, xmax, ymin, ymax), (u0, v0, w, h))) = figure.zoom_2d {
            let mut inset = InsetAxes::new();
            inset
                .set_indicator_line_color(&figure.zoom_indicator_config.0)
                .set_indicator_alpha(figure.zoom_indicator_config.1)
                .set_indicator_line_width(figure.zoom_indicator_config.2);
            inset.add(&figure.canvas_cells);
            if figure.show_cell_ids {
                inset.add(&figure.canvas_cell_ids);
            }
            if figure.show_point_dots {
                inset.add(&figure.canvas_points);
            }
            if figure.show_point_ids {
                inset.add(&figure.canvas_point_ids);
            }
            inset.set_range(xmin, xmax, ymin, ymax).draw(u0, v0, w, h);
            figure.plot.add(&inset);
        }
        figure.plot.save(filepath)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Figure;
    use crate::mesh::{Mesh, Samples};
    use plotpy::{Canvas, Plot, Text};

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

    #[test]
    fn draw_cells_and_points_work() {
        // lin cells ---------------------------------------------------------------------------
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

        // lin cells in 3d ---------------------------------------------------------------------
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

        // tri cells ---------------------------------------------------------------------------
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

        // qua cells ---------------------------------------------------------------------------
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

        // tet cells ---------------------------------------------------------------------------
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

        // hex cells ---------------------------------------------------------------------------
        let mesh = Samples::hex_cells();
        let mut fig = Figure::new();
        mesh.draw_cells(&mut fig, true).unwrap();
        mesh.draw_point_dots(&mut fig);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut fig.plot, &mut labels, &mesh, 0.03, 0.0, 0.03);
            caption.draw_3d(0.5, -0.25, 0.0, "Hex8");
            caption.draw_3d(2.3, -0.25, 0.0, "Hex20");
            caption.draw_3d(0.3, 1.1, 2.0, "Hex32");
            fig.plot.add(&caption);
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range_3d(0.0, 3.0, 0.0, 3.0, 0.0, 2.5)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_5_hex.svg")
                .unwrap();
        }

        // ring --------------------------------------------------------------------------------
        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        let mut fig = Figure::new();
        mesh.draw_cells(&mut fig, true).unwrap();
        mesh.draw_point_dots(&mut fig);

        if SAVE_FIGURE {
            let mut circle_in = Canvas::new();
            let mut circle_mi = Canvas::new();
            let mut circle_ou = Canvas::new();
            circle_in
                .set_face_color("None")
                .set_edge_color("#bfbfbfbb")
                .set_line_width(7.0)
                .draw_circle(0.0, 0.0, 1.0);
            circle_mi
                .set_face_color("None")
                .set_edge_color("#bfbfbfbb")
                .set_line_width(7.0)
                .draw_circle(0.0, 0.0, 1.5);
            circle_ou
                .set_face_color("None")
                .set_edge_color("#bfbfbfbb")
                .set_line_width(7.0)
                .draw_circle(0.0, 0.0, 2.0);
            fig.plot.add(&circle_in);
            fig.plot.add(&circle_mi);
            fig.plot.add(&circle_ou);
            fig.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range(-0.1, 2.1, -0.1, 2.1)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_6_ring.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_works_qua12() {
        if SAVE_FIGURE {
            let mesh = Samples::block_2d_four_qua12();
            let mut fig = Figure::new();
            fig.show_cell_ids(true)
                .show_point_ids(true)
                .show_point_dots(true)
                .range_2d(-0.5, 6.0, -0.5, 6.0)
                .zoom_2d(-0.05, 1.55, -0.05, 1.55, 0.6, 0.6, 0.3, 0.3);
            mesh.draw(Some(fig), "/tmp/gemlab/test_draw_works_qua12.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    fn draw_works_qua16() {
        if SAVE_FIGURE {
            let mesh = Samples::block_2d_four_qua16();
            let mut fig = Figure::new();
            fig.show_cell_ids(true).show_point_ids(true).show_point_dots(true);
            mesh.draw(Some(fig), "/tmp/gemlab/test_draw_works_qua16.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    fn draw_works_qua17() {
        if SAVE_FIGURE {
            let mesh = Samples::block_2d_four_qua17();
            let mut fig = Figure::new();
            fig.show_cell_ids(true).show_point_ids(true).show_point_dots(true);
            mesh.draw(Some(fig), "/tmp/gemlab/test_draw_works_qua17.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    fn draw_works_mixed_2d() {
        if SAVE_FIGURE {
            let mesh = Samples::mixed_shapes_2d();
            let mut fig = Figure::new();
            fig.show_cell_ids(true).show_point_ids(true).show_point_dots(true);
            mesh.draw(Some(fig), "/tmp/gemlab/test_draw_works_mixed_2d.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    fn draw_works_hex8() {
        if SAVE_FIGURE {
            let mesh = Samples::two_hex8();
            let mut fig = Figure::new();
            fig.show_cell_ids(true)
                .show_point_ids(true)
                .show_point_dots(true)
                .size(600.0, 600.0);
            fig.canvas_point_ids()
                .set_align_horizontal("left")
                .set_align_vertical("bottom")
                .set_color("black")
                .set_fontsize(10.0)
                .set_bbox_facecolor("gold")
                .set_bbox_alpha(0.5);
            mesh.draw(Some(fig), "/tmp/gemlab/test_draw_works_hex8.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    fn draw_works_hex20() {
        if SAVE_FIGURE {
            let mesh = Samples::block_3d_eight_hex20();
            let mut fig = Figure::new();
            fig.show_cell_ids(true)
                .show_point_ids(true)
                .show_point_dots(true)
                .size(600.0, 600.0);
            fig.canvas_point_ids()
                .set_align_horizontal("left")
                .set_align_vertical("bottom")
                .set_color("black")
                .set_fontsize(10.0)
                .set_bbox_facecolor("gold")
                .set_bbox_alpha(0.5);
            mesh.draw(Some(fig), "/tmp/gemlab/test_draw_works_hex20.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    fn draw_works_mixed_3d() {
        if SAVE_FIGURE {
            let mesh = Samples::mixed_shapes_3d();
            let mut fig = Figure::new();
            fig.show_cell_ids(true)
                .show_point_ids(true)
                .show_point_dots(true)
                .size(600.0, 600.0);
            fig.canvas_point_ids()
                .set_align_horizontal("left")
                .set_align_vertical("bottom")
                .set_color("black")
                .set_fontsize(10.0)
                .set_bbox_facecolor("gold")
                .set_bbox_alpha(0.5);
            mesh.draw(Some(fig), "/tmp/gemlab/test_works_mixed_3d.svg", |_, _| {})
                .unwrap();
        }
    }

    #[test]
    fn test_figure_setters() {
        let mut fig = Figure::new();
        
        // Test show_cell_ids
        assert!(!fig.show_cell_ids);
        fig.show_cell_ids(true);
        assert!(fig.show_cell_ids);
        
        // Test show_cell_att
        assert!(fig.show_cell_att);
        fig.show_cell_att(false);
        assert!(!fig.show_cell_att);
        
        // Test show_point_ids
        assert!(!fig.show_point_ids);
        fig.show_point_ids(true);
        assert!(fig.show_point_ids);
        
        // Test show_point_marker
        assert!(!fig.show_point_marker);
        fig.show_point_marker(true);
        assert!(fig.show_point_marker);
        
        // Test show_point_dots
        assert!(!fig.show_point_dots);
        fig.show_point_dots(true);
        assert!(fig.show_point_dots);
        
        // Test unequal_axes
        assert!(!fig.unequal_exes);
        fig.unequal_exes(true);
        assert!(fig.unequal_exes);
        
        // Test range_2d
        assert!(fig.range_2d.is_none());
        fig.range_2d(-1.0, 1.0, -2.0, 2.0);
        assert_eq!(fig.range_2d, Some((-1.0, 1.0, -2.0, 2.0)));
        
        // Test range_3d
        assert!(fig.range_3d.is_none());
        fig.range_3d(-1.0, 1.0, -2.0, 2.0, -3.0, 3.0);
        assert_eq!(fig.range_3d, Some((-1.0, 1.0, -2.0, 2.0, -3.0, 3.0)));
        
        // Test size
        assert!(fig.figure_size.is_none());
        fig.size(800.0, 600.0);
        assert_eq!(fig.figure_size, Some((800.0, 600.0)));
        
        // Test zoom_2d
        assert!(fig.zoom_2d.is_none());
        fig.zoom_2d(-0.5, 0.5, -1.0, 1.0, 0.1, 0.2, 0.3, 0.4);
        assert_eq!(fig.zoom_2d, Some(((-0.5, 0.5, -1.0, 1.0), (0.1, 0.2, 0.3, 0.4))));
        
        // Test disable_zoom_2d
        fig.disable_zoom_2d();
        assert!(fig.zoom_2d.is_none());
        
        // Test zoom_indicator
        fig.zoom_indicator("#ff0000", 0.8, 2.5);
        assert_eq!(fig.zoom_indicator_config, ("#ff0000".to_string(), 0.8, 2.5));
        
        // Test method chaining
        let mut fig = Figure::new();
        fig.show_cell_ids(true)
           .show_point_ids(true)
           .show_point_dots(true)
           .size(800.0, 600.0)
           .range_2d(-1.0, 1.0, -1.0, 1.0);
        assert!(fig.show_cell_ids);
        assert!(fig.show_point_ids);
        assert!(fig.show_point_dots);
        assert_eq!(fig.figure_size, Some((800.0, 600.0)));
        assert_eq!(fig.range_2d, Some((-1.0, 1.0, -1.0, 1.0)));
    }
}
