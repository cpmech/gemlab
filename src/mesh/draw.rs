use super::{Features, Mesh};
use crate::mesh::Triangulation;
use crate::shapes::{GeoClass, GeoKind, Scratchpad};
use crate::StrError;
use plotpy::{Canvas, Curve, InsetAxes, Plot, Text};
use russell_lab::math::PI;
use russell_lab::{sort2, sort4, Vector};
use std::ffi::OsStr;

/// Implements functions to draw cells, edges, and faces or the whole mesh
pub struct Draw<'a> {
    /// The plotpy structure to draw figures (plots)
    plot: Plot,

    /// Multiplier used to set the drawing area range
    m_range: f64,

    /// Tolerance to decide if an edge is horizontal or vertical
    tol_edge_marker: f64,

    /// Multiplier to scale the length of the normal vectors
    m_normal_vector: f64,

    /// Multiplier to scale the length of the normal vectors for markers
    m_normal_vector_marker: f64,

    /// Canvas to draw edges
    canvas_edges: Canvas,

    /// Canvas to draw points (dot indicators)
    canvas_point_dots: Curve,

    /// Canvas to draw point ids
    canvas_point_ids: Text,

    /// Canvas to draw cell ids
    canvas_cell_ids: Text,

    /// Canvas to draw cells
    canvas_cells: Canvas,

    /// Canvas to draw lin cells
    canvas_lin_cells: Canvas,

    /// Canvas to draw shells (triangles and quadrilaterals in 3D)
    canvas_shells: Canvas,

    /// Canvas to draw edge markers
    canvas_edge_markers: Text,

    /// Canvas to draw face markers lines
    canvas_edge_markers_lines: Canvas,

    /// Canvas to draw face markers
    canvas_face_markers: Text,

    /// Canvas to draw face markers lines
    canvas_face_markers_lines: Canvas,

    /// Canvas to draw normal vectors (2D)
    canvas_normals_2d: Canvas,

    /// Canvas to draw normal vectors (3D)
    canvas_normals_3d: Canvas,

    /// Canvas to draw boundary edges in 3D
    canvas_boundary_edges_3d: Canvas,

    /// Canvas to draw boundary faces (3D)
    canvas_boundary_faces: Canvas,

    /// Canvas to draw the glyph indicating the X-Y-Z directions in 3D plots
    canvas_glyph_3d: Canvas,

    /// Shows cells (wireframe in 3D)
    show_cells: bool,

    /// Shows cell ids
    show_cell_ids: bool,

    /// Shows cell marker within parenthesis
    show_cell_marker: bool,

    /// Shows point ids
    show_point_ids: bool,

    /// Shows point marker within parenthesis (if not zero)
    show_point_marker: bool,

    /// Shows point dots
    show_point_dots: bool,

    /// Shows edge markers
    show_edge_markers: bool,

    /// Shows face markers
    show_face_markers: bool,

    /// Shows normal vectors on boundaries
    show_normal_vectors: bool,

    /// Draws boundary edges (3D)
    show_boundary_edges_3d: bool,

    /// Draws boundary faces (3D)
    show_boundary_faces: bool,

    /// Shows the glyph indicating the X-Y-Z directions in 3D plots
    show_glyph_3d: bool,

    /// View the figure interactively in a window
    view: bool,

    /// Generates the plot without equal axes
    unequal_exes: bool,

    /// Specifies the camera elevation for 3D plots (default: 30 deg)
    camera_elevation: f64,

    /// Specifies the camera azimuth for 3D plots (default: 30 deg)
    camera_azimuth: f64,

    /// Hides the 3D grid and panes (make them transparent)
    hide_3d_grid: bool,

    /// Hides the axes (2D or 3D)
    hide_axes: bool,

    /// Origin point of the 3D glyph
    glyph_3d_origin: [f64; 3],

    /// Size of the 3D glyph
    glyph_3d_size: f64,

    /// Specifies the plot range (xmin, xmax, ymin, ymax)
    range_2d: Option<(f64, f64, f64, f64)>,

    /// Specifies the plot range (xmin, xmax, ymin, ymax, zmin, zmax)
    range_3d: Option<(f64, f64, f64, f64, f64, f64)>,

    /// Specifies the figure size in points
    size: Option<(f64, f64)>,

    /// Holds a (callback) function to perform some extra drawing on the plot area
    ///
    /// The function is `|plot, before| {}` where:
    ///
    /// * `plot` -- the `plot` reference that can be used perform some extra drawings.
    /// * `before` -- **true** indicates that the callback function is being called before all other drawing functions.
    extra: Option<Box<dyn Fn(&mut Plot, bool) + 'a>>,

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

    /// Holds a (callback) function to perform some extra drawing on the zoom area
    ///
    /// The function is `|inset| {}` where:
    ///
    /// * `inset` -- Is the InsetAxes reference for further configuration
    zoom_extra: Option<Box<dyn Fn(&mut InsetAxes) + 'a>>,
}

impl<'a> Draw<'a> {
    /// Allocates a new instance
    pub fn new() -> Self {
        let mut canvas_edges = Canvas::new();
        let mut canvas_point_dots = Curve::new();
        let mut canvas_point_ids = Text::new();
        let mut canvas_cell_ids = Text::new();
        let mut canvas_cells = Canvas::new();
        let mut canvas_lin_cells = Canvas::new();
        let mut canvas_shells = Canvas::new();
        let mut canvas_edge_markers = Text::new();
        let mut canvas_edge_markers_lines = Canvas::new();
        let mut canvas_face_markers = Text::new();
        let mut canvas_face_markers_lines = Canvas::new();
        let mut canvas_normals_2d = Canvas::new();
        let mut canvas_normals_3d = Canvas::new();
        let mut canvas_boundary_edges_3d = Canvas::new();
        let mut canvas_boundary_faces = Canvas::new();
        let canvas_glyph_3d = Canvas::new();
        canvas_edges
            .set_face_color("None")
            .set_line_width(1.0)
            .set_edge_color("#2440cd");
        canvas_point_dots
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
        canvas_shells
            .set_edge_color("None")
            .set_face_color("#32b31880")
            .set_line_width(1.0);
        canvas_edge_markers
            .set_color("#000000ff")
            .set_fontsize(9.0)
            .set_align_horizontal("center")
            .set_align_vertical("center")
            .set_bbox(true)
            .set_bbox_facecolor("#fff0a3ff")
            .set_bbox_edgecolor("black")
            .set_bbox_style("square,pad=0.15");
        canvas_edge_markers_lines
            .set_edge_color("#000000ff")
            .set_line_width(1.0);
        canvas_face_markers
            .set_color("#000000ff")
            .set_fontsize(9.0)
            .set_align_horizontal("center")
            .set_align_vertical("center")
            .set_bbox(true)
            .set_bbox_facecolor("#aaffb5ff")
            .set_bbox_edgecolor("black")
            .set_bbox_style("square,pad=0.15");
        canvas_face_markers_lines
            .set_face_color("None")
            .set_edge_color("#000000ff")
            .set_line_width(1.0);
        canvas_normals_2d
            .set_face_color("None")
            .set_edge_color("#f400f4ff")
            .set_line_width(2.0)
            .set_arrow_scale(5.0)
            .set_arrow_style("->");
        canvas_normals_3d
            .set_face_color("None")
            .set_edge_color("#f400f4ff")
            .set_line_width(2.0);
        canvas_boundary_edges_3d.set_edge_color("#001a9eff").set_line_width(2.0);
        canvas_boundary_faces.set_edge_color("None").set_face_color("#0095ff80");
        Draw {
            plot: Plot::new(),
            m_range: 0.2,
            tol_edge_marker: 1e-3,
            m_normal_vector: 0.05,
            m_normal_vector_marker: 0.1,
            canvas_edges,
            canvas_point_dots,
            canvas_point_ids,
            canvas_cell_ids,
            canvas_cells,
            canvas_lin_cells,
            canvas_shells,
            canvas_edge_markers,
            canvas_edge_markers_lines,
            canvas_face_markers,
            canvas_face_markers_lines,
            canvas_normals_2d,
            canvas_normals_3d,
            canvas_boundary_edges_3d,
            canvas_boundary_faces,
            canvas_glyph_3d,
            show_cells: true,
            show_cell_ids: false,
            show_cell_marker: true,
            show_point_ids: false,
            show_point_marker: false,
            show_point_dots: false,
            show_edge_markers: false,
            show_face_markers: false,
            show_normal_vectors: false,
            show_boundary_edges_3d: true,
            show_boundary_faces: true,
            show_glyph_3d: false,
            view: false,
            unequal_exes: false,
            camera_elevation: 30.0,
            camera_azimuth: 30.0,
            hide_3d_grid: true,
            hide_axes: false,
            glyph_3d_size: 1.0,
            glyph_3d_origin: [0.0, 0.0, 0.0],
            range_2d: None,
            range_3d: None,
            size: None,
            extra: None,
            zoom_2d: None,
            zoom_indicator_config: ("#f9d835".to_string(), 1.0, 2.0),
            zoom_extra: None,
        }
    }

    /// Get a mutable reference to the canvas to draw edges
    pub fn get_canvas_edges(&mut self) -> &mut Canvas {
        &mut self.canvas_edges
    }

    /// Get a mutable reference to the canvas to draw points
    pub fn get_canvas_points(&mut self) -> &mut Curve {
        &mut self.canvas_point_dots
    }

    /// Get a mutable reference to the canvas to draw point ids
    pub fn get_canvas_point_ids(&mut self) -> &mut Text {
        &mut self.canvas_point_ids
    }

    /// Get a mutable reference to the canvas to draw cell ids
    pub fn get_canvas_cell_ids(&mut self) -> &mut Text {
        &mut self.canvas_cell_ids
    }

    /// Get a mutable reference to the canvas to draw cells
    pub fn get_canvas_cells(&mut self) -> &mut Canvas {
        &mut self.canvas_cells
    }

    /// Get a mutable reference to the canvas to draw lin cells
    pub fn get_canvas_lin_cells(&mut self) -> &mut Canvas {
        &mut self.canvas_lin_cells
    }

    /// Get a mutable reference to the canvas to draw edge markers
    pub fn get_canvas_edge_markers(&mut self) -> &mut Text {
        &mut self.canvas_edge_markers
    }

    /// Get a mutable reference to the canvas to draw edge markers lines
    pub fn get_canvas_edge_markers_lines(&mut self) -> &mut Canvas {
        &mut self.canvas_edge_markers_lines
    }

    /// Get a mutable reference to the canvas to draw face markers
    pub fn get_canvas_face_markers(&mut self) -> &mut Text {
        &mut self.canvas_face_markers
    }

    /// Get a mutable reference to the canvas to draw face markers lines
    pub fn get_canvas_face_markers_lines(&mut self) -> &mut Canvas {
        &mut self.canvas_face_markers_lines
    }

    /// Get a mutable reference to the canvas to draw normals 2D
    pub fn get_canvas_normals_2d(&mut self) -> &mut Canvas {
        &mut self.canvas_normals_2d
    }

    /// Get a mutable reference to the canvas to draw normals 3D
    pub fn get_canvas_normals_3d(&mut self) -> &mut Canvas {
        &mut self.canvas_normals_3d
    }

    /// Get a mutable reference to the canvas tro draw boundary faces (3D)
    pub fn get_canvas_boundary_faces(&mut self) -> &mut Canvas {
        &mut self.canvas_boundary_faces
    }

    /// Sets the multiplier used to set the drawing area range
    ///
    /// Default: `0.2`
    pub fn set_m_range(&mut self, value: f64) -> &mut Self {
        self.m_range = value;
        self
    }

    /// Sets the tolerance to decide if an edge is horizontal or vertical
    ///
    /// Default: `1e-3`
    pub fn set_tol_edge_marker(&mut self, value: f64) -> &mut Self {
        self.tol_edge_marker = value;
        self
    }

    /// Sets the multiplier to scale the length of the normal vectors
    ///
    /// Default: `0.05`
    pub fn set_m_normal_vector(&mut self, value: f64) -> &mut Self {
        self.m_normal_vector = value;
        self
    }

    /// Sets the multiplier to scale the length of the normal vectors for markers
    ///
    /// Default: `0.1`
    pub fn set_m_normal_vector_marker(&mut self, value: f64) -> &mut Self {
        self.m_normal_vector_marker = value;
        self
    }

    /// Shows cells (wireframe in 3D)
    ///
    /// Default: `true`
    pub fn show_cells(&mut self, value: bool) -> &mut Self {
        self.show_cells = value;
        self
    }

    /// Shows cell ids
    ///
    /// Default: `false`
    pub fn show_cell_ids(&mut self, value: bool) -> &mut Self {
        self.show_cell_ids = value;
        self
    }

    /// Shows cell marker within parenthesis
    ///
    /// Default: `true`
    pub fn show_cell_marker(&mut self, value: bool) -> &mut Self {
        self.show_cell_marker = value;
        self
    }

    /// Shows point ids
    ///
    /// Default: `false`
    pub fn show_point_ids(&mut self, value: bool) -> &mut Self {
        self.show_point_ids = value;
        self
    }

    /// Shows point marker within parenthesis (if not zero)
    ///
    /// Default: `false`
    pub fn show_point_marker(&mut self, value: bool) -> &mut Self {
        self.show_point_marker = value;
        self
    }

    /// Shows point dots
    ///
    /// Default: `false`
    pub fn show_point_dots(&mut self, value: bool) -> &mut Self {
        self.show_point_dots = value;
        self
    }

    /// Shows edge markers
    ///
    /// Default: `false`
    pub fn show_edge_markers(&mut self, value: bool) -> &mut Self {
        self.show_edge_markers = value;
        self
    }

    /// Shows face markers
    ///
    /// Default: `false`
    pub fn show_face_markers(&mut self, value: bool) -> &mut Self {
        self.show_face_markers = value;
        self
    }

    /// Shows normal vectors on boundaries
    ///
    /// Default: `false`
    pub fn show_normal_vectors(&mut self, value: bool) -> &mut Self {
        self.show_normal_vectors = value;
        self
    }

    /// Shows boundary edges (3D)
    ///
    /// Default: `true`
    pub fn show_boundary_edges_3d(&mut self, value: bool) -> &mut Self {
        self.show_boundary_edges_3d = value;
        self
    }

    /// Shows boundary faces (3D)
    ///
    /// Default: `true`
    pub fn show_boundary_faces(&mut self, value: bool) -> &mut Self {
        self.show_boundary_faces = value;
        self
    }

    /// Shows the glyph indicating the X-Y-Z directions in 3D plots
    ///
    /// Default: `false`
    pub fn show_glyph_3d(&mut self, value: bool) -> &mut Self {
        self.show_glyph_3d = value;
        self
    }

    /// View the figure interactively in a window
    pub fn set_view_flag(&mut self, value: bool) -> &mut Self {
        self.view = value;
        self
    }

    /// Generates the plot without equal axes
    ///
    /// Default: `false`
    pub fn set_unequal_exes(&mut self, value: bool) -> &mut Self {
        self.unequal_exes = value;
        self
    }

    /// Sets the camera position for 3D plots (azimuth angles in degrees)
    pub fn set_camera(&mut self, elevation: f64, azimuth: f64) -> &mut Self {
        self.camera_elevation = elevation;
        self.camera_azimuth = azimuth;
        self
    }

    /// Sets an option to hide the 3D grid and panes (make them transparent)
    ///
    /// Default: `true`
    pub fn set_hide_3d_grid(&mut self, value: bool) -> &mut Self {
        self.hide_3d_grid = value;
        self
    }

    /// Sets whether to hide the axes (2D or 3D)
    pub fn set_hide_axes(&mut self, value: bool) -> &mut Self {
        self.hide_axes = value;
        self
    }

    /// Specifies the origin point of the 3D glyph and its size
    pub fn set_glyph_3d(&mut self, x: f64, y: f64, z: f64, size: f64) -> &mut Self {
        self.glyph_3d_origin = [x, y, z];
        self.glyph_3d_size = size;
        self
    }

    /// Specifies the plot range in 2D
    ///
    /// Default: `None` (automatic range)
    pub fn set_range_2d(&mut self, xmin: f64, xmax: f64, ymin: f64, ymax: f64) -> &mut Self {
        self.range_2d = Some((xmin, xmax, ymin, ymax));
        self
    }

    /// Specifies the plot range in 3D
    ///
    /// Default: `None` (automatic range)
    pub fn set_range_3d(&mut self, xmin: f64, xmax: f64, ymin: f64, ymax: f64, zmin: f64, zmax: f64) -> &mut Self {
        self.range_3d = Some((xmin, xmax, ymin, ymax, zmin, zmax));
        self
    }

    /// Specifies the figure size in points
    ///
    /// Default: `None` (automatic size)
    pub fn set_size(&mut self, width: f64, height: f64) -> &mut Self {
        self.size = Some((width, height));
        self
    }

    /// Sets a callback function to perform extra drawing on the plot area
    ///
    /// Default: `None`
    ///
    /// The function is `|plot, before| { ... }` where:
    ///
    /// * `plot` -- is the `plot` mutable reference that can be used perform some extra drawings.
    /// * `before` -- indicates if the function is being called before all other drawing functions.
    ///   The value `false` means that the function is being called `after` all other functions,
    ///   and just before the call to `plot.save()`.
    pub fn extra(&mut self, callback: impl Fn(&mut Plot, bool) + 'a) -> &mut Self {
        self.extra = Some(Box::new(callback));
        self
    }

    /// Enables zooming a region of the plot (2D meshes only)
    ///
    /// # Input
    ///
    /// * `xmin, xmax, ymin, ymax` - the region to be zoomed
    /// * `u0, v0, w, h` - the position (normalized coordinates) and size (width,height) of the zoomed region
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::mesh::{Draw, Samples};
    /// use gemlab::StrError;
    /// use plotpy::Text;
    ///
    /// const SAVE_FIGURE: bool = false;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     if SAVE_FIGURE {
    ///         let mesh = Samples::block_2d_four_qua12();
    ///         let mut draw = Draw::new();
    ///         draw.show_cell_ids(true)
    ///             .show_point_ids(true)
    ///             .show_point_dots(true)
    ///             .set_range_2d(-0.5, 6.0, -0.5, 6.0)
    ///             .zoom_2d(-0.05, 1.55, -0.05, 1.55, 0.6, 0.6, 0.3, 0.3)
    ///             .zoom_extra(|inset| {
    ///                 let mut text = Text::new();
    ///                 text.draw(0.3, 1.0, "HELLO");
    ///                 inset.add(&text);
    ///             })
    ///             .all(&mesh, "/tmp/gemlab/doc_draw_works_qua12.svg")?;
    ///     }
    ///     Ok(())
    /// }
    /// ```
    ///
    /// ![doc_draw_works_qua12](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/doc_draw_works_qua12.svg)
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
    ///
    /// Default: `None` (no zoom)
    pub fn disable_zoom_2d(&mut self) -> &mut Self {
        self.zoom_2d = None;
        self
    }

    /// Configures some characteristics of the zoom indicator
    ///
    /// Default: `("#f9d835", 1.0, 2.0)`
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

    /// Sets a callback function to perform extra drawing on the zoom area
    ///
    /// Default: `None`
    ///
    /// The function is `|inset| {}` where:
    ///
    /// * `inset` -- Is the InsetAxes reference for further configuration
    ///
    /// # Examples
    ///
    /// See [zoom_2d](Self::zoom_2d) method for an example.
    pub fn zoom_extra(&mut self, callback: impl Fn(&mut InsetAxes) + 'a) -> &mut Self {
        self.zoom_extra = Some(Box::new(callback));
        self
    }

    /// Draws all points (dots)
    pub fn point_dots(&mut self, mesh: &Mesh) {
        if mesh.ndim == 2 {
            self.canvas_point_dots.points_begin();
            mesh.points.iter().for_each(|point| {
                self.canvas_point_dots.points_add(point.coords[0], point.coords[1]);
            });
            self.canvas_point_dots.points_end();
        } else {
            self.canvas_point_dots.points_3d_begin();
            mesh.points.iter().for_each(|point| {
                self.canvas_point_dots
                    .points_3d_add(point.coords[0], point.coords[1], point.coords[2]);
            });
            self.canvas_point_dots.points_3d_end();
        }
        self.plot.add(&self.canvas_point_dots);
    }

    /// Draws all point ids (labels)
    ///
    /// **Note:** Non-zero point markers are shown within parentheses.
    pub fn point_ids(&mut self, mesh: &Mesh) {
        if mesh.ndim == 2 {
            mesh.points.iter().for_each(|point| {
                let msg = if point.marker != 0 && self.show_point_marker {
                    format!("{}({})", point.id, point.marker)
                } else {
                    format!("{}", point.id)
                };
                self.canvas_point_ids.draw(point.coords[0], point.coords[1], &msg);
            });
        } else {
            mesh.points.iter().for_each(|point| {
                let msg = if point.marker != 0 && self.show_point_marker {
                    format!("{}({})", point.id, point.marker)
                } else {
                    format!("{}", point.id)
                };
                self.canvas_point_ids
                    .draw_3d(point.coords[0], point.coords[1], point.coords[2], &msg);
            });
        }
        self.plot.add(&self.canvas_point_ids);
    }

    /// Draws cells
    pub fn cells(&mut self, mesh: &Mesh) -> Result<(), StrError> {
        // holds the IDs of shell cells
        let mut shell_cell_ids = Vec::new();

        // loop over cells (draws the wireframe if 3D)
        for cell_id in 0..mesh.cells.len() {
            if mesh.ndim == 3 {
                let class = mesh.cells[cell_id].kind.class();
                if class == GeoClass::Tri || class == GeoClass::Qua {
                    shell_cell_ids.push(cell_id);
                }
            };
            let canvas = if mesh.cells[cell_id].kind.is_lin() {
                &mut self.canvas_lin_cells
            } else {
                &mut self.canvas_cells
            };
            mesh.draw_cell(canvas, cell_id)?;
        }

        // loop over shell cells to draw the 3D surface (the wireframe is already drawn)
        if shell_cell_ids.len() > 0 {
            let surface: Vec<_> = shell_cell_ids.iter().map(|&id| &mesh.cells[id]).collect();
            let res = Triangulation::from_surface(mesh, &surface);
            self.canvas_shells
                .draw_triangles_3d(&res.xx, &res.yy, &res.zz, &res.triangles);
            self.plot.add(&self.canvas_shells);
        }

        // add (wireframe) to plot
        self.plot.add(&self.canvas_cells);
        self.plot.add(&self.canvas_lin_cells);
        Ok(())
    }

    /// Draws ids and markers of cells
    ///
    /// **Note:** Cell markers are shown within parentheses.
    pub fn cell_ids(&mut self, mesh: &Mesh) -> Result<(), StrError> {
        // auxiliary
        let mut x = Vector::new(mesh.ndim);

        // loop over all cells
        for cell_id in 0..mesh.cells.len() {
            // compute coordinates of the label
            let cell = &mesh.cells[cell_id];
            x.fill(0.0);
            if cell.kind == GeoKind::Qua9 || cell.kind == GeoKind::Qua17 {
                for i in 0..mesh.ndim {
                    x[i] = (3.0 * mesh.points[cell.points[0]].coords[i] + mesh.points[cell.points[2]].coords[i]) / 4.0;
                }
            } else if cell.kind == GeoKind::Tri10 {
                for i in 0..mesh.ndim {
                    x[i] = (mesh.points[cell.points[0]].coords[i] + mesh.points[cell.points[9]].coords[i]) / 2.0;
                }
            } else {
                for m in 0..cell.points.len() {
                    for i in 0..mesh.ndim {
                        x[i] += mesh.points[cell.points[m]].coords[i];
                    }
                }
                for i in 0..mesh.ndim {
                    x[i] /= cell.points.len() as f64;
                }
            }

            // add label
            let msg = if self.show_cell_marker {
                format!("{}({})", cell.id, cell.marker)
            } else {
                format!("{}", cell.id)
            };
            if mesh.ndim == 2 {
                self.canvas_cell_ids.draw(x[0], x[1], msg.as_str());
            } else {
                self.canvas_cell_ids.draw_3d(x[0], x[1], x[2], msg.as_str());
            }
        }

        // add to plot
        self.plot.add(&self.canvas_cell_ids);
        Ok(())
    }

    /// Draws edge markers
    pub fn edge_markers(&mut self, features: &Features) -> Result<(), StrError> {
        let mesh = &features.mesh;
        let ndim = mesh.ndim;
        let ksi = &[0.0, 0.0];
        let mut un = Vector::new(ndim);
        let mut x = Vector::new(ndim);
        let mut sum = 0.0;
        for i in 0..ndim {
            sum += f64::abs(features.max[i] - features.min[i]) * f64::abs(features.max[i] - features.min[i]);
        }
        let s = f64::sqrt(sum) * self.m_normal_vector_marker;
        for (marker, p1, p2) in &mesh.marked_edges {
            if ndim == 2 {
                let mut key = (*p1, *p2);
                sort2(&mut key);
                if let Some(edge) = features.edges.get(&key) {
                    let mut pad = Scratchpad::new(ndim, edge.kind)?;
                    mesh.set_pad(&mut pad, &edge.points);
                    pad.calc_normal_vector(&mut un, ksi)?;
                    pad.calc_coords(&mut x, ksi)?;
                    let alpha = if f64::abs(un[1]) < self.tol_edge_marker {
                        0.0
                    } else if f64::abs(un[0]) < self.tol_edge_marker {
                        90.0
                    } else {
                        f64::atan(un[1] / un[0]) * 180.0 / PI
                    };
                    self.canvas_edge_markers.set_rotation(alpha).draw(
                        x[0] + 0.5 * s * un[0],
                        x[1] + 0.5 * s * un[1],
                        &format!("{}", marker),
                    );
                    self.canvas_edge_markers_lines
                        .draw_polyline(&[[x[0], x[1]], [x[0] + s * un[0], x[1] + s * un[1]]], false);
                }
            } else {
                for i in 0..ndim {
                    x[i] = 0.5 * (mesh.points[*p1].coords[i] + mesh.points[*p2].coords[i]);
                }
                self.canvas_edge_markers
                    .draw_3d(x[0], x[1], x[2], &format!("{}", marker));
            }
        }
        self.plot.add(&self.canvas_edge_markers_lines);
        self.plot.add(&self.canvas_edge_markers);
        Ok(())
    }

    /// Draws face markers
    pub fn face_markers(&mut self, features: &Features) -> Result<(), StrError> {
        let mesh = &features.mesh;
        let ndim = mesh.ndim;
        if ndim == 2 {
            return Ok(());
        }
        let ksi3 = &[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
        let ksi4 = &[0.0, 0.0, 0.0];
        let mut un = Vector::new(ndim);
        let mut x = Vector::new(ndim);
        let mut sum = 0.0;
        for i in 0..ndim {
            sum += f64::abs(features.max[i] - features.min[i]) * f64::abs(features.max[i] - features.min[i]);
        }
        let s = f64::sqrt(sum) * self.m_normal_vector_marker;
        for (marker, p1, p2, p3, p4) in &mesh.marked_faces {
            let mut key = (*p1, *p2, *p3, *p4);
            sort4(&mut key);
            if let Some(face) = features.faces.get(&key) {
                let mut pad = Scratchpad::new(ndim, face.kind)?;
                mesh.set_pad(&mut pad, &face.points);
                if face.kind.class() == GeoClass::Tri {
                    pad.calc_normal_vector(&mut un, ksi3)?;
                    pad.calc_coords(&mut x, ksi3)?;
                } else {
                    pad.calc_normal_vector(&mut un, ksi4)?;
                    pad.calc_coords(&mut x, ksi4)?;
                }
                self.canvas_face_markers.draw_3d(
                    x[0] + 0.5 * s * un[0],
                    x[1] + 0.5 * s * un[1],
                    x[2] + 0.5 * s * un[2],
                    &format!("{}", marker),
                );
                self.canvas_face_markers_lines.polyline_3d_begin();
                self.canvas_face_markers_lines.polyline_3d_add(x[0], x[1], x[2]);
                self.canvas_face_markers_lines
                    .polyline_3d_add(x[0] + s * un[0], x[1] + s * un[1], x[2] + s * un[2]);
                self.canvas_face_markers_lines.polyline_3d_end();
            }
        }
        self.plot.add(&self.canvas_face_markers_lines);
        self.plot.add(&self.canvas_face_markers);
        Ok(())
    }

    /// Draws normal vectors on boundary edges, faces, and shells
    pub fn normal_vectors(&mut self, features: &Features) -> Result<(), StrError> {
        let ndim = features.mesh.ndim;
        let mut sum = 0.0;
        for i in 0..ndim {
            sum += f64::abs(features.max[i] - features.min[i]) * f64::abs(features.max[i] - features.min[i]);
        }
        let s = f64::sqrt(sum) * self.m_normal_vector;
        let ksi3 = &[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
        let ksi4 = &[0.0, 0.0, 0.0];
        let mut un = Vector::new(ndim);
        let mut x = Vector::new(ndim);
        if ndim == 2 {
            for (edge_key, edge) in &features.edges {
                let ncell = features.all_2d_edges.get(edge_key).unwrap().len();
                if ncell == 1 {
                    // only boundary edges
                    let mut pad = Scratchpad::new(ndim, edge.kind)?;
                    features.mesh.set_pad(&mut pad, &edge.points);
                    pad.calc_normal_vector(&mut un, ksi4)?;
                    pad.calc_coords(&mut x, ksi4)?;
                    self.canvas_normals_2d
                        .draw_arrow(x[0], x[1], x[0] + s * un[0], x[1] + s * un[1]);
                }
            }
            self.plot.add(&self.canvas_normals_2d);
        } else {
            for (face_key, face) in &features.faces {
                let ncell = features.all_faces.get(face_key).unwrap().len();
                if ncell == 1 {
                    // boundary faces only
                    let mut pad = Scratchpad::new(ndim, face.kind)?;
                    features.mesh.set_pad(&mut pad, &face.points);
                    if face.kind.class() == GeoClass::Tri {
                        pad.calc_normal_vector(&mut un, ksi3)?;
                        pad.calc_coords(&mut x, ksi3)?;
                    } else {
                        pad.calc_normal_vector(&mut un, ksi4)?;
                        pad.calc_coords(&mut x, ksi4)?;
                    }
                    self.canvas_normals_3d.polyline_3d_begin();
                    self.canvas_normals_3d.polyline_3d_add(x[0], x[1], x[2]);
                    self.canvas_normals_3d
                        .polyline_3d_add(x[0] + s * un[0], x[1] + s * un[1], x[2] + s * un[2]);
                    self.canvas_normals_3d.polyline_3d_end();
                }
            }
            for cell_id in &features.shells {
                let cell = &features.mesh.cells[*cell_id];
                let mut pad = Scratchpad::new(ndim, cell.kind)?;
                features.mesh.set_pad(&mut pad, &cell.points);
                if cell.kind.class() == GeoClass::Tri {
                    pad.calc_normal_vector(&mut un, ksi3)?;
                    pad.calc_coords(&mut x, ksi3)?;
                } else {
                    pad.calc_normal_vector(&mut un, ksi4)?;
                    pad.calc_coords(&mut x, ksi4)?;
                }
                self.canvas_normals_3d.polyline_3d_begin();
                self.canvas_normals_3d.polyline_3d_add(x[0], x[1], x[2]);
                self.canvas_normals_3d
                    .polyline_3d_add(x[0] + s * un[0], x[1] + s * un[1], x[2] + s * un[2]);
                self.canvas_normals_3d.polyline_3d_end();
            }
            self.plot.add(&self.canvas_normals_3d);
        }
        Ok(())
    }

    /// Draw boundary edges (3D)
    pub fn boundary_edges_3d(&mut self, features: &Features) {
        if features.mesh.ndim == 3 {
            let edge_keys = features.get_boundary_edges();
            if edge_keys.len() > 0 {
                for (i, j) in &edge_keys {
                    let a = &features.mesh.points[*i].coords;
                    let b = &features.mesh.points[*j].coords;
                    self.canvas_boundary_edges_3d.polyline_3d_begin();
                    self.canvas_boundary_edges_3d.polyline_3d_add(a[0], a[1], a[2]);
                    self.canvas_boundary_edges_3d.polyline_3d_add(b[0], b[1], b[2]);
                    self.canvas_boundary_edges_3d.polyline_3d_end();
                }
                self.plot.add(&self.canvas_boundary_edges_3d);
            }
        }
    }

    /// Draw boundary faces (3D)
    pub fn boundary_faces(&mut self, features: &Features) {
        if features.mesh.ndim == 3 {
            let res = features.triangulate_3d_boundary();
            if res.triangles.len() > 0 {
                self.canvas_boundary_faces
                    .draw_triangles_3d(&res.xx, &res.yy, &res.zz, &res.triangles);
                self.plot.add(&self.canvas_boundary_faces);
            }
        }
    }

    /// Configures the plot range
    fn config_plot_range(&mut self, mesh: &Mesh) {
        let (xmin, xmax) = mesh.get_limits();
        let dx = xmax[0] - xmin[0];
        let dy = xmax[1] - xmin[1];
        if dx > 0.0 && dy > 0.0 {
            let gx = self.m_range * dx;
            let gy = self.m_range * dy;
            self.plot
                .set_range(xmin[0] - gx, xmax[0] + gx, xmin[1] - gy, xmax[1] + gy);
        }
    }

    /// Draws the mesh
    ///
    /// # Input
    ///
    /// * `mesh` -- the mesh to be drawn
    /// * `filepath` -- may be a String, &str, or Path
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::mesh::{Draw, Samples};
    /// use gemlab::StrError;
    /// use plotpy::Canvas;
    ///
    /// const SAVE_FIGURE: bool = false;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     if SAVE_FIGURE {
    ///         // use sample mesh
    ///         let mesh = Samples::ring_eight_qua8_rad1_thick1();
    ///
    ///         // draw circles
    ///         let mut circle_in = Canvas::new();
    ///         let mut circle_out = Canvas::new();
    ///         circle_in
    ///             .set_face_color("None")
    ///             .set_edge_color("#ff000080")
    ///             .set_line_width(7.0)
    ///             .draw_circle(0.0, 0.0, 1.0);
    ///         circle_out
    ///             .set_face_color("None")
    ///             .set_edge_color("#0000ff80")
    ///             .set_line_width(7.0)
    ///             .draw_circle(0.0, 0.0, 2.0);
    ///
    ///         // draw mesh
    ///         let mut draw = Draw::new();
    ///         draw.extra(|plot, before| {
    ///             if !before {
    ///                 plot.add(&circle_in);
    ///                 plot.add(&circle_out);
    ///             }
    ///         })
    ///         .show_cell_ids(true)
    ///         .all(&mesh, "/tmp/gemlab/doc_example_mesh_draw.svg")?;
    ///     }
    ///     Ok(())
    /// }
    /// ```
    ///
    /// ![doc_example_mesh_draw](https://raw.githubusercontent.com/cpmech/gemlab/main/data/figures/doc_example_mesh_draw.svg)
    pub fn all<P>(&mut self, mesh: &Mesh, filepath: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        if let Some(callback) = self.extra.as_ref() {
            callback(&mut self.plot, true);
        }
        if self.show_cells {
            self.cells(mesh)?;
        }
        if self.show_cell_ids {
            self.cell_ids(mesh)?;
        }
        if self.show_point_dots {
            self.point_dots(mesh);
        }
        if self.show_point_ids {
            self.point_ids(mesh);
        }
        if self.show_edge_markers
            || self.show_face_markers
            || self.show_normal_vectors
            || self.show_boundary_edges_3d
            || self.show_boundary_faces
        {
            let features = Features::new(mesh, false);
            if self.show_edge_markers {
                self.edge_markers(&features)?;
            }
            if self.show_face_markers {
                self.face_markers(&features)?;
            }
            if self.show_normal_vectors {
                self.normal_vectors(&features)?;
            }
            if self.show_boundary_faces {
                self.boundary_faces(&features);
            }
            if self.show_boundary_edges_3d {
                self.boundary_edges_3d(&features);
            }
        }
        if mesh.ndim == 2 {
            self.plot.grid_and_labels("x", "y");
        }
        if mesh.ndim == 2 {
            if let Some((xmin, xmax, ymin, ymax)) = self.range_2d {
                self.plot.set_range(xmin, xmax, ymin, ymax);
            } else {
                self.config_plot_range(mesh);
            }
        } else {
            if let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = self.range_3d {
                self.plot.set_range_3d(xmin, xmax, ymin, ymax, zmin, zmax);
            } else {
                self.config_plot_range(mesh);
            }
        }
        if !self.unequal_exes {
            self.plot.set_equal_axes(true);
        }
        if let Some((width, height)) = self.size {
            self.plot.set_figure_size_points(width, height);
        }
        if let Some(callback) = self.extra.as_ref() {
            callback(&mut self.plot, false);
        }
        if let Some(((xmin, xmax, ymin, ymax), (u0, v0, w, h))) = self.zoom_2d {
            let mut inset = InsetAxes::new();
            inset
                .set_indicator_line_color(&self.zoom_indicator_config.0)
                .set_indicator_alpha(self.zoom_indicator_config.1)
                .set_indicator_line_width(self.zoom_indicator_config.2);
            inset.add(&self.canvas_cells);
            if self.show_cell_ids {
                inset.add(&self.canvas_cell_ids);
            }
            if self.show_point_dots {
                inset.add(&self.canvas_point_dots);
            }
            if self.show_point_ids {
                inset.add(&self.canvas_point_ids);
            }
            inset.set_range(xmin, xmax, ymin, ymax).draw(u0, v0, w, h);
            if let Some(callback) = self.zoom_extra.as_ref() {
                callback(&mut inset);
            }
            self.plot.add(&inset);
        }
        if mesh.ndim == 3 {
            self.plot.set_camera(self.camera_elevation, self.camera_azimuth);
        }
        if self.hide_3d_grid && mesh.ndim == 3 {
            self.plot.set_hide_3d_grid(true);
        }
        if self.hide_axes {
            self.plot.set_hide_axes(true);
        }
        if self.show_glyph_3d && mesh.ndim == 3 {
            self.canvas_glyph_3d.set_glyph_size(self.glyph_3d_size);
            self.canvas_glyph_3d.draw_glyph_3d(
                self.glyph_3d_origin[0],
                self.glyph_3d_origin[1],
                self.glyph_3d_origin[2],
            );
            self.plot.add(&self.canvas_glyph_3d);
        }
        if self.view {
            self.plot.show(filepath)
        } else {
            self.plot.save(filepath)
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Draw;
    use crate::mesh::{Features, Mesh, Samples};
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
        let mut draw = Draw::new();
        draw.cells(&mesh).unwrap();
        draw.point_dots(&mesh);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut draw.plot, &mut labels, &mesh, 0.0, 0.05, 0.0);
            caption.draw(0.6, -0.05, "Lin2");
            caption.draw(1.4 + 0.6, -0.05, "Lin3");
            caption.draw(2.8 + 0.6, -0.05, "Lin4");
            caption.draw(4.2 + 0.6, -0.05, "Lin5");
            draw.plot.add(&caption);
            draw.plot
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
        let mut draw = Draw::new();
        draw.cells(&mesh).unwrap();
        draw.point_dots(&mesh);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut draw.plot, &mut labels, &mesh, 0.0, 0.05, 0.0);
            caption.draw_3d(1.2, 1.33, 1.33, "Lin2");
            caption.draw_3d(1.4 + 1.2, 1.33, 1.33, "Lin3");
            caption.draw_3d(2.8 + 1.2, 1.33, 1.33, "Lin4");
            caption.draw_3d(4.2 + 1.2, 1.33, 1.33, "Lin5");
            draw.plot.add(&caption);
            draw.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range_3d(0.0, 5.3, 0.0, 1.2, 0.0, 1.2)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_1_lin_3d.svg")
                .unwrap();
        }

        // tri cells ---------------------------------------------------------------------------
        let mesh = Samples::tri_cells();
        let mut draw = Draw::new();
        draw.cells(&mesh).unwrap();
        draw.point_dots(&mesh);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut draw.plot, &mut labels, &mesh, 0.0, 0.02, 0.0);
            caption.draw(0.5, -0.1, "Tri3");
            caption.draw(1.7, -0.1, "Tri6");
            caption.draw(0.5, 1.1, "Tri10");
            caption.draw(1.7, 1.1, "Tri15");
            draw.plot.add(&caption);
            draw.plot
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
        let mut draw = Draw::new();
        draw.cells(&mesh).unwrap();
        draw.point_dots(&mesh);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut draw.plot, &mut labels, &mesh, -0.02, 0.02, 0.0);
            caption.draw(0.4, -0.06, "Qua4");
            caption.draw(1.5, -0.06, "Qua8");
            caption.draw(2.6, -0.06, "Qua9");
            caption.draw(0.4, 1.1, "Qua12");
            caption.draw(1.5, 1.1, "Qua16");
            caption.draw(2.6, 1.1, "Qua17");
            draw.plot.add(&caption);
            draw.plot
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
        let mut draw = Draw::new();
        draw.cells(&mesh).unwrap();
        draw.point_dots(&mesh);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut draw.plot, &mut labels, &mesh, 0.03, 0.0, 0.03);
            caption.draw_3d(0.5, -0.15, 0.0, "Tet4");
            caption.draw_3d(1.7, -0.15, 0.0, "Tet10");
            caption.draw_3d(1.1, 1.05, 0.8, "Tet20");
            draw.plot.add(&caption);
            draw.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range_3d(0.0, 2.2, 0.0, 2.2, 0.0, 2.2)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_4_tet.svg")
                .unwrap();
        }

        // hex cells ---------------------------------------------------------------------------
        let mesh = Samples::hex_cells();
        let mut draw = Draw::new();
        draw.cells(&mesh).unwrap();
        draw.point_dots(&mesh);

        if SAVE_FIGURE {
            let (mut labels, mut caption) = labels_and_caption();
            draw_cell_local_ids(&mut draw.plot, &mut labels, &mesh, 0.03, 0.0, 0.03);
            caption.draw_3d(0.5, -0.25, 0.0, "Hex8");
            caption.draw_3d(2.3, -0.25, 0.0, "Hex20");
            caption.draw_3d(0.3, 1.1, 2.0, "Hex32");
            draw.plot.add(&caption);
            draw.plot
                .set_figure_size_points(600.0, 600.0)
                .set_equal_axes(true)
                .set_range_3d(0.0, 3.0, 0.0, 3.0, 0.0, 2.5)
                .save("/tmp/gemlab/test_draw_cells_and_points_work_5_hex.svg")
                .unwrap();
        }

        // ring --------------------------------------------------------------------------------
        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        let mut draw = Draw::new();
        draw.cells(&mesh).unwrap();
        draw.point_dots(&mesh);

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
            draw.plot.add(&circle_in);
            draw.plot.add(&circle_mi);
            draw.plot.add(&circle_ou);
            draw.plot
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
            let mut draw = Draw::new();
            draw.show_cell_ids(true)
                .show_point_ids(true)
                .show_point_dots(true)
                .set_range_2d(-0.5, 6.0, -0.5, 6.0)
                .zoom_2d(-0.05, 1.55, -0.05, 1.55, 0.6, 0.6, 0.3, 0.3)
                .zoom_extra(|inset| {
                    let mut text = Text::new();
                    text.draw(0.3, 1.0, "HELLO");
                    inset.add(&text);
                })
                .all(&mesh, "/tmp/gemlab/test_draw_works_qua12.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_works_qua16() {
        if SAVE_FIGURE {
            let mesh = Samples::block_2d_four_qua16();
            let mut draw = Draw::new();
            draw.show_cell_ids(true)
                .show_point_ids(true)
                .show_point_dots(true)
                .all(&mesh, "/tmp/gemlab/test_draw_works_qua16.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_works_qua17() {
        if SAVE_FIGURE {
            let mesh = Samples::block_2d_four_qua17();
            let mut draw = Draw::new();
            draw.show_cell_ids(true)
                .show_point_ids(true)
                .show_point_dots(true)
                .all(&mesh, "/tmp/gemlab/test_draw_works_qua17.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_works_mixed_2d() {
        if SAVE_FIGURE {
            let mesh = Samples::mixed_shapes_2d();
            let mut draw = Draw::new();
            draw.show_cell_ids(true)
                .show_point_ids(true)
                .show_point_dots(true)
                .all(&mesh, "/tmp/gemlab/test_draw_works_mixed_2d.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_works_hex8() {
        if SAVE_FIGURE {
            let mesh = Samples::two_hex8();
            let mut draw = Draw::new();
            draw.show_cell_ids(true)
                .show_point_ids(true)
                .show_point_dots(true)
                .set_size(600.0, 600.0);
            draw.get_canvas_point_ids()
                .set_align_horizontal("left")
                .set_align_vertical("bottom")
                .set_color("black")
                .set_fontsize(10.0)
                .set_bbox_facecolor("gold")
                .set_bbox_alpha(0.5);
            draw.all(&mesh, "/tmp/gemlab/test_draw_works_hex8.svg").unwrap();
        }
    }

    #[test]
    fn draw_works_hex20() {
        if SAVE_FIGURE {
            let mesh = Samples::block_3d_eight_hex20();
            let mut draw = Draw::new();
            draw.show_cell_ids(true)
                .show_point_ids(true)
                .show_point_dots(true)
                .set_size(600.0, 600.0);
            draw.get_canvas_point_ids()
                .set_align_horizontal("left")
                .set_align_vertical("bottom")
                .set_color("black")
                .set_fontsize(10.0)
                .set_bbox_facecolor("gold")
                .set_bbox_alpha(0.5);
            draw.all(&mesh, "/tmp/gemlab/test_draw_works_hex20.svg").unwrap();
        }
    }

    #[test]
    fn draw_works_mixed_3d() {
        if SAVE_FIGURE {
            let mesh = Samples::mixed_shapes_3d();
            let mut draw = Draw::new();
            draw.show_cell_ids(true)
                .show_point_ids(true)
                .show_point_dots(true)
                .show_normal_vectors(true)
                .set_m_normal_vector(0.08)
                .set_m_range(0.0)
                .set_glyph_3d(1.8, 0.2, 0.0, 0.25)
                .set_hide_axes(true)
                .set_hide_3d_grid(true)
                .set_range_3d(0.0, 1.0, -0.5, 2.0, 0.0, 1.0)
                .set_size(800.0, 800.0)
                .set_camera(30.0, 30.0)
                .set_view_flag(false);
            draw.get_canvas_point_ids()
                .set_align_horizontal("left")
                .set_align_vertical("bottom")
                .set_color("black")
                .set_fontsize(10.0)
                .set_bbox_facecolor("gold")
                .set_bbox_alpha(0.5);
            draw.all(&mesh, "/tmp/gemlab/test_draw_works_mixed_3d.svg").unwrap();
        }
    }

    #[test]
    fn test_figure_setters() {
        let mut draw = Draw::new();

        // Test show_cell_ids
        assert!(!draw.show_cell_ids);
        draw.show_cell_ids(true);
        assert!(draw.show_cell_ids);

        // Test show_cell_att
        assert!(draw.show_cell_marker);
        draw.show_cell_marker(false);
        assert!(!draw.show_cell_marker);

        // Test show_point_ids
        assert!(!draw.show_point_ids);
        draw.show_point_ids(true);
        assert!(draw.show_point_ids);

        // Test show_point_marker
        assert!(!draw.show_point_marker);
        draw.show_point_marker(true);
        assert!(draw.show_point_marker);

        // Test show_point_dots
        assert!(!draw.show_point_dots);
        draw.show_point_dots(true);
        assert!(draw.show_point_dots);

        // Test unequal_axes
        assert!(!draw.unequal_exes);
        draw.set_unequal_exes(true);
        assert!(draw.unequal_exes);

        // Test range_2d
        assert!(draw.range_2d.is_none());
        draw.set_range_2d(-1.0, 1.0, -2.0, 2.0);
        assert_eq!(draw.range_2d, Some((-1.0, 1.0, -2.0, 2.0)));

        // Test range_3d
        assert!(draw.range_3d.is_none());
        draw.set_range_3d(-1.0, 1.0, -2.0, 2.0, -3.0, 3.0);
        assert_eq!(draw.range_3d, Some((-1.0, 1.0, -2.0, 2.0, -3.0, 3.0)));

        // Test size
        assert!(draw.size.is_none());
        draw.set_size(800.0, 600.0);
        assert_eq!(draw.size, Some((800.0, 600.0)));

        // Test extra
        assert!(draw.extra.is_none());
        draw.extra(|_, _| {});
        assert!(draw.extra.is_some());

        // Test zoom_2d
        assert!(draw.zoom_2d.is_none());
        draw.zoom_2d(-0.5, 0.5, -1.0, 1.0, 0.1, 0.2, 0.3, 0.4);
        assert_eq!(draw.zoom_2d, Some(((-0.5, 0.5, -1.0, 1.0), (0.1, 0.2, 0.3, 0.4))));

        // Test disable_zoom_2d
        draw.disable_zoom_2d();
        assert!(draw.zoom_2d.is_none());

        // Test zoom_indicator
        draw.zoom_indicator("#ff0000", 0.8, 2.5);
        assert_eq!(draw.zoom_indicator_config, ("#ff0000".to_string(), 0.8, 2.5));

        // Test zoom_extra
        assert!(draw.zoom_extra.is_none());
        draw.zoom_extra(|_| {});
        assert!(draw.zoom_extra.is_some());

        // Test method chaining
        let mut draw = Draw::new();
        draw.show_cell_ids(true)
            .show_point_ids(true)
            .show_point_dots(true)
            .set_size(800.0, 600.0)
            .set_range_2d(-1.0, 1.0, -1.0, 1.0);
        assert!(draw.show_cell_ids);
        assert!(draw.show_point_ids);
        assert!(draw.show_point_dots);
        assert_eq!(draw.size, Some((800.0, 600.0)));
        assert_eq!(draw.range_2d, Some((-1.0, 1.0, -1.0, 1.0)));
    }

    #[test]
    fn draw_edge_markers_works() {
        if SAVE_FIGURE {
            let mesh = Samples::two_qua4();
            let features = Features::new(&mesh, false);
            let mut draw = Draw::new();
            draw.cells(&mesh).unwrap();
            draw.point_dots(&mesh);
            draw.edge_markers(&features).unwrap();
            draw.plot
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_edge_markers_works_1.svg")
                .unwrap();

            let mesh = Samples::ring_eight_qua8_rad1_thick1();
            let features = Features::new(&mesh, false);
            let mut draw = Draw::new();
            draw.cells(&mesh).unwrap();
            draw.point_dots(&mesh);
            draw.edge_markers(&features).unwrap();
            draw.plot
                .set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .save("/tmp/gemlab/test_draw_edge_markers_works_2.svg")
                .unwrap();

            let mesh = Samples::two_hex8();
            let features = Features::new(&mesh, true);
            let mut draw = Draw::new();
            draw.cells(&mesh).unwrap();
            draw.point_dots(&mesh);
            draw.edge_markers(&features).unwrap();
            draw.plot
                .set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .save("/tmp/gemlab/test_draw_edge_markers_works_3.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_normals_works() {
        if SAVE_FIGURE {
            let mesh = Samples::two_qua4();
            let features = Features::new(&mesh, false);
            let mut draw = Draw::new();
            draw.cells(&mesh).unwrap();
            draw.point_dots(&mesh);
            draw.normal_vectors(&features).unwrap();
            draw.plot
                .set_equal_axes(true)
                .save("/tmp/gemlab/test_draw_normals_works_1.svg")
                .unwrap();

            let mesh = Samples::ring_eight_qua8_rad1_thick1();
            let features = Features::new(&mesh, false);
            let mut draw = Draw::new();
            draw.cells(&mesh).unwrap();
            draw.point_dots(&mesh);
            draw.normal_vectors(&features).unwrap();
            draw.plot
                .set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .save("/tmp/gemlab/test_draw_normals_works_2.svg")
                .unwrap();

            let mesh = Samples::two_hex8();
            let features = Features::new(&mesh, true);
            let mut draw = Draw::new();
            draw.cells(&mesh).unwrap();
            draw.point_dots(&mesh);
            draw.normal_vectors(&features).unwrap();
            draw.plot
                .set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                // .show("/tmp/gemlab/test_draw_normals_works_3.svg")
                .save("/tmp/gemlab/test_draw_normals_works_3.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_face_markers_works() {
        if SAVE_FIGURE {
            let mesh = Samples::two_hex8();
            let features = Features::new(&mesh, true);
            let mut draw = Draw::new();
            draw.cells(&mesh).unwrap();
            draw.point_dots(&mesh);
            draw.face_markers(&features).unwrap();
            draw.plot
                .set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                // .show("/tmp/gemlab/test_draw_face_markers_works.svg")
                .save("/tmp/gemlab/test_draw_face_markers_works.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_boundary_edges_works() {
        if SAVE_FIGURE {
            // let mesh = Samples::two_hex8();
            let mesh = Samples::block_3d_eight_hex20();
            let features = Features::new(&mesh, true);
            let mut draw = Draw::new();
            draw.boundary_edges_3d(&features);
            draw.point_dots(&mesh);
            draw.plot
                .set_equal_axes(true)
                .set_figure_size_points(800.0, 800.0)
                // .show("/tmp/gemlab/test_draw_boundary_edges_works.svg")
                .save("/tmp/gemlab/test_draw_boundary_edges_works.svg")
                .unwrap();
        }
    }

    #[test]
    fn draw_boundary_faces_works() {
        if SAVE_FIGURE {
            // let mesh = Samples::two_hex8();
            let mesh = Samples::block_3d_eight_hex20();
            let features = Features::new(&mesh, true);
            let mut draw = Draw::new();
            // draw.cells(&mesh).unwrap();
            // draw.point_dots(&mesh);
            draw.boundary_faces(&features);
            draw.boundary_edges_3d(&features);
            draw.plot
                .set_equal_axes(true)
                .set_figure_size_points(800.0, 800.0)
                // .show("/tmp/gemlab/test_draw_boundary_faces_works.svg")
                .save("/tmp/gemlab/test_draw_boundary_faces_works.svg")
                .unwrap();
        }
    }
}
