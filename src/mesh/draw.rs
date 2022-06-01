use super::{Features, Mesh, Region};
use crate::mesh::allocate_state;
use crate::shapes::{geo_class_and_kind, GeoClass, GeoKind, Shape, StateOfShape};
use crate::util::ONE_BY_3;
use crate::StrError;
use plotpy::{Canvas, Curve, Plot, PolyCode, Text};
use russell_lab::Vector;
use std::collections::HashMap;

/// Returns the first point on GeoClass:Lin; Bezier t=0(ξ=-1)
macro_rules! pa {
    ($state:expr,$i:expr) => {
        $state.coords_transp[$i][0]
    };
}

/// Returns the last point on GeoClass:Lin; Bezier t=1(ξ=1)
macro_rules! pb {
    ($state:expr,$i:expr) => {
        $state.coords_transp[$i][1]
    };
}

/// Implements functions to draw edges and faces
pub struct Draw {
    /// Canvas to draw edges
    pub canvas_edges: Canvas,

    /// Canvas to draw nodes (markers)
    pub canvas_nodes: Curve,

    /// Canvas to draw labels of nodes
    pub canvas_nodes_labels: Text,

    /// Canvas to draw labels of cells
    pub canvas_cells_labels: Text,

    /// With drawing of nodes (markers); not available if with_labels = true
    pub with_nodes: bool,

    /// With drawing of labels of nodes; will turn off with_nodes
    pub with_labels_nodes: bool,

    /// Will set the coordinates range (min/max)
    pub with_set_range: bool,
}

impl Draw {
    /// Allocates a new instance
    pub fn new() -> Self {
        let mut canvas_edges = Canvas::new();
        let mut canvas_nodes = Curve::new();
        let mut canvas_nodes_labels = Text::new();
        let mut canvas_cells_labels = Text::new();
        canvas_edges
            .set_stop_clip(true)
            .set_face_color("None")
            .set_line_width(2.0)
            .set_edge_color("#2440cd");
        canvas_nodes
            .set_stop_clip(true)
            .set_marker_color("black")
            .set_marker_line_color("white")
            .set_marker_style("o")
            .set_line_style("None");
        canvas_nodes_labels
            .set_color("red")
            .set_fontsize(8.0)
            .set_align_horizontal("center")
            .set_align_vertical("center")
            .set_bbox(true)
            .set_bbox_facecolor("white")
            .set_bbox_edgecolor("None");
        canvas_cells_labels
            .set_color("#22971f")
            .set_fontsize(9.0)
            .set_align_horizontal("center")
            .set_align_vertical("center")
            .set_bbox(true)
            .set_bbox_facecolor("white")
            .set_bbox_edgecolor("#b7b7b7");
        Draw {
            canvas_edges,
            canvas_nodes,
            canvas_nodes_labels,
            canvas_cells_labels,
            with_nodes: true,
            with_labels_nodes: false,
            with_set_range: true,
        }
    }

    /// Draws ids and attributes of cells
    pub fn cell_ids(&mut self, plot: &mut Plot, region: &Region) -> Result<(), StrError> {
        // auxiliary
        let mesh = &region.mesh;
        let space_ndim = mesh.space_ndim;
        let mut states_memo: HashMap<GeoKind, StateOfShape> = HashMap::new();
        let ksi_hex = &[0.0, 0.0, 0.0];
        let ksi_tet = &[ONE_BY_3, ONE_BY_3, ONE_BY_3];
        let mut x = Vector::new(space_ndim);

        // loop over all cells
        for cell_id in 0..mesh.cells.len() {
            // shape and state
            let (cell, shape) = (&mesh.cells[cell_id], &region.shapes[cell_id]);
            let (_, kind) = geo_class_and_kind(cell.geo_ndim, cell.points.len())?;
            let mut state = states_memo
                .entry(kind)
                .or_insert(allocate_state(&mesh, shape.geo_ndim, &cell.points)?);

            // set coordinates
            for m in 0..cell.points.len() {
                for i in 0..space_ndim {
                    state.coords_transp[i][m] = mesh.points[cell.points[m]].coords[i];
                }
            }

            // calc center of cell by using center of the reference space
            let ksi = match shape.class {
                GeoClass::Lin => ksi_hex,
                GeoClass::Tri => ksi_tet,
                GeoClass::Qua => ksi_hex,
                GeoClass::Tet => ksi_tet,
                GeoClass::Hex => ksi_hex,
            };
            shape.calc_coords(&mut x, &mut state, ksi)?;

            if cell.id == 4 {
                println!("{}: {:?} => {:?}", cell.id, ksi, x.as_data());
            }

            // add label
            if space_ndim == 2 {
                self.canvas_cells_labels
                    .draw(x[0], x[1], format!("{}({})", cell.id, cell.attribute_id).as_str());
            } else {
                self.canvas_cells_labels.draw_3d(
                    x[0],
                    x[1],
                    x[2],
                    format!("{}({})", cell.id, cell.attribute_id).as_str(),
                );
            }
        }

        // add to plot
        plot.add(&self.canvas_cells_labels);
        Ok(())
    }

    /// Draws edges
    pub fn edges(&mut self, plot: &mut Plot, region: &Region) -> Result<(), StrError> {
        if region.mesh.space_ndim == 2 {
            self.edges_2d(plot, &region.mesh, &region.features)
        } else {
            self.edges_3d(plot, &region.mesh, &region.features)
        }
    }

    /// Draws 2D edges
    fn edges_2d(&mut self, plot: &mut Plot, mesh: &Mesh, features: &Features) -> Result<(), StrError> {
        // space dimension
        let space_ndim = mesh.space_ndim;
        assert_eq!(space_ndim, 2);

        // middle points (p) on Bezier curve and control points (q)
        let mut pc = Vector::new(space_ndim); // first middle point on Bezier curve with t=⅓(ξ=-⅓) or t=½(ξ=0)
        let mut pd = Vector::new(space_ndim); // second middle point on Bezier curve with t=⅔(ξ=+⅓)
        let mut qc = Vector::new(space_ndim); // control point corresponding to pc
        let mut qd = Vector::new(space_ndim); // control point corresponding to pd

        // reference coordinate
        let mut ksi = vec![0.0; space_ndim];

        // memoization of shapes and states
        let mut shapes_memo: HashMap<GeoKind, Shape> = HashMap::new();
        let mut states_memo: HashMap<GeoKind, StateOfShape> = HashMap::new();

        // begin poly-curve
        self.canvas_edges.polycurve_begin();

        // loop over edges (GeoKind::Lin)
        const EDGE_GEO_NDIM: usize = 1;
        for (_, edge) in &features.edges {
            // number of nodes on edge and kind
            let nnode = edge.points.len();
            let (_, kind) = geo_class_and_kind(EDGE_GEO_NDIM, nnode)?;

            // retrieve shape or allocate new
            let shape = shapes_memo
                .entry(kind)
                .or_insert(Shape::new(space_ndim, EDGE_GEO_NDIM, nnode)?);

            // retrieve state or allocate new
            let mut state = states_memo
                .entry(kind)
                .or_insert(allocate_state(&mesh, EDGE_GEO_NDIM, &edge.points)?);

            // set coordinates
            for m in 0..nnode {
                for i in 0..space_ndim {
                    state.coords_transp[i][m] = mesh.points[edge.points[m]].coords[i];
                }
            }

            // add poly-curve points (or control points for quadratic/cubic Bezier)
            match nnode {
                2 => {
                    self.canvas_edges
                        .polycurve_add(pa!(state, 0), pa!(state, 1), PolyCode::MoveTo)
                        .polycurve_add(pb!(state, 0), pb!(state, 1), PolyCode::LineTo);
                }
                3 => {
                    ksi[0] = 0.0; // middle
                    shape.calc_coords(&mut pc, &mut state, &ksi)?;
                    qc[0] = (-pa!(state, 0) - pb!(state, 0) + 4.0 * pc[0]) / 2.0;
                    qc[1] = (-pa!(state, 1) - pb!(state, 1) + 4.0 * pc[1]) / 2.0;
                    self.canvas_edges
                        .polycurve_add(pa!(state, 0), pa!(state, 1), PolyCode::MoveTo)
                        .polycurve_add(/*   */ qc[0], /*   */ qc[1], PolyCode::Curve3)
                        .polycurve_add(pb!(state, 0), pb!(state, 1), PolyCode::Curve3);
                }
                _ => {
                    ksi[0] = -1.0 / 3.0; // => t=⅓(ξ=-⅓)
                    shape.calc_coords(&mut pc, &mut state, &ksi)?;
                    ksi[0] = 1.0 / 3.0; // => t=⅔(ξ=+⅓)
                    shape.calc_coords(&mut pd, &mut state, &ksi)?;
                    qc[0] = (-5.0 * pa!(state, 0) + 2.0 * pb!(state, 0) + 18.0 * pc[0] - 9.0 * pd[0]) / 6.0;
                    qc[1] = (-5.0 * pa!(state, 1) + 2.0 * pb!(state, 1) + 18.0 * pc[1] - 9.0 * pd[1]) / 6.0;
                    qd[0] = (2.0 * pa!(state, 0) - 5.0 * pb!(state, 0) - 9.0 * pc[0] + 18.0 * pd[0]) / 6.0;
                    qd[1] = (2.0 * pa!(state, 1) - 5.0 * pb!(state, 1) - 9.0 * pc[1] + 18.0 * pd[1]) / 6.0;
                    self.canvas_edges
                        .polycurve_add(pa!(state, 0), pa!(state, 1), PolyCode::MoveTo)
                        .polycurve_add(/*   */ qc[0], /*   */ qc[1], PolyCode::Curve4)
                        .polycurve_add(/*   */ qd[0], /*   */ qd[1], PolyCode::Curve4)
                        .polycurve_add(pb!(state, 0), pb!(state, 1), PolyCode::Curve4);
                }
            }
        }

        // end polycurve and add to plot
        self.canvas_edges.polycurve_end(false);
        plot.add(&self.canvas_edges);
        if self.with_set_range {
            plot.set_range(features.min[0], features.max[0], features.min[1], features.max[1]);
        }

        // add nodes
        if self.with_nodes && !self.with_labels_nodes {
            self.canvas_nodes.points_begin();
            for p in &features.points {
                self.canvas_nodes
                    .points_add(mesh.points[*p].coords[0], mesh.points[*p].coords[1]);
            }
            self.canvas_nodes.points_end();
            plot.add(&self.canvas_nodes);
        }

        // add nodes labels
        if self.with_labels_nodes {
            for p in &features.points {
                self.canvas_nodes_labels.draw(
                    mesh.points[*p].coords[0],
                    mesh.points[*p].coords[1],
                    format!("{}", *p).as_str(),
                );
            }
            plot.add(&self.canvas_nodes_labels);
        }
        Ok(())
    }

    /// Draws 3D edges
    fn edges_3d(&mut self, plot: &mut Plot, mesh: &Mesh, features: &Features) -> Result<(), StrError> {
        // space dimension
        let space_ndim = mesh.space_ndim;
        assert_eq!(space_ndim, 3);

        // middle points (p) on edge
        let mut pc = Vector::new(space_ndim);
        let mut pd = Vector::new(space_ndim);
        let mut pe = Vector::new(space_ndim);

        // reference coordinate
        let mut ksi = vec![0.0; space_ndim];

        // memoization of shapes and states
        let mut shapes_memo: HashMap<GeoKind, Shape> = HashMap::new();
        let mut states_memo: HashMap<GeoKind, StateOfShape> = HashMap::new();

        // begin poly-line

        // loop over edges (GeoKind::Lin)
        const EDGE_GEO_NDIM: usize = 1;
        for (_, edge) in &features.edges {
            // number of nodes on edge and kind
            let nnode = edge.points.len();
            let (_, kind) = geo_class_and_kind(EDGE_GEO_NDIM, nnode)?;

            // retrieve shape or allocate new
            let shape = shapes_memo
                .entry(kind)
                .or_insert(Shape::new(space_ndim, EDGE_GEO_NDIM, nnode)?);

            // retrieve state or allocate new
            let mut state = states_memo
                .entry(kind)
                .or_insert(allocate_state(&mesh, EDGE_GEO_NDIM, &edge.points)?);

            // set coordinates of GeoKind::Lin
            for m in 0..nnode {
                for i in 0..space_ndim {
                    state.coords_transp[i][m] = mesh.points[edge.points[m]].coords[i];
                }
            }

            // add poly-line points
            self.canvas_edges.polyline_3d_begin();
            match nnode {
                2 => {
                    self.canvas_edges
                        .polyline_3d_add(pa!(state, 0), pa!(state, 1), pa!(state, 2))
                        .polyline_3d_add(pb!(state, 0), pb!(state, 1), pb!(state, 2));
                }
                3 => {
                    ksi[0] = 0.0; // middle
                    shape.calc_coords(&mut pc, &mut state, &ksi)?;
                    self.canvas_edges
                        .polyline_3d_add(pa!(state, 0), pa!(state, 1), pa!(state, 2))
                        .polyline_3d_add(/*   */ pc[0], /*   */ pc[1], /*   */ pc[2])
                        .polyline_3d_add(pb!(state, 0), pb!(state, 1), pb!(state, 2));
                }
                4 => {
                    ksi[0] = -1.0 / 3.0; // middle-left
                    shape.calc_coords(&mut pc, &mut state, &ksi)?;
                    ksi[0] = 1.0 / 3.0; // middle-right
                    shape.calc_coords(&mut pd, &mut state, &ksi)?;
                    self.canvas_edges
                        .polyline_3d_add(pa!(state, 0), pa!(state, 1), pa!(state, 2))
                        .polyline_3d_add(/*   */ pc[0], /*   */ pc[1], /*   */ pc[2])
                        .polyline_3d_add(/*   */ pd[0], /*   */ pd[1], /*   */ pd[2])
                        .polyline_3d_add(pb!(state, 0), pb!(state, 1), pb!(state, 2));
                }
                _ => {
                    ksi[0] = -0.5; // middle-left
                    shape.calc_coords(&mut pc, &mut state, &ksi)?;
                    ksi[0] = 0.0; // middle
                    shape.calc_coords(&mut pd, &mut state, &ksi)?;
                    ksi[0] = 0.5; // middle-right
                    shape.calc_coords(&mut pe, &mut state, &ksi)?;
                    self.canvas_edges
                        .polyline_3d_add(pa!(state, 0), pa!(state, 1), pa!(state, 2))
                        .polyline_3d_add(/*   */ pc[0], /*   */ pc[1], /*   */ pc[2])
                        .polyline_3d_add(/*   */ pd[0], /*   */ pd[1], /*   */ pd[2])
                        .polyline_3d_add(/*   */ pe[0], /*   */ pe[1], /*   */ pe[2])
                        .polyline_3d_add(pb!(state, 0), pb!(state, 1), pb!(state, 2));
                }
            }
            self.canvas_edges.polyline_3d_end();
        }

        // add to plot
        plot.add(&self.canvas_edges);
        if self.with_set_range {
            plot.set_range_3d(
                features.min[0],
                features.max[0],
                features.min[1],
                features.max[1],
                features.min[2],
                features.max[2],
            );
        }

        // add nodes
        if self.with_nodes && !self.with_labels_nodes {
            self.canvas_nodes.points_3d_begin();
            for p in &features.points {
                self.canvas_nodes.points_3d_add(
                    mesh.points[*p].coords[0],
                    mesh.points[*p].coords[1],
                    mesh.points[*p].coords[2],
                );
            }
            self.canvas_nodes.points_3d_end();
            plot.add(&self.canvas_nodes);
        }

        // add nodes labels
        if self.with_labels_nodes {
            for p in &features.points {
                self.canvas_nodes_labels.draw_3d(
                    mesh.points[*p].coords[0],
                    mesh.points[*p].coords[1],
                    mesh.points[*p].coords[2],
                    format!("{}", *p).as_str(),
                );
            }
            plot.add(&self.canvas_nodes_labels);
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Draw;
    use crate::mesh::{Extract, Region, Samples};
    use crate::StrError;
    use plotpy::{Canvas, Plot};

    #[test]
    fn draw_works_2d() -> Result<(), StrError> {
        // draw reference circles
        let mut plot = Plot::new();
        let mut circle_in = Canvas::new();
        let mut circle_mi = Canvas::new();
        let mut circle_ou = Canvas::new();
        circle_in
            .set_face_color("None")
            .set_edge_color("magenta")
            .set_line_width(7.0)
            .draw_circle(0.0, 0.0, 1.0);
        circle_mi
            .set_face_color("None")
            .set_edge_color("cyan")
            .set_line_width(7.0)
            .draw_circle(0.0, 0.0, 1.5);
        circle_ou
            .set_face_color("None")
            .set_edge_color("orange")
            .set_line_width(7.0)
            .draw_circle(0.0, 0.0, 2.0);
        plot.add(&circle_in);
        plot.add(&circle_mi);
        plot.add(&circle_ou);

        // draw edges
        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        let region = Region::with(mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.with_labels_nodes = true;
        draw.edges(&mut plot, &region)?;

        // draw cell ids
        draw.cell_ids(&mut plot, &region)?;

        // save figure
        plot.set_figure_size_points(800.0, 800.0)
            .set_equal_axes(true)
            .save("/tmp/gemlab/draw_works_2d.svg")?;
        Ok(())
    }

    #[test]
    fn draw_works_3d() -> Result<(), StrError> {
        // draw edges
        let mut plot = Plot::new();
        let mesh = Samples::two_cubes_vertical();
        let region = Region::with(mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.with_labels_nodes = true;
        draw.edges(&mut plot, &region)?;

        // draw cell ids
        draw.cell_ids(&mut plot, &region)?;

        // save figure
        plot.set_figure_size_points(800.0, 800.0)
            .set_equal_axes(true)
            .save("/tmp/gemlab/draw_works_2d.svg")?;
        Ok(())
    }

    #[test]
    fn draw_works_3d_wall() -> Result<(), StrError> {
        // draw edges
        let mut plot = Plot::new();
        let mesh = Samples::four_cubes_wall();
        let region = Region::with(mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.with_labels_nodes = true;
        draw.edges(&mut plot, &region)?;

        // draw cell ids
        draw.cell_ids(&mut plot, &region)?;

        // save figure
        plot.set_figure_size_points(800.0, 800.0)
            .set_equal_axes(true)
            .save("/tmp/gemlab/draw_works_3d_wall.svg")?;
        Ok(())
    }

    #[test]
    fn draw_works_3d_mixed() -> Result<(), StrError> {
        // draw edges
        let mut plot = Plot::new();
        let mesh = Samples::mixed_shapes_3d();
        let region = Region::with(mesh, Extract::All)?;
        let mut draw = Draw::new();
        draw.with_labels_nodes = true;
        draw.edges(&mut plot, &region)?;

        // draw cell ids
        draw.cell_ids(&mut plot, &region)?;

        // save figure
        plot.set_figure_size_points(800.0, 800.0)
            .set_equal_axes(true)
            // .save_and_show("/tmp/gemlab/draw_works_3d_mixed.svg")?;
            .save("/tmp/gemlab/draw_works_3d_mixed.svg")?;
        Ok(())
    }
}
