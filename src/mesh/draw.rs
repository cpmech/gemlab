use super::{Features, Mesh, Region};
use crate::mesh::allocate_state;
use crate::shapes::{geo_class_and_kind, GeoKind, Shape, StateOfShape};
use crate::StrError;
use plotpy::{Canvas, Curve, Plot, PolyCode, Text};
use russell_lab::Vector;
use std::collections::HashMap;

// Returns the first point on GeoClass:Lin; Bezier t=0(ξ=-1)
macro_rules! pa {
    ($state:expr,$i:expr) => {
        $state.coords_transp[$i][0]
    };
}

// Returns the last point on GeoClass:Lin; Bezier t=1(ξ=1)
macro_rules! pb {
    ($state:expr,$i:expr) => {
        $state.coords_transp[$i][1]
    };
}

pub struct Draw {
    pub canvas_edges: Canvas,
    pub canvas_nodes: Curve,
    pub canvas_nodes_labels: Text,
    pub nodes: bool,
    pub nodes_labels: bool,
    pub set_range: bool,
}

impl Draw {
    pub fn new() -> Self {
        let mut canvas_edges = Canvas::new();
        let mut canvas_nodes = Curve::new();
        let mut canvas_nodes_labels = Text::new();
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
        canvas_nodes_labels.set_color("red").set_fontsize(6.0);
        Draw {
            canvas_edges,
            canvas_nodes,
            canvas_nodes_labels,
            nodes: true,
            nodes_labels: false,
            set_range: true,
        }
    }

    pub fn edges(&mut self, plot: &mut Plot, region: &Region) -> Result<(), StrError> {
        if region.mesh.space_ndim == 2 {
            self.edges_2d(plot, &region.mesh, &region.features)
        } else {
            self.edges_3d(plot, &region.mesh, &region.features)
        }
    }

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

            // set coordinates of GeoKind::Lin
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
        if self.set_range {
            plot.set_range(features.min[0], features.max[0], features.min[1], features.max[1]);
        }

        // add nodes
        if self.nodes {
            // plt.nodes.begin();
            // for p in &region.features.points {
            //     plt.nodes
            //         .add(region.mesh.points[*p].coords[0], region.mesh.points[*p].coords[1])
            // }
            // plt.nodes.end();
            // plot.add(&plt.nodes);
        }
        Ok(())
    }

    fn edges_3d(&mut self, _plot: &mut Plot, _mesh: &Mesh, _features: &Features) -> Result<(), StrError> {
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
    fn draw_edges_works() -> Result<(), StrError> {
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
        draw.edges(&mut plot, &region)?;

        // save figure
        plot.set_figure_size_points(800.0, 800.0)
            .set_equal_axes(true)
            .save("/tmp/gemlab/draw_edges_works.svg")?;
        Ok(())
    }
}
