use super::calc_coords;
use crate::shapes::{GeoClass, Scratchpad};
use crate::StrError;
use plotpy::{Canvas, Plot, PolyCode, Text};
use russell_lab::Vector;

/// Draws edge
fn draw_edge(canvas: &mut Canvas, edge_pad: &mut Scratchpad, edge_color: &str) -> Result<(), StrError> {
    const N: usize = 11; // number of points along edge (to handle nonlinear edges)
    let (ksi_min, ksi_del) = (-1.0, 2.0); // all Lin shapes go from -1 to +1
    let space_ndim = edge_pad.xmax.len();
    let mut x = Vector::new(space_ndim);
    canvas.set_face_color("None").set_line_width(3.0);
    if edge_color != "" {
        canvas.set_edge_color(edge_color);
    }
    if space_ndim == 2 {
        canvas.polycurve_begin();
        for i in 0..N {
            let kval = ksi_min + ksi_del * (i as f64) / ((N - 1) as f64);
            calc_coords(&mut x, edge_pad, &[kval])?;
            let code = if i == 0 { PolyCode::MoveTo } else { PolyCode::LineTo };
            canvas.polycurve_add(x[0], x[1], code);
        }
        canvas.polycurve_end(false);
    } else {
        canvas.polyline_3d_begin();
        for i in 0..N {
            let kval = ksi_min + ksi_del * (i as f64) / ((N - 1) as f64);
            calc_coords(&mut x, edge_pad, &[kval])?;
            canvas.polyline_3d_add(x[0], x[1], x[2]);
        }
        canvas.polyline_3d_end();
    }
    Ok(())
}

/// Draws the geometric shape
///
/// # Input
///
/// * `plot` -- plot struct
/// * `pad` -- the scratchpad with coordinates set already
/// * `edge_color` -- color to draw edges; an empty string "" means use default color
/// * `with_ids` -- draw IDs of nodes
/// * `set_range` -- sets axis range
pub fn draw_shape(
    plot: &mut Plot,
    pad: &Scratchpad,
    edge_color: &str,
    with_ids: bool,
    set_range: bool,
) -> Result<(), StrError> {
    if !pad.ok_xxt {
        return Err("all components of the coordinates matrix must be set first");
    }
    let space_ndim = pad.xmin.len();
    let mut canvas = Canvas::new();
    if pad.kind.ndim() == 1 {
        let mut lin_pad = pad.clone();
        draw_edge(&mut canvas, &mut lin_pad, edge_color)?;
    } else {
        let mut edge_pad = Scratchpad::new(space_ndim, pad.kind.edge_kind().unwrap())?;
        for e in 0..pad.kind.nedge() {
            for i in 0..pad.kind.edge_nnode() {
                let m = pad.kind.edge_node_id(e, i);
                for j in 0..space_ndim {
                    edge_pad.set_xx(i, j, pad.xxt[j][m]);
                }
            }
            draw_edge(&mut canvas, &mut edge_pad, edge_color)?;
        }
    }
    plot.add(&canvas);
    if with_ids {
        let mut labels_corner = Text::new();
        let mut labels_middle = Text::new();
        labels_corner
            .set_align_horizontal("center")
            .set_align_vertical("center")
            .set_bbox(true)
            .set_bbox_style("circle,pad=0.3")
            .set_bbox_facecolor("white")
            .set_bbox_edgecolor("red");
        labels_middle
            .set_align_horizontal("center")
            .set_align_vertical("center")
            .set_fontsize(8.0)
            .set_bbox(true)
            .set_bbox_style("circle,pad=0.2")
            .set_bbox_facecolor("white")
            .set_bbox_edgecolor("gold");
        if space_ndim == 2 {
            let m_corner_max = match pad.kind.class() {
                GeoClass::Lin => 2,
                GeoClass::Tri => 3,
                GeoClass::Qua => 4,
                GeoClass::Tet => 4,
                GeoClass::Hex => 8,
            };
            for m in 0..pad.kind.nnode() {
                if m < m_corner_max {
                    labels_corner.draw(pad.xxt[0][m], pad.xxt[1][m], format!("{}", m).as_str());
                } else {
                    labels_middle.draw(pad.xxt[0][m], pad.xxt[1][m], format!("{}", m).as_str());
                }
            }
        } else {
            let m_corner_max = if pad.kind.is_tri_or_tet() { 4 } else { 8 };
            for m in 0..pad.kind.nnode() {
                if m < m_corner_max {
                    labels_corner.draw_3d(pad.xxt[0][m], pad.xxt[1][m], pad.xxt[2][m], format!("{}", m).as_str());
                } else {
                    labels_middle.draw_3d(pad.xxt[0][m], pad.xxt[1][m], pad.xxt[2][m], format!("{}", m).as_str());
                }
            }
        }
        plot.add(&labels_corner).add(&labels_middle);
    }
    if set_range {
        let dx = pad.xmax[0] - pad.xmin[0];
        let dy = pad.xmax[1] - pad.xmin[1];
        const PCT: f64 = 2.0 / 100.0;
        if space_ndim == 2 {
            plot.set_range(
                pad.xmin[0] - dx * PCT,
                pad.xmax[0] + dx * PCT,
                pad.xmin[1] - dy * PCT,
                pad.xmax[1] + dy * PCT,
            );
        } else {
            let dz = pad.xmax[2] - pad.xmin[2];
            plot.set_range_3d(
                pad.xmin[0] - dx * PCT,
                pad.xmax[0] + dx * PCT,
                pad.xmin[1] - dy * PCT,
                pad.xmax[1] + dy * PCT,
                pad.xmin[2] - dz * PCT,
                pad.xmax[2] + dz * PCT,
            );
        }
    }
    Ok(())
}

/// Draws shape (simple version)
pub fn draw_shape_simple(pad: &Scratchpad, filename: &str) -> Result<(), StrError> {
    let mut plot = Plot::new();
    draw_shape(&mut plot, &pad, "", true, true)?;
    plot.set_equal_axes(true)
        .set_figure_size_points(600.0, 600.0)
        .save(filename)?;
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::draw_shape;
    use crate::shapes::op::testing::aux::{gen_canvas_mapping, gen_scratchpad_with_coords};
    use crate::shapes::GeoKind;
    use crate::StrError;
    use plotpy::Plot;

    #[test]
    fn draw_shape_works() -> Result<(), StrError> {
        // select all kinds
        let kinds = vec![
            // Lin
            GeoKind::Lin2,
            GeoKind::Lin3,
            GeoKind::Lin4,
            GeoKind::Lin5,
            // Tri
            GeoKind::Tri3,
            GeoKind::Tri6,
            GeoKind::Tri10,
            GeoKind::Tri15,
            // Qua
            GeoKind::Qua4,
            GeoKind::Qua8,
            GeoKind::Qua9,
            GeoKind::Qua12,
            GeoKind::Qua16,
            GeoKind::Qua17,
            // Tet
            GeoKind::Tet4,
            GeoKind::Tet10,
            GeoKind::Tet20,
            // Hex
            GeoKind::Hex8,
            GeoKind::Hex20,
            GeoKind::Hex32,
        ];
        for kind in kinds {
            let space_ndim = usize::max(2, kind.ndim());
            let pad = gen_scratchpad_with_coords(space_ndim, kind);
            let mut plot = Plot::new();
            if space_ndim == 2 {
                let canvas = gen_canvas_mapping();
                plot.add(&canvas);
                plot.set_range(1.0, 10.0, 1.0, 10.0);
            }
            draw_shape(&mut plot, &pad, "", true, true)?;
            if space_ndim == 2 {
                plot.set_range(1.0, 10.0, 1.0, 10.0);
            }
            if false {
                plot.set_figure_size_points(600.0, 600.0)
                    .set_equal_axes(true)
                    .save(format!("/tmp/gemlab/test_draw_shape_{}.svg", kind.to_string()).as_str())?;
            }
        }
        Ok(())
    }
}
