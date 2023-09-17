use super::{Mesh, PointId};
use crate::shapes::{GeoKind, Scratchpad};
use crate::StrError;
use plotpy::{Canvas, PolyCode};
use russell_lab::math::ONE_BY_3;
use russell_lab::Vector;
use std::collections::HashMap;

impl Mesh {
    /// Draws a single cell specified by an array of point ids
    ///
    /// # Input
    ///
    /// * `canvas` -- where to draw the cell
    /// * `cell_kind` -- the GeoKind of the cell
    /// * `cell_points` -- the connectivity list (IDs of points)
    /// * `pads` -- an access to a map of Scratchpads (exploiting memoization)
    #[rustfmt::skip]
    pub(crate) fn draw_cell(
        &self,
        canvas: &mut Canvas,
        cell_kind: GeoKind,
        cell_points: &Vec<PointId>,
        pads: &mut HashMap<GeoKind, Scratchpad>,
    ) -> Result<(), StrError> {
        let ii = cell_points;
        match cell_kind {
            // Lin
            GeoKind::Lin2 => {
                self.add_curve(canvas, cell_kind, cell_points, true, true, pads)?;
            }
            GeoKind::Lin3 => {
                self.add_curve(canvas, cell_kind, cell_points, true, true, pads)?;
            }
            GeoKind::Lin4 => {
                self.add_curve(canvas, cell_kind, cell_points, true, true, pads)?;
            }
            GeoKind::Lin5 => {
                self.add_curve(canvas, cell_kind, cell_points, true, true, pads)?;
            }
            // Tri
            GeoKind::Tri3 => {
                self.add_curve(canvas, GeoKind::Lin2, &[ii[0], ii[1]], true, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin2, &[ii[1], ii[2]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin2, &[ii[2], ii[0]], false, true, pads)?;
            }
            GeoKind::Tri6 => {
                self.add_curve(canvas, GeoKind::Lin3, &[ii[0], ii[1], ii[3]], true, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin3, &[ii[1], ii[2], ii[4]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin3, &[ii[2], ii[0], ii[5]], false, true, pads)?;
            }
            GeoKind::Tri10 => {
                self.add_curve(canvas, GeoKind::Lin4, &[ii[0], ii[1], ii[3], ii[6]], true, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin4, &[ii[1], ii[2], ii[4], ii[7]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin4, &[ii[2], ii[0], ii[5], ii[8]], false, true, pads)?;
            }
            GeoKind::Tri15 => {
                self.add_curve(canvas, GeoKind::Lin5, &[ii[0], ii[1], ii[3], ii[6], ii[7]], true, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin5, &[ii[1], ii[2], ii[4], ii[8], ii[9]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin5, &[ii[2], ii[0], ii[5], ii[10], ii[11]], false, true, pads)?;
            }
            // Qua
            GeoKind::Qua4 => {
                self.add_curve(canvas, GeoKind::Lin2, &[ii[0], ii[1]], true, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin2, &[ii[1], ii[2]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin2, &[ii[2], ii[3]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin2, &[ii[3], ii[0]], false, true, pads)?;
            }
            GeoKind::Qua8 => {
                self.add_curve(canvas, GeoKind::Lin3, &[ii[0], ii[1], ii[4]], true, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin3, &[ii[1], ii[2], ii[5]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin3, &[ii[2], ii[3], ii[6]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin3, &[ii[3], ii[0], ii[7]], false, true, pads)?;
            }
            GeoKind::Qua9 => {
                self.add_curve(canvas, GeoKind::Lin3, &[ii[0], ii[1], ii[4]], true, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin3, &[ii[1], ii[2], ii[5]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin3, &[ii[2], ii[3], ii[6]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin3, &[ii[3], ii[0], ii[7]], false, true, pads)?;
            }
            GeoKind::Qua12 => {
                self.add_curve(canvas, GeoKind::Lin4, &[ii[0], ii[1], ii[4], ii[8]], true, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin4, &[ii[1], ii[2], ii[5], ii[9]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin4, &[ii[2], ii[3], ii[6], ii[10]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin4, &[ii[3], ii[0], ii[7], ii[11]], false, true, pads)?;
            }
            GeoKind::Qua16 => {
                self.add_curve(canvas, GeoKind::Lin4, &[ii[0], ii[1], ii[4], ii[8]], true, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin4, &[ii[1], ii[2], ii[5], ii[9]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin4, &[ii[2], ii[3], ii[6], ii[10]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin4, &[ii[3], ii[0], ii[7], ii[11]], false, true, pads)?;
            }
            GeoKind::Qua17 => {
                self.add_curve(canvas, GeoKind::Lin5, &[ii[0], ii[1], ii[4], ii[9], ii[10]], true, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin5, &[ii[1], ii[2], ii[5], ii[11], ii[12]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin5, &[ii[2], ii[3], ii[6], ii[13], ii[14]], false, false, pads)?;
                self.add_curve(canvas, GeoKind::Lin5, &[ii[3], ii[0], ii[7], ii[15], ii[16]], false, true, pads)?;
            }
            // Tet
            GeoKind::Tet4 => {
                for e in 0..cell_kind.nedge() {
                    let a = cell_kind.edge_node_id(e, 0);
                    let b = cell_kind.edge_node_id(e, 1);
                    self.add_curve(canvas, GeoKind::Lin2, &[ii[a], ii[b]], true, true, pads)?;
                }
            }
            GeoKind::Tet10 => {
                for e in 0..cell_kind.nedge() {
                    let a = cell_kind.edge_node_id(e, 0);
                    let b = cell_kind.edge_node_id(e, 1);
                    let c = cell_kind.edge_node_id(e, 2);
                    self.add_curve(canvas, GeoKind::Lin2, &[ii[a], ii[b], ii[c]], true, true, pads)?;
                }
            }
            GeoKind::Tet20 => {
                for e in 0..cell_kind.nedge() {
                    let a = cell_kind.edge_node_id(e, 0);
                    let b = cell_kind.edge_node_id(e, 1);
                    let c = cell_kind.edge_node_id(e, 2);
                    let d = cell_kind.edge_node_id(e, 3);
                    self.add_curve(canvas, GeoKind::Lin3, &[ii[a], ii[b], ii[c], ii[d]], true, true, pads)?;
                }
            }
            // Hex
            GeoKind::Hex8 => {
                for e in 0..cell_kind.nedge() {
                    let a = cell_kind.edge_node_id(e, 0);
                    let b = cell_kind.edge_node_id(e, 1);
                    self.add_curve(canvas, GeoKind::Lin2, &[ii[a], ii[b]], true, true, pads)?;
                }
            }
            GeoKind::Hex20 => {
                for e in 0..cell_kind.nedge() {
                    let a = cell_kind.edge_node_id(e, 0);
                    let b = cell_kind.edge_node_id(e, 1);
                    let c = cell_kind.edge_node_id(e, 2);
                    self.add_curve(canvas, GeoKind::Lin2, &[ii[a], ii[b], ii[c]], true, true, pads)?;
                }
            }
            GeoKind::Hex32 => {
                for e in 0..cell_kind.nedge() {
                    let a = cell_kind.edge_node_id(e, 0);
                    let b = cell_kind.edge_node_id(e, 1);
                    let c = cell_kind.edge_node_id(e, 2);
                    let d = cell_kind.edge_node_id(e, 3);
                    self.add_curve(canvas, GeoKind::Lin3, &[ii[a], ii[b], ii[c], ii[d]], true, true, pads)?;
                }
            }
        }
        Ok(())
    }

    /// Adds a Polyline or Bezier curve to the canvas
    ///
    /// # Input
    ///
    /// * `canvas` -- the canvas
    /// * `lin_kind` -- the type of line
    /// * `begin` -- a flag telling to begin the polyline/Bezier
    /// * `end` -- a flag telling to stop the polyline/Bezier
    /// * `pads` -- an access to a map of Scratchpads (exploiting memoization)
    fn add_curve(
        &self,
        canvas: &mut Canvas,
        lin_kind: GeoKind,
        lin_points: &[PointId],
        begin: bool,
        end: bool,
        pads: &mut HashMap<GeoKind, Scratchpad>,
    ) -> Result<(), StrError> {
        let pp = &self.points;
        let ii = lin_points;
        let mut x = Vector::new(self.ndim);
        match lin_kind {
            GeoKind::Lin2 => {
                if self.ndim == 2 {
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
                if self.ndim == 2 {
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
                if self.ndim == 2 {
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
                if self.ndim == 2 {
                    // add poly-curve points or control points (q..) for cubic Bezier
                    // (see file bezier-curves-math.pdf under data/derivations)
                    // in this case, we have to interpolate the Lin5 to make it a Lin4
                    let pad = pads.entry(lin_kind).or_insert(Scratchpad::new(self.ndim, lin_kind)?);
                    for m in 0..lin_points.len() {
                        for j in 0..self.ndim {
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
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::mesh::Samples;
    use crate::shapes::GeoKind;
    use plotpy::Canvas;
    use std::collections::HashMap;

    #[test]
    fn add_curve_catches_errors() {
        let mesh = Samples::one_lin2();
        let mut pads = HashMap::new();
        let mut canvas = Canvas::new();
        assert_eq!(
            mesh.add_curve(&mut canvas, GeoKind::Tri3, &[], true, true, &mut pads)
                .err(),
            Some("lin_kind is not Lin")
        );
    }
}
