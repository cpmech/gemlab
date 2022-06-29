#![allow(unused)]

use super::{GS_DEFAULT_BORDER_TOL, GS_DEFAULT_TOLERANCE, SQRT_2, SQRT_3};
use crate::StrError;
use plotpy::{Canvas, Plot, PolyCode, Text};
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Specifies the key of containers (or bins in the grid)
type ContainerKey = usize;

/// Specifies the identification number of items
type ItemID = usize;

/// Defines the container type: ID to Coordinates
type Container = HashSet<ItemID>;

/// Defines the containers type: Key to Container
type Containers = HashMap<ContainerKey, Container>;

/// Defines the bounding box of a cell
const N_MIN_MAX: usize = 2; // 2 means {min,max}
const I_MIN: usize = 0;
const I_MAX: usize = 1;
type BboxMinMax = Vec<Vec<f64>>; // [ndim][N_MIN_MAX]
type BoundingBoxes = HashMap<ItemID, BboxMinMax>;

pub struct GridSearchCell {
    ndim: usize,                   // space dimension
    ndiv: Vec<usize>,              // (ndim) number of divisions along each direction
    xmin: Vec<f64>,                // (ndim) min values
    xmax: Vec<f64>,                // (ndim) max values
    bbox_large: Vec<f64>,          // largest bounding box
    side_length: f64,              // side length of a container
    tol_dist: f64,                 // tolerance to find points using the distance between points
    bounding_boxes: BoundingBoxes, // bounding boxes
    containers: Containers,        // structure to hold all items
}

impl GridSearchCell {
    pub fn new<'a, F, G>(
        ndim: usize,
        ncell: usize,
        get_nnode: F,
        get_x: G,
        tolerance: Option<f64>,
        border_tol: Option<f64>,
    ) -> Result<Self, StrError>
    where
        F: Fn(usize) -> Result<usize, StrError>,
        G: Fn(usize, usize) -> Result<&'a [f64], StrError>,
    {
        // check input
        if ndim < 2 || ndim > 3 {
            return Err("ndim must be 2 or 3");
        }

        // tolerance
        let tolerance = match tolerance {
            Some(v) => v,
            None => GS_DEFAULT_TOLERANCE,
        };
        if tolerance <= 0.0 {
            return Err("tolerance must be > 0.0");
        }

        // border tolerance
        let border_tol = match border_tol {
            Some(v) => v,
            None => GS_DEFAULT_BORDER_TOL,
        };
        if border_tol < 0.0 {
            return Err("border_tol must be ≥ 0.0");
        }

        // find limits, bounding boxes, and largest cell
        let mut xmin = vec![f64::MAX; ndim];
        let mut xmax = vec![f64::MIN; ndim];
        let mut x_min_max = vec![vec![0.0; N_MIN_MAX]; ndim];
        let mut bbox_large = vec![f64::MIN; ndim];
        let mut bounding_boxes = HashMap::new();
        for c in 0..ncell {
            let nnode = get_nnode(c)?;
            for m in 0..nnode {
                let x = get_x(c, m)?;
                if x.len() != ndim {
                    return Err("x.len() must be equal to ndim");
                }
                for i in 0..ndim {
                    // limits
                    xmin[i] = f64::min(xmin[i], x[i]);
                    xmax[i] = f64::max(xmax[i], x[i]);
                    // bounding box
                    if m == 0 {
                        for i in 0..ndim {
                            x_min_max[i][I_MIN] = x[i];
                            x_min_max[i][I_MAX] = x[i];
                        }
                    } else {
                        x_min_max[i][I_MIN] = f64::min(x_min_max[i][I_MIN], x[i]);
                        x_min_max[i][I_MAX] = f64::max(x_min_max[i][I_MAX], x[i]);
                    }
                }
            }
            // largest cell
            for i in 0..ndim {
                bbox_large[i] = f64::max(bbox_large[i], x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
            }
            // add to bounding box maps
            bounding_boxes.insert(c, x_min_max.clone());
        }

        // make the side_length equal to the largest bounding box dimension
        let mut side_length = bbox_large[0];
        for i in 1..ndim {
            side_length = f64::max(side_length, bbox_large[i]);
        }

        // expand side_length by two times the tolerance (left and right)
        side_length += 2.0 * tolerance;

        // expand borders
        if border_tol > 0.0 {
            for i in 0..ndim {
                xmin[i] -= border_tol;
                xmax[i] += border_tol;
            }
        }

        // number of divisions
        let mut ndiv = vec![0; ndim];
        for i in 0..ndim {
            assert!(xmax[i] > xmin[i]);
            ndiv[i] = f64::ceil((xmax[i] - xmin[i]) / side_length) as usize;
        }

        // update xmax after deciding on the side_length and number of divisions
        for i in 0..ndim {
            xmax[i] = xmin[i] + side_length * (ndiv[i] as f64);
        }

        // tolerance to find points
        let tol_dist = if ndim == 2 {
            SQRT_2 * tolerance
        } else {
            SQRT_3 * tolerance
        };

        // insert cells
        let coefficient = vec![1, ndiv[0], ndiv[0] * ndiv[1]];
        let mut containers = HashMap::new();
        let mut ratio = vec![0; ndim]; // ratio = trunc(δx[i]/Δx[i]) (Eq. 8)
        let mut x = vec![0.0; ndim];
        for c in 0..ncell {
            let x_min_max = bounding_boxes.get(&c).unwrap();
            for r in 0..2 {
                for s in 0..2 {
                    for t in 0..(ndim - 1) {
                        x[0] = x_min_max[0][r];
                        x[1] = x_min_max[1][s];
                        if ndim == 3 {
                            x[2] = x_min_max[2][t];
                        }
                        let mut key = 0;
                        for i in 0..ndim {
                            ratio[i] = ((x[i] - xmin[i]) / side_length) as usize;
                            if ratio[i] == ndiv[i] {
                                // the point is exactly on the max edge, thus select inner container
                                ratio[i] -= 1; // move to the inside
                            }
                            key += ratio[i] * coefficient[i];
                        }
                        let container = containers.entry(key).or_insert(HashSet::new());
                        container.insert(c);
                    }
                }
            }
        }

        // allocate grid
        Ok(GridSearchCell {
            ndim,
            ndiv,
            xmin,
            xmax,
            bbox_large,
            side_length,
            tol_dist,
            bounding_boxes,
            containers,
        })
    }

    /// Draws grid and items
    pub fn draw(&self, plot: &mut Plot) -> Result<(), StrError> {
        // draw grid
        let mut xmin = vec![0.0; self.ndim];
        let mut xmax = vec![0.0; self.ndim];
        let mut ndiv = vec![0; self.ndim];
        for i in 0..self.ndim {
            xmin[i] = self.xmin[i];
            xmax[i] = self.xmax[i];
            ndiv[i] = self.ndiv[i];
        }
        let mut canvas = Canvas::new();
        canvas
            .set_alt_text_color("#5d5d5d")
            .draw_grid(&xmin, &xmax, &ndiv, false, true)?;
        plot.add(&canvas);

        // draw bounding boxes
        let mut bbox = Canvas::new();
        bbox.set_edge_color("#3aff79")
            .set_face_color("None")
            .set_line_width(0.75);
        for x_min_max in self.bounding_boxes.values() {
            if self.ndim == 2 {
                bbox.polycurve_begin();
                bbox.polycurve_add(x_min_max[0][I_MIN], x_min_max[1][I_MIN], PolyCode::MoveTo);
                bbox.polycurve_add(x_min_max[0][I_MAX], x_min_max[1][I_MIN], PolyCode::LineTo);
                bbox.polycurve_add(x_min_max[0][I_MAX], x_min_max[1][I_MAX], PolyCode::LineTo);
                bbox.polycurve_add(x_min_max[0][I_MIN], x_min_max[1][I_MAX], PolyCode::LineTo);
                bbox.polycurve_end(true);
            } else {
                bbox.polyline_3d_begin();
                bbox.polyline_3d_add(x_min_max[0][I_MIN], x_min_max[1][I_MIN], x_min_max[2][I_MIN]);
                bbox.polyline_3d_add(x_min_max[0][I_MAX], x_min_max[1][I_MIN], x_min_max[2][I_MIN]);
                bbox.polyline_3d_add(x_min_max[0][I_MAX], x_min_max[1][I_MAX], x_min_max[2][I_MIN]);
                bbox.polyline_3d_add(x_min_max[0][I_MIN], x_min_max[1][I_MAX], x_min_max[2][I_MIN]);
                bbox.polyline_3d_add(x_min_max[0][I_MIN], x_min_max[1][I_MIN], x_min_max[2][I_MIN]);
                bbox.polyline_3d_end();
                bbox.polyline_3d_begin();
                bbox.polyline_3d_add(x_min_max[0][I_MIN], x_min_max[1][I_MIN], x_min_max[2][I_MAX]);
                bbox.polyline_3d_add(x_min_max[0][I_MAX], x_min_max[1][I_MIN], x_min_max[2][I_MAX]);
                bbox.polyline_3d_add(x_min_max[0][I_MAX], x_min_max[1][I_MAX], x_min_max[2][I_MAX]);
                bbox.polyline_3d_add(x_min_max[0][I_MIN], x_min_max[1][I_MAX], x_min_max[2][I_MAX]);
                bbox.polyline_3d_add(x_min_max[0][I_MIN], x_min_max[1][I_MIN], x_min_max[2][I_MAX]);
                bbox.polyline_3d_end();
                bbox.polyline_3d_begin();
                bbox.polyline_3d_add(x_min_max[0][I_MIN], x_min_max[1][I_MIN], x_min_max[2][I_MIN]);
                bbox.polyline_3d_add(x_min_max[0][I_MIN], x_min_max[1][I_MIN], x_min_max[2][I_MAX]);
                bbox.polyline_3d_end();
                bbox.polyline_3d_begin();
                bbox.polyline_3d_add(x_min_max[0][I_MAX], x_min_max[1][I_MIN], x_min_max[2][I_MIN]);
                bbox.polyline_3d_add(x_min_max[0][I_MAX], x_min_max[1][I_MIN], x_min_max[2][I_MAX]);
                bbox.polyline_3d_end();
                bbox.polyline_3d_begin();
                bbox.polyline_3d_add(x_min_max[0][I_MAX], x_min_max[1][I_MAX], x_min_max[2][I_MIN]);
                bbox.polyline_3d_add(x_min_max[0][I_MAX], x_min_max[1][I_MAX], x_min_max[2][I_MAX]);
                bbox.polyline_3d_end();
                bbox.polyline_3d_begin();
                bbox.polyline_3d_add(x_min_max[0][I_MIN], x_min_max[1][I_MAX], x_min_max[2][I_MIN]);
                bbox.polyline_3d_add(x_min_max[0][I_MIN], x_min_max[1][I_MAX], x_min_max[2][I_MAX]);
                bbox.polyline_3d_end();
            }
        }

        // draw items
        let mut text = Text::new();
        text.set_color("#cd0000");
        let mut xcen = vec![0.0; self.ndim];
        for container in self.containers.values() {
            for id in container {
                let x_min_max = self
                    .bounding_boxes
                    .get(&id)
                    .ok_or("INTERNAL ERROR: bounding box is missing")?;
                let txt = format!("{}", id);
                for i in 0..self.ndim {
                    xcen[i] = (x_min_max[i][I_MIN] + x_min_max[i][I_MAX]) / 2.0;
                }
                if self.ndim == 2 {
                    text.draw(xcen[0], xcen[1], &txt);
                } else {
                    text.draw_3d(xcen[0], xcen[1], xcen[2], &txt);
                }
            }
        }
        plot.add(&bbox).add(&text);
        Ok(())
    }
}

impl fmt::Display for GridSearchCell {
    /// Shows info about the items in the grid containers
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // items
        let mut unique_items: HashMap<usize, bool> = HashMap::new();
        let mut indices: Vec<_> = self.containers.keys().collect();
        indices.sort();
        for index in indices {
            let container = self.containers.get(index).unwrap();
            let mut ids: Vec<_> = container.iter().map(|id| *id).collect();
            ids.sort();
            write!(f, "{}: {:?}\n", index, ids).unwrap();
            for id in ids {
                unique_items.insert(id, true);
            }
        }
        // summary
        let mut ids: Vec<_> = unique_items.keys().collect();
        ids.sort();
        write!(f, "ids = {:?}\n", ids).unwrap();
        write!(f, "nitem = {}\n", unique_items.len()).unwrap();
        write!(f, "ncontainer = {}\n", self.containers.len()).unwrap();
        write!(f, "ndiv = {:?}\n", self.ndiv).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::GridSearchCell;
    use crate::util::{SQRT_2, SQRT_3};
    use crate::StrError;
    use plotpy::{Canvas, Curve, Plot, PolyCode, RayEndpoint, Surface};
    use russell_chk::{assert_approx_eq, assert_vec_approx_eq};

    fn draw_triangles(plot: &mut Plot, triangles: &[[[f64; 2]; 3]]) {
        let mut canvas = Canvas::new();
        canvas.set_face_color("#fefddc").set_edge_color("#fcb827");
        for t in 0..triangles.len() {
            canvas.polycurve_begin();
            for m in 0..3 {
                let code = if m == 0 { PolyCode::MoveTo } else { PolyCode::LineTo };
                canvas.polycurve_add(&triangles[t][m][0], &triangles[t][m][1], code);
            }
            canvas.polycurve_end(true);
        }
        plot.add(&canvas);
    }

    fn draw_tetrahedra(plot: &mut Plot, tets: &[[[f64; 3]; 4]]) {
        const EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        let mut canvas = Canvas::new();
        canvas.set_face_color("None").set_edge_color("#fcb827");
        for t in 0..tets.len() {
            for (ma, mb) in &EDGES {
                let xa = &tets[t][*ma];
                let xb = &tets[t][*mb];
                canvas.polyline_3d_begin();
                canvas.polyline_3d_add(xa[0], xa[1], xa[2]);
                canvas.polyline_3d_add(xb[0], xb[1], xb[2]);
                canvas.polyline_3d_end();
            }
        }
        plot.add(&canvas);
    }

    #[test]
    fn new_works_1() -> Result<(), StrError> {
        const TRIANGLES: [[[f64; 2]; 3]; 2] = [
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        ];
        let tolerance = 1e-3;
        let border_tol = 0.1;
        let get_nnode = |_| Ok(3);
        let get_x = |t: usize, m: usize| Ok(&TRIANGLES[t][m][..]);
        let mut grid = GridSearchCell::new(2, TRIANGLES.len(), get_nnode, get_x, Some(tolerance), Some(border_tol))?;
        let max_len = 1.0;
        let side_len = max_len + 2.0 * tolerance; // because the bbox is expanded
        assert_eq!(grid.ndim, 2);
        assert_eq!(grid.side_length, side_len);
        assert_eq!(grid.ndiv, &[2, 2]);
        assert_eq!(grid.xmin, &[-0.1, -0.1]);
        assert_eq!(grid.xmax, &[-0.1 + side_len * 2.0, -0.1 + side_len * 2.0]);
        assert_eq!(grid.bbox_large, &[1.0, 1.0]);
        assert_eq!(grid.tol_dist, SQRT_2 * tolerance);
        assert_eq!(grid.bounding_boxes.len(), 2);
        assert_eq!(grid.containers.len(), 4);
        let bbox_0 = grid.bounding_boxes.get(&0).unwrap();
        let bbox_1 = grid.bounding_boxes.get(&1).unwrap();
        assert_eq!(bbox_0, &[[0.0, 1.0], [0.0, 1.0]]);
        assert_eq!(bbox_1, &[[0.0, 1.0], [0.0, 1.0]]);
        if false {
            let mut plot = Plot::new();
            draw_triangles(&mut plot, &TRIANGLES);
            grid.draw(&mut plot)?;
            plot.set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .grid_and_labels("x", "y")
                .set_ticks_x(0.2, 0.0, "")
                .set_ticks_y(0.2, 0.0, "")
                .save("/tmp/gemlab/test_grid_search_cell_new_1.svg")?;
        }
        Ok(())
    }

    #[test]
    fn new_works_2() -> Result<(), StrError> {
        const TRIANGLES: [[[f64; 2]; 3]; 2] = [
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            [[1.0, 0.0], [1.2, 1.5], [0.0, 1.0]],
        ];
        let tolerance = 1e-3;
        let border_tol = 0.1;
        let get_nnode = |_| Ok(3);
        let get_x = |t: usize, m: usize| Ok(&TRIANGLES[t][m][..]);
        let mut grid = GridSearchCell::new(2, TRIANGLES.len(), get_nnode, get_x, Some(tolerance), Some(border_tol))?;
        let max_len = 1.5;
        let side_len = max_len + 2.0 * tolerance; // because the bbox is expanded
        assert_eq!(grid.ndim, 2);
        assert_eq!(grid.side_length, side_len);
        assert_eq!(grid.ndiv, &[1, 2]);
        assert_eq!(grid.xmin, &[-0.1, -0.1]);
        assert_eq!(grid.xmax, &[-0.1 + side_len, -0.1 + side_len * 2.0]);
        assert_eq!(grid.bbox_large, &[1.2, 1.5]);
        assert_eq!(grid.tol_dist, SQRT_2 * tolerance);
        assert_eq!(grid.bounding_boxes.len(), 2);
        assert_eq!(grid.containers.len(), 2);
        let bbox_0 = grid.bounding_boxes.get(&0).unwrap();
        let bbox_1 = grid.bounding_boxes.get(&1).unwrap();
        assert_eq!(bbox_0, &[[0.0, 1.0], [0.0, 1.0]]);
        assert_eq!(bbox_1, &[[0.0, 1.2], [0.0, 1.5]]);
        if false {
            let mut plot = Plot::new();
            draw_triangles(&mut plot, &TRIANGLES);
            grid.draw(&mut plot)?;
            plot.set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .grid_and_labels("x", "y")
                .set_ticks_x(0.2, 0.0, "")
                .set_ticks_y(0.2, 0.0, "")
                .save("/tmp/gemlab/test_grid_search_cell_new_2.svg")?;
        }
        Ok(())
    }

    #[test]
    fn new_works_3() -> Result<(), StrError> {
        const TETS: [[[f64; 3]; 4]; 1] = [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]];
        let tolerance = 1e-3;
        let border_tol = 0.1;
        let get_nnode = |_| Ok(4);
        let get_x = |t: usize, m: usize| Ok(&TETS[t][m][..]);
        let mut grid = GridSearchCell::new(3, TETS.len(), get_nnode, get_x, Some(tolerance), Some(border_tol))?;
        let max_len = 1.0;
        let side_len = max_len + 2.0 * tolerance; // because the bbox is expanded
        assert_eq!(grid.ndim, 3);
        assert_eq!(grid.side_length, side_len);
        assert_eq!(grid.ndiv, &[2, 2, 2]);
        assert_eq!(grid.xmin, &[-0.1, -0.1, -0.1]);
        assert_eq!(
            grid.xmax,
            &[-0.1 + side_len * 2.0, -0.1 + side_len * 2.0, -0.1 + side_len * 2.0]
        );
        assert_eq!(grid.bbox_large, &[1.0, 1.0, 1.0]);
        assert_eq!(grid.tol_dist, SQRT_3 * tolerance);
        assert_eq!(grid.bounding_boxes.len(), 1);
        assert_eq!(grid.containers.len(), 8);
        let bbox_0 = grid.bounding_boxes.get(&0).unwrap();
        assert_eq!(bbox_0, &[[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]);
        if false {
            let mut plot = Plot::new();
            draw_tetrahedra(&mut plot, &TETS);
            grid.draw(&mut plot)?;
            plot.set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .save("/tmp/gemlab/test_grid_search_cell_new_3.svg")?;
        }
        Ok(())
    }

    #[test]
    fn gs_cell_works_2d() -> Result<(), StrError> {
        // [num_triangle][nnode=3][ndim=2]
        const TRIANGLES: [[[f64; 2]; 3]; 12] = [
            [[0.230951, 0.558482], [0.133721, 0.348832], [0.540745, 0.331184]],
            [[0.13928, 0.180603], [0.133721, 0.348832], [0.0307942, 0.459123]],
            [[0.0307942, 0.459123], [0.230951, 0.558482], [0.0980015, 0.981755]],
            [[0.230951, 0.558482], [0.0307942, 0.459123], [0.133721, 0.348832]],
            [[0.0980015, 0.981755], [0.230951, 0.558482], [0.578587, 0.760349]],
            [[0.133721, 0.348832], [0.13928, 0.180603], [0.540745, 0.331184]],
            [[0.540745, 0.331184], [0.578587, 0.760349], [0.230951, 0.558482]],
            [[0.540745, 0.331184], [0.478554, 0.00869692], [0.648071, 0.369534]],
            [[0.578587, 0.760349], [0.648071, 0.369534], [0.903726, 0.975904]],
            [[0.648071, 0.369534], [0.578587, 0.760349], [0.540745, 0.331184]],
            [[0.578587, 0.760349], [0.903726, 0.975904], [0.0980015, 0.981755]],
            [[0.540745, 0.331184], [0.13928, 0.180603], [0.478554, 0.00869692]],
        ];
        let get_nnode = |_| Ok(3);
        let get_x = |t: usize, m: usize| Ok(&TRIANGLES[t][m][..]);
        let mut grid = GridSearchCell::new(2, TRIANGLES.len(), get_nnode, get_x, None, None)?;
        if false {
            let mut plot = Plot::new();
            draw_triangles(&mut plot, &TRIANGLES);
            grid.draw(&mut plot)?;
            plot.set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .grid_and_labels("x", "y")
                .set_ticks_x(0.2, 0.0, "")
                .set_ticks_y(0.2, 0.0, "")
                .save("/tmp/gemlab/test_grid_search_cell_works_2d.svg")?;
        }
        Ok(())
    }
}
