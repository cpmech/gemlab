use super::{draw_cell, Cell, Mesh};
use crate::util::calc_container_key;
use crate::StrError;
use plotpy::{Canvas, Plot, PolyCode, Text};
use russell_stat::Histogram;
use std::collections::{HashMap, HashSet};
use std::fmt::{self, Write};

/// Default tolerance for all directions
pub const GRID_CELLS_TOLERANCE: f64 = 1e-4;

/// Default border tolerance to handle imprecision near the borders
pub const GRID_CELLS_BORDER_TOL: f64 = 1e-2;

/// Defines a bounding-box ratio in [0,1] to mark cells as "large" compared others
///
/// Example: If a cell has a max(bounding_box_length) / largest >= 0.75,
///          Then this cell is put into the large_cells list
///
/// This constant enables the strategy of putting some cells in a selected list of large cells.
/// The other constant [GRID_CELLS_LARGE_PCT] also influences the algorithm by enabling or disabling this strategy.
const GRID_CELLS_LARGE_CUTOFF: f64 = 0.75;

/// Controls the count percentage in [0,1] of large cells allowed.
///
/// Example: If 20% of cells are too large, then we put them in a separated list of large cells;
///          otherwise, the mesh is homogeneous enough so put all cells in the containers.
///
/// If the actual ratio is greater than this constant (i.e., many cells are large such as
/// in a homogeneous mesh), the "strategy" of selecting large cells is abandoned.
/// The other constant [GRID_CELLS_LARGE_CUTOFF] also influences the algorithm.
const GRID_CELLS_LARGE_PCT: f64 = 0.2;

/// Specifies the key of containers (or bins in the grid)
type ContainerKey = usize;

/// Specifies the identification number of cells (must be sequential from 0 to ncell - 1)
type CellId = usize;

/// Defines the container type
type Container = HashSet<CellId>;

/// Defines the containers type: Key to Container
type Containers = HashMap<ContainerKey, Container>;

/// Defines the bounding box of a cell
const N_MIN_MAX: usize = 2; // 2 means {min,max}
const I_MIN: usize = 0;
const I_MAX: usize = 1;
type BboxMinMax = Vec<Vec<f64>>; // [ndim][N_MIN_MAX]

pub struct GridCells<'a> {
    mesh: &'a Mesh,                  // the mesh
    ndiv: Vec<usize>,                // (ndim) number of divisions along each direction
    xmin: Vec<f64>,                  // (ndim) min values
    xmax: Vec<f64>,                  // (ndim) max values
    side_length: f64,                // side length of a container
    bounding_boxes: Vec<BboxMinMax>, // (ncell) bounding boxes
    containers: Containers,          // structure to hold all items
    large_cells: Container,          // holds the CellId of large cells
}

impl<'a> GridCells<'a> {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `tolerance` -- is a tolerance to expand the bounding box of triangles and compare points; e.g. 1e-4
    ///     - If None, [GRID_CELLS_TOLERANCE] is used
    /// * `border_tol` -- is a tolerance used to expand the border a little bit and then
    ///   accommodate eventual imprecision near the borders; e.g. 1e-2
    ///     - If None, [GRID_CELLS_BORDER_TOL] is used
    /// * `large_cutoff` -- is the ratio to define large cells compared to the largest bounding box length; e.g. 0.75
    ///     - If None, [GRID_CELLS_LARGE_CUTOFF] is used
    /// * `large_pct` -- is the ratio (percentage) of allowed large cells; e.g. 0.2
    ///     - If None, [GRID_CELLS_LARGE_PCT] is used
    pub fn new(
        mesh: &'a Mesh,
        tolerance: Option<f64>,
        border_tol: Option<f64>,
        large_cutoff: Option<f64>,
        large_pct: Option<f64>,
    ) -> Result<Self, StrError> {
        // constants
        let ndim = mesh.ndim;
        let ncell = mesh.cells.len();

        // tolerance
        let tolerance = match tolerance {
            Some(v) => v,
            None => GRID_CELLS_TOLERANCE,
        };
        if tolerance <= 0.0 {
            return Err("tolerance must be > 0.0");
        }

        // border tolerance
        let border_tol = match border_tol {
            Some(v) => v,
            None => GRID_CELLS_BORDER_TOL,
        };
        if border_tol < 0.0 {
            return Err("border_tol must be ≥ 0.0");
        }

        // large cells cutoff
        let large_cells_cutoff = match large_cutoff {
            Some(v) => v,
            None => GRID_CELLS_LARGE_CUTOFF,
        };
        if large_cells_cutoff < 0.0 || large_cells_cutoff > 1.0 {
            return Err("large_cutoff must be in [0.0, 1.0]");
        }

        // large cells pct
        let large_cells_pct = match large_pct {
            Some(v) => v,
            None => GRID_CELLS_LARGE_PCT,
        };
        if large_cells_pct < 0.0 || large_cells_pct > 1.0 {
            return Err("large_pct must be in [0.0, 1.0]");
        }

        // find limits, bounding boxes, and largest bounding box
        let mut x = vec![0.0; ndim];
        let mut xmin = vec![f64::MAX; ndim];
        let mut xmax = vec![f64::MIN; ndim];
        let mut x_min_max = vec![vec![0.0; N_MIN_MAX]; ndim];
        let mut bbox_largest = vec![f64::MIN; ndim];
        let mut bounding_boxes = Vec::new();
        for cell in &mesh.cells {
            for m in 0..cell.points.len() {
                let point_id = cell.points[m];
                for i in 0..ndim {
                    // coords
                    x[i] = mesh.points[point_id].coords[i];
                    // limits
                    xmin[i] = f64::min(xmin[i], x[i]);
                    xmax[i] = f64::max(xmax[i], x[i]);
                    // bounding box
                    if m == 0 {
                        x_min_max[i][I_MIN] = x[i];
                        x_min_max[i][I_MAX] = x[i];
                    } else {
                        x_min_max[i][I_MIN] = f64::min(x_min_max[i][I_MIN], x[i]);
                        x_min_max[i][I_MAX] = f64::max(x_min_max[i][I_MAX], x[i]);
                    }
                }
            }
            // largest bounding box
            for i in 0..ndim {
                bbox_largest[i] = f64::max(bbox_largest[i], x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
            }
            // add to bounding box maps
            bounding_boxes.push(x_min_max.clone());
        }

        // compute the largest length of the largest bounding box
        let mut bbox_largest_length = bbox_largest[0];
        for i in 1..ndim {
            bbox_largest_length = f64::max(bbox_largest_length, bbox_largest[i]);
        }

        // handle large cells
        let mut bbox_large_length = f64::MIN; // largest length lower than the cutoff
        let mut large_cells = HashSet::new();
        for cell_id in 0..ncell {
            let x_min_max = &bounding_boxes[cell_id];
            let mut largest_length = x_min_max[0][I_MAX] - x_min_max[0][I_MIN];
            for i in 1..ndim {
                largest_length = f64::max(largest_length, x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
            }
            if largest_length >= large_cells_cutoff * bbox_largest_length {
                large_cells.insert(cell_id);
            } else {
                for i in 0..ndim {
                    bbox_large_length = f64::max(bbox_large_length, x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
                }
            }
        }

        // make the side_length equal to the largest bounding box dimension
        let pct = (large_cells.len() as f64) / (ncell as f64);
        let mut side_length = if pct <= large_cells_pct {
            bbox_large_length // use the largest length lower than the cutoff
        } else {
            large_cells.clear(); // abandon the large cells strategy
            bbox_largest_length // use the largest length among them all
        };

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

        // insert (not large) cells into containers
        let mut containers = HashMap::new();
        for cell_id in 0..ncell {
            if large_cells.contains(&cell_id) {
                continue; // skip large cells
            }
            let x_min_max = &bounding_boxes[cell_id];
            for r in 0..N_MIN_MAX {
                for s in 0..N_MIN_MAX {
                    for t in 0..(ndim - 1) {
                        x[0] = x_min_max[0][r];
                        x[1] = x_min_max[1][s];
                        if ndim == 3 {
                            x[2] = x_min_max[2][t];
                        }
                        let key = calc_container_key(ndim, side_length, &ndiv, &xmin, &x);
                        let container = containers.entry(key).or_insert(HashSet::new());
                        container.insert(cell_id);
                    }
                }
            }
        }

        // allocate grid
        Ok(GridCells {
            mesh,
            ndiv,
            xmin,
            xmax,
            side_length,
            bounding_boxes,
            containers,
            large_cells,
        })
    }

    /// Find the cell (e.g., triangle or tetrahedron) where the given coordinate falls in
    ///
    /// # Output
    ///
    /// * Returns the CellId or None if no cell contains the point
    ///
    /// # Input
    ///
    /// * `in_cell` is a function(cell_id, x) that tells whether the point is in the cell or not.
    ///     - If None and the cells are Triangle of Tetrahedron, the proper function is selected automatically
    pub fn find_cell<F>(&self, _x: &[f64], _in_cell: Option<F>) -> Result<Option<CellId>, StrError>
    where
        F: Fn(usize, &[f64]) -> Result<bool, StrError>,
    {
        Err("")
    }

    /// Returns the grid limits (xmin,xmax)
    pub fn limits(&self) -> (Vec<f64>, Vec<f64>) {
        (self.xmin.clone(), self.xmax.clone())
    }

    /// Draws grid and items
    pub fn draw(&self, plot: &mut Plot, with_ids: bool) -> Result<(), StrError> {
        // define function to draw ids
        let mut ids = Text::new();
        ids.set_color("#f16100");
        let ndim = self.mesh.ndim;
        let mut xcen = vec![0.0; ndim];
        let mut draw_ids = |cell: &Cell| {
            let txt = format!("{}", cell.id);
            for m in 0..cell.points.len() {
                let x = &self.mesh.points[cell.points[m]].coords;
                for i in 0..ndim {
                    if m == 0 {
                        xcen[i] = x[i];
                    } else {
                        xcen[i] += x[i];
                    }
                }
            }
            for i in 0..ndim {
                xcen[i] /= cell.points.len() as f64;
            }
            if ndim == 2 {
                ids.draw(xcen[0], xcen[1], &txt);
            } else {
                ids.draw_3d(xcen[0], xcen[1], xcen[2], &txt);
            }
        };

        // draw cells and ids
        let mut canvas_cells = Canvas::new();
        canvas_cells.set_edge_color("#fcb827").set_face_color("#fefddc");
        let mut pads = HashMap::new();
        for container in self.containers.values() {
            for cell_id in container {
                let cell = &self.mesh.cells[*cell_id];
                draw_cell(&mut canvas_cells, self.mesh, cell.kind, &cell.points, &mut pads)?;
                if with_ids {
                    draw_ids(cell);
                }
            }
        }
        canvas_cells.set_face_color("#fffca3");
        for cell_id in &self.large_cells {
            let cell = &self.mesh.cells[*cell_id];
            draw_cell(&mut canvas_cells, self.mesh, cell.kind, &cell.points, &mut pads)?;
            if with_ids {
                draw_ids(cell);
            }
        }
        plot.add(&canvas_cells);
        if with_ids {
            plot.add(&ids);
        }

        // draw grid
        let ndim = self.mesh.ndim;
        let mut xmin = vec![0.0; ndim];
        let mut xmax = vec![0.0; ndim];
        let mut ndiv = vec![0; ndim];
        for i in 0..ndim {
            xmin[i] = self.xmin[i];
            xmax[i] = self.xmax[i];
            ndiv[i] = self.ndiv[i];
        }
        let mut canvas_grid = Canvas::new();
        canvas_grid
            .set_alt_text_color("#5d5d5d")
            .draw_grid(&xmin, &xmax, &ndiv, false, with_ids)?;
        plot.add(&canvas_grid);

        // draw bounding boxes
        let mut bbox = Canvas::new();
        bbox.set_edge_color("#3aff79")
            .set_face_color("None")
            .set_line_width(0.75);
        for x_min_max in &self.bounding_boxes {
            if ndim == 2 {
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
        plot.add(&bbox);
        Ok(())
    }

    /// Compute some statistics
    pub fn stat(&self) -> String {
        let mut unique_items = HashSet::new();
        let mut max_container_num_items = 0;
        let mut container_num_items = Vec::new();
        for container in self.containers.values() {
            let n = container.len();
            max_container_num_items = usize::max(max_container_num_items, n);
            container_num_items.push(n);
            for id in container {
                unique_items.insert(id);
            }
        }
        let mut b = String::new();
        let ncell = unique_items.len() + self.large_cells.len();
        write!(&mut b, "Summary\n").unwrap();
        write!(&mut b, "=======\n").unwrap();
        write!(&mut b, "ncell total = {:?}\n", ncell).unwrap();
        write!(&mut b, "xmin = {:?}\n", self.xmin).unwrap();
        write!(&mut b, "xmax = {:?}\n", self.xmax).unwrap();
        write!(&mut b, "side_length = {:?}\n", self.side_length).unwrap();
        write!(&mut b, "num of non-empty containers = {}\n", self.containers.len()).unwrap();
        write!(&mut b, "max container num items = {}\n", max_container_num_items).unwrap();
        write!(&mut b, "\nHistogram of container num items\n").unwrap();
        write!(&mut b, "================================\n").unwrap();
        let stations: Vec<_> = (0..max_container_num_items + 2).collect();
        let mut hist = Histogram::new(&stations).unwrap();
        hist.count(&container_num_items);
        hist.set_bar_char('*');
        hist.set_bar_max_len(60);
        write!(&mut b, "{}", hist).unwrap();
        b
    }
}

impl<'a> fmt::Display for GridCells<'a> {
    /// Shows info about GridCells
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // items
        let mut unique_items = HashSet::new();
        let mut indices: Vec<_> = self.containers.keys().collect();
        indices.sort();
        for index in indices {
            let container = self.containers.get(index).unwrap();
            let mut ids: Vec<_> = container.iter().map(|id| *id).collect();
            ids.sort();
            write!(f, "{}: {:?}\n", index, ids).unwrap();
            for id in ids {
                unique_items.insert(id);
            }
        }
        // summary
        let mut large_ids: Vec<_> = self.large_cells.iter().collect();
        let mut ids: Vec<_> = unique_items.iter().collect();
        large_ids.sort();
        ids.sort();
        write!(f, "large_cells = {:?}\n", large_ids).unwrap();
        write!(f, "ncell = {}\n", unique_items.len() + self.large_cells.len()).unwrap();
        write!(f, "ncontainer = {}\n", self.containers.len()).unwrap();
        write!(f, "ndiv = {:?}\n", self.ndiv).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::GridCells;
    use crate::geometry::{in_triangle, triangle_coords};
    use crate::mesh::{Cell, Mesh, Point, Samples, GRID_CELLS_BORDER_TOL, GRID_CELLS_TOLERANCE};
    use crate::shapes::GeoKind;

    #[allow(unused_imports)]
    use plotpy::Plot;

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        assert_eq!(
            GridCells::new(&mesh, Some(-1.0), None, None, None).err(),
            Some("tolerance must be > 0.0")
        );
        assert_eq!(
            GridCells::new(&mesh, None, Some(-1.0), None, None).err(),
            Some("border_tol must be ≥ 0.0")
        );
        assert_eq!(
            GridCells::new(&mesh, None, None, Some(-1.0), None).err(),
            Some("large_cutoff must be in [0.0, 1.0]")
        );
        assert_eq!(
            GridCells::new(&mesh, None, None, Some(2.0), None).err(),
            Some("large_cutoff must be in [0.0, 1.0]")
        );
        assert_eq!(
            GridCells::new(&mesh, None, None, None, Some(-1.0)).err(),
            Some("large_pct must be in [0.0, 1.0]")
        );
        assert_eq!(
            GridCells::new(&mesh, None, None, None, Some(2.0)).err(),
            Some("large_pct must be in [0.0, 1.0]")
        );
    }

    #[test]
    fn new_works_1() {
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0] },
                Point { id: 3, coords: vec![0.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri3, points: vec![0, 1, 3] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Tri3, points: vec![1, 2, 3] },
            ],
        };
        let tolerance = 1e-3;
        let border_tol = 0.1;
        let grid = GridCells::new(&mesh, Some(tolerance), Some(border_tol), None, None).unwrap();
        let max_len = 1.0;
        let sl = max_len + 2.0 * tolerance; // because the bbox is expanded
        let xmin = &[-0.1, -0.1];
        let xmax = &[-0.1 + sl * 2.0, -0.1 + sl * 2.0];
        assert_eq!(grid.side_length, sl);
        assert_eq!(grid.ndiv, &[2, 2]);
        assert_eq!(grid.xmin, xmin);
        assert_eq!(grid.xmax, xmax);
        assert_eq!(grid.bounding_boxes.len(), 2);
        assert_eq!(grid.containers.len(), 4);
        let bbox_0 = &grid.bounding_boxes[0];
        let bbox_1 = &grid.bounding_boxes[1];
        assert_eq!(bbox_0, &[[0.0, 1.0], [0.0, 1.0]]);
        assert_eq!(bbox_1, &[[0.0, 1.0], [0.0, 1.0]]);
        let container_0 = grid.containers.get(&0).unwrap();
        let container_1 = grid.containers.get(&1).unwrap();
        let container_2 = grid.containers.get(&2).unwrap();
        let container_3 = grid.containers.get(&3).unwrap();
        assert_eq!(container_0.len(), 2);
        assert_eq!(container_1.len(), 2);
        assert_eq!(container_2.len(), 2);
        assert_eq!(container_3.len(), 2);
        assert_eq!(
            format!("{}", grid),
            "0: [0, 1]\n\
             1: [0, 1]\n\
             2: [0, 1]\n\
             3: [0, 1]\n\
             large_cells = []\n\
             ncell = 2\n\
             ncontainer = 4\n\
             ndiv = [2, 2]\n"
        );
        assert_eq!(
            format!("{}", grid.stat()),
            "Summary\n\
             =======\n\
             ncell total = 2\n\
             xmin = [-0.1, -0.1]\n\
             xmax = [1.904, 1.904]\n\
             side_length = 1.002\n\
             num of non-empty containers = 4\n\
             max container num items = 2\n\
             \n\
             Histogram of container num items\n\
             ================================\n\
             [0,1) | 0 \n\
             [1,2) | 0 \n\
             [2,3) | 4 ************************************************************\n\
             \x20\x20sum = 4\n"
        );
        let (mi, ma) = grid.limits();
        assert_eq!(mi, xmin);
        assert_eq!(ma, xmax);

        let mut plot = Plot::new();
        grid.draw(&mut plot, true).unwrap();
        // if true {
        //     plot.set_equal_axes(true)
        //         .set_figure_size_points(600.0, 600.0)
        //         .grid_and_labels("x", "y")
        //         .set_ticks_x(0.2, 0.0, "")
        //         .set_ticks_y(0.2, 0.0, "")
        //         .save("/tmp/gemlab/test_grid_cells_new_works_1.svg")
        //         .unwrap();
        // }

        // with zero border
        let grid = GridCells::new(&mesh, None, Some(0.0), None, None).unwrap();
        let sl = max_len + 2.0 * GRID_CELLS_TOLERANCE;
        assert_eq!(grid.xmin, &[0.0, 0.0]);
        assert_eq!(grid.ndiv, &[1, 1]);
        assert_eq!(grid.xmax, &[sl, sl]);
    }

    #[test]
    fn new_works_2() {
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.2, 1.5] },
                Point { id: 3, coords: vec![0.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri3, points: vec![0, 1, 3] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Tri3, points: vec![1, 2, 3] },
            ],
        };
        let tolerance = 1e-3;
        let border_tol = 0.1;
        let grid = GridCells::new(&mesh, Some(tolerance), Some(border_tol), None, None).unwrap();
        let max_len = 1.5;
        let sl = max_len + 2.0 * tolerance; // because the bbox is expanded
        assert_eq!(grid.side_length, sl);
        assert_eq!(grid.ndiv, &[1, 2]);
        assert_eq!(grid.xmin, &[-0.1, -0.1]);
        assert_eq!(grid.xmax, &[-0.1 + sl, -0.1 + sl * 2.0]);
        assert_eq!(grid.bounding_boxes.len(), 2);
        assert_eq!(grid.containers.len(), 2);
        let bbox_0 = &grid.bounding_boxes[0];
        let bbox_1 = &grid.bounding_boxes[1];
        assert_eq!(bbox_0, &[[0.0, 1.0], [0.0, 1.0]]);
        assert_eq!(bbox_1, &[[0.0, 1.2], [0.0, 1.5]]);
        let container_0 = grid.containers.get(&0).unwrap();
        let container_1 = grid.containers.get(&1).unwrap();
        assert_eq!(container_0.len(), 2);
        assert_eq!(container_1.len(), 1);
        assert_eq!(
            format!("{}", grid),
            "0: [0, 1]\n\
             1: [1]\n\
             large_cells = []\n\
             ncell = 2\n\
             ncontainer = 2\n\
             ndiv = [1, 2]\n"
        );
        let mut plot = Plot::new();
        grid.draw(&mut plot, true).unwrap();
        // if true {
        //     plot.set_equal_axes(true)
        //         .set_figure_size_points(600.0, 600.0)
        //         .grid_and_labels("x", "y")
        //         .set_ticks_x(0.2, 0.0, "")
        //         .set_ticks_y(0.2, 0.0, "")
        //         .save("/tmp/gemlab/test_grid_cells_new_works_2.svg")
        //         .unwrap();
        // }
    }

    #[test]
    fn new_works_3() {
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 3,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0, 0.0] },
                Point { id: 2, coords: vec![0.0, 1.0, 0.0] },
                Point { id: 3, coords: vec![0.0, 0.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Tet4, points: vec![0, 1, 2, 3] },
            ],
        };
        let tolerance = 1e-3;
        let border_tol = 0.1;
        let grid = GridCells::new(&mesh, Some(tolerance), Some(border_tol), None, None).unwrap();
        let max_len = 1.0;
        let sl = max_len + 2.0 * tolerance; // because the bbox is expanded
        assert_eq!(grid.side_length, sl);
        assert_eq!(grid.ndiv, &[2, 2, 2]);
        assert_eq!(grid.xmin, &[-0.1, -0.1, -0.1]);
        assert_eq!(grid.xmax, &[-0.1 + sl * 2.0, -0.1 + sl * 2.0, -0.1 + sl * 2.0]);
        assert_eq!(grid.bounding_boxes.len(), 1);
        assert_eq!(grid.containers.len(), 8);
        let bbox_0 = &grid.bounding_boxes[0];
        assert_eq!(bbox_0, &[[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]);
        assert_eq!(
            format!("{}", grid),
            "0: [0]\n\
             1: [0]\n\
             2: [0]\n\
             3: [0]\n\
             4: [0]\n\
             5: [0]\n\
             6: [0]\n\
             7: [0]\n\
             large_cells = []\n\
             ncell = 1\n\
             ncontainer = 8\n\
             ndiv = [2, 2, 2]\n"
        );
        let mut plot = Plot::new();
        grid.draw(&mut plot, true).unwrap();
        // if true {
        //     plot.set_equal_axes(true)
        //         .set_figure_size_points(600.0, 600.0)
        //         .save("/tmp/gemlab/test_grid_cells_new_works_3.svg")
        //         .unwrap();
        // }
    }

    #[test]
    fn new_works_large_cells() {
        let mesh = Samples::tri3_from_delaunay();
        let grid = GridCells::new(&mesh, None, None, None, None).unwrap();
        let max_len = mesh.points[1].coords[1] - mesh.points[0].coords[1];
        let sl = max_len + 2.0 * GRID_CELLS_TOLERANCE; // because the bbox is expanded
        let g = GRID_CELLS_BORDER_TOL;
        assert_eq!(grid.side_length, sl);
        assert_eq!(grid.ndiv, &[2, 2]);
        assert_eq!(grid.xmin, &[mesh.points[0].coords[0] - g, mesh.points[5].coords[1] - g]);
        assert_eq!(
            grid.xmax,
            &[
                mesh.points[0].coords[0] - g + sl * 2.0,
                mesh.points[5].coords[1] - g + sl * 2.0
            ]
        );
        assert_eq!(grid.bounding_boxes.len(), mesh.cells.len());
        assert_eq!(grid.containers.len(), 4);
        let bbox_0 = &grid.bounding_boxes[0];
        assert_eq!(
            bbox_0,
            &[
                [mesh.points[2].coords[0], mesh.points[6].coords[0]],
                [mesh.points[6].coords[1], mesh.points[4].coords[1]]
            ]
        );
        assert_eq!(
            format!("{}", grid),
            "0: [0, 1, 2, 3, 5, 6, 7, 9, 11]\n\
             1: [6, 7, 9]\n\
             2: [0, 2, 3, 4, 6, 9]\n\
             3: [4, 6, 9]\n\
             large_cells = [8, 10]\n\
             ncell = 12\n\
             ncontainer = 4\n\
             ndiv = [2, 2]\n"
        );
        let mut plot = Plot::new();
        grid.draw(&mut plot, true).unwrap();
        // if true {
        //     plot.set_equal_axes(true)
        //         .grid_and_labels("x", "y")
        //         .set_ticks_x(0.1, 0.0, "")
        //         .set_ticks_y(0.1, 0.0, "")
        //         .set_figure_size_points(600.0, 600.0)
        //         .save("/tmp/gemlab/test_grid_cells_new_works_large_cells.svg")
        //         .unwrap();
        // }
    }
}
