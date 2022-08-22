#![allow(unused)]

use super::{calc_container_key, AsArray2D, GS_DEFAULT_BORDER_TOL, GS_DEFAULT_TOLERANCE};
use crate::StrError;
use plotpy::{Canvas, Plot, PolyCode, Text};
use russell_stat::Histogram;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Defines a bounding-box ratio coefficient in [0,1] to mark cells as "large" cells compared the others
///
/// Example: If a cell has a max(bounding box length) / largest >= 0.75,
///          Then this cell is put into the large_cells list
///
/// This constant enables the strategy of putting some cells in a selected list of large cells.
/// The other constant [GS_LARGE_CELLS_PCT] also influences the algorithm by enabling or disabling this strategy.
const GS_LARGE_CELLS_CUTOFF: f64 = 0.75;

/// Controls the count percentage in [0,1] of large cells allowed in the selected array of large cells.
///
/// Example: If 20% of cells are too large, put them in a separated list of large cells;
///          otherwise, the mesh is homogeneous enough so put all cells in the containers.
///
/// If the actual ratio is greater than this constant (i.e., many cells are large such as in a homogeneous mesh),
/// the "strategy" of selecting large cells is abandoned. The other constant [GS_LARGE_CELLS_CUTOFF] also influences the algorithm.
const GS_LARGE_CELLS_PCT: f64 = 0.2;

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

pub struct GridSearchCellNew {
    ndim: usize,                     // space dimension
    ndiv: Vec<usize>,                // (ndim) number of divisions along each direction
    xmin: Vec<f64>,                  // (ndim) min values
    xmax: Vec<f64>,                  // (ndim) max values
    side_length: f64,                // side length of a container
    bounding_boxes: Vec<BboxMinMax>, // (ncell) bounding boxes
    containers: Containers,          // structure to hold all items
    large_cells: Container,          // holds the CellId of large cells
}

impl GridSearchCellNew {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `ndim` -- is the space dimension (2 or 3)
    /// * `points` -- (npoint) is a list of coordinates such as `[[x0,y0,z0,...], [x1,y1,z1,...], [x2,y2,z2,...], [x3,y3,z3,...]]`
    /// * `cells` -- (ncell) is a list of connectivity (topology) such as `[[0,2,1,3], [2,1,0,3]]`
    /// * `tolerance` -- is a tolerance to expand the bounding box of triangles and compare points; e.g. 1e-4
    ///     - If None, [GS_DEFAULT_TOLERANCE] is used
    /// * `border_tol` -- is a tolerance used to expand the border a little bit and then
    ///   accommodate eventual imprecision near the borders; e.g. 1e-2
    ///     - If None, [GS_DEFAULT_BORDER_TOL] is used
    /// * `large_cells_cutoff` -- is the ratio to define large cells compared to the largest bounding box length; e.g. 0.75
    ///     - If None, [GS_LARGE_CELLS_CUTOFF] is used
    /// * `large_cells_pct` -- is the ratio (percentage) of allowed large cells; e.g. 0.2
    ///     - If None, [GS_LARGE_CELLS_PCT] is used
    ///
    /// # Panics
    ///
    /// If the list of points or cells contain more entries at the first row, a panic will occur.
    /// The first row is used to compute the number of points/cells.
    pub fn new<'a, P, C>(
        ndim: usize,
        points: &'a P,
        cells: &'a C,
        tolerance: Option<f64>,
        border_tol: Option<f64>,
        large_cells_cutoff: Option<f64>,
        large_cells_pct: Option<f64>,
    ) -> Result<Self, StrError>
    where
        P: AsArray2D<'a, f64>,
        C: AsArray2D<'a, usize>,
    {
        // check input
        if ndim < 2 || ndim > 3 {
            return Err("ndim must be 2 or 3");
        }
        let (npoint, ncol) = points.size();
        let (ncell, nnode) = cells.size();
        if npoint < 3 {
            return Err("number of points must be ≥ 3");
        }
        if ncol < ndim {
            return Err("points.ncol must be ≥ ndim");
        }
        if nnode < ndim + 1 {
            return Err("number of cell nodes (nnode) must be ≥ ndim + 1");
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

        // large cells cutoff
        let large_cells_cutoff = match large_cells_cutoff {
            Some(v) => v,
            None => GS_LARGE_CELLS_CUTOFF,
        };
        if large_cells_cutoff < 0.0 || large_cells_cutoff > 1.0 {
            return Err("large_cells_cutoff must be in [0.0, 1.0]");
        }

        // large cells pct
        let large_cells_pct = match large_cells_pct {
            Some(v) => v,
            None => GS_LARGE_CELLS_PCT,
        };
        if large_cells_pct < 0.0 || large_cells_pct > 1.0 {
            return Err("large_cells_pct must be in [0.0, 1.0]");
        }

        // find limits, bounding boxes, and largest bounding box
        let mut x = vec![0.0; ndim];
        let mut xmin = vec![f64::MAX; ndim];
        let mut xmax = vec![f64::MIN; ndim];
        let mut x_min_max = vec![vec![0.0; N_MIN_MAX]; ndim];
        let mut bbox_largest = vec![f64::MIN; ndim];
        let mut bounding_boxes = Vec::new();
        for cell_id in 0..ncell {
            for m in 0..nnode {
                let point_id = cells.at(cell_id, m);
                for i in 0..ndim {
                    // coords
                    x[i] = points.at(point_id, i);
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
        Ok(GridSearchCellNew {
            ndim,
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
    pub fn find_cell<F>(&self, x: &[f64], in_cell: Option<F>) -> Result<Option<CellId>, StrError>
    where
        F: Fn(usize, &[f64]) -> Result<bool, StrError>,
    {
        Err("")
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::GridSearchCellNew;
    use crate::geometry::{in_triangle, triangle_coords};
    use crate::util::{GS_DEFAULT_BORDER_TOL, GS_DEFAULT_TOLERANCE};
    use plotpy::Plot;

    #[allow(unused_imports)]
    use plotpy::{Canvas, PolyCode};

    // DO NOT DELETE the lines below
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

    // DO NOT DELETE the lines below
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
    fn new_handles_errors() {
        let points = vec![vec![0.0]];
        let cells = vec![vec![]];
        let x = vec![0.0, 0.0, 0.0];
        assert_eq!(
            GridSearchCellNew::new(1, &points, &cells, None, None, None, None).err(),
            Some("ndim must be 2 or 3")
        );
        assert_eq!(
            GridSearchCellNew::new(2, &points, &cells, None, None, None, None).err(),
            Some("number of points must be ≥ 3")
        );
        let points = vec![vec![0.0, 0.0], vec![1.0, 1.0], vec![2.0, 2.0]];
        assert_eq!(
            GridSearchCellNew::new(3, &points, &cells, None, None, None, None).err(),
            Some("points.ncol must be ≥ ndim")
        );
        let points = vec![vec![0.0], vec![1.0], vec![2.0]];
        assert_eq!(
            GridSearchCellNew::new(2, &points, &cells, None, None, None, None).err(),
            Some("points.ncol must be ≥ ndim")
        );
        let points = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0]];
        let cells = vec![vec![1, 2], vec![0, 1, 2], vec![0, 1, 2]];
        assert_eq!(
            GridSearchCellNew::new(2, &points, &cells, None, None, None, None).err(),
            Some("number of cell nodes (nnode) must be ≥ ndim + 1")
        );
        let points = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0]];
        let cells = vec![vec![0, 1, 2], vec![0, 1, 2], vec![0, 1, 2]];
        assert_eq!(
            GridSearchCellNew::new(2, &points, &cells, Some(-1.0), None, None, None).err(),
            Some("tolerance must be > 0.0")
        );
        assert_eq!(
            GridSearchCellNew::new(2, &points, &cells, None, Some(-1.0), None, None).err(),
            Some("border_tol must be ≥ 0.0")
        );
        assert_eq!(
            GridSearchCellNew::new(2, &points, &cells, None, None, Some(-1.0), None).err(),
            Some("large_cells_cutoff must be in [0.0, 1.0]")
        );
        assert_eq!(
            GridSearchCellNew::new(2, &points, &cells, None, None, Some(2.0), None).err(),
            Some("large_cells_cutoff must be in [0.0, 1.0]")
        );
        assert_eq!(
            GridSearchCellNew::new(2, &points, &cells, None, None, None, Some(-1.0)).err(),
            Some("large_cells_pct must be in [0.0, 1.0]")
        );
        assert_eq!(
            GridSearchCellNew::new(2, &points, &cells, None, None, None, Some(2.0)).err(),
            Some("large_cells_pct must be in [0.0, 1.0]")
        );
    }
}
