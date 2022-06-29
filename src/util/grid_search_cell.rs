#![allow(unused)]

use super::{GS_DEFAULT_BORDER_TOL, GS_DEFAULT_TOLERANCE, SQRT_2, SQRT_3};
use crate::StrError;
use plotpy::{Canvas, Curve, Plot, Text};
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
    coefficient: Vec<usize>,       // (3) coefficients [1, ndiv[0], ndiv[0]*ndiv[1]] (Eq. 8)
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

        // coefficient
        let coefficient = vec![1, ndiv[0], ndiv[0] * ndiv[1]];

        // tolerance to find points
        let tol_dist = if ndim == 2 {
            SQRT_2 * tolerance
        } else {
            SQRT_3 * tolerance
        };

        Ok(GridSearchCell {
            ndim,
            ndiv,
            xmin,
            xmax,
            bbox_large,
            side_length,
            coefficient,
            tol_dist,
            bounding_boxes,
            containers: HashMap::new(),
        })
    }

    /// Inserts a cell (e.g., triangle/tetrahedron) to the grid
    ///
    /// # Input
    ///
    /// * `id` -- is the id of the cell
    /// * `nnode` -- is the number of nodes of the cell
    /// * `get_x` -- is a function that returns the coordinates of a cell node
    pub fn insert_cell<'a, F>(&mut self, id: usize, nnode: usize, get_x: F) -> Result<(), StrError>
    where
        F: Fn(usize) -> Result<&'a [f64], StrError>,
    {
        // TODO
        Ok(())
    }

    /// Calculates the key of the container where the point should fall in
    ///
    /// **Note:** Returns None if the point is out-of-range
    #[inline]
    fn calc_container_key(&self, x: &[f64]) -> Option<usize> {
        let mut ratio = vec![0; self.ndim]; // ratio = trunc(δx[i]/Δx[i]) (Eq. 8)
        let mut key = 0;
        for i in 0..self.ndim {
            if x[i] < self.xmin[i] || x[i] > self.xmax[i] {
                return None;
            }
            ratio[i] = ((x[i] - self.xmin[i]) / self.side_length) as usize;
            if ratio[i] == self.ndiv[i] {
                // the point is exactly on the max edge, thus select inner container
                ratio[i] -= 1; // move to the inside
            }
            key += ratio[i] * self.coefficient[i];
        }
        Some(key)
    }

    /// Updates a container or inserts a point into an existing container
    #[inline]
    fn update_or_insert(&mut self, key: ContainerKey, id: ItemID) {
        let container = self.containers.entry(key).or_insert(HashSet::new());
        container.insert(id);
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

    #[test]
    fn new_works() -> Result<(), StrError> {
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
        assert_eq!(grid.coefficient, &[1, 2, 2 * 2]);
        assert_eq!(grid.tol_dist, SQRT_2 * tolerance);
        assert_eq!(grid.bounding_boxes.len(), 2);
        assert_eq!(grid.containers.len(), 0);
        let bbox_0 = grid.bounding_boxes.get(&0).unwrap();
        let bbox_1 = grid.bounding_boxes.get(&0).unwrap();
        assert_eq!(bbox_0, &[[0.0, 1.0], [0.0, 1.0]]);
        assert_eq!(bbox_1, &[[0.0, 1.0], [0.0, 1.0]]);
        Ok(())
    }

    #[test]
    fn insert_cell_works_2d() -> Result<(), StrError> {
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
        println!("{:?}", grid.bbox_large);
        /*
        println!("{}", grid);
        let mut id = 0;
        for t in 0..TRIANGLES.len() {
            grid.insert_cell(id, 3, |m| Ok(&TRIANGLES[t][m]))?;
            id += 1;
        }
        assert_eq!(
            format!("{}", grid),
            "0: [100, 101]\n\
             ids = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]\n\
             nitem = 12\n\
             ncontainer = 21\n\
             ndiv = [5, 5]\n"
        );
        if false {
            let mut plot = Plot::new();
            let mut canvas = Canvas::new();
            canvas.set_face_color("#fefddc").set_edge_color("#fcb827");
            for t in 0..TRIANGLES.len() {
                canvas.polycurve_begin();
                for m in 0..3 {
                    let code = if m == 0 { PolyCode::MoveTo } else { PolyCode::LineTo };
                    canvas.polycurve_add(&TRIANGLES[t][m][0], &TRIANGLES[t][m][1], code);
                }
                canvas.polycurve_end(true);
            }
            plot.add(&canvas);
            grid.draw(&mut plot)?;
            plot.set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .grid_and_labels("x", "y")
                .set_ticks_x(0.1, 0.01, "")
                .set_ticks_y(0.1, 0.01, "")
                .save("/tmp/gemlab/test_insert_cell_2d.svg")?;
        }
        */
        Ok(())
    }
}
