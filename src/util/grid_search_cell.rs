use super::{GS_DEFAULT_BORDER_TOL, GS_DEFAULT_TOLERANCE};
use crate::StrError;
use plotpy::{Canvas, Plot, PolyCode, Text};
use std::collections::{HashMap, HashSet};
use std::fmt;

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

/// Defines a tool to search the cell where a point is located within a mesh
///
/// # Concept
///
/// A grid is composed of containers. The container (aka "bin") is uniquely identified by its "left-most" corner (aka **key**).
/// Any coordinate in 2D or 3D can be quickly compared with the "left-most" corner; therefore, we can find coordinates efficiently.
///
/// The main idea revolves around the following expression that calculates the index of a container
/// where the coordinates `x` fall in (Eq. 8 of the Reference):
///
/// ```text
/// ratio[i] = truncate((x[i] - xmin[i]) / Δx[i])
/// key = ratio[0] + ratio[1] * ndiv[0]
/// ```
///
/// Below is an illustration of a grid with two triangles:
///
/// ![test_grid_search_cell_new_2](https://github.com/cpmech/gemlab/raw/main/data/figures/test_grid_search_cell_new_2.svg)
///
/// Only the containers touched by the bounding box of the triangle are saved. Thus, we use a `HashMap` connecting the container
/// key to a `HashSet` of cells belonging to this container. Therefore, all data is stored in the following data structure;
///
/// ```text
/// type Container = HashSet<CellId>;
/// type Containers = HashMap<ContainerKey, Container>;
/// ```
///
/// This data structure avoids repetition, saves space, and is somewhat efficient.
///
/// # Limitations
///
/// The `GridSearchCell` only works if the minimum container size (edge/side length) is greater than
/// the maximum dimension of the largest triangle. Therefore, if one triangle is much larger that the other ones,
/// the algorithm won't perform as well as it could. On the other hand, if the triangles have similar "sizes," then the search should be fast.
///
/// In summary, we need to make sure that:
///
/// * The container's `side_length` must be greater than the maximum dimension of the largest triangle
///
/// # Examples
///
/// ```
/// use gemlab::geometry::is_point_inside_triangle;
/// use gemlab::util::GridSearchCell;
/// use gemlab::StrError;
///
/// fn main() -> Result<(), StrError> {
///     // [num_triangle][nnode=3][ndim=2]
///     #[rustfmt::skip]
///     const TRIS: [[[f64; 2]; 3]; 8] = [
///         [[0.0, 0.0],  [1.0, 0.0],  [0.5, 0.85]],
///         [[1.0, 0.0],  [2.0, 0.0],  [1.5, 0.85]],
///         [[0.5, 0.85], [1.0, 0.0],  [1.5, 0.85]],
///         [[1.5, 0.85], [2.0, 0.0],  [2.5, 0.85]],
///         [[0.5, 0.85], [1.5, 0.85], [1.0, 1.7]],
///         [[1.5, 0.85], [2.5, 0.85], [2.0, 1.7]],
///         [[1.0, 1.7],  [1.5, 0.85], [2.0, 1.7]],
///         [[2.0, 1.7],  [2.5, 0.85], [3.0, 1.7]],
///     ];
///
///     // closure that returns the number of nodes of a cell `t`
///     let get_nnode = |_t| Ok(3);
///
///     // closure that returns the coordinates of point `m` of cell `t`
///     let get_x = |t: usize, m: usize| Ok(&TRIS[t][m][..]);
///
///     // allocate grid search tool
///     let ndim = 2;
///     let grid = GridSearchCell::new(ndim, TRIS.len(), get_nnode, get_x, None, None)?;
///
///     // closure that tells whether the point is in the cell or not
///     let is_in_cell = |t: usize, x: &[f64]| Ok(is_point_inside_triangle(&TRIS[t][0], &TRIS[t][1], &TRIS[t][2], x));
///
///     // find triangle given coords
///     assert_eq!(grid.find_cell(&[1.0, 0.5], is_in_cell)?, Some(2));
///     assert_eq!(grid.find_cell(&[2.9, 1.6], is_in_cell)?, Some(7));
///     assert_eq!(grid.find_cell(&[3.0, 1.0], is_in_cell)?, None);
///     Ok(())
/// }
/// ```
///
/// ![example_grid_search_triangles](https://github.com/cpmech/gemlab/raw/main/data/figures/example_grid_search_triangles.svg)
pub struct GridSearchCell {
    ndim: usize,                     // space dimension
    ndiv: Vec<usize>,                // (ndim) number of divisions along each direction
    xmin: Vec<f64>,                  // (ndim) min values
    xmax: Vec<f64>,                  // (ndim) max values
    side_length: f64,                // side length of a container
    bounding_boxes: Vec<BboxMinMax>, // (ncell) bounding boxes
    containers: Containers,          // structure to hold all items
}

impl GridSearchCell {
    /// Calculates the key of the container where the point should fall in
    #[inline]
    fn calc_container_key(ndim: usize, side_length: f64, ndiv: &[usize], xmin: &[f64], x: &[f64]) -> ContainerKey {
        let mut ix = ((x[0] - xmin[0]) / side_length) as usize; // (Eq. 8)
        let mut iy = ((x[1] - xmin[1]) / side_length) as usize;
        if ix == ndiv[0] {
            ix -= 1; // point is on max edge => move to inner container
        }
        if iy == ndiv[1] {
            iy -= 1; // point is on max edge => move to inner container
        }
        if ndim == 2 {
            return ix + iy * ndiv[0];
        }
        let mut iz = ((x[2] - xmin[2]) / side_length) as usize;
        if iz == ndiv[2] {
            iz -= 1; // point is on max edge => move to inner container
        }
        return ix + iy * ndiv[0] + iz * ndiv[0] * ndiv[1];
    }

    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `ndim` -- is the space dimension (2 or 3)
    /// * `ncell` -- is the number of cells (e.g., triangle/tetrahedron) in the mesh.
    ///     - All cells are numbered from `0` to `ncell - 1`
    ///     - The index of the cell in a mesh is also called CellId (`cell_id`)
    /// * `get_nnode` -- is a function of the `cell_id` that returns the number of nodes `nnode` of the cell
    /// * `get_x` -- is a function of the `cell_id` and the local index of the node/point `m`.
    ///    This function returns the coordinates `x` of the point.
    /// * `tolerance` -- is a tolerance to expand the bounding box of cells and compare points; e.g. 1e-4
    ///     - If None, [GS_DEFAULT_TOLERANCE] is used
    /// * `border_tol` -- is a tolerance used to expand the border a little bit and then
    ///   accommodate eventual imprecision near the borders; e.g. 1e-2
    ///     - If None, [GS_DEFAULT_BORDER_TOL] is used
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
        let mut bounding_boxes = Vec::new();
        for cell_id in 0..ncell {
            let nnode = get_nnode(cell_id)?;
            for m in 0..nnode {
                let x = get_x(cell_id, m)?;
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
            bounding_boxes.push(x_min_max.clone());
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

        // insert cells
        let mut containers = HashMap::new();
        let mut x = vec![0.0; ndim];
        for cell_id in 0..ncell {
            let x_min_max = &bounding_boxes[cell_id];
            for r in 0..2 {
                for s in 0..2 {
                    for t in 0..(ndim - 1) {
                        x[0] = x_min_max[0][r];
                        x[1] = x_min_max[1][s];
                        if ndim == 3 {
                            x[2] = x_min_max[2][t];
                        }
                        let key = GridSearchCell::calc_container_key(ndim, side_length, &ndiv, &xmin, &x);
                        let container = containers.entry(key).or_insert(HashSet::new());
                        container.insert(cell_id);
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
            side_length,
            bounding_boxes,
            containers,
        })
    }

    /// Find the cell (e.g., triangle or tetrahedron) where the given coordinate falls in
    ///
    /// * Returns the CellId or None if no cell contains the point
    /// * `is_in_cell` is a function of cell_id and x that tells whether the point si in the cell or not
    pub fn find_cell<F>(&self, x: &[f64], is_in_cell: F) -> Result<Option<CellId>, StrError>
    where
        F: Fn(usize, &[f64]) -> Result<bool, StrError>,
    {
        // check if the point is out-of-bounds
        for i in 0..self.ndim {
            if x[i] < self.xmin[i] || x[i] > self.xmax[i] {
                return Err("given point coordinates are outside the grid");
            }
        }

        // get the container where `x` falls in
        let key = GridSearchCell::calc_container_key(self.ndim, self.side_length, &self.ndiv, &self.xmin, x);
        let container = match self.containers.get(&key) {
            Some(c) => c,
            None => return Ok(None), // no container with cells in it
        };

        // find the cell where the point falls in
        for cell_id in container {
            let x_min_max = &self.bounding_boxes[*cell_id];
            for i in 0..self.ndim {
                if x[i] < x_min_max[i][I_MIN] || x[i] > x_min_max[i][I_MAX] {
                    continue; // outside the bounding box
                }
            }
            if (is_in_cell)(*cell_id, x)? {
                return Ok(Some(*cell_id));
            }
        }

        // not found
        Ok(None)
    }

    /// Draws grid and items
    pub fn draw(&self, plot: &mut Plot, with_ids: bool) -> Result<(), StrError> {
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
        for x_min_max in &self.bounding_boxes {
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
            for cell_id in container {
                let x_min_max = &self.bounding_boxes[*cell_id];
                let txt = format!("{}", cell_id);
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
        plot.add(&bbox);
        if with_ids {
            plot.add(&text);
        }
        Ok(())
    }

    /// Returns the grid limits (xmin,xmax)
    pub fn limits(&self) -> (Vec<f64>, Vec<f64>) {
        (self.xmin.clone(), self.xmax.clone())
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
        write!(f, "ncell = {}\n", unique_items.len()).unwrap();
        write!(f, "ncontainer = {}\n", self.containers.len()).unwrap();
        write!(f, "ndiv = {:?}\n", self.ndiv).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::GridSearchCell;
    use crate::geometry::is_point_inside_triangle;
    use crate::util::{GS_DEFAULT_BORDER_TOL, GS_DEFAULT_TOLERANCE};
    use crate::StrError;
    use plotpy::{Canvas, Plot, PolyCode};

    #[no_coverage]
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

    #[no_coverage]
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
    fn calc_container_key_works() {
        // 2d
        let ndim = 2;
        let side_length = 0.30000000000000004;
        let ndiv = &[4, 8];
        let xmin = &[-0.30000000000000004, -0.30000000000000004];
        let x = &[0.1, 0.5];
        assert_eq!(GridSearchCell::calc_container_key(ndim, side_length, ndiv, xmin, x), 9);
        let x = &[0.7, 0.8];
        assert_eq!(GridSearchCell::calc_container_key(ndim, side_length, ndiv, xmin, x), 15);
        let x = &[-0.2, 1.8];
        assert_eq!(GridSearchCell::calc_container_key(ndim, side_length, ndiv, xmin, x), 24);
        let x = &[0.8, 1.8];
        assert_eq!(GridSearchCell::calc_container_key(ndim, side_length, ndiv, xmin, x), 27);
        let xmax = &[
            xmin[0] + (ndiv[0] as f64) * side_length,
            xmin[1] + (ndiv[1] as f64) * side_length,
        ];
        let x = &[xmax[0], xmax[1]];
        assert_eq!(
            GridSearchCell::calc_container_key(ndim, side_length, ndiv, xmin, x),
            (ndiv[0] * ndiv[1] - 1)
        );

        // 3d
        let ndim = 3;
        let side_length = 1.1;
        let ndiv = &[2, 2, 2];
        let xmin = &[-1.1, -1.1, -1.1];
        let x = &[-1.0, -1.0, -1.0];
        assert_eq!(GridSearchCell::calc_container_key(ndim, side_length, ndiv, xmin, x), 0);
        let x = &[1.0, 1.0, 1.0];
        assert_eq!(GridSearchCell::calc_container_key(ndim, side_length, ndiv, xmin, x), 7);
        let xmax = &[
            xmin[0] + (ndiv[0] as f64) * side_length,
            xmin[1] + (ndiv[1] as f64) * side_length,
            xmin[2] + (ndiv[2] as f64) * side_length,
        ];
        let x = &[xmax[0], xmax[1], xmax[2]];
        assert_eq!(
            GridSearchCell::calc_container_key(ndim, side_length, ndiv, xmin, x),
            (ndiv[0] * ndiv[1] * ndiv[2] - 1)
        );
    }

    #[test]
    fn new_handles_errors() {
        let x = vec![0.0, 0.0, 0.0];
        assert_eq!(
            GridSearchCell::new(1, 1, |_| Ok(3), |_, _| Ok(&x), None, None).err(),
            Some("ndim must be 2 or 3")
        );
        assert_eq!(
            GridSearchCell::new(2, 1, |_| Ok(3), |_, _| Ok(&x), Some(-1.0), None).err(),
            Some("tolerance must be > 0.0")
        );
        assert_eq!(
            GridSearchCell::new(2, 1, |_| Ok(3), |_, _| Ok(&x), None, Some(-1.0)).err(),
            Some("border_tol must be ≥ 0.0")
        );
        assert_eq!(
            GridSearchCell::new(2, 1, |_| Ok(3), |_, _| Ok(&x), None, None).err(),
            Some("x.len() must be equal to ndim")
        );
        let x = vec![0.0, 0.0];
        let get_nnode = |_| Err("get_nnode returns error");
        assert_eq!(
            GridSearchCell::new(2, 1, get_nnode, |_, _| Ok(&x), None, None).err(),
            Some("get_nnode returns error")
        );
        let get_x = |_: usize, _: usize| Err("get_x returns error");
        assert_eq!(
            GridSearchCell::new(2, 1, |_| Ok(3), get_x, None, None).err(),
            Some("get_x returns error")
        );
    }

    #[test]
    fn new_works_1() -> Result<(), StrError> {
        const TRIS: [[[f64; 2]; 3]; 2] = [
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        ];
        let tolerance = 1e-3;
        let border_tol = 0.1;
        let get_nnode = |_| Ok(3);
        let get_x = |t: usize, m: usize| Ok(&TRIS[t][m][..]);
        let grid = GridSearchCell::new(2, TRIS.len(), get_nnode, get_x, Some(tolerance), Some(border_tol))?;
        let max_len = 1.0;
        let sl = max_len + 2.0 * tolerance; // because the bbox is expanded
        assert_eq!(grid.ndim, 2);
        assert_eq!(grid.side_length, sl);
        assert_eq!(grid.ndiv, &[2, 2]);
        assert_eq!(grid.xmin, &[-0.1, -0.1]);
        assert_eq!(grid.xmax, &[-0.1 + sl * 2.0, -0.1 + sl * 2.0]);
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
             ncell = 2\n\
             ncontainer = 4\n\
             ndiv = [2, 2]\n"
        );
        if false {
            let mut plot = Plot::new();
            draw_triangles(&mut plot, &TRIS);
            grid.draw(&mut plot, true)?;
            plot.set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .grid_and_labels("x", "y")
                .set_ticks_x(0.2, 0.0, "")
                .set_ticks_y(0.2, 0.0, "")
                .save("/tmp/gemlab/test_grid_search_cell_new_1.svg")?;
        }
        // with zero border
        let grid = GridSearchCell::new(2, TRIS.len(), get_nnode, get_x, None, Some(0.0))?;
        let sl = max_len + 2.0 * GS_DEFAULT_TOLERANCE;
        assert_eq!(grid.xmin, &[0.0, 0.0]);
        assert_eq!(grid.ndiv, &[1, 1]);
        assert_eq!(grid.xmax, &[sl, sl]);
        Ok(())
    }

    #[test]
    fn draw_works() -> Result<(), StrError> {
        const TRIS: [[[f64; 2]; 3]; 2] = [
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            [[1.0, 0.0], [1.2, 1.5], [0.0, 1.0]],
        ];
        let get_x = |t: usize, m: usize| Ok(&TRIS[t][m][..]);
        let grid = GridSearchCell::new(2, 1, |_| Ok(3), get_x, None, None)?;
        let mut plot = Plot::new();
        grid.draw(&mut plot, false)?;

        const TETS: [[[f64; 3]; 4]; 1] = [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]];
        let get_x = |t: usize, m: usize| Ok(&TETS[t][m][..]);
        let grid = GridSearchCell::new(3, 1, |_| Ok(4), get_x, None, None)?;
        let mut plot = Plot::new();
        grid.draw(&mut plot, true)?;
        Ok(())
    }

    #[test]
    fn new_works_2() -> Result<(), StrError> {
        const TRIS: [[[f64; 2]; 3]; 2] = [
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            [[1.0, 0.0], [1.2, 1.5], [0.0, 1.0]],
        ];
        let tolerance = 1e-3;
        let border_tol = 0.1;
        let get_nnode = |_| Ok(3);
        let get_x = |t: usize, m: usize| Ok(&TRIS[t][m][..]);
        let grid = GridSearchCell::new(2, TRIS.len(), get_nnode, get_x, Some(tolerance), Some(border_tol))?;
        let max_len = 1.5;
        let sl = max_len + 2.0 * tolerance; // because the bbox is expanded
        assert_eq!(grid.ndim, 2);
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
             ncell = 2\n\
             ncontainer = 2\n\
             ndiv = [1, 2]\n"
        );
        if false {
            let mut plot = Plot::new();
            draw_triangles(&mut plot, &TRIS);
            grid.draw(&mut plot, true)?;
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
        let grid = GridSearchCell::new(3, TETS.len(), get_nnode, get_x, Some(tolerance), Some(border_tol))?;
        let max_len = 1.0;
        let sl = max_len + 2.0 * tolerance; // because the bbox is expanded
        assert_eq!(grid.ndim, 3);
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
             ncell = 1\n\
             ncontainer = 8\n\
             ndiv = [2, 2, 2]\n"
        );
        if false {
            let mut plot = Plot::new();
            draw_tetrahedra(&mut plot, &TETS);
            grid.draw(&mut plot, true)?;
            plot.set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .save("/tmp/gemlab/test_grid_search_cell_new_3.svg")?;
        }
        Ok(())
    }

    #[test]
    fn find_cell_handles_errors() -> Result<(), StrError> {
        let x = vec![0.0, 0.0];
        let grid = GridSearchCell::new(2, 1, |_| Ok(3), |_, _| Ok(&x), None, None)?;
        let y = vec![10.0, 0.0];
        assert_eq!(
            grid.find_cell(&y, |_, _| Ok(true)).err(),
            Some("given point coordinates are outside the grid")
        );
        Ok(())
    }

    #[test]
    fn find_cell_works_2d() -> Result<(), StrError> {
        // [num_triangle][nnode=3][ndim=2]
        #[rustfmt::skip]
        const TRIS: [[[f64; 2]; 3]; 12] = [
            [[0.230951,  0.558482], [0.133721,  0.348832],   [0.540745,  0.331184]],   //  0
            [[0.13928,   0.180603], [0.133721,  0.348832],   [0.0307942, 0.459123]],   //  1
            [[0.0307942, 0.459123], [0.230951,  0.558482],   [0.0980015, 0.981755]],   //  2
            [[0.230951,  0.558482], [0.0307942, 0.459123],   [0.133721,  0.348832]],   //  3
            [[0.0980015, 0.981755], [0.230951,  0.558482],   [0.578587,  0.760349]],   //  4
            [[0.133721,  0.348832], [0.13928,   0.180603],   [0.540745,  0.331184]],   //  5
            [[0.540745,  0.331184], [0.578587,  0.760349],   [0.230951,  0.558482]],   //  6
            [[0.540745,  0.331184], [0.478554,  0.00869692], [0.648071,  0.369534]],   //  7
            [[0.578587,  0.760349], [0.648071,  0.369534],   [0.903726,  0.975904]],   //  8
            [[0.648071,  0.369534], [0.578587,  0.760349],   [0.540745,  0.331184]],   //  9
            [[0.578587,  0.760349], [0.903726,  0.975904],   [0.0980015, 0.981755]],   // 10
            [[0.540745,  0.331184], [0.13928,   0.180603],   [0.478554,  0.00869692]], // 11
        ];
        let get_nnode = |_| Ok(3);
        let get_x = |t: usize, m: usize| Ok(&TRIS[t][m][..]);
        let grid = GridSearchCell::new(2, TRIS.len(), get_nnode, get_x, None, None)?;
        let max_len = TRIS[10][1][0] - TRIS[10][2][0];
        let sl = max_len + 2.0 * GS_DEFAULT_TOLERANCE; // because the bbox is expanded
        let g = GS_DEFAULT_BORDER_TOL;
        assert_eq!(grid.ndim, 2);
        assert_eq!(grid.side_length, sl);
        assert_eq!(grid.ndiv, &[2, 2]);
        assert_eq!(grid.xmin, &[TRIS[1][2][0] - g, TRIS[11][2][1] - g]);
        assert_eq!(
            grid.xmax,
            &[TRIS[1][2][0] - g + sl * 2.0, TRIS[11][2][1] - g + sl * 2.0]
        );
        assert_eq!(grid.bounding_boxes.len(), TRIS.len());
        assert_eq!(grid.containers.len(), 4);
        let bbox_0 = &grid.bounding_boxes[0];
        assert_eq!(
            bbox_0,
            &[[TRIS[0][1][0], TRIS[0][2][0]], [TRIS[0][2][1], TRIS[0][0][1]]]
        );
        assert_eq!(
            format!("{}", grid),
            "0: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]\n\
             1: [8, 10]\n\
             2: [2, 4, 8, 10]\n\
             3: [8, 10]\n\
             ncell = 12\n\
             ncontainer = 4\n\
             ndiv = [2, 2]\n"
        );
        let is_in_cell = |t: usize, x: &[f64]| Ok(is_point_inside_triangle(&TRIS[t][0], &TRIS[t][1], &TRIS[t][2], x));
        assert_eq!(grid.find_cell(&[0.4, 0.2], is_in_cell)?, Some(11));
        assert_eq!(grid.find_cell(&[0.6, 0.3], is_in_cell)?, Some(7));
        assert_eq!(grid.find_cell(&[0.1, 0.7], is_in_cell)?, Some(2));
        assert_eq!(grid.find_cell(&[0.8, 0.8], is_in_cell)?, Some(8));
        let res = grid.find_cell(&TRIS[7][1], is_in_cell)?;
        if res != Some(7) {
            assert_eq!(res, Some(11));
        }
        assert_eq!(grid.find_cell(&[0.1, 0.1], is_in_cell)?, None);
        assert_eq!(grid.find_cell(&[0.6, 0.2], is_in_cell)?, None);
        assert_eq!(grid.find_cell(&[0.4, 1.0], is_in_cell)?, None);
        assert_eq!(
            grid.find_cell(&[10.0, 1.0], is_in_cell).err(),
            Some("given point coordinates are outside the grid")
        );
        if false {
            let mut plot = Plot::new();
            draw_triangles(&mut plot, &TRIS);
            grid.draw(&mut plot, true)?;
            plot.set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .grid_and_labels("x", "y")
                .set_ticks_x(0.2, 0.0, "")
                .set_ticks_y(0.2, 0.0, "")
                .save("/tmp/gemlab/test_grid_search_cell_find_works_2d.svg")?;
        }
        Ok(())
    }
}
