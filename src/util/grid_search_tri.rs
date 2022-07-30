use super::{calc_container_key, AsArray2D, GS_DEFAULT_BORDER_TOL, GS_DEFAULT_TOLERANCE};
use crate::geometry::{is_point_inside_triangle, triangle_interpolation};
use crate::StrError;
use plotpy::{Canvas, Plot, PolyCode, Text};
use russell_stat::Histogram;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// 2D
const NDIM: usize = 2;

/// Specifies the key of containers (or bins in the grid)
type ContainerKey = usize;

/// Specifies the identification number of triangles (must be sequential from 0 to ntriangle - 1)
type TriangleId = usize;

/// Defines the container type
type Container = HashSet<TriangleId>;

/// Defines the containers type: Key to Container
type Containers = HashMap<ContainerKey, Container>;

/// Defines the bounding box of a triangle
const N_MIN_MAX: usize = 2; // 2 means {min,max}
const I_MIN: usize = 0;
const I_MAX: usize = 1;
type BboxMinMax = Vec<Vec<f64>>; // [ndim][N_MIN_MAX]

/// Defines a tool to search the triangle where a point is located within a mesh
///
/// **Note::* This is a specialization of [super::GridSearchCell] for the 2D case using triangles.
///
/// See [super::GridSearchCell] for more details.
///
pub struct GridSearchTri {
    ndiv: Vec<usize>,                // (NDIM) number of divisions along each direction
    xmin: Vec<f64>,                  // (NDIM) min values
    xmax: Vec<f64>,                  // (NDIM) max values
    side_length: f64,                // side length of a container
    bounding_boxes: Vec<BboxMinMax>, // (ntriangle) bounding boxes
    containers: Containers,          // structure to hold all items
}

impl GridSearchTri {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `coordinates` -- is a list of coordinates such as `[[x0,y0], [x1,y1], [x2,y2], [x3,y3]]`
    /// * `triangles` -- is a list of connectivity (topology) such as `[[0,2,1], [2,1,0]]`
    /// * `tolerance` -- is a tolerance to expand the bounding box of triangles and compare points; e.g. 1e-4
    ///     - If None, [GS_DEFAULT_TOLERANCE] is used
    /// * `border_tol` -- is a tolerance used to expand the border a little bit and then
    ///   accommodate eventual imprecision near the borders; e.g. 1e-2
    ///     - If None, [GS_DEFAULT_BORDER_TOL] is used
    pub fn new<'a, C, T>(
        coordinates: &'a C,
        triangles: &'a T,
        tolerance: Option<f64>,
        border_tol: Option<f64>,
    ) -> Result<Self, StrError>
    where
        C: AsArray2D<'a, f64>,
        T: AsArray2D<'a, usize>,
    {
        // constants
        let (npoint, ncol) = coordinates.size();
        let (ntriangle, nnode) = triangles.size();
        if npoint < 3 {
            return Err("number of points must be ≥ 3");
        }
        if ncol < 2 {
            return Err("coordinates.ncol must be ≥ 2");
        }
        if nnode != 3 {
            return Err("number of triangle nodes (nnode) must be = 3");
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

        // find limits, bounding boxes, and largest triangle
        let mut xmin = vec![f64::MAX; NDIM];
        let mut xmax = vec![f64::MIN; NDIM];
        let mut x_min_max = vec![vec![0.0; N_MIN_MAX]; NDIM];
        let mut bbox_large = vec![f64::MIN; NDIM];
        let mut bounding_boxes = Vec::new();
        let mut x = vec![0.0; NDIM];
        for cell_id in 0..ntriangle {
            for m in 0..nnode {
                let p = triangles.at(cell_id, m);
                x[0] = coordinates.at(p, 0);
                x[1] = coordinates.at(p, 1);
                for i in 0..NDIM {
                    // limits
                    xmin[i] = f64::min(xmin[i], x[i]);
                    xmax[i] = f64::max(xmax[i], x[i]);
                    // bounding box
                    if m == 0 {
                        for i in 0..NDIM {
                            x_min_max[i][I_MIN] = x[i];
                            x_min_max[i][I_MAX] = x[i];
                        }
                    } else {
                        x_min_max[i][I_MIN] = f64::min(x_min_max[i][I_MIN], x[i]);
                        x_min_max[i][I_MAX] = f64::max(x_min_max[i][I_MAX], x[i]);
                    }
                }
            }
            // largest triangle
            for i in 0..NDIM {
                bbox_large[i] = f64::max(bbox_large[i], x_min_max[i][I_MAX] - x_min_max[i][I_MIN]);
            }
            // add to bounding box maps
            bounding_boxes.push(x_min_max.clone());
        }

        // make the side_length equal to the largest bounding box dimension
        let mut side_length = bbox_large[0];
        for i in 1..NDIM {
            side_length = f64::max(side_length, bbox_large[i]);
        }

        // expand side_length by two times the tolerance (left and right)
        side_length += 2.0 * tolerance;

        // expand borders
        if border_tol > 0.0 {
            for i in 0..NDIM {
                xmin[i] -= border_tol;
                xmax[i] += border_tol;
            }
        }

        // number of divisions
        let mut ndiv = vec![0; NDIM];
        for i in 0..NDIM {
            assert!(xmax[i] > xmin[i]);
            ndiv[i] = f64::ceil((xmax[i] - xmin[i]) / side_length) as usize;
        }

        // update xmax after deciding on the side_length and number of divisions
        for i in 0..NDIM {
            xmax[i] = xmin[i] + side_length * (ndiv[i] as f64);
        }

        // insert triangles
        let mut containers = HashMap::new();
        let mut x = vec![0.0; NDIM];
        for cell_id in 0..ntriangle {
            let x_min_max = &bounding_boxes[cell_id];
            for r in 0..N_MIN_MAX {
                for s in 0..N_MIN_MAX {
                    x[0] = x_min_max[0][r];
                    x[1] = x_min_max[1][s];
                    let key = calc_container_key(NDIM, side_length, &ndiv, &xmin, &x);
                    let container = containers.entry(key).or_insert(HashSet::new());
                    container.insert(cell_id);
                }
            }
        }

        // allocate grid
        Ok(GridSearchTri {
            ndiv,
            xmin,
            xmax,
            side_length,
            bounding_boxes,
            containers,
        })
    }

    /// Finds the triangle where the given coordinate falls in
    ///
    /// # Input
    ///
    /// * `coordinates` -- is a list of coordinates such as `[[x0,y0], [x1,y1], [x2,y2], [x3,y3]]`
    /// * `triangles` -- is a list of connectivity (topology) such as `[[0,2,1], [2,1,0]]`
    ///
    /// # Output
    ///
    /// Returns the index of the triangle in `triangles` or None if no triangle contains the point
    ///
    /// # Warning (Panics)
    ///
    /// The pair `coordinates` and `triangles` must be the same as the ones used in the `new` function,
    /// otherwise **panics** may occur or, even worse, you may get **incorrect results**.
    pub fn find_triangle<'a, C, T>(
        &self,
        x: &[f64],
        coordinates: &'a C,
        triangles: &'a T,
    ) -> Result<Option<TriangleId>, StrError>
    where
        C: AsArray2D<'a, f64>,
        T: AsArray2D<'a, usize>,
    {
        // check if the point is out-of-bounds
        for i in 0..NDIM {
            if x[i] < self.xmin[i] || x[i] > self.xmax[i] {
                return Err("given point coordinates are outside the grid");
            }
        }

        // get the container where `x` falls in
        let key = calc_container_key(NDIM, self.side_length, &self.ndiv, &self.xmin, x);
        let container = match self.containers.get(&key) {
            Some(c) => c,
            None => return Ok(None), // no container with triangles in it
        };

        // find the triangle where the point falls in
        let mut xa = vec![0.0; NDIM];
        let mut xb = vec![0.0; NDIM];
        let mut xc = vec![0.0; NDIM];
        for cell_id in container {
            let x_min_max = &self.bounding_boxes[*cell_id];
            let a = triangles.at(*cell_id, 0);
            let b = triangles.at(*cell_id, 1);
            let c = triangles.at(*cell_id, 2);
            let mut outside = false;
            for i in 0..NDIM {
                if x[i] < x_min_max[i][I_MIN] || x[i] > x_min_max[i][I_MAX] {
                    outside = true; // outside the bounding box
                    break;
                }
                xa[i] = coordinates.at(a, i);
                xb[i] = coordinates.at(b, i);
                xc[i] = coordinates.at(c, i);
            }
            if outside {
                continue;
            }
            if is_point_inside_triangle(&xa, &xb, &xc, x) {
                return Ok(Some(*cell_id));
            }
        }

        // not found
        Ok(None)
    }

    /// Finds the triangle where the given coordinate falls in and performs an interpolation of coordinates
    ///
    /// # Input
    ///
    /// * `coordinates` -- is a list of coordinates such as `[[x0,y0,T0], [x1,y1,T1], [x2,y2,T2], [x3,y3,T3]]`
    ///   where `T[i]` are the values (e.g., temperatures) at that coordinate and will be used for interpolation
    /// * `triangles` -- is a list of connectivity (topology) such as `[[0,2,1], [2,1,0]]`
    ///
    /// # Output
    ///
    /// Returns the value (e.g., temperature) at the target point (xp) inside the triangle.
    /// Returns None if no triangle contains the point.
    ///
    /// # Warning (Panics)
    ///
    /// The pair `coordinates` and `triangles` must be the same as the ones used in the `new` function,
    /// otherwise **panics** may occur or, even worse, you may get **incorrect results**.
    pub fn find_triangle_and_interpolate<'a, C, T>(
        &self,
        x: &[f64],
        coordinates: &'a C,
        triangles: &'a T,
    ) -> Result<Option<f64>, StrError>
    where
        C: AsArray2D<'a, f64>,
        T: AsArray2D<'a, usize>,
    {
        // check if the point is out-of-bounds
        for i in 0..NDIM {
            if x[i] < self.xmin[i] || x[i] > self.xmax[i] {
                return Err("given point coordinates are outside the grid");
            }
        }

        // check if the temperature is present in the coordinates list
        let (_, ncol) = coordinates.size();
        if ncol < 3 {
            return Err("coordinates must contain a third column with the temperature values");
        }

        // get the container where `x` falls in
        let key = calc_container_key(NDIM, self.side_length, &self.ndiv, &self.xmin, x);
        let container = match self.containers.get(&key) {
            Some(c) => c,
            None => return Ok(None), // no container with triangles in it
        };

        // find the triangle where the point falls in
        let mut xa = vec![0.0; NDIM];
        let mut xb = vec![0.0; NDIM];
        let mut xc = vec![0.0; NDIM];
        let mut temp = vec![0.0; 3]; // 3 nodes
        for cell_id in container {
            let x_min_max = &self.bounding_boxes[*cell_id];
            let a = triangles.at(*cell_id, 0);
            let b = triangles.at(*cell_id, 1);
            let c = triangles.at(*cell_id, 2);
            let mut outside = false;
            for i in 0..NDIM {
                if x[i] < x_min_max[i][I_MIN] || x[i] > x_min_max[i][I_MAX] {
                    outside = true; // outside the bounding box
                    break;
                }
                xa[i] = coordinates.at(a, i);
                xb[i] = coordinates.at(b, i);
                xc[i] = coordinates.at(c, i);
            }
            if outside {
                continue;
            }
            if is_point_inside_triangle(&xa, &xb, &xc, x) {
                temp[0] = coordinates.at(a, 2);
                temp[1] = coordinates.at(b, 2);
                temp[2] = coordinates.at(c, 2);
                let res = triangle_interpolation(&xa, &xb, &xc, &temp, x);
                return Ok(Some(res));
            }
        }

        // not found
        Ok(None)
    }

    /// Draws grid and triangles
    ///
    /// # Warning
    ///
    /// The pair `coordinates` and `triangles` must be the same as the ones used in the `new` function,
    /// otherwise **panics** may occur or, even worse, you may get **incorrect results**.
    pub fn draw<'a, C, T>(
        &self,
        plot: &mut Plot,
        coordinates: &'a C,
        triangles: &'a T,
        with_ids: bool,
    ) -> Result<(), StrError>
    where
        C: AsArray2D<'a, f64>,
        T: AsArray2D<'a, usize>,
    {
        // draw triangles
        let mut tri_canvas = Canvas::new();
        tri_canvas.set_face_color("#fefddc").set_edge_color("#fcb827");
        let ntriangle = triangles.size().0;
        for t in 0..ntriangle {
            tri_canvas.polycurve_begin();
            for m in 0..3 {
                let code = if m == 0 { PolyCode::MoveTo } else { PolyCode::LineTo };
                let p = triangles.at(t, m);
                tri_canvas.polycurve_add(coordinates.at(p, 0), coordinates.at(p, 1), code);
            }
            tri_canvas.polycurve_end(true);
        }
        plot.add(&tri_canvas);

        // draw grid
        let mut xmin = vec![0.0; NDIM];
        let mut xmax = vec![0.0; NDIM];
        let mut ndiv = vec![0; NDIM];
        for i in 0..NDIM {
            xmin[i] = self.xmin[i];
            xmax[i] = self.xmax[i];
            ndiv[i] = self.ndiv[i];
        }
        let mut canvas = Canvas::new();
        canvas
            .set_alt_text_color("#5d5d5d")
            .draw_grid(&xmin, &xmax, &ndiv, false, with_ids)?;
        plot.add(&canvas);

        // draw bounding boxes
        let mut bbox = Canvas::new();
        bbox.set_edge_color("#3aff79")
            .set_face_color("None")
            .set_line_width(0.75);
        for x_min_max in &self.bounding_boxes {
            bbox.polycurve_begin();
            bbox.polycurve_add(x_min_max[0][I_MIN], x_min_max[1][I_MIN], PolyCode::MoveTo);
            bbox.polycurve_add(x_min_max[0][I_MAX], x_min_max[1][I_MIN], PolyCode::LineTo);
            bbox.polycurve_add(x_min_max[0][I_MAX], x_min_max[1][I_MAX], PolyCode::LineTo);
            bbox.polycurve_add(x_min_max[0][I_MIN], x_min_max[1][I_MAX], PolyCode::LineTo);
            bbox.polycurve_end(true);
        }
        plot.add(&bbox);

        // draw ids
        if with_ids {
            let mut ids = Text::new();
            ids.set_color("#cd0000");
            let mut xcen = vec![0.0; NDIM];
            for container in self.containers.values() {
                for t in container {
                    let x_min_max = &self.bounding_boxes[*t];
                    let txt = format!("{}", t);
                    for i in 0..NDIM {
                        xcen[i] = (x_min_max[i][I_MIN] + x_min_max[i][I_MAX]) / 2.0;
                    }
                    ids.draw(xcen[0], xcen[1], &txt);
                }
            }
            plot.add(&ids);
        }
        Ok(())
    }

    /// Returns the grid limits (xmin,xmax)
    pub fn limits(&self) -> (Vec<f64>, Vec<f64>) {
        (self.xmin.clone(), self.xmax.clone())
    }

    /// Print some statistics
    pub fn print_stat(&self) {
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
        println!("Summary");
        println!("=======");
        println!("ntriangle = {:?}", unique_items.len());
        println!("xmin = {:?}", self.xmin);
        println!("xmax = {:?}", self.xmax);
        println!("side_length = {:?}", self.side_length);
        println!("num of non-empty containers = {}", self.containers.len());
        println!("max container num items = {}", max_container_num_items);
        println!("\nHistogram of container num items");
        println!("================================");
        let stations: Vec<_> = (0..max_container_num_items + 2).collect();
        let mut hist = Histogram::new(&stations).unwrap();
        hist.count(&container_num_items);
        hist.set_bar_char('*');
        hist.set_bar_max_len(60);
        println!("{}", hist);
    }
}

impl fmt::Display for GridSearchTri {
    /// Shows info about the items in the grid containers
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
        let mut ids: Vec<_> = unique_items.iter().collect();
        ids.sort();
        write!(f, "ntriangle = {}\n", unique_items.len()).unwrap();
        write!(f, "ncontainer = {}\n", self.containers.len()).unwrap();
        write!(f, "ndiv = {:?}\n", self.ndiv).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::GridSearchTri;
    use crate::util::{GS_DEFAULT_BORDER_TOL, GS_DEFAULT_TOLERANCE};
    use plotpy::Plot;
    use russell_chk::assert_approx_eq;

    #[test]
    fn new_handles_errors() {
        let coordinates = vec![vec![0.0]];
        let triangles = vec![vec![]];
        assert_eq!(
            GridSearchTri::new(&coordinates, &triangles, None, None).err(),
            Some("number of points must be ≥ 3")
        );
        let coordinates = vec![vec![0.0], vec![1.0], vec![2.0]];
        assert_eq!(
            GridSearchTri::new(&coordinates, &triangles, None, None).err(),
            Some("coordinates.ncol must be ≥ 2")
        );
        let coordinates = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0]];
        let triangles = vec![vec![0, 1, 2, 3], vec![0, 1, 2], vec![0, 1, 2]];
        assert_eq!(
            GridSearchTri::new(&coordinates, &triangles, None, None).err(),
            Some("number of triangle nodes (nnode) must be = 3")
        );
        let triangles = vec![vec![0, 1, 2], vec![0, 1, 2], vec![0, 1, 2]];
        assert_eq!(
            GridSearchTri::new(&coordinates, &triangles, Some(-1.0), None).err(),
            Some("tolerance must be > 0.0")
        );
        assert_eq!(
            GridSearchTri::new(&coordinates, &triangles, None, Some(-1.0)).err(),
            Some("border_tol must be ≥ 0.0")
        );
    }

    #[test]
    fn new_works_1() {
        let coordinates = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let triangles = &[[0, 1, 3], [1, 2, 3]];
        let tolerance = 1e-3;
        let border_tol = 0.1;
        let grid = GridSearchTri::new(coordinates, triangles, Some(tolerance), Some(border_tol)).unwrap();
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
             ntriangle = 2\n\
             ncontainer = 4\n\
             ndiv = [2, 2]\n"
        );
        grid.print_stat();
        let (mi, ma) = grid.limits();
        assert_eq!(mi, xmin);
        assert_eq!(ma, xmax);

        // let mut plot = Plot::new();
        // grid.draw(&mut plot, coordinates, triangles, true).unwrap();
        // plot.set_equal_axes(true)
        //     .set_figure_size_points(600.0, 600.0)
        //     .grid_and_labels("x", "y")
        //     .set_ticks_x(0.2, 0.0, "")
        //     .set_ticks_y(0.2, 0.0, "")
        //     .save("/tmp/gemlab/test_grid_search_tri_new_1.svg")
        //     .unwrap();

        // with zero border
        let grid = GridSearchTri::new(coordinates, triangles, None, Some(0.0)).unwrap();
        let sl = max_len + 2.0 * GS_DEFAULT_TOLERANCE;
        assert_eq!(grid.xmin, &[0.0, 0.0]);
        assert_eq!(grid.ndiv, &[1, 1]);
        assert_eq!(grid.xmax, &[sl, sl]);
    }

    #[test]
    fn draw_works() {
        let coordinates = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let triangles = &[[0, 1, 3], [1, 2, 3]];
        let grid = GridSearchTri::new(coordinates, triangles, None, None).unwrap();
        let mut plot = Plot::new();
        grid.draw(&mut plot, coordinates, triangles, true).unwrap();
    }

    #[test]
    fn find_triangle_handles_errors() {
        let coordinates = &[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let triangles = &[[0, 1, 2]];
        let grid = GridSearchTri::new(coordinates, triangles, None, None).unwrap();
        let y = vec![10.0, 0.0];
        assert_eq!(
            grid.find_triangle(&y, coordinates, triangles).err(),
            Some("given point coordinates are outside the grid")
        );
        assert_eq!(
            grid.find_triangle_and_interpolate(&y, coordinates, triangles).err(),
            Some("given point coordinates are outside the grid")
        );
        let y = vec![0.0, 0.0];
        assert_eq!(
            grid.find_triangle_and_interpolate(&y, coordinates, triangles).err(),
            Some("coordinates must contain a third column with the temperature values")
        );
    }

    #[test]
    fn find_triangle_works() {
        #[rustfmt::skip]
        let coordinates = &[
            [0.0307942, 0.459123  ], // 0
            [0.0980015, 0.981755  ], // 1
            [0.133721,  0.348832  ], // 2
            [0.13928,   0.180603  ], // 3
            [0.230951,  0.558482  ], // 4
            [0.478554,  0.00869692], // 5
            [0.540745,  0.331184  ], // 6
            [0.578587,  0.760349  ], // 7
            [0.648071,  0.369534  ], // 8
            [0.903726,  0.975904  ], // 9
        ];
        let triangles = &[
            [4, 2, 6], //  0
            [3, 2, 0], //  1
            [0, 4, 1], //  2
            [4, 0, 2], //  3
            [1, 4, 7], //  4
            [2, 3, 6], //  5
            [6, 7, 4], //  6
            [6, 5, 8], //  7
            [7, 8, 9], //  8
            [8, 7, 6], //  9
            [7, 9, 1], // 10
            [6, 3, 5], // 11
        ];
        let grid = GridSearchTri::new(coordinates, triangles, None, None).unwrap();
        let max_len = coordinates[9][0] - coordinates[1][0];
        let sl = max_len + 2.0 * GS_DEFAULT_TOLERANCE; // because the bbox is expanded
        let g = GS_DEFAULT_BORDER_TOL;
        assert_eq!(grid.side_length, sl);
        assert_eq!(grid.ndiv, &[2, 2]);
        assert_eq!(grid.xmin, &[coordinates[0][0] - g, coordinates[5][1] - g]);
        assert_eq!(
            grid.xmax,
            &[coordinates[0][0] - g + sl * 2.0, coordinates[5][1] - g + sl * 2.0]
        );
        assert_eq!(grid.bounding_boxes.len(), triangles.len());
        assert_eq!(grid.containers.len(), 4);
        let bbox_0 = &grid.bounding_boxes[0];
        assert_eq!(
            bbox_0,
            &[
                [coordinates[2][0], coordinates[6][0]],
                [coordinates[6][1], coordinates[4][1]]
            ]
        );
        assert_eq!(
            format!("{}", grid),
            "0: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]\n\
             1: [8, 10]\n\
             2: [2, 4, 8, 10]\n\
             3: [8, 10]\n\
             ntriangle = 12\n\
             ncontainer = 4\n\
             ndiv = [2, 2]\n"
        );
        assert_eq!(
            grid.find_triangle(&[0.4, 0.2], coordinates, triangles).unwrap(),
            Some(11)
        );
        assert_eq!(
            grid.find_triangle(&[0.6, 0.3], coordinates, triangles).unwrap(),
            Some(7)
        );
        assert_eq!(
            grid.find_triangle(&[0.1, 0.7], coordinates, triangles).unwrap(),
            Some(2)
        );
        assert_eq!(
            grid.find_triangle(&[0.8, 0.8], coordinates, triangles).unwrap(),
            Some(8)
        );
        let res = grid.find_triangle(&coordinates[5], coordinates, triangles).unwrap();
        if res != Some(7) {
            assert_eq!(res, Some(11));
        }
        assert_eq!(grid.find_triangle(&[0.1, 0.1], coordinates, triangles).unwrap(), None);
        assert_eq!(grid.find_triangle(&[0.6, 0.2], coordinates, triangles).unwrap(), None);
        assert_eq!(grid.find_triangle(&[0.4, 1.0], coordinates, triangles).unwrap(), None);
        assert_eq!(
            grid.find_triangle(&[10.0, 1.0], coordinates, triangles).err(),
            Some("given point coordinates are outside the grid")
        );
        // let mut plot = Plot::new();
        // grid.draw(&mut plot, coordinates, triangles, true).unwrap();
        // plot.set_equal_axes(true)
        //     .set_figure_size_points(600.0, 600.0)
        //     .grid_and_labels("x", "y")
        //     .set_ticks_x(0.2, 0.0, "")
        //     .set_ticks_y(0.2, 0.0, "")
        //     .save("/tmp/gemlab/test_grid_search_tri_find_works.svg")
        //     .unwrap();
    }

    #[test]
    fn find_triangle_and_interpolate_works() {
        let coordinates = &[
            [0.0, 0.0, 0.0], // last column is the temperature
            [0.5, 0.85, 0.986154146165801],
            [1.0, 0.0, 1.0],
            [1.0, 1.7, 1.972308292331602],
            [1.5, 0.85, 1.724093964956667],
            [2.0, 0.0, 2.0],
            [2.0, 1.7, 2.6248809496813372],
            [2.5, 0.85, 2.640549185302179],
            [3.0, 1.7, 3.448187929913334],
        ];
        let triangles = &[
            [0, 2, 1],
            [2, 5, 4],
            [1, 2, 4],
            [4, 5, 7],
            [1, 4, 3],
            [4, 7, 6],
            [3, 4, 6],
            [6, 7, 8],
        ];
        let grid = GridSearchTri::new(coordinates, triangles, None, None).unwrap();
        for x_y_tt in coordinates {
            assert_eq!(
                grid.find_triangle_and_interpolate(x_y_tt, coordinates, triangles)
                    .unwrap(),
                Some(x_y_tt[2])
            );
            let tt = f64::sqrt(x_y_tt[0] * x_y_tt[0] + x_y_tt[1] * x_y_tt[1]);
            assert_eq!(
                grid.find_triangle_and_interpolate(x_y_tt, coordinates, triangles)
                    .unwrap(),
                Some(tt)
            );
        }
        let temp = grid
            .find_triangle_and_interpolate(&[1.5, 1.0], coordinates, triangles)
            .unwrap()
            .unwrap();
        assert_approx_eq!(temp, f64::sqrt(1.5 * 1.5 + 1.0 * 1.0), 0.025);
        // let mut plot = Plot::new();
        // grid.draw(&mut plot, coordinates, triangles, true).unwrap();
        // plot.set_equal_axes(true)
        //     .set_figure_size_points(600.0, 600.0)
        //     .grid_and_labels("x", "y")
        //     .set_ticks_x(0.2, 0.0, "")
        //     .set_ticks_y(0.2, 0.0, "")
        //     .save("/tmp/gemlab/test_grid_search_tri_interpolate_works.svg")
        //     .unwrap();
    }

    #[test]
    fn find_triangle_works_2d_particles() {
        #[rustfmt::skip]
        let coordinates = &[
            [0.0, 0.0, 10.0], [1.0, 0.0, 10.0], [0.5, 0.85, 10.0],
            [2.0, 0.0, 10.0], [3.0, 0.0, 10.0], [2.5, 0.85, 10.0],
            [1.0, 1.0, 10.0], [2.0, 1.0, 10.0], [1.5, 1.85, 10.0],
        ];
        let triangles = &[[0, 1, 2], [3, 4, 5], [6, 7, 8]];
        let grid = GridSearchTri::new(coordinates, triangles, None, None).unwrap();
        assert_eq!(
            grid.find_triangle(&[0.4, 0.2], coordinates, triangles).unwrap(),
            Some(0)
        );
        assert_eq!(
            grid.find_triangle(&[3.0, 0.0], coordinates, triangles).unwrap(),
            Some(1)
        );
        assert_eq!(
            grid.find_triangle(&[1.5, 1.4], coordinates, triangles).unwrap(),
            Some(2)
        );
        assert_eq!(grid.find_triangle(&[0.5, 1.6], coordinates, triangles).unwrap(), None);
        assert_eq!(grid.find_triangle(&[1.5, 0.5], coordinates, triangles).unwrap(), None);
        assert_eq!(grid.find_triangle(&[2.5, 1.6], coordinates, triangles).unwrap(), None);
        assert_eq!(grid.find_triangle(&[2.5, 1.6], coordinates, triangles).unwrap(), None);
        let (c, t) = (coordinates, triangles);
        assert_eq!(
            grid.find_triangle_and_interpolate(&[0.4, 0.2], c, t).unwrap(),
            Some(10.0)
        );
        assert_eq!(
            grid.find_triangle_and_interpolate(&[3.0, 0.0], c, t).unwrap(),
            Some(10.0)
        );
        assert_eq!(
            grid.find_triangle_and_interpolate(&[1.5, 1.4], c, t).unwrap(),
            Some(10.0)
        );
        assert_eq!(grid.find_triangle_and_interpolate(&[0.5, 1.6], c, t).unwrap(), None);
        assert_eq!(grid.find_triangle_and_interpolate(&[1.5, 0.5], c, t).unwrap(), None);
        assert_eq!(grid.find_triangle_and_interpolate(&[2.5, 1.6], c, t).unwrap(), None);
        assert_eq!(
            grid.find_triangle_and_interpolate(&[2.5, 1.6], coordinates, triangles)
                .unwrap(),
            None
        );
    }
}
