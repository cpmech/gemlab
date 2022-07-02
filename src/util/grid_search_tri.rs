use super::{calc_container_key, AsArray2D, GS_DEFAULT_BORDER_TOL, GS_DEFAULT_TOLERANCE};
use crate::geometry::is_point_inside_triangle;
use crate::StrError;
use plotpy::{Canvas, Plot, PolyCode, Text};
use russell_stat::Histogram;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// 2D
const NDIM: usize = 2;

/// Specifies the key of containers (or bins in the grid)
type ContainerKey = usize;

/// Specifies the identification number of cells (must be sequential from 0 to ncell - 1)
type TriangleId = usize;

/// Defines the container type
type Container = HashSet<TriangleId>;

/// Defines the containers type: Key to Container
type Containers = HashMap<ContainerKey, Container>;

/// Defines the bounding box of a cell
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
    ndiv: Vec<usize>,                // (ndim) number of divisions along each direction
    xmin: Vec<f64>,                  // (ndim) min values
    xmax: Vec<f64>,                  // (ndim) max values
    side_length: f64,                // side length of a container
    bounding_boxes: Vec<BboxMinMax>, // (ncell) bounding boxes
    containers: Containers,          // structure to hold all items
}

impl GridSearchTri {
    /// Allocates a new instance
    ///
    /// # Input
    ///
    /// * `coordinates` -- is a list of coordinates such as `[[x0,y0], [x1,y1], [x2,y2], [x3,y3]]`
    /// * `triangles` -- is a list of connectivity (topology) such as `[[0,2,1], [2,1,0]]`
    /// * `tolerance` -- is a tolerance to expand the bounding box of cells and compare points; e.g. 1e-4
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
        let (npoint, ndim) = coordinates.size();
        let (ntriangle, nnode) = triangles.size();
        if npoint < 3 {
            return Err("number of points must be ≥ 3");
        }
        if ndim != 2 {
            return Err("ndim must be = 2");
        }
        if ntriangle < 1 {
            return Err("number of triangles must be ≥ 1");
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

        // find limits, bounding boxes, and largest cell
        let mut xmin = vec![f64::MAX; ndim];
        let mut xmax = vec![f64::MIN; ndim];
        let mut x_min_max = vec![vec![0.0; N_MIN_MAX]; ndim];
        let mut bbox_large = vec![f64::MIN; ndim];
        let mut bounding_boxes = Vec::new();
        let mut x = vec![0.0; ndim];
        for t in 0..ntriangle {
            for m in 0..nnode {
                let p = triangles.at(t, m);
                x[0] = coordinates.at(p, 0);
                x[1] = coordinates.at(p, 1);
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
        for t in 0..ntriangle {
            let x_min_max = &bounding_boxes[t];
            for r in 0..N_MIN_MAX {
                for s in 0..N_MIN_MAX {
                    for t in 0..(ndim - 1) {
                        x[0] = x_min_max[0][r];
                        x[1] = x_min_max[1][s];
                        let key = calc_container_key(ndim, side_length, &ndiv, &xmin, &x);
                        let container = containers.entry(key).or_insert(HashSet::new());
                        container.insert(t);
                    }
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

    /// Find the triangle where the given coordinate falls in
    ///
    /// * Returns the TriangleId or None if no triangle contains the point
    ///
    /// # Warning
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
        for t in container {
            let x_min_max = &self.bounding_boxes[*t];
            for i in 0..NDIM {
                if x[i] < x_min_max[i][I_MIN] || x[i] > x_min_max[i][I_MAX] {
                    continue; // outside the bounding box
                }
                xa[i] = coordinates.at(triangles.at(*t, 0), i);
                xb[i] = coordinates.at(triangles.at(*t, 1), i);
                xc[i] = coordinates.at(triangles.at(*t, 2), i);
            }
            if is_point_inside_triangle(&xa, &xb, &xc, x) {
                return Ok(Some(*t));
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
                for cell_id in container {
                    let x_min_max = &self.bounding_boxes[*cell_id];
                    let txt = format!("{}", cell_id);
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
        println!("ncell = {:?}", unique_items.len());
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
        write!(f, "ncell = {}\n", unique_items.len()).unwrap();
        write!(f, "ncontainer = {}\n", self.containers.len()).unwrap();
        write!(f, "ndiv = {:?}\n", self.ndiv).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::GridSearchTri;
    use crate::util::GS_DEFAULT_TOLERANCE;
    use crate::StrError;
    use plotpy::Plot;

    #[test]
    fn new_handles_errors() {
        let coordinates = &[[0.0]];
        let triangles = &[[0, 1]];
        assert_eq!(
            GridSearchTri::new(coordinates, triangles, None, None).err(),
            Some("number of points must be ≥ 3")
        );
        assert_eq!(
            GridSearchTri::new(coordinates, triangles, None, None).err(),
            Some("ndim must be = 2")
        );
        assert_eq!(
            GridSearchTri::new(coordinates, triangles, None, None).err(),
            Some("number of triangles must be ≥ 1")
        );
        assert_eq!(
            GridSearchTri::new(coordinates, triangles, None, None).err(),
            Some("number of triangle nodes (nnode) must be = 3")
        );
        assert_eq!(
            GridSearchTri::new(coordinates, triangles, Some(-1.0), None).err(),
            Some("tolerance must be > 0.0")
        );
        assert_eq!(
            GridSearchTri::new(coordinates, triangles, None, Some(-1.0)).err(),
            Some("border_tol must be ≥ 0.0")
        );
    }

    #[test]
    fn new_works_1() -> Result<(), StrError> {
        let coordinates = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let triangles = &[[0, 1, 3], [1, 2, 3]];
        let tolerance = 1e-3;
        let border_tol = 0.1;
        let grid = GridSearchTri::new(coordinates, triangles, Some(tolerance), Some(border_tol))?;
        let max_len = 1.0;
        let sl = max_len + 2.0 * tolerance; // because the bbox is expanded
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
            grid.draw(&mut plot, coordinates, triangles, true)?;
            plot.set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .grid_and_labels("x", "y")
                .set_ticks_x(0.2, 0.0, "")
                .set_ticks_y(0.2, 0.0, "")
                .save("/tmp/gemlab/test_grid_search_tri_new_1.svg")?;
        }
        // with zero border
        let grid = GridSearchTri::new(coordinates, triangles, None, Some(0.0))?;
        let sl = max_len + 2.0 * GS_DEFAULT_TOLERANCE;
        assert_eq!(grid.xmin, &[0.0, 0.0]);
        assert_eq!(grid.ndiv, &[1, 1]);
        assert_eq!(grid.xmax, &[sl, sl]);
        Ok(())
    }

    #[test]
    fn draw_works() -> Result<(), StrError> {
        let coordinates = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let triangles = &[[0, 1, 3], [1, 2, 3]];
        let grid = GridSearchTri::new(coordinates, triangles, None, None)?;
        let mut plot = Plot::new();
        grid.draw(&mut plot, coordinates, triangles, false)?;
        Ok(())
    }

    #[test]
    fn find_triangle_handles_errors() -> Result<(), StrError> {
        let coordinates = &[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let triangles = &[[0, 1, 2]];
        let grid = GridSearchTri::new(coordinates, triangles, None, None)?;
        let y = vec![10.0, 0.0];
        assert_eq!(
            grid.find_triangle(&y, coordinates, triangles).err(),
            Some("given point coordinates are outside the grid")
        );
        Ok(())
    }

    #[test]
    fn find_triangle_works() -> Result<(), StrError> {
        let coordinates = &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let triangles = &[[0, 1, 3], [1, 2, 3]];
        let grid = GridSearchTri::new(coordinates, triangles, None, None)?;
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
        assert_eq!(grid.find_triangle(&[0.4, 0.2], coordinates, triangles)?, Some(11));
        assert_eq!(grid.find_triangle(&[0.6, 0.3], coordinates, triangles)?, Some(7));
        assert_eq!(grid.find_triangle(&[0.1, 0.7], coordinates, triangles)?, Some(2));
        assert_eq!(grid.find_triangle(&[0.8, 0.8], coordinates, triangles)?, Some(8));
        assert_eq!(grid.find_triangle(&[0.1, 0.1], coordinates, triangles)?, None);
        assert_eq!(grid.find_triangle(&[0.6, 0.2], coordinates, triangles)?, None);
        assert_eq!(grid.find_triangle(&[0.4, 1.0], coordinates, triangles)?, None);
        if false {
            let mut plot = Plot::new();
            grid.draw(&mut plot, coordinates, triangles, true)?;
            plot.set_equal_axes(true)
                .set_figure_size_points(600.0, 600.0)
                .grid_and_labels("x", "y")
                .set_ticks_x(0.2, 0.0, "")
                .set_ticks_y(0.2, 0.0, "")
                .save("/tmp/gemlab/test_grid_search_tri_find_works.svg")?;
        }
        Ok(())
    }
}
