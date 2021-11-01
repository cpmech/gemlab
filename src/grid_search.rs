use plotpy::{Curve, Plot, Shapes, Text};
use std::collections::HashMap;
use std::fmt;

/// Holds the id and coordinates of an item
struct Item {
    id: usize,   // identification number
    x: Vec<f64>, // (ndim) coordinates
}

/// Holds items
struct Container {
    items: Vec<Item>,
}

/// Implements a grid for fast searching entries by coordinates
///
/// # Reference
///
/// * Durand, Farias, and Pedroso (2015) Computing intersections between
///   non-compatible curves and finite elements, Computational Mechanics;
///   DOI=10.1007/s00466-015-1181-y
pub struct GridSearch {
    // constants
    ndim: usize,      // space dimension
    ndiv: [usize; 3], // number of divisions along each direction
    min: [f64; 3],    // min values
    max: [f64; 3],    // max values
    delta: [f64; 3],  // difference between max and min
    size: [f64; 3],   // size of each container
    cf: [usize; 3],   // coefficients [1, ndiv[0], ndiv[0]*ndiv[1]] (Eq. 8)

    // auxiliary variable
    ratio: [usize; 3], // ratio = trunc(δx[i]/Δx[i]) (Eq. 8)

    // square/cubic halo: bounding box corners around point, including the point
    tol: [f64; 3],       // tolerances to compare coordinates and define the halo
    halo: [[f64; 3]; 8], // 4 in 2D or 8 in 3D (each contains 3 coords)
    ncorner: usize,      // number of halo corners 4 in 2D or 8 in 3D

    // holds non-empty containers. maps container.index to container.data
    // a point may be located in more than one container (e.g., when at boundaries)
    containers: HashMap<usize, Container>,
}

impl GridSearch {
    /// Creates a new GridSearch
    ///
    /// # Input
    ///
    /// * `ndiv` -- number of divisions along each dimension (len = 2 or 3 == ndim)
    /// * `min` -- min coordinates (len = ndim)
    /// * `max` -- max coordinates (len = ndim)
    /// * `tol` -- tolerances to compare coordinates (len = ndim)
    pub fn new(ndiv: &[usize], min: &[f64], max: &[f64], tol: &[f64]) -> Result<Self, &'static str> {
        // check input
        let ndim = ndiv.len();
        if ndim < 2 || ndim > 3 {
            return Err("ndiv.len() must be 2 or 3 = ndim");
        }
        if min.len() != ndim {
            return Err("xmin.len() must equal ndiv.len()");
        }
        if max.len() != ndim {
            return Err("xmax.len() must equal ndiv.len()");
        }
        if tol.len() != ndim {
            return Err("tol.len() must equal ndiv.len()");
        }

        // allocate gs
        let mut grid = GridSearch {
            ndim,
            ndiv: [0; 3],
            min: [0.0; 3],
            max: [0.0; 3],
            delta: [0.0; 3],
            size: [0.0; 3],
            cf: [0; 3],
            ratio: [0; 3],
            tol: [0.0; 3],
            halo: [[0.0; 3]; 8],
            ncorner: usize::pow(2, ndim as u32),
            containers: HashMap::new(),
        };

        // check and compute sizes
        for i in 0..ndim {
            grid.ndiv[i] = ndiv[i];
            if ndiv[i] < 1 {
                return Err("ndiv must be at least equal to 1");
            }
            grid.min[i] = min[i];
            grid.max[i] = max[i];
            grid.delta[i] = max[i] - min[i];
            if grid.delta[i] <= 0.0 {
                return Err("xmax must be greater than xmin");
            }
            grid.size[i] = grid.delta[i] / (ndiv[i] as f64);
            grid.tol[i] = tol[i];
        }

        // coefficient
        grid.cf[0] = 1;
        grid.cf[1] = grid.ndiv[0];
        grid.cf[2] = grid.ndiv[0] * grid.ndiv[1];

        // done
        Ok(grid)
    }

    /// Inserts a new item to the right container in the grid
    ///
    /// # Input
    ///
    /// * `id` -- identification number for the item
    /// * `x` -- coordinates (ndim) of the item
    pub fn insert(&mut self, id: usize, x: &[f64]) -> Result<(), &'static str> {
        // check
        if x.len() != self.ndim {
            return Err("x.len() must equal ndim");
        }

        // add point to container
        let index = match self.container_index(x) {
            Some(i) => i,
            None => return Err("point is outside the grid"),
        };
        self.update_or_insert(index, id, x);

        // add point to containers touched by halo corners
        self.set_halo(x);
        let mut tmp = vec![0.0; self.ndim];
        for c in 0..self.ncorner {
            tmp.copy_from_slice(&self.halo[c][0..self.ndim]);
            if let Some(index_corner) = self.container_index(&tmp) {
                if index_corner != index {
                    self.update_or_insert(index_corner, id, x); // make sure to use original `x`
                }
            }
        }
        Ok(())
    }

    /// Find previously inserted item to the grid
    ///
    /// # Input
    ///
    /// * `x` -- coordinates (ndim) of the item
    ///
    /// # Output
    ///
    /// * `id` -- if found, returns the identification number of the item
    pub fn find(&mut self, x: &[f64]) -> Result<Option<usize>, &'static str> {
        // check
        if x.len() != self.ndim {
            return Err("x.len() must equal ndim");
        }

        // find index of container where x should be
        let index = match self.container_index(x) {
            Some(i) => i,
            None => return Err("point is outside the grid"),
        };

        // find container that should have a point close to `x`
        let container = match self.containers.get(&index) {
            Some(c) => c,
            None => return Ok(None), // no container has a point close to `x`
        };

        // find closest point to `x` in the container
        for item in &container.items {
            if self.ndim == 2 {
                if f64::abs(x[0] - item.x[0]) <= self.tol[0] && f64::abs(x[1] - item.x[1]) <= self.tol[1] {
                    return Ok(Some(item.id));
                }
            } else {
                if f64::abs(x[0] - item.x[0]) <= self.tol[0]
                    && f64::abs(x[1] - item.x[1]) <= self.tol[1]
                    && f64::abs(x[2] - item.x[2]) <= self.tol[2]
                {
                    return Ok(Some(item.id));
                }
            }
        }
        Ok(None)
    }

    /// Returns a drawing of this object
    pub fn plot(&self) -> Result<Plot, &'static str> {
        // create plot
        let mut plot = Plot::new();

        // draw grid
        let mut xmin = vec![0.0; self.ndim];
        let mut xmax = vec![0.0; self.ndim];
        let mut ndiv = vec![0; self.ndim];
        for i in 0..self.ndim {
            xmin[i] = self.min[i];
            xmax[i] = self.max[i];
            ndiv[i] = self.ndiv[i];
        }
        let mut shapes = Shapes::new();
        shapes
            .set_alt_text_color("#5d5d5d")
            .draw_grid(&xmin, &xmax, &ndiv, false, true)?;
        plot.add(&shapes);

        // draw items
        let mut curve = Curve::new();
        let mut text = Text::new();
        curve
            .set_marker_style("o")
            .set_marker_color("#fab32faa")
            .set_marker_line_color("black")
            .set_marker_line_width(0.5);
        text.set_color("#cd0000");
        for container in self.containers.values() {
            for item in &container.items {
                let txt = format!("{}", item.id);
                if self.ndim == 2 {
                    curve.draw(&[item.x[0]], &[item.x[1]]);
                    text.draw(item.x[0], item.x[1], &txt);
                } else {
                    curve.draw_3d(&[item.x[0]], &[item.x[1]], &[item.x[2]]);
                    text.draw_3d(item.x[0], item.x[1], item.x[2], &txt);
                }
            }
        }
        plot.add(&curve);
        plot.add(&text);

        // done
        Ok(plot)
    }

    /// Calculates the container index where the point x should be located
    ///
    /// # Output
    ///
    /// * returns the index of the container or None if the point is out-of-range
    fn container_index(&mut self, x: &[f64]) -> Option<usize> {
        let mut index = 0;
        for i in 0..self.ndim {
            if x[i] < self.min[i] || x[i] > self.max[i] {
                return None;
            }
            self.ratio[i] = ((x[i] - self.min[i]) / self.size[i]) as usize;
            if self.ratio[i] == self.ndiv[i] {
                // the point is exactly on the max edge, thus select inner container
                self.ratio[i] -= 1; // move to the inside
            }
            index += self.ratio[i] * self.cf[i];
        }
        Some(index)
    }

    /// Update container or insert point in container
    fn update_or_insert(&mut self, index: usize, id: usize, x: &[f64]) {
        let item = Item { id, x: Vec::from(x) };
        if self.containers.contains_key(&index) {
            let container = self.containers.get_mut(&index).unwrap();
            container.items.push(item);
        } else {
            self.containers.insert(index, Container { items: vec![item] });
        }
    }

    /// Sets square/cubic halo around point
    fn set_halo(&mut self, x: &[f64]) {
        if self.ndim == 2 {
            self.halo[0][0] = x[0] - self.tol[0];
            self.halo[0][1] = x[1] - self.tol[1];

            self.halo[1][0] = x[0] + self.tol[0];
            self.halo[1][1] = x[1] - self.tol[1];

            self.halo[2][0] = x[0] + self.tol[0];
            self.halo[2][1] = x[1] + self.tol[1];

            self.halo[3][0] = x[0] - self.tol[0];
            self.halo[3][1] = x[1] + self.tol[1];
        } else {
            self.halo[0][0] = x[0] - self.tol[0];
            self.halo[0][1] = x[1] - self.tol[1];
            self.halo[0][2] = x[2] - self.tol[2];

            self.halo[1][0] = x[0] + self.tol[0];
            self.halo[1][1] = x[1] - self.tol[1];
            self.halo[1][2] = x[2] - self.tol[2];

            self.halo[2][0] = x[0] + self.tol[0];
            self.halo[2][1] = x[1] + self.tol[1];
            self.halo[2][2] = x[2] - self.tol[2];

            self.halo[3][0] = x[0] - self.tol[0];
            self.halo[3][1] = x[1] + self.tol[1];
            self.halo[3][2] = x[2] - self.tol[2];

            self.halo[4][0] = x[0] - self.tol[0];
            self.halo[4][1] = x[1] - self.tol[1];
            self.halo[4][2] = x[2] + self.tol[2];

            self.halo[5][0] = x[0] + self.tol[0];
            self.halo[5][1] = x[1] - self.tol[1];
            self.halo[5][2] = x[2] + self.tol[2];

            self.halo[6][0] = x[0] + self.tol[0];
            self.halo[6][1] = x[1] + self.tol[1];
            self.halo[6][2] = x[2] + self.tol[2];

            self.halo[7][0] = x[0] - self.tol[0];
            self.halo[7][1] = x[1] + self.tol[1];
            self.halo[7][2] = x[2] + self.tol[2];
        }
    }
}

impl fmt::Display for GridSearch {
    /// Shows info about the items in the grid containers
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // items
        let mut unique_items: HashMap<usize, bool> = HashMap::new();
        let mut indices: Vec<_> = self.containers.keys().collect();
        indices.sort();
        for index in indices {
            let container = self.containers.get(index).unwrap();
            let mut ids: Vec<_> = container.items.iter().map(|item| item.id).collect();
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
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    struct TestData<'a> {
        id: usize,
        x: &'a [f64],
        container: usize,        // container index calc with x only
        containers: &'a [usize], // container index calc with halo and x
    }

    fn get_test_grid_2d() -> GridSearch {
        GridSearch::new(&[5, 5], &[-0.2, -0.2], &[0.8, 1.8], &[1e-2, 1e-2]).unwrap()
    }

    fn get_test_grid_3d() -> GridSearch {
        GridSearch::new(&[3, 3, 3], &[-1.0, -1.0, -1.0], &[1.0, 1.0, 1.0], &[1e-2, 1e-2, 1e-2]).unwrap()
    }

    fn get_test_data_2d<'a>() -> Vec<TestData<'a>> {
        vec![
            TestData {
                id: 100,
                x: &[0.0, 0.2],
                container: 6,
                containers: &[0, 1, 5, 6],
            },
            TestData {
                id: 200,
                x: &[0.2, 0.6000000000000001],
                container: 12,
                containers: &[6, 7, 11, 12],
            },
            TestData {
                id: 300,
                x: &[0.4, 1.0],
                // in this case, the ratio is:
                // 3.0000000000000004, 2.9999999999999996
                // thus, the point falls in #13 instead of #18
                container: 13,
                containers: &[12, 13, 17, 18],
            },
            TestData {
                id: 400,
                x: &[0.6, 1.4],
                // in this case, the ratio is:
                // 4, 3.9999999999999996
                // thus, the point falls in #19 instead of #24
                container: 19,
                containers: &[18, 19, 23, 24],
            },
            TestData {
                id: 500,
                x: &[0.8, 1.8],
                container: 24,
                containers: &[24],
            },
        ]
    }

    fn get_test_data_3d<'a>() -> Vec<TestData<'a>> {
        const L: f64 = 2.0;
        vec![
            TestData {
                id: 100,
                x: &[-1.0, -1.0, -1.0],
                container: 0,
                containers: &[0],
            },
            TestData {
                id: 200,
                x: &[0.0, 0.0, 0.0],
                container: 13,
                containers: &[13],
            },
            TestData {
                id: 300,
                x: &[-1.0 + L * 2.0 / 3.0, -1.0 + L * 2.0 / 3.0, -1.0 + L * 2.0 / 3.0],
                container: 26,
                containers: &[13, 14, 16, 17, 22, 23, 25, 26],
            },
        ]
    }

    #[test]
    fn new_fails_on_wrong_input() {
        let grid = GridSearch::new(&[1], &[0.0, 0.0], &[1.0, 1.0], &[0.1, 0.1]);
        assert_eq!(grid.err(), Some("ndiv.len() must be 2 or 3 = ndim"));

        let grid = GridSearch::new(&[1, 1], &[0.0], &[1.0, 1.0], &[0.1, 0.1]);
        assert_eq!(grid.err(), Some("xmin.len() must equal ndiv.len()"));

        let grid = GridSearch::new(&[1, 1], &[0.0, 0.0], &[1.0], &[0.1, 0.1]);
        assert_eq!(grid.err(), Some("xmax.len() must equal ndiv.len()"));

        let grid = GridSearch::new(&[0, 1], &[0.0, 0.0], &[1.0, 1.0], &[0.1, 0.1]);
        assert_eq!(grid.err(), Some("ndiv must be at least equal to 1"));

        let grid = GridSearch::new(&[0, 1], &[0.0, 0.0], &[1.0, 1.0], &[0.1]);
        assert_eq!(grid.err(), Some("tol.len() must equal ndiv.len()"));

        let grid = GridSearch::new(&[1, 1], &[0.0, 0.0], &[0.0, 1.0], &[0.1, 0.1]);
        assert_eq!(grid.err(), Some("xmax must be greater than xmin"));
    }

    #[test]
    fn new_works() {
        let g2d = get_test_grid_2d();
        assert_eq!(g2d.ndim, 2);
        assert_eq!(g2d.ndiv, [5, 5, 0]);
        assert_eq!(g2d.min, [-0.2, -0.2, 0.0]);
        assert_eq!(g2d.max, [0.8, 1.8, 0.0]);
        assert_eq!(g2d.delta, [1.0, 2.0, 0.0]);
        assert_eq!(g2d.size, [0.2, 0.4, 0.0]);
        assert_eq!(g2d.cf, [1, 5, 25]);
        assert_eq!(g2d.ratio, [0, 0, 0]);
        assert_eq!(g2d.tol, [1e-2, 1e-2, 0.0]);
        assert_eq!(g2d.halo.len(), 8);
        assert_eq!(g2d.ncorner, 4);
        assert_eq!(g2d.containers.len(), 0);

        let g3d = get_test_grid_3d();
        assert_eq!(g3d.ndim, 3);
        assert_eq!(g3d.ndiv, [3, 3, 3]);
        assert_eq!(g3d.min, [-1.0, -1.0, -1.0]);
        assert_eq!(g3d.max, [1.0, 1.0, 1.0]);
        assert_eq!(g3d.delta, [2.0, 2.0, 2.0]);
        assert_eq!(g3d.size, [2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]);
        assert_eq!(g3d.cf, [1, 3, 9]);
        assert_eq!(g3d.ratio, [0, 0, 0]);
        assert_eq!(g3d.tol, [1e-2, 1e-2, 1e-2]);
        assert_eq!(g3d.halo.len(), 8);
        assert_eq!(g3d.ncorner, 8);
        assert_eq!(g3d.containers.len(), 0);
    }

    #[test]
    fn display_trait_works() {
        let g2d = get_test_grid_2d();
        assert_eq!(
            format!("{}", g2d),
            "ids = []\n\
             nitem = 0\n\
             ncontainer = 0\n"
        );

        let g3d = get_test_grid_3d();
        assert_eq!(
            format!("{}", g3d),
            "ids = []\n\
             nitem = 0\n\
             ncontainer = 0\n"
        );
    }

    #[test]
    fn set_halo_works() {
        let mut g2d = get_test_grid_2d();
        g2d.set_halo(&[0.5, 0.5]);
        assert_eq!(g2d.halo[0], [0.49, 0.49, 0.0]);
        assert_eq!(g2d.halo[1], [0.51, 0.49, 0.0]);
        assert_eq!(g2d.halo[2], [0.51, 0.51, 0.0]);
        assert_eq!(g2d.halo[3], [0.49, 0.51, 0.0]);

        let mut g3d = get_test_grid_3d();
        g3d.set_halo(&[0.5, 0.5, 0.5]);
        assert_eq!(g3d.halo[0], [0.49, 0.49, 0.49]);
        assert_eq!(g3d.halo[1], [0.51, 0.49, 0.49]);
        assert_eq!(g3d.halo[2], [0.51, 0.51, 0.49]);
        assert_eq!(g3d.halo[3], [0.49, 0.51, 0.49]);
        assert_eq!(g3d.halo[4], [0.49, 0.49, 0.51]);
        assert_eq!(g3d.halo[5], [0.51, 0.49, 0.51]);
        assert_eq!(g3d.halo[6], [0.51, 0.51, 0.51]);
        assert_eq!(g3d.halo[7], [0.49, 0.51, 0.51]);
    }

    #[test]
    fn container_index_works() -> Result<(), &'static str> {
        let mut g2d = get_test_grid_2d();
        for data in get_test_data_2d() {
            let index = g2d.container_index(data.x).unwrap();
            assert_eq!(index, data.container);
        }
        let index = g2d.container_index(&[0.80001, 0.0]);
        assert_eq!(index, None); // outside

        let mut g3d = get_test_grid_3d();
        for data in get_test_data_3d() {
            let index = g3d.container_index(data.x).unwrap();
            assert_eq!(index, data.container);
        }
        let index = g3d.container_index(&[1.00001, 0.0, 0.0]);
        assert_eq!(index, None); // outside
        Ok(())
    }

    #[test]
    fn insert_fails_on_wrong_input() {
        let mut g2d = get_test_grid_2d();
        let res = g2d.insert(0, &[0.0, 0.0, 0.0]);
        assert_eq!(res, Err("x.len() must equal ndim"));
        let res = g2d.insert(1000, &[0.80001, 0.0]);
        assert_eq!(res, Err("point is outside the grid"));

        let mut g3d = get_test_grid_3d();
        let res = g3d.insert(0, &[0.0, 0.0]);
        assert_eq!(res, Err("x.len() must equal ndim"));
        let res = g3d.insert(1000, &[1.00001, 0.0, 0.0]);
        assert_eq!(res, Err("point is outside the grid"));
    }

    #[test]
    fn insert_2d_works() -> Result<(), &'static str> {
        let mut grid = get_test_grid_2d();
        for data in get_test_data_2d() {
            grid.insert(data.id, data.x)?;
            for index in data.containers {
                let container = grid.containers.get(index).unwrap();
                container.items.iter().find(|item| item.id == data.id).unwrap();
            }
        }
        assert_eq!(
            format!("{}", grid),
            "0: [100]\n\
             1: [100]\n\
             5: [100]\n\
             6: [100, 200]\n\
             7: [200]\n\
             11: [200]\n\
             12: [200, 300]\n\
             13: [300]\n\
             17: [300]\n\
             18: [300, 400]\n\
             19: [400]\n\
             23: [400]\n\
             24: [400, 500]\n\
             ids = [100, 200, 300, 400, 500]\n\
             nitem = 5\n\
             ncontainer = 13\n"
        );
        let mut indices: Vec<_> = grid.containers.into_keys().collect();
        indices.sort();
        assert_eq!(indices, &[0, 1, 5, 6, 7, 11, 12, 13, 17, 18, 19, 23, 24]);
        Ok(())
    }

    #[test]
    fn insert_3d_works() -> Result<(), &'static str> {
        let mut grid = get_test_grid_3d();
        for data in get_test_data_3d() {
            grid.insert(data.id, data.x)?;
            for index in data.containers {
                let container = grid.containers.get(index).unwrap();
                container.items.iter().find(|item| item.id == data.id).unwrap();
            }
        }
        assert_eq!(
            format!("{}", grid),
            "0: [100]\n\
             13: [200, 300]\n\
             14: [300]\n\
             16: [300]\n\
             17: [300]\n\
             22: [300]\n\
             23: [300]\n\
             25: [300]\n\
             26: [300]\n\
             ids = [100, 200, 300]\n\
             nitem = 3\n\
             ncontainer = 9\n"
        );
        let mut indices: Vec<_> = grid.containers.into_keys().collect();
        indices.sort();
        assert_eq!(indices, &[0, 13, 14, 16, 17, 22, 23, 25, 26]);
        Ok(())
    }

    #[test]
    fn find_fails_on_wrong_input() {
        let mut g2d = get_test_grid_2d();
        let res = g2d.find(&[0.0, 0.0, 0.0]);
        assert_eq!(res, Err("x.len() must equal ndim"));
        let res = g2d.find(&[0.80001, 0.0]);
        assert_eq!(res, Err("point is outside the grid"));

        let mut g3d = get_test_grid_3d();
        let res = g3d.find(&[0.0, 0.0]);
        assert_eq!(res, Err("x.len() must equal ndim"));
        let res = g3d.find(&[1.00001, 0.0, 0.0]);
        assert_eq!(res, Err("point is outside the grid"));
    }

    #[test]
    fn find_works() -> Result<(), &'static str> {
        let mut g2d = get_test_grid_2d();
        for data in get_test_data_2d() {
            g2d.insert(data.id, data.x)?;
            let id = g2d.find(data.x)?;
            assert_eq!(id, Some(data.id));
        }
        let id = g2d.find(&[0.5, 0.5])?;
        assert_eq!(id, None);

        let mut g3d = get_test_grid_3d();
        for data in get_test_data_3d() {
            g3d.insert(data.id, data.x)?;
            let id = g3d.find(data.x)?;
            assert_eq!(id, Some(data.id));
        }
        let id = g3d.find(&[0.5, 0.5, 0.5])?;
        assert_eq!(id, None);
        Ok(())
    }

    #[test]
    fn plot_works() -> Result<(), &'static str> {
        let mut g2d = get_test_grid_2d();
        for data in get_test_data_2d() {
            g2d.insert(data.id, data.x)?;
        }
        let plot = g2d.plot()?;
        plot.save("/tmp/gemlab/search_grid_plot_2d_works.svg")?;

        let mut g3d = get_test_grid_3d();
        for data in get_test_data_3d() {
            g3d.insert(data.id, data.x)?;
        }
        let plot = g3d.plot()?;
        plot.save("/tmp/gemlab/search_grid_plot_3d_works.svg")?;
        Ok(())
    }
}
