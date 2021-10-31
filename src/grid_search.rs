use plotpy::{Curve, Plot, Shapes, Text};
use std::collections::HashMap;

/// Holds the id and coordinates of an item
struct Item {
    id: usize,   // identification number
    x: [f64; 3], // coordinates
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
    xmin: [f64; 3],   // min values
    xmax: [f64; 3],   // max values
    xdelta: [f64; 3], // difference between max and min
    size: [f64; 3],   // size of each container
    cf: [usize; 3],   // coefficients [1, ndiv[0], ndiv[0]*ndiv[1]] (Eq. 8)

    // auxiliary variable
    ratio: [usize; 3], // ratio = trunc(δx[i]/Δx[i]) (Eq. 8)

    // holds non-empty containers. maps container.index to container.data
    containers: HashMap<usize, Container>,
}

impl GridSearch {
    /// Creates a new GridSearch
    ///
    /// # Input
    ///
    /// * `xmin, xmax` -- min and max coordinates (len = 2 or 3 == ndim)
    /// * `ndiv` -- number of divisions along each dimension (len = 2 or 3 == ndim)
    pub fn new(xmin: &[f64], xmax: &[f64], ndiv: &[usize]) -> Result<Self, &'static str> {
        // check input
        let ndim = ndiv.len();
        if ndim < 2 || ndim > 3 {
            return Err("len(ndiv) == ndim must be 2 or 3");
        }
        if xmin.len() != ndim {
            return Err("len(xmin) must equal ndim == len(ndiv)");
        }
        if xmax.len() != ndim {
            return Err("len(xmax) must equal ndim == len(ndiv)");
        }

        // allocate gs
        let mut grid = GridSearch {
            ndim,
            ndiv: [0; 3],
            xmin: [0.0; 3],
            xmax: [0.0; 3],
            xdelta: [0.0; 3],
            size: [0.0; 3],
            ratio: [0; 3],
            cf: [0; 3],
            containers: HashMap::new(),
        };

        // check and compute sizes
        for i in 0..ndim {
            grid.ndiv[i] = ndiv[i];
            if ndiv[i] < 1 {
                return Err("ndiv must be at least equal to 1");
            }
            grid.xmin[i] = xmin[i];
            grid.xmax[i] = xmax[i];
            grid.xdelta[i] = xmax[i] - xmin[i];
            if grid.xdelta[i] <= 0.0 {
                return Err("xmax must be greater than xmin");
            }
            grid.size[i] = grid.xdelta[i] / (ndiv[i] as f64);
        }

        // coefficient
        grid.cf[0] = 1;
        grid.cf[1] = grid.ndiv[0];
        grid.cf[2] = grid.ndiv[0] * grid.ndiv[1];

        // done
        Ok(grid)
    }

    /// Inserts a new item to some container in the grid
    ///
    /// # Input
    ///
    /// * `id` -- identification number for the item
    /// * `x`, `y`, `z` -- coordinates of the item (`z` is ignored if ndim == 2)
    pub fn insert(&mut self, id: usize, x: f64, y: f64, z: f64) -> Result<(), &'static str> {
        match self.container_index(&[x, y, z]) {
            Some(index) => {
                // new item data
                let item = Item { id, x: [x, y, z] };

                // update container
                if self.containers.contains_key(&index) {
                    let container = self.containers.get_mut(&index).unwrap();
                    container.items.push(item);

                // initialize container
                } else {
                    self.containers.insert(index, Container { items: vec![item] });
                }
            }
            None => return Err("point is outside the grid"),
        }
        Ok(())
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
            xmin[i] = self.xmin[i];
            xmax[i] = self.xmax[i];
            ndiv[i] = self.ndiv[i];
        }
        let mut shapes = Shapes::new();
        shapes
            .set_alt_text_color("#5d5d5d")
            .draw_grid(&xmin, &xmax, &ndiv, false, true)?;
        plot.add(&shapes);

        // draw items
        if self.ndim == 2 {
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
                    curve.draw(&[item.x[0]], &[item.x[1]]);
                    text.draw(item.x[0], item.x[1], &txt);
                }
            }
            plot.add(&curve);
            plot.add(&text);
        }
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
            if x[i] < self.xmin[i] || x[i] > self.xmax[i] {
                return None;
            }
            self.ratio[i] = ((x[i] - self.xmin[i]) / self.size[i]) as usize;
            if self.ratio[i] == self.ndiv[i] {
                // the point is exactly on the max edge, thus select inner container
                self.ratio[i] -= 1; // move to the inside
            }
            index += self.ratio[i] * self.cf[i];
        }
        Some(index)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    struct TestData2D {
        id: usize,
        x: f64,
        y: f64,
        container_index: usize,
    }

    struct TestData3D {
        id: usize,
        x: f64,
        y: f64,
        z: f64,
        container_index: usize,
    }

    fn get_test_data_2d() -> Vec<TestData2D> {
        vec![
            TestData2D {
                id: 100,
                x: 0.0,
                y: 0.2,
                container_index: 6,
            },
            TestData2D {
                id: 200,
                x: 0.2,
                y: 0.6000000000000001,
                container_index: 12,
            },
            TestData2D {
                id: 300,
                x: 0.4,
                y: 1.0,
                // in this case, the ratio is:
                // 3.0000000000000004, 2.9999999999999996
                // thus, the point falls in #13 instead of #18
                container_index: 13,
            },
            TestData2D {
                id: 400,
                x: 0.6,
                y: 1.4,
                // in this case, the ratio is:
                // 4, 3.9999999999999996
                // thus, the point falls in #19 instead of #24
                container_index: 19,
            },
            TestData2D {
                id: 500,
                x: 0.8,
                y: 1.8,
                container_index: 24,
            },
        ]
    }

    fn get_test_data_3d() -> Vec<TestData3D> {
        vec![
            TestData3D {
                id: 100,
                x: -1.0,
                y: -1.0,
                z: -1.0,
                container_index: 0,
            },
            TestData3D {
                id: 200,
                x: 0.0,
                y: 0.0,
                z: 0.0,
                container_index: 7,
            },
            TestData3D {
                id: 300,
                x: 0.25,
                y: 0.25,
                z: 0.25,
                container_index: 7,
            },
        ]
    }

    #[test]
    fn new_fails_on_wrong_input() {
        let grid = GridSearch::new(&[0.0, 0.0], &[1.0, 1.0], &[1]);
        assert_eq!(grid.err(), Some("len(ndiv) == ndim must be 2 or 3"));
        let grid = GridSearch::new(&[0.0], &[1.0, 1.0], &[1, 1]);
        assert_eq!(grid.err(), Some("len(xmin) must equal ndim == len(ndiv)"));
        let grid = GridSearch::new(&[0.0, 0.0], &[1.0], &[1, 1]);
        assert_eq!(grid.err(), Some("len(xmax) must equal ndim == len(ndiv)"));
        let grid = GridSearch::new(&[0.0, 0.0], &[1.0, 1.0], &[0, 1]);
        assert_eq!(grid.err(), Some("ndiv must be at least equal to 1"));
        let grid = GridSearch::new(&[0.0, 0.0], &[0.0, 1.0], &[1, 1]);
        assert_eq!(grid.err(), Some("xmax must be greater than xmin"));
    }

    #[test]
    fn new_2d_works() -> Result<(), &'static str> {
        let grid = GridSearch::new(&[-0.2, -0.2], &[0.8, 1.8], &[5, 5])?;
        assert_eq!(grid.ndim, 2);
        assert_eq!(grid.ndiv, [5, 5, 0]);
        assert_eq!(grid.xmin, [-0.2, -0.2, 0.0]);
        assert_eq!(grid.xmax, [0.8, 1.8, 0.0]);
        assert_eq!(grid.xdelta, [1.0, 2.0, 0.0]);
        assert_eq!(grid.size, [0.2, 0.4, 0.0]);
        assert_eq!(grid.cf, [1, 5, 25]);
        assert_eq!(grid.ratio, [0, 0, 0]);
        assert_eq!(grid.containers.len(), 0);
        println!("{:.20}", grid.size[0]);
        assert_eq!(grid.size[0], 0.20000000000000001110); // TODO: check if this is a problem
        let b0 = 0.20000000000000001110_f64.to_bits();
        let b1 = 0.2_f64.to_bits();
        println!("{} =? {}", b0, b1);
        assert_eq!(b0, b1);
        Ok(())
    }

    #[test]
    fn new_3d_works() -> Result<(), &'static str> {
        let grid = GridSearch::new(&[-1.0, -1.0, -1.0], &[1.0, 1.0, 1.0], &[2, 2, 2])?;
        assert_eq!(grid.ndim, 3);
        assert_eq!(grid.ndiv, [2, 2, 2]);
        assert_eq!(grid.xmin, [-1.0, -1.0, -1.0]);
        assert_eq!(grid.xmax, [1.0, 1.0, 1.0]);
        assert_eq!(grid.xdelta, [2.0, 2.0, 2.0]);
        assert_eq!(grid.size, [1.0, 1.0, 1.0]);
        assert_eq!(grid.cf, [1, 2, 4]);
        assert_eq!(grid.ratio, [0, 0, 0]);
        assert_eq!(grid.containers.len(), 0);
        Ok(())
    }

    #[test]
    fn container_index_2d_works() -> Result<(), &'static str> {
        let mut grid = GridSearch::new(&[-0.2, -0.2], &[0.8, 1.8], &[5, 5])?;
        for data in get_test_data_2d() {
            let index = grid.container_index(&[data.x, data.y]).unwrap();
            assert_eq!(index, data.container_index);
        }
        let index = grid.container_index(&[0.80001, 0.0]);
        assert_eq!(index, None); // outside
        Ok(())
    }

    #[test]
    fn container_index_3d_works() -> Result<(), &'static str> {
        let mut grid = GridSearch::new(&[-1.0, -1.0, -1.0], &[1.0, 1.0, 1.0], &[2, 2, 2])?;
        for data in get_test_data_3d() {
            let index = grid.container_index(&[data.x, data.y, data.z]).unwrap();
            assert_eq!(index, data.container_index);
        }
        let index = grid.container_index(&[1.00001, 0.0, 0.0]);
        assert_eq!(index, None); // outside
        Ok(())
    }

    #[test]
    fn insert_2d_works() -> Result<(), &'static str> {
        let mut grid = GridSearch::new(&[-0.2, -0.2], &[0.8, 1.8], &[5, 5])?;
        for data in get_test_data_2d() {
            grid.insert(data.id, data.x, data.y, 0.0)?;
            let container = grid.containers.get(&data.container_index).unwrap();
            let item = container.items.iter().find(|item| item.id == data.id).unwrap();
            assert_eq!(item.id, data.id);
        }
        assert_eq!(grid.containers.len(), 5);
        let res = grid.insert(1000, 0.80001, 0.0, 0.0);
        assert_eq!(res, Err("point is outside the grid"));
        Ok(())
    }

    #[test]
    fn insert_3d_works() -> Result<(), &'static str> {
        let mut grid = GridSearch::new(&[-1.0, -1.0, -1.0], &[1.0, 1.0, 1.0], &[2, 2, 2])?;
        for data in get_test_data_3d() {
            grid.insert(data.id, data.x, data.y, data.z)?;
            let container = grid.containers.get(&data.container_index).unwrap();
            let item = container.items.iter().find(|item| item.id == data.id).unwrap();
            assert_eq!(item.id, data.id);
        }
        assert_eq!(grid.containers.len(), 2);
        let res = grid.insert(1000, 1.00001, 0.0, 0.0);
        assert_eq!(res, Err("point is outside the grid"));
        Ok(())
    }

    #[test]
    fn plot_2d_works() -> Result<(), &'static str> {
        let mut grid = GridSearch::new(&[-0.2, -0.2], &[0.8, 1.8], &[5, 5])?;
        for data in get_test_data_2d() {
            grid.insert(data.id, data.x, data.y, 0.0)?;
        }
        let plot = grid.plot()?;
        plot.save("/tmp/gemlab/search_grid_plot_2d_works.svg")?;
        Ok(())
    }

    #[test]
    fn plot_3d_works() -> Result<(), &'static str> {
        let mut grid = GridSearch::new(&[-1.0, -1.0, -1.0], &[1.0, 1.0, 1.0], &[2, 2, 2])?;
        for data in get_test_data_3d() {
            grid.insert(data.id, data.x, data.y, 0.0)?;
        }
        let plot = grid.plot()?;
        plot.save("/tmp/gemlab/search_grid_plot_3d_works.svg")?;
        Ok(())
    }
}
