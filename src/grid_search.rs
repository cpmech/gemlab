use crate::AsArray1D;
use std::collections::HashMap;

/// Holds the id and coordinates of an item
#[derive(Clone)]
struct Item {
    id: usize,   // identification number
    x: [f64; 3], // coordinates
}

/// Holds items
#[derive(Clone)]
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
    size: [f64; 3],   // size of each cell
    cf: [usize; 3],   // coefficients [1, ndiv[0], ndiv[0]*ndiv[1]] (Eq. 8)

    // auxiliary variable
    ratio: [usize; 3], // ratio = trunc(δx[i]/Δx[i]) (Eq. 8)

    // holds non-empty cells. maps cell.id to container
    cells: HashMap<usize, Container>,
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
            return Err("ndim must be 2 or 3");
        }
        if xmin.len() != ndim {
            return Err("size of xmin must equal ndim");
        }
        if xmax.len() != ndim {
            return Err("size of xmax must equal ndim");
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
            cells: HashMap::new(),
        };

        // check and compute sizes
        for i in 0..ndim {
            grid.ndiv[i] = ndiv[i];
            if ndiv[i] < 1 {
                return Err("ndiv must be at least 1");
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
    /// * `x` -- coordinates of the item (an 1D array with len = ndim)
    pub fn insert<'a, T, U>(&mut self, id: usize, x: &'a T) -> Result<(), &'static str>
    where
        T: AsArray1D<'a, U>,
        U: 'a + Into<f64>,
    {
        if x.size() != self.ndim {
            return Err("x must have len equal to ndim");
        }
        match self.container_index(x) {
            Some(index) => {
                // new item data
                let mut item = Item {
                    id,
                    x: [x.at(0).into(), x.at(1).into(), 0.0],
                };
                if self.ndim == 3 {
                    item.x[2] = x.at(2).into();
                }

                // update container
                if self.cells.contains_key(&index) {
                    let cell = self.cells.get_mut(&index).unwrap();
                    cell.items.push(item);

                // initialize container
                } else {
                    self.cells.insert(index, Container { items: vec![item] });
                }
            }
            None => return Err("point is outside the grid"),
        }
        Ok(())
    }

    /// Draws grid
    pub fn draw(&self) {}

    /// Calculates the container index where the point x should be located
    ///
    /// # Output
    ///
    /// * returns the index of the container or None if the point is out-of-range
    fn container_index<'a, T, U>(&mut self, x: &'a T) -> Option<usize>
    where
        T: AsArray1D<'a, U>,
        U: 'a + Into<f64>,
    {
        let mut index = 0;
        for i in 0..self.ndim {
            if x.at(i).into() < self.xmin[i] || x.at(i).into() > self.xmax[i] {
                return None;
            }
            self.ratio[i] = ((x.at(i).into() - self.xmin[i]) / self.size[i]) as usize;
            if self.ratio[i] == self.ndiv[i] {
                // the point is exactly on the max edge, thus select inner bin
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

    #[test]
    fn new_works() -> Result<(), &'static str> {
        let grid = GridSearch::new(&[-0.2, -0.2], &[0.8, 1.8], &[5, 5])?;
        assert_eq!(grid.ndim, 2);
        assert_eq!(grid.ndiv, [5, 5, 0]);
        assert_eq!(grid.xmin, [-0.2, -0.2, 0.0]);
        assert_eq!(grid.xmax, [0.8, 1.8, 0.0]);
        assert_eq!(grid.xdelta, [1.0, 2.0, 0.0]);
        assert_eq!(grid.size, [0.2, 0.4, 0.0]);
        assert_eq!(grid.cf, [1, 5, 25]);
        assert_eq!(grid.ratio, [0, 0, 0]);
        assert_eq!(grid.cells.len(), 0);
        println!("{:.20}", grid.size[0]);
        Ok(())
    }

    #[test]
    fn container_index_works() -> Result<(), &'static str> {
        let mut grid = GridSearch::new(&[-0.2, -0.2], &[0.8, 1.8], &[5, 5])?;
        let index = grid.container_index(&[0.4, 1.0]).unwrap();
        assert_eq!(index, 13);
        Ok(())
    }

    #[test]
    fn insert_works() -> Result<(), &'static str> {
        let mut grid = GridSearch::new(&[-0.2, -0.2], &[0.8, 1.8], &[5, 5])?;
        grid.insert(33, &[0.4, 1.0])?;
        Ok(())
    }
}
