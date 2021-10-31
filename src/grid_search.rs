use std::collections::HashMap;

/// Holds the id and coordinates of an item
#[derive(Clone)]
struct Item {
    id: usize,       // identification number
    coord: Vec<f64>, // coordinates
}

/// Holds items
#[derive(Clone)]
struct Container {
    items: Vec<Item>,
}

/// Implements a grid for fast searching entries by coordinates
pub struct GridSearch {
    ndim: usize,      // space dimension
    ndiv: [usize; 3], // number of divisions along each direction
    xmin: [f64; 3],   // min values
    xmax: [f64; 3],   // max values
    xdelta: [f64; 3], // difference between max and min
    size: [f64; 3],   // size of each cell

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
            if grid.xdelta[i] < 0.0 {
                return Err("xmax must be greater than xmin");
            }
            grid.size[i] = grid.xdelta[i] / (ndiv[i] as f64);
        }

        // done
        Ok(grid)
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
        assert_eq!(grid.cells.len(), 0);
        println!("{:.20}", grid.size[0]);
        Ok(())
    }
}
