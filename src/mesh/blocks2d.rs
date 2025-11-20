use crate::mesh::Constraint2d;
use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds a collection of 2D blocks for structured mesh generation.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Blocks2d {
    /// List of points
    ///
    /// Each point is defined by x and y coordinates.
    ///
    /// (length = number of points)
    pub points: Vec<(f64, f64)>,

    /// List of regions defining each block
    ///
    /// Each block is defined by `(marker, p1, p2, p3, p4)` where
    /// `marker` is an integer to identify a group of blocks and
    /// `p1`, `p2`, `p3`, `p4` are the indices of the corner points of the block.
    ///
    /// (length = number of blocks)
    pub regions: Vec<(i32, usize, usize, usize, usize)>,

    /// List of division weights for each block
    ///
    /// Each entry is defined as `(weights_x, weights_y)` where
    /// `weights_x` is a vector of weights for divisions in the x direction and
    /// `weights_y` is a vector of weights for divisions in the y direction.
    /// The number of divisions along each direction is determined by the length of the respective weight vector.
    /// Only vectors with length > 0 are considered, otherwise default uniform divisions are used.
    ///
    /// Note: Both `weights_x` and `weights_y` are required to be > 0 when setting the number of divisions/weights.
    ///
    /// (length = number of blocks)
    pub div_weights: Vec<(Vec<f64>, Vec<f64>)>,

    /// List of edge constraints
    ///
    /// Defined as `(local_edge_index, constraint)` where
    /// `local_edge_index` is the index of the edge in the block (0 to 3),
    /// and `constraint` is the `Constraint2D` applied to that edge
    ///
    /// (length = number of blocks)
    pub edge_constraints: Vec<Option<(usize, Constraint2d)>>,

    /// Holds all marked edges
    ///
    /// Each entry contains `(marker, p1, p2)`, where `marker` is the edge marker,
    /// and `p1` and `p2` are the point ids. The point ids may be unsorted.
    pub marked_edges: Vec<(i32, usize, usize)>,
}

impl Blocks2d {
    /// Reads a JSON file containing the data
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn read_json<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let input = File::open(path).map_err(|_| "cannot open file")?;
        let buffered = BufReader::new(input);
        let mesh = serde_json::from_reader(buffered).map_err(|_| "cannot parse JSON file")?;
        Ok(mesh)
    }

    /// Writes a JSON file with the data
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn write_json<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }
        let mut file = File::create(&path).map_err(|_| "cannot create file")?;
        serde_json::to_writer(&mut file, &self).map_err(|_| "cannot write file")?;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Blocks2d;
    use crate::mesh::Constraint2d;

    #[test]
    fn test_blocks2d_write_works() {
        //           -100       -200
        // 1.0  3-----------2-----------5
        //      |           |           |
        //      |           |           |
        //      |           |           |
        //      |           |           |
        // 0.0  0-----------1-----------4
        //     0.0         1.0         2.0
        let blocks = Blocks2d {
            points: vec![
                (0.0, 0.0), // 0
                (1.0, 0.0), // 1
                (1.0, 1.0), // 2
                (0.0, 1.0), // 3
                (2.0, 0.0), // 4
                (2.0, 1.0), // 5
            ],
            regions: vec![
                (1, 0, 1, 2, 3), // marker, p1, p2, p3, p4
                (2, 1, 4, 5, 2),
            ],
            div_weights: vec![
                (Vec::new(), Vec::new()), // default number of divisions and weighting
                (vec![1.0], vec![1.0]),   // one division with weight 1.0 along each direction
            ],
            edge_constraints: vec![
                None,                                           // block 0
                Some((0, Constraint2d::Circle(0.0, 0.0, 1.0))), // block 1
            ],
            marked_edges: vec![(-100, 3, 2), (-200, 2, 5)],
        };

        assert_eq!(blocks.points.len(), 6);
        assert_eq!(blocks.regions.len(), 2);
        assert_eq!(blocks.div_weights.len(), 2);
        assert_eq!(blocks.edge_constraints.len(), 2);

        blocks.write_json("/tmp/gemlab/test_blocks2d_write_works.json").unwrap();

        let loaded_blocks = Blocks2d::read_json("/tmp/gemlab/test_blocks2d_write_works.json").unwrap();
        assert_eq!(&loaded_blocks.points, &blocks.points);
        assert_eq!(&loaded_blocks.regions, &blocks.regions);
        assert_eq!(&loaded_blocks.div_weights, &blocks.div_weights);
        assert_eq!(&loaded_blocks.edge_constraints.len(), &blocks.edge_constraints.len());
        for i in 0..blocks.edge_constraints.len() {
            if let Some(c) = &blocks.edge_constraints[i] {
                assert_eq!(c.0, loaded_blocks.edge_constraints[i].as_ref().unwrap().0);
            }
        }
    }
}
