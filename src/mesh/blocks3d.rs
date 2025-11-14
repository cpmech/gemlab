use crate::mesh::Constraint3d;
use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds a collection of 3d blocks for structured mesh generation.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Blocks3d {
    /// List of points
    ///
    /// Each point is defined by x, y, and z coordinates.
    ///
    /// (length = number of points)
    pub points: Vec<(f64, f64, f64)>,

    /// List of regions defining each block
    ///
    /// Each block is defined by `(attribute, p1, p2, p3, p4, p5, p6, p7, p8)` where
    /// `attribute` is an integer attribute for the block and
    /// `p1`, `p2`, `p3`, `p4`, `p5`, `p6`, `p7`, `p8` are the indices of the corner points of the block.
    ///
    /// (length = number of blocks)
    pub regions: Vec<(i32, usize, usize, usize, usize, usize, usize, usize, usize)>,

    /// List of division weights for each block
    ///
    /// Each entry is defined as `(weights_x, weights_y, weights_z)` where
    /// `weights_x` is a vector of weights for divisions in the x direction,
    /// `weights_y` is a vector of weights for divisions in the y direction, and
    /// `weights_z` is a vector of weights for divisions in the z direction.
    /// The number of divisions along each direction is determined by the length of the respective weight vector.
    /// Only vectors with length > 0 are considered, otherwise default uniform divisions are used.
    ///
    /// Note: All `weights_x`, `weights_y`, and `weights_z` are required to be > 0 when setting the number of divisions/weights.
    ///
    /// (length = number of blocks)
    pub div_weights: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)>,

    /// List of face constraints
    ///
    /// Defined as `(local_face_index, constraint)` where
    /// `local_face_index` is the index of the face in the block (0 to 5),
    /// and `constraint` is the `Constraint3D` applied to that face
    ///
    /// (length = number of blocks)
    pub face_constraints: Vec<Option<(usize, Constraint3d)>>,

    /// Holds all marked edges
    ///
    /// Each entry contains `(marker, p1, p2)`, where `marker` is the edge marker,
    /// and `p1` and `p2` are the point ids. The point ids may be unsorted.
    pub marked_edges: Vec<(i32, usize, usize)>,

    /// Holds all marked faces
    ///
    /// Each entry contains `(marker, p1, p2, p3, p4)`, where `marker` is the face marker,
    /// and `p1`, `p2`, `p3`, and `p4` are the point ids. The point ids may be unsorted.
    /// For triangular faces, the fourth point (p4) must be set to [usize::MAX].
    pub marked_faces: Vec<(i32, usize, usize, usize, usize)>,
}

impl Blocks3d {
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
    use super::Blocks3d;
    use crate::mesh::Constraint3d;

    #[test]
    fn test_blocks3d_write_works() {
        let r = 0.866;
        let w = 2.0;
        let h = 1.0;
        let l = h / 2.0;
        let m = f64::sqrt(r * r - l * l);
        let blocks = Blocks3d {
            points: vec![
                (0.0, 0.0, 0.0), //  0
                (1.0, 0.0, 0.0), //  1
                (1.0, 1.0, 0.0), //  2
                (0.0, 1.0, 0.0), //  3
                (2.0, 0.0, 0.0), //  4
                (2.0, 1.0, 0.0), //  5
                (0.0, 0.0, 1.0), //  6
                (1.0, 0.0, 1.0), //  7
                (1.0, 1.0, 1.0), //  8
                (0.0, 1.0, 1.0), //  9
                (2.0, 0.0, 1.0), // 10
                (2.0, 1.0, 1.0), // 11
            ],
            regions: vec![
                (1, 0, 1, 2, 3, 6, 7, 8, 9), // attribute, p1, p2, p3, p4, p5, p6, p7, p8
                (2, 1, 4, 5, 2, 7, 10, 11, 8),
            ],
            div_weights: vec![
                (Vec::new(), Vec::new(), Vec::new()), // default number of divisions and weighting
                (vec![1.0], vec![1.0], vec![1.0]),    // one division with weight 1.0 along each direction
            ],
            face_constraints: vec![
                None,                                            // block 0
                Some((0, Constraint3d::CylinderZ(w + m, l, r))), // block 1
            ],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };

        assert_eq!(blocks.points.len(), 12);
        assert_eq!(blocks.regions.len(), 2);
        assert_eq!(blocks.div_weights.len(), 2);
        assert_eq!(blocks.face_constraints.len(), 2);

        blocks.write_json("/tmp/gemlab/test_blocks3d_write_works.json").unwrap();

        let loaded_blocks = Blocks3d::read_json("/tmp/gemlab/test_blocks3d_write_works.json").unwrap();
        assert_eq!(&loaded_blocks.points, &blocks.points);
        assert_eq!(&loaded_blocks.regions, &blocks.regions);
        assert_eq!(&loaded_blocks.div_weights, &blocks.div_weights);
        assert_eq!(&loaded_blocks.face_constraints.len(), &blocks.face_constraints.len());
        for i in 0..blocks.face_constraints.len() {
            if let Some(c) = &blocks.face_constraints[i] {
                assert_eq!(c.0, loaded_blocks.face_constraints[i].as_ref().unwrap().0);
            }
        }
    }
}
