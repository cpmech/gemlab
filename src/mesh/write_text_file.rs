use super::Mesh;
use crate::StrError;
use std::ffi::OsStr;
use std::fmt::Write;
use std::fs::{self, File};
use std::io::Write as IoWrite;
use std::path::Path;

impl Mesh {
    /// Writes a text file with the mesh description
    ///
    /// # File format
    ///
    /// The text file format includes three sections:
    ///
    /// 1. The header with the space dimension (`ndim`), number of points (`npoint`), and number of cells (`ncell`);
    /// 2. The points list where each line contains the `id` of the point, which must be **equal to the position** in the list,
    ///    followed by the `x` and `y` (and `z`) coordinates;
    /// 3. The cells list where each line contains the `id` of the cell, which must be **equal to the position** in the list,
    ///    the attribute (`att`) of the cell, the `kind` of the cell, followed by the IDs of the points that define the cell (connectivity).
    ///
    /// The text file looks like this (the hash tag indicates a comment/the mesh below is just an example which won't work):
    ///
    /// ```text
    /// # header
    /// # ndim npoint ncell
    ///      2      8     5
    ///
    /// # points
    /// # id marker x y
    ///    0 0 0.0 0.0
    ///    1 0 0.5 0.0
    ///    2 0 1.0 0.0
    /// # ... more points should follow
    ///
    /// # cells
    /// # id attribute kind point_ids...
    ///    0 1 tri3 0 1 3
    ///    1 1 qua4 1 4 6 3
    /// ```
    ///
    /// where we can see that different cell (shape) kinds can be present in the same mesh.
    /// However, this function does not check for element compatibility as required by finite element analyses.
    ///
    /// See [crate::shapes::GeoKind::from] for the keys used to identify the cell kind (the keys are **lowercase**).
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn write_text_file<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let ncell = self.cells.len();
        if ncell < 1 {
            return Err("there are no cells to write");
        }

        // output buffer
        let mut buffer = String::new();
        write!(&mut buffer, "{}", self).map_err(|_| "cannot write to buffer")?;

        // create directory
        let path = Path::new(full_path);
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }

        // write file
        let mut file = File::create(path).map_err(|_| "cannot create file")?;
        file.write_all(buffer.as_bytes()).map_err(|_| "cannot write file")?;

        // force sync
        file.sync_all().map_err(|_| "cannot sync file")?;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::{mesh::Samples, StrError};
    use std::fs;

    #[test]
    fn write_text_file_handles_unavailable_types() {
        assert_eq!(
            Samples::lin_cells().write_vtu("/tmp/gemlab/nothing.vtu").err(),
            Some("cannot generate VTU file because VTK cell type is not available")
        );
        assert_eq!(
            Samples::tri_cells().write_vtu("/tmp/gemlab/nothing.vtu").err(),
            Some("cannot generate VTU file because VTK cell type is not available")
        );
        assert_eq!(
            Samples::qua_cells().write_vtu("/tmp/gemlab/nothing.vtu").err(),
            Some("cannot generate VTU file because VTK cell type is not available")
        );
    }

    #[test]
    fn write_text_file_works_qua8_tri6_lin2() -> Result<(), StrError> {
        let mesh = Samples::qua8_tri6_lin2();
        let file_path = "/tmp/gemlab/test_qua8_tri6_lin2.msh";
        mesh.write_text_file(file_path)?;
        let contents = fs::read_to_string(file_path).map_err(|_| "cannot open file")?;
        assert_eq!(
            contents,
            "# header\n\
             # ndim npoint ncell\n\
             2 11 4\n\
             \n\
             # points\n\
             # id marker x y {z}\n\
             0 -100 0.0 0.0\n\
             1 0 0.5 0.0\n\
             2 -200 1.0 0.0\n\
             3 0 1.433 0.25\n\
             4 -300 1.866 0.5\n\
             5 0 1.433 0.75\n\
             6 -400 1.0 1.0\n\
             7 0 0.5 1.0\n\
             8 -500 0.0 1.0\n\
             9 0 0.0 0.5\n\
             10 0 1.0 0.5\n\
             \n\
             # cells\n\
             # id attribute kind points\n\
             0 1 qua8 0 2 6 8 1 10 7 9\n\
             1 2 tri6 2 4 6 3 5 10\n\
             2 3 lin2 2 10\n\
             3 3 lin2 10 6\n"
        );
        Ok(())
    }
}
