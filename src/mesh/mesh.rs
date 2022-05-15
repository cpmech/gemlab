use super::read_text_mesh::{parse_text_mesh, read_text_mesh};
use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::Path;

/// Aliases usize as Point ID
pub type PointId = usize;

/// Aliases usize as Cell ID
pub type CellId = usize;

/// Aliases usize as Cell's attribute ID
pub type CellAttributeId = usize;

/// Holds point data
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Point {
    /// Identification number which equals the index of the point in the mesh
    pub id: PointId,

    /// Point coordinates (2D or 3D)
    pub coords: Vec<f64>,
}

/// Holds cell (aka geometric shape, polygon, polyhedra) data
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Cell {
    /// Identification number which equals the index of the cell in the mesh
    pub id: CellId,

    /// Attribute identification number
    pub attribute_id: CellAttributeId,

    /// Space dimension of this cell
    ///
    /// The cell's ndim may be different than the space dimension of the mesh.
    /// For example, a 1D line in the 2D or 3D space or a 2D triangle in the 3D space.
    pub geo_ndim: usize,

    /// List of points defining this cell (nodes); in the right order (unsorted)
    ///
    /// Note: The list of nodes must follow a **counter-clockwise order**.
    pub points: Vec<PointId>,
}

/// Holds mesh data
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Mesh {
    /// Space dimension of the mesh
    ///
    /// The mesh's ndim may be different that an cell's ndim.
    /// For example, a 3D mesh may contain 1D lines or 2D triangles.
    pub space_ndim: usize,

    /// All points in the mesh
    pub points: Vec<Point>,

    /// All cells (aka geometric shape, polygon, polyhedra) in the mesh
    pub cells: Vec<Cell>,
}

impl Mesh {
    /// Returns a new empty mesh
    pub(super) fn new(space_ndim: usize) -> Result<Self, StrError> {
        if space_ndim < 2 || space_ndim > 3 {
            return Err("space_ndim must be 2 or 3");
        }
        Ok(Mesh {
            space_ndim,
            points: Vec::new(),
            cells: Vec::new(),
        })
    }

    /// Parses raw mesh data from a text string
    ///
    /// # Note
    ///
    /// This function calls `compute_derived_props` already.
    pub fn from_text(raw_mesh_data: &str) -> Result<Self, StrError> {
        parse_text_mesh(raw_mesh_data)
    }

    /// Reads raw mesh data from text file
    ///
    /// # Note
    ///
    /// This function calls `compute_derived_props` already.
    pub fn from_text_file<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        read_text_mesh(full_path)
    }

    /// Reads a binary file containing the mesh data
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn read<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let mut file = File::open(&path).map_err(|_| "file not found")?;
        let metadata = fs::metadata(&path).map_err(|_| "unable to read metadata")?;
        let mut bin = vec![0; metadata.len() as usize];
        file.read(&mut bin).expect("buffer overflow");
        let mut des = rmp_serde::Deserializer::new(&bin[..]);
        let mesh: Mesh = Deserialize::deserialize(&mut des).map_err(|_| "deserialize failed")?;
        Ok(mesh)
    }

    /// Writes a binary file with the mesh data
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn write<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }
        let mut bin = Vec::new();
        let mut ser = rmp_serde::Serializer::new(&mut bin);
        self.serialize(&mut ser).map_err(|_| "serialize failed")?;
        let mut file = File::create(&path).map_err(|_| "cannot create file")?;
        file.write_all(&bin).map_err(|_| "cannot write file")?;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::mesh::Mesh;
    use crate::StrError;

    #[test]
    fn new_fails_on_wrong_input() {
        assert_eq!(Mesh::new(1).err(), Some("space_ndim must be 2 or 3"));
        assert_eq!(Mesh::new(4).err(), Some("space_ndim must be 2 or 3"));
    }

    #[test]
    fn new_works() -> Result<(), StrError> {
        let mesh = Mesh::new(2)?;
        assert_eq!(mesh.space_ndim, 2);
        assert_eq!(mesh.points.len(), 0);
        assert_eq!(mesh.cells.len(), 0);
        Ok(())
    }

    #[test]
    fn from_text_fails_on_invalid_data() -> Result<(), StrError> {
        assert_eq!(
            Mesh::from_text("").err(),
            Some("text string is empty or header is missing")
        );
        Ok(())
    }

    #[test]
    fn from_text_file_fails_on_invalid_data() -> Result<(), StrError> {
        assert_eq!(Mesh::from_text_file("").err(), Some("cannot open file"));
        Ok(())
    }

    #[test]
    fn from_text_fails_on_wrong_jacobian_2d() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let res = Mesh::from_text(
            r"# header
            # space_ndim npoint ncell
                       2      6     2
            
            # points
            # id   x   y
               0 0.0 0.0
               1 1.0 0.0
               2 1.0 1.0
               3 0.0 1.0
               4 2.0 0.0
               5 2.0 1.0
            
            # cells
            # id att geo_ndim nnode  (wrong) point_ids...
               0   1        2     4  0 3 2 1
               1   0        2     4  1 2 5 4",
        );
        assert_eq!(res.err(), Some("a cell has incorrect ordering of nodes"));
        Ok(())
    }

    #[test]
    fn from_text_works() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let mesh = Mesh::from_text(
            r"# header
            # space_ndim npoint ncell
                       2      6     2
            
            # points
            # id   x   y
               0 0.0 0.0
               1 1.0 0.0
               2 1.0 1.0
               3 0.0 1.0
               4 2.0 0.0
               5 2.0 1.0
            
            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        2     4  0 1 2 3
               1   0        2     4  1 4 5 2",
        )?;
        assert_eq!(mesh.space_ndim, 2);
        assert_eq!(mesh.points.len(), 6);
        assert_eq!(mesh.cells.len(), 2);
        Ok(())
    }

    #[test]
    fn from_text_file_fails_on_wrong_jacobian_3d() -> Result<(), StrError> {
        let res = Mesh::from_text_file("./data/meshes/bad_wrong_jacobian.msh");
        assert_eq!(res.err(), Some("a cell has incorrect ordering of nodes"));
        Ok(())
    }

    #[test]
    fn from_text_file_fails_on_wrong_nodes_3d() -> Result<(), StrError> {
        let res = Mesh::from_text_file("./data/meshes/bad_wrong_nodes.msh");
        assert_eq!(res.err(), Some("cannot compute inverse due to zero determinant"));
        Ok(())
    }

    #[test]
    fn from_text_file_works() -> Result<(), StrError> {
        //
        //       8-------------11
        //      /.             /|
        //     / .            / |
        //    /  .           /  |
        //   /   .          /   |
        //  9-------------10    |
        //  |    .         |    |
        //  |    4---------|----7
        //  |   /.         |   /|
        //  |  / .         |  / |
        //  | /  .         | /  |
        //  |/   .         |/   |
        //  5--------------6    |
        //  |    .         |    |
        //  |    0---------|----3
        //  |   /          |   /
        //  |  /           |  /
        //  | /            | /
        //  |/             |/
        //  1--------------2
        //
        let mesh = Mesh::from_text_file("./data/meshes/ok2.msh")?;
        assert_eq!(mesh.space_ndim, 3);
        assert_eq!(mesh.points.len(), 12);
        assert_eq!(mesh.cells.len(), 2);
        Ok(())
    }

    #[test]
    fn from_text_file_works_with_mixed() -> Result<(), StrError> {
        //
        //          4--------3
        //          |        |
        //          |        |
        //          |        |
        //  0-------1--------2
        //
        let mesh = Mesh::from_text_file("./data/meshes/ok_mixed_shapes2.msh")?;

        //
        //                       4------------7-----------10
        //                      /.           /|            |
        //                     / .          / |            |
        //                    /  .         /  |            |
        //                   /   .        /   |            |
        //                  5------------6    |            |
        //                  |    .       |`.  |            |
        //                  |    0-------|--`.3------------9
        //                  |   /        |   /`.          /
        //                  |  /         |  /   `.       /
        //                  | /          | /      `.    /
        //                  |/           |/         `. /
        //  12-----11-------1------------2------------8
        //
        let mesh = Mesh::from_text_file("./data/meshes/ok_mixed_shapes.msh")?;
        Ok(())
    }

    #[test]
    fn read_write_capture_errors() -> Result<(), StrError> {
        assert_eq!(Mesh::read("/tmp/not_found").err(), Some("file not found"));
        assert_eq!(
            Mesh::read("./data/meshes/ok_mixed_shapes2.msh").err(),
            Some("deserialize failed")
        );
        let mesh = Mesh::new(2)?;
        assert_eq!(mesh.write("/tmp/").err(), Some("cannot create file"));
        Ok(())
    }

    #[test]
    fn read_write_work() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let mesh = Mesh::from_text(
            r"# header
            # space_ndim npoint ncell
                       2      6     2
            
            # points
            # id   x   y
               0 0.0 0.0
               1 1.0 0.0
               2 1.0 1.0
               3 0.0 1.0
               4 2.0 0.0
               5 2.0 1.0
            
            # cells
            # id att geo_ndim nnode  point_ids...
               0   1        2     4  0 1 2 3
               1   0        2     4  1 4 5 2",
        )?;
        mesh.write("/tmp/gemlab/test.msh")?;
        let mesh_read = Mesh::read("/tmp/gemlab/test.msh")?;
        assert_eq!(format!("{:?}", mesh), format!("{:?}", mesh_read));
        Ok(())
    }
}
