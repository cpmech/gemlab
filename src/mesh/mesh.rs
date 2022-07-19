use crate::shapes::GeoKind;
use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fmt::{self, Write as FmtWrite};
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

    /// The kind of cell
    pub kind: GeoKind,

    /// List of points defining this cell (nodes); in the right order (unsorted)
    ///
    /// Note: The list of nodes must follow a **counter-clockwise order**.
    pub points: Vec<PointId>,
}

/// Holds mesh data
///
/// # Examples
///
/// ```
/// use gemlab::mesh::{Cell, Mesh, Point};
/// use gemlab::shapes::GeoKind;
///
/// //          [#] indicates id
/// //      y   (#) indicates attribute_id
/// //      ↑
/// // 1.0  3-----------2-----------5
/// //      |           |           |
/// //      |    [0]    |    [1]    |
/// //      |    (1)    |    (2)    |
/// //      |           |           |
/// // 0.0  0-----------1-----------4  → x
/// //     0.0         1.0         2.0
///
/// let mesh = Mesh {
///     ndim: 2,
///     points: vec![
///         Point { id: 0, coords: vec![0.0, 0.0] },
///         Point { id: 1, coords: vec![1.0, 0.0] },
///         Point { id: 2, coords: vec![1.0, 1.0] },
///         Point { id: 3, coords: vec![0.0, 1.0] },
///         Point { id: 4, coords: vec![2.0, 0.0] },
///         Point { id: 5, coords: vec![2.0, 1.0] },
///     ],
///     cells: vec![
///         Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
///         Cell { id: 1, attribute_id: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
///     ],
/// };
/// ```
///
/// See more examples in [super::Samples]
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Mesh {
    /// Space dimension of the mesh
    ///
    /// The mesh's ndim may be different that an cell's ndim.
    /// For example, a 3D mesh may contain 1D lines or 2D triangles.
    pub ndim: usize,

    /// All points in the mesh
    pub points: Vec<Point>,

    /// All cells (aka geometric shape, polygon, polyhedra) in the mesh
    pub cells: Vec<Cell>,
}

impl Mesh {
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

impl fmt::Display for Mesh {
    /// Returns a text representation of the Mesh (can be used with [Mesh::from_text])
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "# header\n").unwrap();
        write!(f, "# ndim npoint ncell\n").unwrap();
        write!(f, "{} {} {}\n", self.ndim, self.points.len(), self.cells.len()).unwrap();
        write!(f, "\n# points\n").unwrap();
        write!(f, "# id x y {{z}}\n").unwrap();
        self.points.iter().for_each(|point| {
            if self.ndim == 2 {
                write!(f, "{} {} {}\n", point.id, point.coords[0], point.coords[1]).unwrap();
            } else {
                write!(
                    f,
                    "{} {} {} {}\n",
                    point.id, point.coords[0], point.coords[1], point.coords[2]
                )
                .unwrap();
            }
        });
        write!(f, "\n# cells\n").unwrap();
        write!(f, "# id att kind  points\n").unwrap();
        self.cells.iter().for_each(|cell| {
            write!(
                f,
                "{} {} {} {}\n",
                cell.id,
                cell.attribute_id,
                cell.kind.to_string(),
                cell.points.iter().fold(&mut String::new(), |acc, cur| {
                    write!(acc, " {}", cur).unwrap();
                    acc
                }),
            )
            .unwrap();
        });
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::mesh::{Mesh, Samples};
    use crate::StrError;
    use serde_json;

    #[test]
    fn read_and_write_capture_errors() -> Result<(), StrError> {
        assert_eq!(Mesh::read("/tmp/not_found").err(), Some("file not found"));
        assert_eq!(
            Mesh::read("./data/meshes/two_quads_horizontal.msh").err(),
            Some("deserialize failed")
        );
        let mesh = Mesh {
            ndim: 2,
            points: Vec::new(),
            cells: Vec::new(),
        };
        assert_eq!(mesh.write("/tmp/").err(), Some("cannot create file"));
        Ok(())
    }

    #[test]
    fn read_and_write_work() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let mesh = Mesh::from_text(
            r"# header
            # ndim npoint ncell
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
            # id att kind point_ids...
               0   1 qua4 0 1 2 3
               1   0 qua4 1 4 5 2",
        )?;
        mesh.write("/tmp/gemlab/test.msh")?;
        let mesh_read = Mesh::read("/tmp/gemlab/test.msh")?;
        assert_eq!(format!("{:?}", mesh), format!("{:?}", mesh_read));
        Ok(())
    }

    #[test]
    fn derive_works() -> Result<(), StrError> {
        let mesh = Samples::two_qua4();
        let mesh_clone = mesh.clone();
        let correct ="Mesh { ndim: 2, points: [Point { id: 0, coords: [0.0, 0.0] }, Point { id: 1, coords: [1.0, 0.0] }, Point { id: 2, coords: [1.0, 1.0] }, Point { id: 3, coords: [0.0, 1.0] }, Point { id: 4, coords: [2.0, 0.0] }, Point { id: 5, coords: [2.0, 1.0] }], cells: [Cell { id: 0, attribute_id: 1, kind: Qua4, points: [0, 1, 2, 3] }, Cell { id: 1, attribute_id: 2, kind: Qua4, points: [1, 4, 5, 2] }] }";
        assert_eq!(format!("{:?}", mesh), correct);
        assert_eq!(mesh_clone.ndim, mesh.ndim);
        assert_eq!(mesh_clone.points.len(), mesh.points.len());
        assert_eq!(mesh_clone.cells.len(), mesh.cells.len());
        // serialize
        let mesh_json = serde_json::to_string(&mesh).unwrap();
        // deserialize
        let mesh_read: Mesh = serde_json::from_str(&mesh_json).unwrap();
        assert_eq!(format!("{:?}", mesh_read), correct);
        Ok(())
    }

    #[test]
    fn display_works_2d() {
        let mesh = Samples::two_qua4();
        let text = format!("{}", mesh);
        assert_eq!(
            text,
            "# header\n\
             # ndim npoint ncell\n\
             2 6 2\n\
             \n\
             # points\n\
             # id x y {z}\n\
             0 0 0\n\
             1 1 0\n\
             2 1 1\n\
             3 0 1\n\
             4 2 0\n\
             5 2 1\n\
             \n\
             # cells\n\
             # id att kind  points\n\
             0 1 qua4  0 1 2 3\n\
             1 2 qua4  1 4 5 2\n"
        );
        let mesh_in = Mesh::from_text(&text).unwrap();
        assert_eq!(format!("{}", mesh_in), text);
    }

    #[test]
    fn display_works_3d() {
        let mesh = Samples::two_hex8();
        let text = format!("{}", mesh);
        assert_eq!(
            text,
            "# header\n\
             # ndim npoint ncell\n\
             3 12 2\n\
             \n\
             # points\n\
             # id x y {z}\n\
             0 0 0 0\n\
             1 1 0 0\n\
             2 1 1 0\n\
             3 0 1 0\n\
             4 0 0 1\n\
             5 1 0 1\n\
             6 1 1 1\n\
             7 0 1 1\n\
             8 0 0 2\n\
             9 1 0 2\n\
             10 1 1 2\n\
             11 0 1 2\n\
             \n\
             # cells\n\
             # id att kind  points\n\
             0 1 hex8  0 1 2 3 4 5 6 7\n\
             1 2 hex8  4 5 6 7 8 9 10 11\n"
        );
        let mesh_in = Mesh::from_text(&text).unwrap();
        assert_eq!(format!("{}", mesh_in), text);
    }

    #[test]
    fn display_works_mixed_2d() {
        let mesh = Samples::mixed_shapes_2d();
        let text = format!("{}", mesh);
        let mesh_in = Mesh::from_text(&text).unwrap();
        assert_eq!(format!("{}", mesh_in), text);
        assert_eq!(mesh_in.cells[0].points.len(), 2);
        assert_eq!(mesh_in.cells[1].points.len(), 4);
    }

    #[test]
    fn display_works_mixed_3d() {
        let mesh = Samples::mixed_shapes_3d();
        let text = format!("{}", mesh);
        let mesh_in = Mesh::from_text(&text).unwrap();
        assert_eq!(format!("{}", mesh_in), text);
        assert_eq!(mesh_in.cells[0].points.len(), 8);
        assert_eq!(mesh_in.cells[4].points.len(), 3);
    }
}
