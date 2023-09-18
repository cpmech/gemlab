use crate::shapes::{GeoKind, Scratchpad};
use crate::StrError;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fmt::{self, Write as FmtWrite};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::Path;

/// Aliases usize as Point ID
pub type PointId = usize;

/// Aliases i32 as Point Marker
pub type PointMarker = i32;

/// Aliases usize as Cell ID
pub type CellId = usize;

/// Aliases usize as Cell's attribute
pub type CellAttribute = usize;

/// Holds point data
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Point {
    /// Identification number which equals the index of the point in the mesh
    pub id: PointId,

    /// Holds a marker that can be used to group points (e.g., on the boundary)
    pub marker: PointMarker,

    /// Point coordinates (2D or 3D)
    pub coords: Vec<f64>,
}

/// Holds cell (aka geometric shape, polygon, polyhedra) data
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Cell {
    /// Identification number which equals the index of the cell in the mesh
    pub id: CellId,

    /// Attribute number
    pub attribute: CellAttribute,

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
/// //      y   (#) indicates attribute
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
///         Point { id: 0, marker: 0, coords: vec![0.0, 0.0] },
///         Point { id: 1, marker: 0, coords: vec![1.0, 0.0] },
///         Point { id: 2, marker: 0, coords: vec![1.0, 1.0] },
///         Point { id: 3, marker: 0, coords: vec![0.0, 1.0] },
///         Point { id: 4, marker: 0, coords: vec![2.0, 0.0] },
///         Point { id: 5, marker: 0, coords: vec![2.0, 1.0] },
///     ],
///     cells: vec![
///         Cell { id: 0, attribute: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
///         Cell { id: 1, attribute: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
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

    /// Searches marked points
    ///
    /// # Input
    ///
    /// * `mark` -- the point marker
    /// * `filter` -- function `fn(x) -> bool` that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true).
    ///   Use `|_| true` or [crate::util::any_x] to allow any point in the resulting array.
    ///   Another example of filter: `|x| x[0] > 1.4 && x[0] < 1.6`
    ///
    /// # Output
    ///
    /// * If at least one point has been found, returns a **sorted** array of point ids.
    /// * Otherwise, returns an error
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::prelude::*;
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     // ```text
    ///     // -500        -400
    ///     // 8------7------6._
    ///     // |       [3](3)|  '-.5
    ///     // |  [0]        |     '-._
    ///     // 9  (1)       10  [1]    '4 -300
    ///     // |             |  (2)  .-'
    ///     // |       [2](3)|   _.3'
    ///     // 0------1------2.-'
    ///     // -100         -200
    ///     // ```
    ///     let h = 0.866; // ~ SQRT_3 / 2
    ///     let m = h / 2.0;
    ///     #[rustfmt::skip]
    ///     let mesh = Mesh {
    ///         ndim: 2,
    ///         points: vec![
    ///             Point { id:  0, marker: -100, coords: vec![0.0,   0.0 ] },
    ///             Point { id:  1, marker:    0, coords: vec![0.5,   0.0 ] },
    ///             Point { id:  2, marker: -200, coords: vec![1.0,   0.0 ] },
    ///             Point { id:  3, marker:    0, coords: vec![1.0+m, 0.25] },
    ///             Point { id:  4, marker: -300, coords: vec![1.0+h, 0.5 ] },
    ///             Point { id:  5, marker:    0, coords: vec![1.0+m, 0.75] },
    ///             Point { id:  6, marker: -400, coords: vec![1.0,   1.0 ] },
    ///             Point { id:  7, marker:    0, coords: vec![0.5,   1.0 ] },
    ///             Point { id:  8, marker: -500, coords: vec![0.0,   1.0 ] },
    ///             Point { id:  9, marker:    0, coords: vec![0.0,   0.5 ] },
    ///             Point { id: 10, marker:    0, coords: vec![1.0,   0.5 ] },
    ///         ],
    ///         cells: vec![
    ///             Cell { id: 0, attribute: 1, kind: GeoKind::Qua8, points: vec![0, 2, 6, 8, 1, 10, 7, 9] },
    ///             Cell { id: 1, attribute: 2, kind: GeoKind::Tri6, points: vec![2, 4, 6, 3, 5, 10] },
    ///             Cell { id: 2, attribute: 3, kind: GeoKind::Lin2, points: vec![2, 10] },
    ///             Cell { id: 3, attribute: 3, kind: GeoKind::Lin2, points: vec![10, 6] },
    ///         ],
    ///     };
    ///     assert_eq!(mesh.search_marked_points(-200, |_| true)?, &[2]);
    ///     assert_eq!(mesh.search_marked_points(-400, |_| true)?, &[6]);
    ///     assert_eq!(mesh.search_marked_points(0, |x| x[1] > 0.49 && x[1] < 0.51)?, &[9, 10]);
    ///     Ok(())
    /// }
    /// ```
    pub fn search_marked_points<F>(&self, marker: PointMarker, mut filter: F) -> Result<Vec<PointId>, StrError>
    where
        F: FnMut(&Vec<f64>) -> bool,
    {
        let mut point_ids: Vec<_> = self
            .points
            .iter()
            .filter_map(|p| {
                if p.marker == marker && filter(&p.coords) {
                    Some(p.id)
                } else {
                    None
                }
            })
            .collect();
        if point_ids.len() == 0 {
            return Err("cannot find at least one point with the given mark (and filter)");
        }
        point_ids.sort();
        Ok(point_ids)
    }

    /// Searches the first marked point for a given mark
    ///
    /// # Input
    ///
    /// * `mark` -- the point marker
    /// * `filter` -- function `fn(x) -> bool` that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true).
    ///   Use `|_| true` or [crate::util::any_x] to allow any point in the resulting array.
    ///   Another example of filter: `|x| x[0] > 1.4 && x[0] < 1.6`
    ///
    /// # Output
    ///
    /// * If **one** point has been found, returns its ID
    /// * Otherwise, returns an error
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::prelude::*;
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     // ```text
    ///     // +3     +3     +3
    ///     //  8------7------6._
    ///     //  |       [3](3)|  '-.5
    ///     //  |  [0]        |     '-._
    ///     //  9  (1)       10  [1]    '4 +2
    ///     //  |             |  (2)  .-'
    ///     //  |       [2](3)|   _.3'
    ///     //  0------1------2.-'
    ///     // +1     +1     +1
    ///     // ```
    ///     let h = 0.866; // ~ SQRT_3 / 2
    ///     let m = h / 2.0;
    ///     #[rustfmt::skip]
    ///     let mesh = Mesh {
    ///         ndim: 2,
    ///         points: vec![
    ///             Point { id:  0, marker: 1, coords: vec![0.0,   0.0 ] },
    ///             Point { id:  1, marker: 1, coords: vec![0.5,   0.0 ] },
    ///             Point { id:  2, marker: 1, coords: vec![1.0,   0.0 ] },
    ///             Point { id:  3, marker: 0, coords: vec![1.0+m, 0.25] },
    ///             Point { id:  4, marker: 2, coords: vec![1.0+h, 0.5 ] },
    ///             Point { id:  5, marker: 0, coords: vec![1.0+m, 0.75] },
    ///             Point { id:  6, marker: 3, coords: vec![1.0,   1.0 ] },
    ///             Point { id:  7, marker: 3, coords: vec![0.5,   1.0 ] },
    ///             Point { id:  8, marker: 3, coords: vec![0.0,   1.0 ] },
    ///             Point { id:  9, marker: 0, coords: vec![0.0,   0.5 ] },
    ///             Point { id: 10, marker: 0, coords: vec![1.0,   0.5 ] },
    ///         ],
    ///         cells: vec![
    ///             Cell { id: 0, attribute: 1, kind: GeoKind::Qua8, points: vec![0, 2, 6, 8, 1, 10, 7, 9] },
    ///             Cell { id: 1, attribute: 2, kind: GeoKind::Tri6, points: vec![2, 4, 6, 3, 5, 10] },
    ///             Cell { id: 2, attribute: 3, kind: GeoKind::Lin2, points: vec![2, 10] },
    ///             Cell { id: 3, attribute: 3, kind: GeoKind::Lin2, points: vec![10, 6] },
    ///         ],
    ///     };
    ///     assert_eq!(mesh.search_first_marked_point(1, |_| true)?, 0);
    ///     assert_eq!(mesh.search_first_marked_point(2, |_| true)?, 4);
    ///     assert_eq!(mesh.search_first_marked_point(3, |_| true)?, 6);
    ///     assert_eq!(mesh.search_first_marked_point(3, |x| x[0] > 0.49 && x[0] < 0.51)?, 7);
    ///     Ok(())
    /// }
    /// ```
    pub fn search_first_marked_point<F>(&self, marker: PointMarker, mut filter: F) -> Result<PointId, StrError>
    where
        F: FnMut(&[f64]) -> bool,
    {
        for i in 0..self.points.len() {
            if self.points[i].marker == marker && filter(&self.points[i].coords) {
                return Ok(self.points[i].id);
            }
        }
        return Err("cannot find at least one point with the given mark (and filter)");
    }

    /// Sets the pad's matrix of coordinates X given a list of point ids
    ///
    /// # Panics
    ///
    /// 1. Make sure `pad.kind.nnode() == points.len()`; otherwise a panic will occur
    /// 2. This function does not check for bounds on point indices and dimensions
    /// 3. Use [Mesh::check_all] to capture (some) errors
    pub fn set_pad(&self, pad: &mut Scratchpad, points: &[PointId]) {
        let nnode = pad.kind.nnode();
        assert_eq!(nnode, points.len());
        for m in 0..nnode {
            for j in 0..self.ndim {
                pad.set_xx(m, j, self.points[points[m]].coords[j]);
            }
        }
    }

    /// Returns the (min,max) point coordinates in a mesh
    pub fn get_limits(&self) -> (Vec<f64>, Vec<f64>) {
        let mut min = vec![f64::MAX; self.ndim];
        let mut max = vec![f64::MIN; self.ndim];
        for point in &self.points {
            for i in 0..self.ndim {
                min[i] = f64::min(min[i], point.coords[i]);
                max[i] = f64::max(max[i], point.coords[i]);
            }
        }
        (min, max)
    }
}

impl fmt::Display for Mesh {
    /// Returns a text representation of the Mesh (can be used with [Mesh::from_text])
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "# header\n").unwrap();
        write!(f, "# ndim npoint ncell\n").unwrap();
        write!(f, "{} {} {}\n", self.ndim, self.points.len(), self.cells.len()).unwrap();
        write!(f, "\n# points\n").unwrap();
        write!(f, "# id marker x y {{z}}\n").unwrap();
        self.points.iter().for_each(|point| {
            if self.ndim == 2 {
                write!(
                    f,
                    "{} {} {:?} {:?}\n",
                    point.id, point.marker, point.coords[0], point.coords[1]
                )
                .unwrap();
            } else {
                write!(
                    f,
                    "{} {} {:?} {:?} {:?}\n",
                    point.id, point.marker, point.coords[0], point.coords[1], point.coords[2]
                )
                .unwrap();
            }
        });
        write!(f, "\n# cells\n").unwrap();
        write!(f, "# id attribute kind points\n").unwrap();
        self.cells.iter().for_each(|cell| {
            write!(
                f,
                "{} {} {}{}\n",
                cell.id,
                cell.attribute,
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
    use crate::shapes::Scratchpad;
    use serde_json;

    #[test]
    fn read_and_write_capture_errors() {
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
    }

    #[test]
    fn read_and_write_work() {
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
            # id marker x y
               0 0 0.0 0.0
               1 0 1.0 0.0
               2 0 1.0 1.0
               3 0 0.0 1.0
               4 0 2.0 0.0
               5 0 2.0 1.0
            
            # cells
            # id attribute kind point_ids...
               0   1 qua4 0 1 2 3
               1   0 qua4 1 4 5 2",
        )
        .unwrap();
        mesh.write("/tmp/gemlab/test.msh").unwrap();
        let mesh_read = Mesh::read("/tmp/gemlab/test.msh").unwrap();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", mesh_read));
    }

    #[test]
    fn derive_works() {
        let mesh = Samples::two_qua4();
        let mesh_clone = mesh.clone();
        let correct ="Mesh { ndim: 2, points: [Point { id: 0, marker: 0, coords: [0.0, 0.0] }, Point { id: 1, marker: 0, coords: [1.0, 0.0] }, Point { id: 2, marker: 0, coords: [1.0, 1.0] }, Point { id: 3, marker: 0, coords: [0.0, 1.0] }, Point { id: 4, marker: 0, coords: [2.0, 0.0] }, Point { id: 5, marker: 0, coords: [2.0, 1.0] }], cells: [Cell { id: 0, attribute: 1, kind: Qua4, points: [0, 1, 2, 3] }, Cell { id: 1, attribute: 2, kind: Qua4, points: [1, 4, 5, 2] }] }";
        assert_eq!(format!("{:?}", mesh), correct);
        assert_eq!(mesh_clone.ndim, mesh.ndim);
        assert_eq!(mesh_clone.points.len(), mesh.points.len());
        assert_eq!(mesh_clone.cells.len(), mesh.cells.len());
        // serialize
        let mesh_json = serde_json::to_string(&mesh).unwrap();
        // deserialize
        let mesh_read: Mesh = serde_json::from_str(&mesh_json).unwrap();
        assert_eq!(format!("{:?}", mesh_read), correct);
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
             # id marker x y {z}\n\
             0 0 0.0 0.0\n\
             1 0 1.0 0.0\n\
             2 0 1.0 1.0\n\
             3 0 0.0 1.0\n\
             4 0 2.0 0.0\n\
             5 0 2.0 1.0\n\
             \n\
             # cells\n\
             # id attribute kind points\n\
             0 1 qua4 0 1 2 3\n\
             1 2 qua4 1 4 5 2\n"
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
             # id marker x y {z}\n\
             0 0 0.0 0.0 0.0\n\
             1 0 1.0 0.0 0.0\n\
             2 0 1.0 1.0 0.0\n\
             3 0 0.0 1.0 0.0\n\
             4 0 0.0 0.0 1.0\n\
             5 0 1.0 0.0 1.0\n\
             6 0 1.0 1.0 1.0\n\
             7 0 0.0 1.0 1.0\n\
             8 0 0.0 0.0 2.0\n\
             9 0 1.0 0.0 2.0\n\
             10 0 1.0 1.0 2.0\n\
             11 0 0.0 1.0 2.0\n\
             \n\
             # cells\n\
             # id attribute kind points\n\
             0 1 hex8 0 1 2 3 4 5 6 7\n\
             1 2 hex8 4 5 6 7 8 9 10 11\n"
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

    #[test]
    fn search_marked_points_works() {
        let mesh = Samples::four_tri3();
        assert_eq!(mesh.search_marked_points(-1, |_| true).unwrap(), &[0]);
        assert_eq!(mesh.search_marked_points(-3, |_| true).unwrap(), &[2]);
        assert_eq!(mesh.search_marked_points(-4, |_| true).unwrap(), &[3]);
        assert_eq!(mesh.search_marked_points(-5, |_| true).unwrap(), &[4]);
        assert_eq!(
            mesh.search_marked_points(-10, |_| true).err(),
            Some("cannot find at least one point with the given mark (and filter)")
        );

        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        assert_eq!(mesh.search_marked_points(-1, |_| true).unwrap(), &[0]);
        assert_eq!(mesh.search_marked_points(-2, |_| true).unwrap(), &[2]);
        assert_eq!(mesh.search_marked_points(-3, |_| true).unwrap(), &[14]);
        assert_eq!(mesh.search_marked_points(-4, |_| true).unwrap(), &[12]);
        assert_eq!(
            mesh.search_marked_points(-10, |_| true).unwrap(),
            &[3, 6, 9, 25, 28, 31, 34]
        );
        assert_eq!(
            mesh.search_marked_points(-20, |_| true).unwrap(),
            &[5, 8, 11, 27, 30, 33, 36]
        );
        assert_eq!(mesh.search_marked_points(-30, |_| true).unwrap(), &[1, 15, 16]);
        assert_eq!(mesh.search_marked_points(-40, |_| true).unwrap(), &[13, 23, 24]);
        assert_eq!(
            mesh.search_marked_points(0, |_| true).unwrap(),
            &[4, 7, 10, 17, 18, 19, 20, 21, 22, 26, 29, 32, 35]
        );
        assert_eq!(
            mesh.search_marked_points(8, |_| true).err(),
            Some("cannot find at least one point with the given mark (and filter)")
        );
        assert_eq!(
            mesh.search_marked_points(-30, |x| x[0] > 10.0).err(),
            Some("cannot find at least one point with the given mark (and filter)")
        );
        assert_eq!(
            mesh.search_marked_points(-30, |x| x[0] > 1.4 && x[0] < 1.6).unwrap(),
            &[1]
        );
    }

    #[test]
    fn search_first_marked_point_works() {
        let mesh = Samples::two_tri3_one_qua4();
        assert_eq!(mesh.search_first_marked_point(1, |_| true).unwrap(), 0);
        assert_eq!(mesh.search_first_marked_point(2, |_| true).unwrap(), 2);
        assert_eq!(mesh.search_first_marked_point(2, |x| x[0] > 1.5).unwrap(), 5);
        assert_eq!(
            mesh.search_first_marked_point(-10, |_| true).err(),
            Some("cannot find at least one point with the given mark (and filter)")
        );

        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        assert_eq!(mesh.search_first_marked_point(-1, |_| true).unwrap(), 0);
        assert_eq!(mesh.search_first_marked_point(-2, |_| true).unwrap(), 2);
        assert_eq!(mesh.search_first_marked_point(-3, |_| true).unwrap(), 14);
        assert_eq!(mesh.search_first_marked_point(-4, |_| true).unwrap(), 12);
        assert_eq!(mesh.search_first_marked_point(-10, |_| true).unwrap(), 3);
        assert_eq!(mesh.search_first_marked_point(-20, |_| true).unwrap(), 5);
        assert_eq!(mesh.search_first_marked_point(-30, |_| true).unwrap(), 1);
        assert_eq!(mesh.search_first_marked_point(-40, |_| true).unwrap(), 13);
        assert_eq!(mesh.search_first_marked_point(0, |_| true).unwrap(), 4);
        assert_eq!(
            mesh.search_first_marked_point(8, |_| true).err(),
            Some("cannot find at least one point with the given mark (and filter)")
        );
        assert_eq!(
            mesh.search_first_marked_point(-30, |x| x[0] > 10.0).err(),
            Some("cannot find at least one point with the given mark (and filter)")
        );
        assert_eq!(
            mesh.search_first_marked_point(-30, |x| x[0] > 1.24 && x[0] < 1.26)
                .unwrap(),
            15
        );
    }

    #[test]
    fn set_pad_works_2d() {
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        let mesh = Samples::two_qua4();
        let cell = &mesh.cells[0];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        mesh.set_pad(&mut pad, &cell.points);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌         ┐\n\
             │ 0 1 1 0 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
        let cell = &mesh.cells[1];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        mesh.set_pad(&mut pad, &cell.points);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌         ┐\n\
             │ 1 2 2 1 │\n\
             │ 0 0 1 1 │\n\
             └         ┘"
        );
    }

    #[test]
    fn set_pad_works_3d() {
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
        let mesh = Samples::two_hex8();
        let cell = &mesh.cells[0];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        mesh.set_pad(&mut pad, &cell.points);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 0 0 0 0 1 1 1 1 │\n\
             └                 ┘"
        );
        let cell = &mesh.cells[1];
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        mesh.set_pad(&mut pad, &cell.points);
        assert_eq!(
            format!("{}", pad.xxt),
            "┌                 ┐\n\
             │ 0 1 1 0 0 1 1 0 │\n\
             │ 0 0 1 1 0 0 1 1 │\n\
             │ 1 1 1 1 2 2 2 2 │\n\
             └                 ┘"
        );
    }

    #[test]
    fn get_mesh_limits_works() {
        let mesh = &Samples::two_qua4();
        let (min, max) = mesh.get_limits();
        assert_eq!(min, &[0.0, 0.0]);
        assert_eq!(max, &[2.0, 1.0]);

        let mesh = &Samples::two_hex8();
        let (min, max) = mesh.get_limits();
        assert_eq!(min, &[0.0, 0.0, 0.0]);
        assert_eq!(max, &[1.0, 1.0, 2.0]);
    }
}
