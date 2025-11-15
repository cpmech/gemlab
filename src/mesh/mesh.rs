use crate::shapes::{GeoKind, Scratchpad};
use crate::StrError;
use russell_lab::{argsort2_f64, argsort3_f64};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fmt::{self, Write as FmtWrite};
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Defines a tolerance to compare points in [Mesh::get_sorted_points()]
pub const TOL_COMPARE_POINTS: f64 = 1e-6;

/// Aliases usize as Point ID
pub type PointId = usize;

/// Aliases i32 as Point Marker
pub type PointMarker = i32;

/// Aliases usize as Cell ID
pub type CellId = usize;

/// Aliases usize as Cell's attribute
pub type CellAttribute = i32;

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
///     marked_edges: Vec::new(),
///     marked_faces: Vec::new(),
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

impl Mesh {
    /// Allocates a new Mesh with zeroed coordinates and connectivities and homogeneous cells (all cells are of the same kind)
    ///
    /// Use this function to allocate memory and later fill the coordinates and cell connectivities as appropriate.
    ///
    /// # Arguments
    ///
    /// * `space_ndim` -- space dimension (2 or 3). Also, `space_ndim` must be greater than or equal to
    ///   `geo_ndim`; e.g., you cannot have a 3D shape in a 2D space.
    /// * `npoint` -- number of points (≥ 2)
    /// * `ncell` -- number of cells (≥ 1)
    /// * `kind` -- cell kind, such that `geo_kind` ≤ `space_ndim` with `geo_ndim = kind.ndim()`
    ///
    /// **Warning:** This function does not check for validity of the mesh. Use [Mesh::check_all()] to capture (some) errors.
    ///
    /// Geometry cases regarding the number of dimensions (geo vs space)
    ///
    /// 1. Case `CABLE` -- `geo_ndim = 1` and `space_ndim = 2 or 3`; e.g., line in 2D or 3D (cables and rods)
    /// 2. Case `SHELL` -- `geo_ndim = 2` and `space_ndim = 3`; e.g. Tri or Qua in 3D (shells and surfaces)
    /// 3. Case `SOLID` -- `geo_ndim = space_ndim`; e.g., Tri and Qua in 2D or Tet and Hex in 3D
    ///
    /// | `geo_ndim` | `space_ndim = 2` | `space_ndim = 3` |
    /// |:----------:|:----------------:|:----------------:|
    /// |     1      |     `CABLE`      |     `CABLE`      |
    /// |     2      |     `SOLID`      |     `SHELL`      |
    /// |     3      |    impossible    |     `SOLID`      |
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::prelude::*;
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     let mesh = Mesh::new_zero_homogeneous(2, 4, 2, GeoKind::Tri3).unwrap();
    ///     assert_eq!(mesh.ndim, 2);
    ///     assert_eq!(mesh.points.len(), 4);
    ///     assert_eq!(mesh.cells.len(), 2);
    ///     for i in 0..4 {
    ///         let point = &mesh.points[i];
    ///         assert_eq!(point.id, i);
    ///         assert_eq!(point.coords, vec![0.0, 0.0]);
    ///     }
    ///     for e in 0..2 {
    ///         let cell = &mesh.cells[e];
    ///         assert_eq!(cell.id, e);
    ///         assert_eq!(cell.attribute, 1);
    ///         assert_eq!(cell.kind, GeoKind::Tri3);
    ///         assert_eq!(cell.points, vec![0, 0, 0]);
    ///     }
    ///     assert_eq!(
    ///         mesh.check_all().err(),
    ///         Some("cannot compute inverse due to zero determinant")
    ///     );
    ///     Ok(())
    /// }
    /// ```
    pub fn new_zero_homogeneous(
        space_ndim: usize,
        npoint: usize,
        ncell: usize,
        kind: GeoKind,
    ) -> Result<Self, StrError> {
        if space_ndim < 2 || space_ndim > 3 {
            return Err("space_ndim must be 2 or 3");
        }
        if npoint < 2 {
            return Err("npoint must be ≥ 2");
        }
        if ncell < 1 {
            return Err("ncell must be ≥ 1");
        }
        let geo_ndim = kind.ndim();
        if geo_ndim > space_ndim {
            return Err("geo_ndim cannot be greater than space_ndim");
        }
        let nnode = kind.nnode();
        let mut points: Vec<Point> = Vec::with_capacity(npoint);
        let mut cells: Vec<Cell> = Vec::with_capacity(ncell);
        for i in 0..npoint {
            points.push(Point {
                id: i,
                marker: 0,
                coords: vec![0.0; space_ndim],
            });
        }
        for i in 0..ncell {
            cells.push(Cell {
                id: i,
                attribute: 1,
                kind,
                points: vec![0; nnode],
            });
        }
        Ok(Mesh {
            ndim: space_ndim,
            points,
            cells,
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        })
    }

    /// Reads a JSON file containing the mesh data
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

    /// Writes a JSON file with the mesh data
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

    /// Searches points with a given marker
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
    ///         marked_edges: Vec::new(),
    ///         marked_faces: Vec::new(),
    ///     };
    ///     assert_eq!(mesh.search_marked_points(-200, |_| true)?, &[2]);
    ///     assert_eq!(mesh.search_marked_points(-400, |_| true)?, &[6]);
    ///     assert_eq!(mesh.search_marked_points(0, |x| x[1] > 0.49 && x[1] < 0.51)?, &[9, 10]);
    ///     Ok(())
    /// }
    /// ```
    pub fn search_marked_points<F>(&self, marker: PointMarker, mut filter: F) -> Result<Vec<PointId>, StrError>
    where
        F: FnMut(&[f64]) -> bool,
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
    ///         marked_edges: Vec::new(),
    ///         marked_faces: Vec::new(),
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

    /// Allocates a Scratchpad for numerical integration and sets its coordinates
    ///
    /// # Panics
    ///
    /// 1. This function does not check for bounds on point indices and dimensions
    /// 2. Use [Mesh::check_all] to capture (some) errors
    pub fn get_pad(&self, cell_id: CellId) -> Scratchpad {
        let cell = &self.cells[cell_id];
        let mut pad = Scratchpad::new(self.ndim, cell.kind).unwrap();
        self.set_pad(&mut pad, &cell.points);
        pad
    }

    /// Returns the min and max coordinates of all points in the mesh
    ///
    /// Returns `(min, max)`, each being an `ndim` vector.
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

    /// Returns the bounding box of a cell
    ///
    /// Returns `(min, max)`, each being an `ndim` vector.
    pub fn get_cell_bounding_box(&self, cell_id: CellId) -> (Vec<f64>, Vec<f64>) {
        let cell = &self.cells[cell_id];
        let mut min = vec![f64::MAX; self.ndim];
        let mut max = vec![f64::MIN; self.ndim];
        for m in 0..cell.points.len() {
            let p = cell.points[m];
            for i in 0..self.ndim {
                min[i] = f64::min(min[i], self.points[p].coords[i]);
                max[i] = f64::max(max[i], self.points[p].coords[i]);
            }
        }
        (min, max)
    }

    /// Returns the IDs of points such that their coordinates are sorted by z → y → x (ascending)
    ///
    /// This function sorts the points and returns the IDs of the points. In 2D, the sorting
    /// is done by y → x (ascending). In 3D, the sorting is done by z → y → x (ascending).
    ///
    /// # Input
    ///
    /// * `point_ids` -- list of point ids to be sorted
    /// * `filter` -- function `fn(x, y, z) -> bool` that returns true to **keep** the coordinate just found
    ///
    /// # Notes
    ///
    /// 1. The filter is applied before sorting the points.
    /// 2. The tolerance to compare points is [TOL_COMPARE_POINTS] times the range of the coordinates; i.e.:
    ///
    /// ```text
    /// tol_x = TOL_COMPARE_POINTS * (xmax - xmin)
    /// tol_y = TOL_COMPARE_POINTS * (ymax - ymin)
    /// tol_z = TOL_COMPARE_POINTS * (zmax - zmin)
    /// ```
    ///
    /// # Output
    ///
    /// Returns the IDs of the points sorted by z → y → x (ascending)
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::mesh::Samples;
    ///
    /// fn main() {
    ///     // 8------7------6._
    ///     // |             |  '-.5
    ///     // |             |     '-._
    ///     // 9            10         '4
    ///     // |             |       .-'
    ///     // |             |   _.3'
    ///     // 0------1------2.-'
    ///     let mesh = Samples::qua8_tri6_lin2();
    ///
    ///     // define the point IDs to be sorted
    ///     let point_ids = vec![8, 7, 6, 5, 4, 3, 2, 1, 0, 10, 9];
    ///
    ///     // use the get_sorted_points method with a filter
    ///     let sorted_points = mesh.get_sorted_points(&point_ids, |x, y, _| {
    ///         x > 0.0 && y > 0.0
    ///     });
    ///
    ///     // check
    ///     assert_eq!(sorted_points, &[3, 10, 4, 5, 7, 6]);
    /// }
    /// ```
    pub fn get_sorted_points<F>(&self, point_ids: &[PointId], filter: F) -> Vec<PointId>
    where
        F: Fn(f64, f64, f64) -> bool,
    {
        // build arrays of coordinates
        let d3 = self.ndim == 3;
        let np = point_ids.len();
        let mut k2id = HashMap::new();
        let mut xx = Vec::with_capacity(np);
        let mut yy = Vec::with_capacity(np);
        let mut zz = if d3 { Vec::with_capacity(np) } else { Vec::new() };
        for nid in point_ids {
            let x = self.points[*nid].coords[0];
            let y = self.points[*nid].coords[1];
            let z = if d3 { self.points[*nid].coords[2] } else { 0.0 };
            if filter(x, y, z) {
                k2id.insert(xx.len(), *nid);
                xx.push(x);
                yy.push(y);
                if d3 {
                    zz.push(z);
                }
            }
        }
        // calculate tolerances to compare point coordinates
        let (min, max) = self.get_limits();
        let mut tol = vec![TOL_COMPARE_POINTS; self.ndim];
        for i in 0..self.ndim {
            tol[i] *= max[i] - min[i];
        }
        // sort nodes by x → y → z
        let sorted_indices = if d3 {
            argsort3_f64(&zz, &yy, &xx, &tol)
        } else {
            argsort2_f64(&yy, &xx, &tol)
        };
        sorted_indices.iter().map(|&i| k2id[&i]).collect()
    }

    /// Renumber all points (reordering)
    ///
    /// # Input
    ///
    /// * `old_to_new` -- maps old point ids to new point ids. Must have `len = npoint`
    ///
    /// # Panics
    ///
    /// This function will panic if old_to_new.len() != points.len()
    pub fn renumber_points(&mut self, old_to_new: &[PointId]) -> Result<(), StrError> {
        let npoint = self.points.len();
        if old_to_new.len() != npoint {
            return Err("old_to_new.len() must be equal to the number of points");
        }
        let mut updated = vec![false; npoint];
        let points = self.points.clone();
        for old in 0..npoint {
            let new = old_to_new[old];
            if new >= npoint {
                return Err("new point id is out of range");
            }
            self.points[new].id = new;
            self.points[new].marker = points[old].marker;
            for k in 0..self.ndim {
                self.points[new].coords[k] = points[old].coords[k];
            }
            updated[new] = true;
        }
        for ok in &updated {
            if !ok {
                return Err("not all points have been updated");
            }
        }
        for cell in &mut self.cells {
            for i in 0..cell.points.len() {
                let old = cell.points[i];
                cell.points[i] = old_to_new[old];
            }
        }
        Ok(())
    }
}

impl fmt::Display for Mesh {
    /// Returns a text representation of the Mesh (can be used with [Mesh::from_text])
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // write header
        write!(f, "# header\n").unwrap();
        write!(f, "# ndim npoint ncell nmarked_edge nmarked_face\n").unwrap();
        write!(
            f,
            "{} {} {} {} {}\n",
            self.ndim,
            self.points.len(),
            self.cells.len(),
            self.marked_edges.len(),
            self.marked_faces.len()
        )
        .unwrap();

        // write points
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

        // write cells
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

        // write marked edges
        if self.marked_edges.len() > 0 {
            write!(f, "\n# marked edges\n").unwrap();
            write!(f, "# marker p1 p2\n").unwrap();
            for entry in &self.marked_edges {
                write!(f, "{} {} {}\n", entry.0, entry.1, entry.2).unwrap();
            }
        }

        // write marked faces
        if self.marked_faces.len() > 0 {
            write!(f, "\n# marked faces\n").unwrap();
            write!(f, "# marker p1 p2 p3 {{p4}}\n").unwrap();
            for entry in &self.marked_faces {
                if entry.4 == usize::MAX {
                    write!(f, "{} {} {} {}\n", entry.0, entry.1, entry.2, entry.3).unwrap();
                } else {
                    write!(f, "{} {} {} {} {}\n", entry.0, entry.1, entry.2, entry.3, entry.4).unwrap();
                }
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::mesh::{Mesh, Point, Samples};
    use crate::shapes::{GeoKind, Scratchpad};

    #[test]
    fn new_zero_homogeneous_captures_errors() {
        assert_eq!(
            Mesh::new_zero_homogeneous(1, 2, 1, GeoKind::Lin2).err(),
            Some("space_ndim must be 2 or 3")
        );
        assert_eq!(
            Mesh::new_zero_homogeneous(2, 1, 1, GeoKind::Lin2).err(),
            Some("npoint must be ≥ 2")
        );
        assert_eq!(
            Mesh::new_zero_homogeneous(2, 2, 0, GeoKind::Lin2).err(),
            Some("ncell must be ≥ 1")
        );
        assert_eq!(
            Mesh::new_zero_homogeneous(2, 2, 1, GeoKind::Hex8).err(),
            Some("geo_ndim cannot be greater than space_ndim")
        );
    }

    #[test]
    fn new_zero_homogeneous_works() {
        let mesh = Mesh::new_zero_homogeneous(2, 4, 2, GeoKind::Tri3).unwrap();
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 4);
        assert_eq!(mesh.cells.len(), 2);
        for i in 0..4 {
            let point = &mesh.points[i];
            assert_eq!(point.id, i);
            assert_eq!(point.coords, vec![0.0, 0.0]);
        }
        for i in 0..2 {
            let cell = &mesh.cells[i];
            assert_eq!(cell.id, i);
            assert_eq!(cell.attribute, 1);
            assert_eq!(cell.kind, GeoKind::Tri3);
            assert_eq!(cell.points, vec![0, 0, 0]);
        }
        assert_eq!(
            mesh.check_all().err(),
            Some("cannot compute inverse due to zero determinant")
        );
    }

    #[test]
    fn read_and_write_capture_errors() {
        assert_eq!(Mesh::read_json("/tmp/not_found").err(), Some("cannot open file"));
        assert_eq!(
            Mesh::read_json("./data/meshes/two_quads_horizontal.msh").err(),
            Some("cannot parse JSON file")
        );
        let mesh = Mesh {
            ndim: 2,
            points: Vec::new(),
            cells: Vec::new(),
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };
        assert_eq!(mesh.write_json("/tmp/").err(), Some("cannot create file"));
    }

    #[test]
    fn read_and_write_work() {
        //     -100
        //  3--------2--------5
        //  |(-4)    |(-3)    |(-6)
        //  |        |        |
        //  |(-1)    |(-2)    |(-5)
        //  0--------1--------4
        //
        let mesh = Mesh::from_text(
            r"# header
            # ndim npoint ncell nmarked_edge nmarked_face
                 2      6     2            1            0
            
            # points
            # id marker x y
               0 -1 0.0 0.0
               1 -2 1.0 0.0
               2 -3 1.0 1.0
               3 -4 0.0 1.0
               4 -5 2.0 0.0
               5 -6 2.0 1.0
            
            # cells
            # id attribute kind point_ids...
               0   1 qua4 0 1 2 3
               1   0 qua4 1 4 5 2
               
            # marked edges
            # marker p1 p2
                -100  3  2",
        )
        .unwrap();
        mesh.write_json("/tmp/gemlab/test.json").unwrap();
        let mesh_read = Mesh::read_json("/tmp/gemlab/test.json").unwrap();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", mesh_read));
    }

    #[test]
    fn derive_works() {
        let mesh = Samples::two_qua4();
        let mesh_clone = mesh.clone();
        let correct ="Mesh { ndim: 2, points: [Point { id: 0, marker: -1, coords: [0.0, 0.0] }, Point { id: 1, marker: -2, coords: [1.0, 0.0] }, Point { id: 2, marker: -3, coords: [1.0, 1.0] }, Point { id: 3, marker: -4, coords: [0.0, 1.0] }, Point { id: 4, marker: -5, coords: [2.0, 0.0] }, Point { id: 5, marker: -6, coords: [2.0, 1.0] }], cells: [Cell { id: 0, attribute: 1, kind: Qua4, points: [0, 1, 2, 3] }, Cell { id: 1, attribute: 2, kind: Qua4, points: [1, 4, 5, 2] }], marked_edges: [(-100, 3, 2), (-200, 2, 5)], marked_faces: [] }";
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
             # ndim npoint ncell nmarked_edge nmarked_face\n\
             2 6 2 2 0\n\
             \n\
             # points\n\
             # id marker x y {z}\n\
             0 -1 0.0 0.0\n\
             1 -2 1.0 0.0\n\
             2 -3 1.0 1.0\n\
             3 -4 0.0 1.0\n\
             4 -5 2.0 0.0\n\
             5 -6 2.0 1.0\n\
             \n\
             # cells\n\
             # id attribute kind points\n\
             0 1 qua4 0 1 2 3\n\
             1 2 qua4 1 4 5 2\n\
             \n\
             # marked edges\n\
             # marker p1 p2\n\
             -100 3 2\n\
             -200 2 5\n"
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
             # ndim npoint ncell nmarked_edge nmarked_face\n\
             3 12 2 4 2\n\
             \n\
             # points\n\
             # id marker x y {z}\n\
             0 0 0.0 0.0 0.0\n\
             1 -1 1.0 0.0 0.0\n\
             2 -1 1.0 1.0 0.0\n\
             3 0 0.0 1.0 0.0\n\
             4 0 0.0 0.0 1.0\n\
             5 0 1.0 0.0 1.0\n\
             6 0 1.0 1.0 1.0\n\
             7 -1 0.0 1.0 1.0\n\
             8 0 0.0 0.0 2.0\n\
             9 0 1.0 0.0 2.0\n\
             10 0 1.0 1.0 2.0\n\
             11 0 0.0 1.0 2.0\n\
             \n\
             # cells\n\
             # id attribute kind points\n\
             0 1 hex8 0 1 2 3 4 5 6 7\n\
             1 2 hex8 4 5 6 7 8 9 10 11\n\
             \n\
             # marked edges\n\
             # marker p1 p2\n\
             123 7 11\n\
             -5 11 10\n\
             -4 7 3\n\
             -5 8 9\n\
             \n\
             # marked faces\n\
             # marker p1 p2 p3 {p4}\n\
             -8 3 2 7 6\n\
             -9 8 10 9 11\n"
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
    fn set_and_get_pad_works_2d() {
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        let mesh = Samples::two_qua4();
        let cell = &mesh.cells[0];
        let mut pad_0a = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        mesh.set_pad(&mut pad_0a, &cell.points);
        let pad_0b = mesh.get_pad(0);
        let correct = "┌         ┐\n\
                       │ 0 1 1 0 │\n\
                       │ 0 0 1 1 │\n\
                       └         ┘";
        assert_eq!(format!("{}", pad_0a.xxt), correct);
        assert_eq!(format!("{}", pad_0b.xxt), correct);
        let cell = &mesh.cells[1];
        let mut pad_1a = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        mesh.set_pad(&mut pad_1a, &cell.points);
        let pad_1b = mesh.get_pad(1);
        let correct = "┌         ┐\n\
                       │ 1 2 2 1 │\n\
                       │ 0 0 1 1 │\n\
                       └         ┘";
        assert_eq!(format!("{}", pad_1a.xxt), correct);
        assert_eq!(format!("{}", pad_1b.xxt), correct);
    }

    #[test]
    fn set_and_get_pad_works_3d() {
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
        let mut pad_0a = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        mesh.set_pad(&mut pad_0a, &cell.points);
        let pad_0b = mesh.get_pad(0);
        let correct = "┌                 ┐\n\
                       │ 0 1 1 0 0 1 1 0 │\n\
                       │ 0 0 1 1 0 0 1 1 │\n\
                       │ 0 0 0 0 1 1 1 1 │\n\
                       └                 ┘";
        assert_eq!(format!("{}", pad_0a.xxt), correct);
        assert_eq!(format!("{}", pad_0b.xxt), correct);
        let cell = &mesh.cells[1];
        let mut pad_1a = Scratchpad::new(mesh.ndim, cell.kind).unwrap();
        mesh.set_pad(&mut pad_1a, &cell.points);
        let pad_1b = mesh.get_pad(1);
        let correct = "┌                 ┐\n\
                       │ 0 1 1 0 0 1 1 0 │\n\
                       │ 0 0 1 1 0 0 1 1 │\n\
                       │ 1 1 1 1 2 2 2 2 │\n\
                       └                 ┘";
        assert_eq!(format!("{}", pad_1a.xxt), correct);
        assert_eq!(format!("{}", pad_1b.xxt), correct);
    }

    #[test]
    fn get_limits_works() {
        let mesh = &Samples::two_qua4();
        let (min, max) = mesh.get_limits();
        assert_eq!(min, &[0.0, 0.0]);
        assert_eq!(max, &[2.0, 1.0]);

        let mesh = &Samples::two_hex8();
        let (min, max) = mesh.get_limits();
        assert_eq!(min, &[0.0, 0.0, 0.0]);
        assert_eq!(max, &[1.0, 1.0, 2.0]);
    }

    #[test]
    fn get_cell_bounding_box_works() {
        let mesh = &Samples::two_qua4();
        let (min, max) = mesh.get_cell_bounding_box(0);
        assert_eq!(min, &[0.0, 0.0]);
        assert_eq!(max, &[1.0, 1.0]);
        let (min, max) = mesh.get_cell_bounding_box(1);
        assert_eq!(min, &[1.0, 0.0]);
        assert_eq!(max, &[2.0, 1.0]);

        let mesh = &Samples::two_hex8();
        let (min, max) = mesh.get_cell_bounding_box(0);
        assert_eq!(min, &[0.0, 0.0, 0.0]);
        assert_eq!(max, &[1.0, 1.0, 1.0]);
        let (min, max) = mesh.get_cell_bounding_box(1);
        assert_eq!(min, &[0.0, 0.0, 1.0]);
        assert_eq!(max, &[1.0, 1.0, 2.0]);
    }

    #[test]
    fn reorder_vertices_works() {
        //  0        1        2 (new)
        //  ↑        ↑        ↑
        //  3--------2--------5
        //  |(-4)    |(-3)    |(-6)
        //  |        |        |
        //  |(-1)    |(-2)    |(-5)
        //  0--------1--------4
        //  ↓        ↓        ↓
        //  3        4        5 (new)
        let mut mesh = Samples::two_qua4();
        //                 0  1  2  3  4  5   // old
        let old_to_new = &[3, 4, 1, 0, 5, 2]; // new
        mesh.renumber_points(old_to_new).unwrap();

        // check
        assert_eq!(mesh.cells[0].points, &[3, 4, 1, 0]);
        assert_eq!(mesh.cells[1].points, &[4, 5, 2, 1]);
        assert_eq!(mesh.points[0].id, 0);
        assert_eq!(mesh.points[0].marker, -4);
        assert_eq!(mesh.points[0].coords, &[0.0, 1.0]);
        assert_eq!(mesh.points[1].id, 1);
        assert_eq!(mesh.points[1].marker, -3);
        assert_eq!(mesh.points[1].coords, &[1.0, 1.0]);
        assert_eq!(mesh.points[2].id, 2);
        assert_eq!(mesh.points[2].marker, -6);
        assert_eq!(mesh.points[2].coords, &[2.0, 1.0]);
        assert_eq!(mesh.points[3].id, 3);
        assert_eq!(mesh.points[3].marker, -1);
        assert_eq!(mesh.points[3].coords, &[0.0, 0.0]);
        assert_eq!(mesh.points[4].id, 4);
        assert_eq!(mesh.points[4].marker, -2);
        assert_eq!(mesh.points[4].coords, &[1.0, 0.0]);
        assert_eq!(mesh.points[5].id, 5);
        assert_eq!(mesh.points[5].marker, -5);
        assert_eq!(mesh.points[5].coords, &[2.0, 0.0]);

        // catch errors
        assert_eq!(
            mesh.renumber_points(&[1, 2]).err(),
            Some("old_to_new.len() must be equal to the number of points")
        );
        assert_eq!(
            mesh.renumber_points(&[0, 0, 6, 0, 0, 0]).err(),
            Some("new point id is out of range")
        );
        assert_eq!(
            mesh.renumber_points(&[0, 0, 0, 0, 0, 0]).err(),
            Some("not all points have been updated")
        );
    }

    #[test]
    fn test_get_sorted_points_2d() {
        #[rustfmt::skip]
        let points = vec![
            Point { id: 0, marker: 0, coords: vec![1.0, 2.0] },
            Point { id: 1, marker: 0, coords: vec![3.0, 4.0] },
            Point { id: 2, marker: 0, coords: vec![1.0, 1.0] },
            Point { id: 3, marker: 0, coords: vec![2.0, 3.0] },
            Point { id: 4, marker: 0, coords: vec![0.0, 1.0] },
        ];
        let mesh = Mesh {
            points,
            ndim: 2,
            cells: vec![],
            marked_edges: Vec::new(),
            marked_faces: Vec::new(),
        };
        let point_ids = vec![0, 1, 2, 3, 4];
        let sorted_points = mesh.get_sorted_points(&point_ids, |_, _, _| true);
        assert_eq!(sorted_points, vec![4, 2, 0, 3, 1]);
    }

    #[test]
    fn test_get_sorted_points_3d() {
        //       8-------------11  2.0
        //      /.             /|
        //     / .            / |
        //    /  .           /  |
        //   /   .          /   |
        //  9-------------10    |
        //  |    .         |    |
        //  |    4---------|----7  1.0
        //  |   /. [1]     |   /|
        //  |  / . (2)     |  / |
        //  | /  .         | /  |
        //  |/   .         |/   |
        //  5--------------6    |          z
        //  |    .         |    |          ↑
        //  |    0---------|----3  0.0     o → y
        //  |   /  [0]     |   /          ↙
        //  |  /   (1)     |  /          x
        //  | /            | /
        //  |/             |/
        //  1--------------2   1.0
        // 0.0            1.0
        let mesh = Samples::two_hex8();
        let point_ids = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
        let sorted_points = mesh.get_sorted_points(&point_ids, |_, _, _| true);
        assert_eq!(sorted_points, vec![0, 1, 3, 2, 4, 5, 7, 6, 8, 9, 11, 10]);
    }

    #[test]
    fn test_get_sorted_points_with_filter() {
        //       8-------------11  2.0
        //      /.             /|
        //     / .            / |
        //    /  .           /  |
        //   /   .          /   |
        //  9-------------10    |
        //  |    .         |    |
        //  |    4---------|----7  1.0
        //  |   /. [1]     |   /|
        //  |  / . (2)     |  / |
        //  | /  .         | /  |
        //  |/   .         |/   |
        //  5--------------6    |          z
        //  |    .         |    |          ↑
        //  |    0---------|----3  0.0     o → y
        //  |   /  [0]     |   /          ↙
        //  |  /   (1)     |  /          x
        //  | /            | /
        //  |/             |/
        //  1--------------2   1.0
        // 0.0            1.0
        let mesh = Samples::two_hex8();
        let point_ids = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
        let sorted_points = mesh.get_sorted_points(&point_ids, |x, y, z| x > 0.0 && y > 0.0 && z > 0.0);
        assert_eq!(sorted_points, vec![6, 10]);
    }
}
