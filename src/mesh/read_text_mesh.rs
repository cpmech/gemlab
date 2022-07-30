use super::{Cell, Mesh, Point, PointId};
use crate::shapes::GeoKind;
use crate::StrError;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

struct DataForReadTextMesh {
    ndim: usize,
    npoint: usize,
    ncell: usize,
    current_npoint: usize,
    current_ncell: usize,
}

impl DataForReadTextMesh {
    fn new() -> Self {
        DataForReadTextMesh {
            ndim: 0,
            npoint: 0,
            ncell: 0,
            current_npoint: 0,
            current_ncell: 0,
        }
    }

    fn parse_sizes(&mut self, line: &str) -> Result<bool, StrError> {
        let maybe_data = line.trim_start().trim_end_matches("\n");
        if maybe_data.starts_with("#") || maybe_data == "" {
            return Ok(false); // ignore comments or empty lines; returns false == not parsed
        }

        let mut data = maybe_data.split_whitespace();

        self.ndim = data
            .next()
            .unwrap() // must panic because no error expected here
            .parse()
            .map_err(|_| "cannot parse ndim")?;

        self.npoint = match data.next() {
            Some(v) => v.parse().map_err(|_| "cannot parse npoint")?,
            None => return Err("cannot read npoint"),
        };

        self.ncell = match data.next() {
            Some(v) => v.parse().map_err(|_| "cannot parse ncell")?,
            None => return Err("cannot read ncell"),
        };

        Ok(true) // returns true == parsed
    }

    fn parse_point(&mut self, mesh: &mut Mesh, line: &str) -> Result<bool, StrError> {
        let maybe_data = line.trim_start().trim_end_matches("\n");
        if maybe_data.starts_with("#") || maybe_data == "" {
            return Ok(false); // ignore comments or empty lines
        }

        let mut data = maybe_data.split_whitespace();

        let id: usize = data
            .next()
            .unwrap() // must panic because no error expected here
            .parse()
            .map_err(|_| "cannot parse point id")?;

        if id != self.current_npoint {
            return Err("the id and index of points must equal each other");
        }

        let mut coords = vec![0.0; self.ndim];

        coords[0] = match data.next() {
            Some(v) => v.parse().map_err(|_| "cannot parse point x coordinate")?,
            None => return Err("cannot read point x coordinate"),
        };

        coords[1] = match data.next() {
            Some(v) => v.parse().map_err(|_| "cannot parse point y coordinate")?,
            None => return Err("cannot read point y coordinate"),
        };

        if self.ndim == 3 {
            coords[2] = match data.next() {
                Some(v) => v.parse().map_err(|_| "cannot parse point z coordinate")?,
                None => return Err("cannot read point z coordinate"),
            };
        }

        if data.next() != None {
            return Err("point data contains extra values");
        }

        mesh.points.push(Point { id, coords });

        self.current_npoint += 1; // next point

        Ok(true) // returns true == parsed
    }

    fn parse_cell(&mut self, mesh: &mut Mesh, line: &str) -> Result<bool, StrError> {
        let maybe_data = line.trim_start().trim_end_matches("\n");
        if maybe_data.starts_with("#") || maybe_data == "" {
            return Ok(false); // ignore comments or empty lines
        }

        let mut data = maybe_data.split_whitespace();

        let id: usize = data
            .next()
            .unwrap() // must panic because no error expected here
            .parse()
            .map_err(|_| "cannot parse cell id")?;

        if id != self.current_ncell {
            return Err("the id and index of cells must equal each other");
        }

        let attribute_id: usize = match data.next() {
            Some(v) => v.parse().map_err(|_| "cannot parse cell attribute id")?,
            None => return Err("cannot read cell attribute id"),
        };

        let str_kind = match data.next() {
            Some(v) => v,
            None => return Err("cannot read cell kind"),
        };

        let kind = GeoKind::from(&str_kind)?;
        let mut points: Vec<PointId> = vec![0; kind.nnode()];

        for m in 0..kind.nnode() {
            match data.next() {
                Some(v) => {
                    let point_id: usize = v.parse().map_err(|_| "cannot parse cell point id")?;
                    points[m] = point_id;
                }
                None => return Err("cannot read cell point id"),
            }
        }

        if data.next() != None {
            return Err("cell data contains extra values");
        }

        mesh.cells.push(Cell {
            id,
            attribute_id,
            kind,
            points,
        });

        self.current_ncell += 1; // next cell

        Ok(true) // returns true == parsed
    }
}

impl Mesh {
    /// Allocates a new Mesh by reading raw mesh data from a text file
    ///
    /// # File format
    ///
    /// The text file format includes three sections:
    ///
    /// 1. The header with the space dimension (`ndim`), number of points (`npoint`), and number of cells (`ncell`);
    /// 2. The points list where each line contains the `id` of the point, which must be **equal to the position** in the list,
    ///    followed by the `x` and `y` (and `z`) coordinates;
    /// 3. The cells list where each line contains the `id` of the cell, which must be **equal to the position** in the list,
    ///    the attribute ID (`att`) of the cell, the `kind` of the cell, followed by the IDs of the points that define the cell (connectivity).
    ///
    /// The text file looks like this (the hash tag indicates a comment/the mesh below is just an example which won't work):
    ///
    /// ```text
    /// ## header
    /// ## ndim npoint ncell
    ///      2      8     5
    ///
    /// ## points
    /// ## id    x   y
    ///    0  0.0 0.0
    ///    1  0.5 0.0
    ///    2  1.0 0.0
    /// ## ... more points should follow
    ///
    /// ## cells
    /// ## id att kind  point_ids...
    ///    0   1 tri3  0 1 3
    ///    1   1 qua4  1 4 6 3
    /// ```
    ///
    /// where we can see that different cell (shape) kinds can be present in the same mesh.
    /// However, this function does not check for element compatibility as required by finite element analyses.
    ///
    /// See [GeoKind::from] for the keys used to identify the cell kind (the keys are **lowercase**).
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn from_text_file<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let input = File::open(path).map_err(|_| "cannot open file")?;
        let buffered = BufReader::new(input);
        let mut lines_iter = buffered.lines();

        // auxiliary data structure
        let mut data = DataForReadTextMesh::new();

        // read and parse sizes
        loop {
            match lines_iter.next() {
                Some(v) => {
                    let line = v.unwrap(); // must panic because no error expected here
                    if data.parse_sizes(&line)? {
                        break;
                    }
                }
                None => return Err("file is empty or header is missing"),
            }
        }

        // allocate mesh
        let mut mesh = Mesh {
            ndim: data.ndim,
            points: Vec::new(),
            cells: Vec::new(),
        };

        // read and parse points
        loop {
            match lines_iter.next() {
                Some(v) => {
                    let line = v.unwrap(); // must panic because no error expected here
                    if data.parse_point(&mut mesh, &line)? {
                        if data.current_npoint == data.npoint {
                            break;
                        }
                    }
                }
                None => break,
            }
        }

        // check data
        if data.current_npoint != data.npoint {
            return Err("not all points have been found");
        }

        // read and parse cells
        loop {
            match lines_iter.next() {
                Some(v) => {
                    let line = v.unwrap(); // must panic because no error expected here
                    if data.parse_cell(&mut mesh, &line)? {
                        if data.current_ncell == data.ncell {
                            break;
                        }
                    }
                }
                None => break,
            }
        }

        // check data
        if data.current_ncell != data.ncell {
            return Err("not all cells have been found");
        }

        // done
        Ok(mesh)
    }

    /// Allocates a new Mesh by parsing raw mesh data from a text string
    ///
    /// # Text format
    ///
    /// The text includes three sections:
    ///
    /// 1. The header with the space dimension (`ndim`), number of points (`npoint`), and number of cells (`ncell`);
    /// 2. The points list where each line contains the `id` of the point, which must be **equal to the position** in the list,
    ///    followed by the `x` and `y` (and `z`) coordinates;
    /// 3. The cells list where each line contains the `id` of the cell, which must be **equal to the position** in the list,
    ///    the attribute ID (`att`) of the cell, the `kind` of the cell, followed by the IDs of the points that define the cell (connectivity).
    ///
    /// The text looks like this (the hash tag indicates a comment/the mesh below is just an example which won't work):
    ///
    /// ```text
    /// ## header
    /// ## ndim npoint ncell
    ///      2      8     5
    ///
    /// ## points
    /// ## id    x   y
    ///    0  0.0 0.0
    ///    1  0.5 0.0
    ///    2  1.0 0.0
    /// ## ... more points should follow
    ///
    /// ## cells
    /// ## id att kind  point_ids...
    ///    0   1 tri3  0 1 3
    ///    1   1 qua4  1 4 6 3
    /// ```
    ///
    /// where we can see that different cell (shape) kinds can be present in the same mesh.
    /// However, this function does not check for element compatibility as required by finite element analyses.
    ///
    /// See [GeoKind::from] for the keys used to identify the cell kind (the keys are **lowercase**).
    ///
    /// # Examples
    ///
    /// See `examples` and `data/meshes` directories for more examples
    /// (with pretty formatted strings).
    ///
    /// ```
    /// use gemlab::mesh::Mesh;
    /// use gemlab::StrError;
    ///
    /// fn main() -> Result<(), StrError> {
    ///     // 1.0  3-------2
    ///     //      |`. [1] |
    ///     //      |  `.   |
    ///     //      | [0]`. |
    ///     // 0.0  0------`1
    ///     //     0.0     1.0
    ///     let mesh = Mesh::from_text(
    ///         "2 4 2\n# points\n0 0.0 0.0\n1 1.0 0.0\n2 1.0 1.0\n3 0.0 1.0\n# cells\n0 1 tri3  0 1 3\n1 1 tri3  2 3 1\n",
    ///     )?;
    ///     assert_eq!(mesh.points.len(), 4);
    ///     assert_eq!(mesh.cells.len(), 2);
    ///     assert_eq!(mesh.cells[0].points, &[0, 1, 3]);
    ///     assert_eq!(mesh.cells[1].points, &[2, 3, 1]);
    ///     Ok(())
    /// }
    /// ```
    pub fn from_text(text: &str) -> Result<Self, StrError> {
        // auxiliary data structure
        let mut data = DataForReadTextMesh::new();

        // read and parse sizes
        let mut lines_iter = text.lines();
        loop {
            match lines_iter.next() {
                Some(line) => {
                    if data.parse_sizes(line)? {
                        break;
                    }
                }
                None => return Err("text string is empty or header is missing"),
            }
        }

        // allocate mesh
        let mut mesh = Mesh {
            ndim: data.ndim,
            points: Vec::new(),
            cells: Vec::new(),
        };

        // read and parse points
        loop {
            match lines_iter.next() {
                Some(line) => {
                    if data.parse_point(&mut mesh, line)? {
                        if data.current_npoint == data.npoint {
                            break;
                        }
                    }
                }
                None => break,
            }
        }

        // check data
        if data.current_npoint != data.npoint {
            return Err("not all points have been found");
        }

        // read and parse cells
        loop {
            match lines_iter.next() {
                Some(line) => {
                    if data.parse_cell(&mut mesh, line)? {
                        if data.current_ncell == data.ncell {
                            break;
                        }
                    }
                }
                None => break,
            }
        }

        // check data
        if data.current_ncell != data.ncell {
            return Err("not all cells have been found");
        }

        // done
        Ok(mesh)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::DataForReadTextMesh;
    use crate::mesh::{Mesh, Samples};

    #[test]
    fn parse_sizes_captures_errors() {
        let mut data = DataForReadTextMesh::new();

        assert_eq!(
            data.parse_sizes(&String::from(" wrong \n")).err(),
            Some("cannot parse ndim")
        );

        assert_eq!(
            data.parse_sizes(&String::from(" 2 \n")).err(),
            Some("cannot read npoint")
        );
        assert_eq!(
            data.parse_sizes(&String::from(" 1 wrong")).err(),
            Some("cannot parse npoint")
        );

        assert_eq!(
            data.parse_sizes(&String::from(" 2 4   \n")).err(),
            Some("cannot read ncell")
        );
        assert_eq!(
            data.parse_sizes(&String::from(" 2 4  wrong")).err(),
            Some("cannot parse ncell")
        );
    }

    #[test]
    fn parse_point_captures_errors() {
        let mut data = DataForReadTextMesh::new();
        data.ndim = 3;
        data.npoint = 2;
        data.ncell = 1;

        let mut mesh = Mesh {
            ndim: data.ndim,
            points: Vec::new(),
            cells: Vec::new(),
        };

        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" wrong \n")).err(),
            Some("cannot parse point id")
        );

        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 2 0.0 0.0 0.0 \n")).err(),
            Some("the id and index of points must equal each other")
        );

        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0    \n")).err(),
            Some("cannot read point x coordinate")
        );
        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0   wrong")).err(),
            Some("cannot parse point x coordinate")
        );

        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0  0.0  \n")).err(),
            Some("cannot read point y coordinate")
        );
        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0  0.0 wrong")).err(),
            Some("cannot parse point y coordinate")
        );

        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0  0.0 0.0 \n")).err(),
            Some("cannot read point z coordinate")
        );
        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0  0.0 0.0 wrong")).err(),
            Some("cannot parse point z coordinate")
        );
    }

    #[test]
    fn parse_cell_captures_errors() {
        let mut data = DataForReadTextMesh::new();
        data.ndim = 3;
        data.npoint = 2;
        data.ncell = 1;

        let mut mesh = Mesh {
            ndim: data.ndim,
            points: Vec::new(),
            cells: Vec::new(),
        };

        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" wrong \n")).err(),
            Some("cannot parse cell id")
        );

        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 2 1 0.0 0.0 0.0 \n")).err(),
            Some("the id and index of cells must equal each other")
        );

        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 \n")).err(),
            Some("cannot read cell attribute id")
        );
        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 wrong")).err(),
            Some("cannot parse cell attribute id")
        );

        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1 ")).err(),
            Some("cannot read cell kind")
        );
        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1  wrong")).err(),
            Some("string representation of GeoKind is incorrect")
        );
        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1 Lin2")).err(), // must be lowercase
            Some("string representation of GeoKind is incorrect")
        );

        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1  tri3  ")).err(),
            Some("cannot read cell point id")
        );
        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1  tri3  0 1 wrong")).err(),
            Some("cannot parse cell point id")
        );

        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1  tri3  0 1 2  extra"))
                .err(),
            Some("cell data contains extra values")
        );
    }

    #[test]
    fn from_text_file_captures_errors() {
        assert_eq!(
            Mesh::from_text_file(&String::from("__wrong__")).err(),
            Some("cannot open file")
        );
        assert_eq!(
            Mesh::from_text_file(&String::from("./data/meshes/bad_empty.msh")).err(),
            Some("file is empty or header is missing")
        );
        assert_eq!(
            Mesh::from_text_file(&String::from("./data/meshes/bad_extra_cell_data.msh")).err(),
            Some("cell data contains extra values")
        );
        assert_eq!(
            Mesh::from_text_file(&String::from("./data/meshes/bad_extra_point_data.msh")).err(),
            Some("point data contains extra values")
        );
        assert_eq!(
            Mesh::from_text_file(&String::from("./data/meshes/bad_missing_header.msh")).err(),
            Some("file is empty or header is missing")
        );
        assert_eq!(
            Mesh::from_text_file(&String::from("./data/meshes/bad_missing_points.msh")).err(),
            Some("not all points have been found")
        );
        assert_eq!(
            Mesh::from_text_file(&String::from("./data/meshes/bad_missing_cells.msh")).err(),
            Some("not all cells have been found")
        );
        assert_eq!(
            Mesh::from_text_file(&String::from("./data/meshes/bad_wrong_cell_kind.msh")).err(),
            Some("string representation of GeoKind is incorrect")
        );
    }

    #[test]
    fn from_text_file_works() {
        let mesh = Mesh::from_text_file("./data/meshes/two_quads_horizontal.msh").unwrap();
        let sample = Samples::two_qua4();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", sample));

        let mesh = Mesh::from_text_file("./data/meshes/two_cubes_vertical.msh").unwrap();
        let sample = Samples::two_hex8();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", sample));

        let mesh = Mesh::from_text_file("./data/meshes/mixed_shapes_2d.msh").unwrap();
        let sample = Samples::mixed_shapes_2d();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", sample));

        let mesh = Mesh::from_text_file("./data/meshes/mixed_shapes_3d.msh").unwrap();
        let sample = Samples::mixed_shapes_3d();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", sample));
    }

    #[test]
    fn from_text_captures_errors() {
        assert_eq!(
            Mesh::from_text(
                "# header\n\
                 # ndim npoint ncell\n"
            )
            .err(),
            Some("text string is empty or header is missing")
        );

        assert_eq!(
            Mesh::from_text(
                "# header\n\
                 # ndim npoint ncell\n\
                      2      4     1\n\
                 \n\
                 # points\n\
                 # id   x   y\n\
                    0 0.0 0.0\n"
            )
            .err(),
            Some("not all points have been found")
        );

        assert_eq!(
            Mesh::from_text(
                "# header\n\
                 # ndim npoint ncell\n\
                      2      6     2\n\
                 \n\
                 # points\n\
                 # id   x   y\n\
                    0 0.0 0.0\n\
                    1 1.0 0.0\n\
                    2 1.0 1.0\n\
                    3 0.0 1.0\n\
                    4 2.0 0.0\n\
                    5 2.0 1.0\n\
                 \n\
                 # cells\n\
                 # id att kind  point_ids...\n\
                    0   1 qua4  0 1 2 3\n"
            )
            .err(),
            Some("not all cells have been found")
        );

        assert_eq!(
            Mesh::from_text(
                "# header\n\
                 # ndim npoint ncell\n\
                      2      4     1\n\
                 # points\n\
                 # id wrong   x   y\n\
                    0     1 0.0 0.0\n\
                    1     1 1.0 0.0\n\
                    2     1 1.0 1.0\n\
                    3     1 0.0 1.0\n\
                 # cells\n\
                 # id att kind  point_ids...\n\
                    0   1 qua4  0 1 2 3\n"
            )
            .err(),
            Some("point data contains extra values")
        );

        assert_eq!(
            Mesh::from_text(
                "# header\n\
                 # ndim npoint ncell\n\
                      2      4     1\n\
                 # points\n\
                 # id   x   y\n\
                    0 0.0 0.0\n\
                    1 1.0 0.0\n\
                    2 1.0 1.0\n\
                    3 0.0 1.0\n\
                 # cells\n\
                 # id att kind  point_ids + wrong...\n\
                    0   1 qua4  0 1 2 3       4\n"
            )
            .err(),
            Some("cell data contains extra values")
        );

        assert_eq!(
            Mesh::from_text(
                "# header\n\
                 # ndim npoint ncell\n\
                      2      4     1\n\
                 # points\n\
                 # id   x   y\n\
                    0 0.0 0.0\n\
                    1 1.0 0.0\n\
                    2 1.0 1.0\n\
                    3 0.0 1.0\n\
                 # cells\n\
                 # id att kind  point_ids...\n\
                    0   1 Qua4  0 1 2 3     \n"
            )
            .err(),
            Some("string representation of GeoKind is incorrect")
        );
    }

    #[test]
    fn from_text_works() {
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
            # id att kind  point_ids...
               0   1 qua4  0 1 2 3
               1   2 qua4  1 4 5 2",
        )
        .unwrap();
        let sample = Samples::two_qua4();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", sample));

        let mesh = Mesh::from_text(
            r"# header
            # ndim npoint ncell
                 3     12     2
            
            # points
            # id    x   y   z
               0  0.0 0.0 0.0
               1  1.0 0.0 0.0
               2  1.0 1.0 0.0
               3  0.0 1.0 0.0
               4  0.0 0.0 1.0
               5  1.0 0.0 1.0
               6  1.0 1.0 1.0
               7  0.0 1.0 1.0
               8  0.0 0.0 2.0
               9  1.0 0.0 2.0
              10  1.0 1.0 2.0
              11  0.0 1.0 2.0
            
            # cells
            # id att kind  point_ids...
               0   1 hex8  0 1 2 3 4 5  6  7
               1   2 hex8  4 5 6 7 8 9 10 11",
        )
        .unwrap();
        let sample = Samples::two_hex8();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", sample));
    }
}
