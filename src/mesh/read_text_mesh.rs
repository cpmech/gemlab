use super::{Cell, Mesh, Point, PointId};
use crate::StrError;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

struct DataForReadTextMesh {
    space_ndim: usize,
    npoint: usize,
    ncell: usize,
    current_npoint: usize,
    current_ncell: usize,
}

impl DataForReadTextMesh {
    fn new() -> Self {
        DataForReadTextMesh {
            space_ndim: 0,
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

        self.space_ndim = data
            .next()
            .unwrap() // must panic because no error expected here
            .parse()
            .map_err(|_| "cannot parse space_ndim")?;

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

        let mut coords = vec![0.0; self.space_ndim];

        coords[0] = match data.next() {
            Some(v) => v.parse().map_err(|_| "cannot parse point x coordinate")?,
            None => return Err("cannot read point x coordinate"),
        };

        coords[1] = match data.next() {
            Some(v) => v.parse().map_err(|_| "cannot parse point y coordinate")?,
            None => return Err("cannot read point y coordinate"),
        };

        if self.space_ndim == 3 {
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

        let geo_ndim: usize = match data.next() {
            Some(v) => v.parse().map_err(|_| "cannot parse cell geo_ndim")?,
            None => return Err("cannot read cell geo_ndim"),
        };

        let nnode: usize = match data.next() {
            Some(v) => v.parse().map_err(|_| "cannot parse cell nnode")?,
            None => return Err("cannot read cell nnode"),
        };

        let mut points: Vec<PointId> = vec![0; nnode];

        for m in 0..nnode {
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
            geo_ndim,
            points,
        });

        self.current_ncell += 1; // next cell

        Ok(true) // returns true == parsed
    }
}

impl Mesh {
    /// Allocates a new Mesh by reading raw mesh data from a text file
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
            space_ndim: data.space_ndim,
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
            space_ndim: data.space_ndim,
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
    use crate::StrError;

    #[test]
    fn parse_sizes_captures_errors() -> Result<(), StrError> {
        let mut data = DataForReadTextMesh::new();

        assert_eq!(
            data.parse_sizes(&String::from(" wrong \n")).err(),
            Some("cannot parse space_ndim")
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
        Ok(())
    }

    #[test]
    fn parse_point_captures_errors() -> Result<(), StrError> {
        let mut data = DataForReadTextMesh::new();
        data.space_ndim = 3;
        data.npoint = 2;
        data.ncell = 1;

        let mut mesh = Mesh {
            space_ndim: data.space_ndim,
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
        Ok(())
    }

    #[test]
    fn parse_cell_captures_errors() -> Result<(), StrError> {
        let mut data = DataForReadTextMesh::new();
        data.space_ndim = 3;
        data.npoint = 2;
        data.ncell = 1;

        let mut mesh = Mesh {
            space_ndim: data.space_ndim,
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
            Some("cannot read cell geo_ndim")
        );
        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1  wrong")).err(),
            Some("cannot parse cell geo_ndim")
        );

        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1 2")).err(),
            Some("cannot read cell nnode")
        );
        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1  2 wrong")).err(),
            Some("cannot parse cell nnode")
        );

        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1  2 4  ")).err(),
            Some("cannot read cell point id")
        );
        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1  2 4  0 1 2 wrong"))
                .err(),
            Some("cannot parse cell point id")
        );
        Ok(())
    }

    #[test]
    fn from_text_file_captures_errors() -> Result<(), StrError> {
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
        Ok(())
    }

    #[test]
    fn from_text_file_works() -> Result<(), StrError> {
        let mesh = Mesh::from_text_file("./data/meshes/two_quads_horizontal.msh")?;
        let sample = Samples::two_quads_horizontal();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", sample));

        let mesh = Mesh::from_text_file("./data/meshes/two_cubes_vertical.msh")?;
        let sample = Samples::two_cubes_vertical();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", sample));

        let mesh = Mesh::from_text_file("./data/meshes/mixed_shapes_2d.msh")?;
        let sample = Samples::mixed_shapes_2d();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", sample));

        let mesh = Mesh::from_text_file("./data/meshes/mixed_shapes_3d.msh")?;
        let sample = Samples::mixed_shapes_3d();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", sample));
        Ok(())
    }

    #[test]
    fn from_text_captures_errors() -> Result<(), StrError> {
        assert_eq!(
            Mesh::from_text(
                "# header\n\
                 # space_ndim npoint ncell\n"
            )
            .err(),
            Some("text string is empty or header is missing")
        );

        assert_eq!(
            Mesh::from_text(
                "# header\n\
                 # space_ndim npoint ncell\n\
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
                 # space_ndim npoint ncell\n\
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
                 # id att geo_ndim nnode  point_ids...\n\
                    0   1        2     4  0 1 2 3\n"
            )
            .err(),
            Some("not all cells have been found")
        );

        assert_eq!(
            Mesh::from_text(
                "# header\n\
                 # space_ndim npoint ncell\n\
                            2      4     1\n\
                 # points\n\
                 # id wrong   x   y\n\
                    0     1 0.0 0.0\n\
                    1     1 1.0 0.0\n\
                    2     1 1.0 1.0\n\
                    3     1 0.0 1.0\n\
                 # cells\n\
                 # id att geo_ndim nnode  point_ids...\n\
                    0   1        2     4  0 1 2 3\n"
            )
            .err(),
            Some("point data contains extra values")
        );

        assert_eq!(
            Mesh::from_text(
                "# header\n\
                 # space_ndim npoint ncell\n\
                            2      4     1\n\
                 # points\n\
                 # id   x   y\n\
                    0 0.0 0.0\n\
                    1 1.0 0.0\n\
                    2 1.0 1.0\n\
                    3 0.0 1.0\n\
                 # cells\n\
                 # id att geo_ndim nnode  point_ids + wrong...\n\
                    0   1        2     4  0 1 2 3       4\n"
            )
            .err(),
            Some("cell data contains extra values")
        );
        Ok(())
    }

    #[test]
    fn from_text_works() -> Result<(), StrError> {
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
               1   2        2     4  1 4 5 2",
        )?;
        let sample = Samples::two_quads_horizontal();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", sample));

        let mesh = Mesh::from_text(
            r"# header
            # space_ndim npoint ncell
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
            # id att geo_ndim nnode  point_ids...
               0   1        3     8  0 1 2 3 4 5  6  7
               1   2        3     8  4 5 6 7 8 9 10 11",
        )?;
        let sample = Samples::two_cubes_vertical();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", sample));
        Ok(())
    }

    /*
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
    */
}
