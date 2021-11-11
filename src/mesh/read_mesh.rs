use super::Mesh;
use crate::StrError;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

struct DataForReadMesh {
    ndim: usize,
    npoint: usize,
    ncell: usize,
    current_npoint: usize,
    current_ncell: usize,
}

impl DataForReadMesh {
    fn new() -> Self {
        DataForReadMesh {
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

        match data.next() {
            Some(v) => self.npoint = v.parse().map_err(|_| "cannot parse npoint")?,
            None => return Err("cannot read npoint"),
        };

        match data.next() {
            Some(v) => self.ncell = v.parse().map_err(|_| "cannot parse ncell")?,
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

        let i: usize = data
            .next()
            .unwrap() // must panic because no error expected here
            .parse()
            .map_err(|_| "cannot parse point id")?;

        if i != self.current_npoint {
            return Err("the id and index of points must equal each other");
        }

        mesh.points[i].id = i;

        match data.next() {
            Some(v) => mesh.points[i].coords[0] = v.parse().map_err(|_| "cannot parse point x coordinate")?,
            None => return Err("cannot read point x coordinate"),
        };

        match data.next() {
            Some(v) => mesh.points[i].coords[1] = v.parse().map_err(|_| "cannot parse point y coordinate")?,
            None => return Err("cannot read point y coordinate"),
        };

        if self.ndim == 3 {
            match data.next() {
                Some(v) => mesh.points[i].coords[2] = v.parse().map_err(|_| "cannot parse point z coordinate")?,
                None => return Err("cannot read point z coordinate"),
            };
        }

        self.current_npoint += 1; // next point

        Ok(true) // returns true == parsed
    }

    fn parse_cell(&mut self, mesh: &mut Mesh, line: &str) -> Result<bool, StrError> {
        let maybe_data = line.trim_start().trim_end_matches("\n");
        if maybe_data.starts_with("#") || maybe_data == "" {
            return Ok(false); // ignore comments or empty lines
        }

        let mut data = maybe_data.split_whitespace();

        let i: usize = data
            .next()
            .unwrap() // must panic because no error expected here
            .parse()
            .map_err(|_| "cannot parse cell id")?;

        if i != self.current_ncell {
            return Err("the id and index of cells must equal each other");
        }

        mesh.cells[i].id = i;

        match data.next() {
            Some(v) => mesh.cells[i].attribute_id = v.parse().map_err(|_| "cannot parse cell attribute id")?,
            None => return Err("cannot read cell attribute id"),
        };

        match data.next() {
            Some(v) => mesh.cells[i].shape_ndim = v.parse().map_err(|_| "cannot parse cell ndim")?,
            None => return Err("cannot read cell ndim"),
        };

        let cell_npoint: usize = match data.next() {
            Some(v) => v.parse().map_err(|_| "cannot parse cell npoint")?,
            None => return Err("cannot read cell npoint"),
        };

        for _ in 0..cell_npoint {
            match data.next() {
                Some(v) => {
                    let point_id: usize = v.parse().map_err(|_| "cannot parse cell point id")?;
                    mesh.cells[i].points.push(point_id);
                }
                None => return Err("cannot read cell point id"),
            }
        }

        self.current_ncell += 1; // next cell

        Ok(true) // returns true == parsed
    }
}

/// Reads raw mesh data from text file
pub(super) fn read_mesh<P>(full_path: &P) -> Result<Mesh, StrError>
where
    P: AsRef<OsStr> + ?Sized,
{
    let path = Path::new(full_path).to_path_buf();
    let input = File::open(path).map_err(|_| "cannot open file")?;
    let buffered = BufReader::new(input);
    let mut lines_iter = buffered.lines();

    // auxiliary data structure
    let mut data = DataForReadMesh::new();

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
    let mut mesh = Mesh::new_sized(data.ndim, data.npoint, data.ncell)?;

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

/// Parses raw mesh data from text string
pub(super) fn parse_mesh(text: &str) -> Result<Mesh, StrError> {
    // auxiliary data structure
    let mut data = DataForReadMesh::new();

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
    let mut mesh = Mesh::new_sized(data.ndim, data.npoint, data.ncell)?;

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{parse_mesh, read_mesh, DataForReadMesh, Mesh, StrError};

    #[test]
    fn parse_sizes_captures_errors() -> Result<(), StrError> {
        let mut data = DataForReadMesh::new();

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
        Ok(())
    }

    #[test]
    fn parse_point_captures_errors() -> Result<(), StrError> {
        let mut data = DataForReadMesh::new();
        data.ndim = 3;
        data.npoint = 2;
        data.ncell = 1;

        let mut mesh = Mesh::new_sized(data.ndim, data.npoint, data.ncell)?;

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
        let mut data = DataForReadMesh::new();
        data.ndim = 3;
        data.npoint = 2;
        data.ncell = 1;

        let mut mesh = Mesh::new_sized(data.ndim, data.npoint, data.ncell)?;

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
            Some("cannot read cell ndim")
        );
        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1  wrong")).err(),
            Some("cannot parse cell ndim")
        );

        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1 2")).err(),
            Some("cannot read cell npoint")
        );
        assert_eq!(
            data.parse_cell(&mut mesh, &String::from(" 0 1  2 wrong")).err(),
            Some("cannot parse cell npoint")
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
    fn read_mesh_handle_wrong_files() -> Result<(), StrError> {
        assert_eq!(read_mesh(&String::from("__wrong__")).err(), Some("cannot open file"));
        assert_eq!(
            read_mesh(&String::from("./data/meshes/bad_empty.msh")).err(),
            Some("file is empty or header is missing")
        );
        assert_eq!(
            read_mesh(&String::from("./data/meshes/bad_missing_header.msh")).err(),
            Some("file is empty or header is missing")
        );
        assert_eq!(
            read_mesh(&String::from("./data/meshes/bad_missing_points.msh")).err(),
            Some("not all points have been found")
        );
        assert_eq!(
            read_mesh(&String::from("./data/meshes/bad_missing_cells.msh")).err(),
            Some("not all cells have been found")
        );
        Ok(())
    }

    #[test]
    fn read_mesh_2d_works() -> Result<(), StrError> {
        let mesh = read_mesh("./data/meshes/ok1.msh")?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 2\n\
             npoint = 6\n\
             ncell = 2\n\
             n_boundary_point = 0\n\
             n_boundary_edge = 0\n\
             n_boundary_face = 0\n\
             \n\
             points\n\
             i:0 x:[0.0, 0.0] e:[] f:[]\n\
             i:1 x:[1.0, 0.0] e:[] f:[]\n\
             i:2 x:[1.0, 1.0] e:[] f:[]\n\
             i:3 x:[0.0, 1.0] e:[] f:[]\n\
             i:4 x:[2.0, 0.0] e:[] f:[]\n\
             i:5 x:[2.0, 1.0] e:[] f:[]\n\
             \n\
             cells\n\
             i:0 a:1 n:2 p:[0, 1, 2, 3]\n\
             i:1 a:0 n:2 p:[1, 4, 5, 2]\n\
             \n\
             boundary_points\n\
             \n\
             boundary_edges\n\
             \n\
             boundary_faces\n"
        );
        Ok(())
    }

    #[test]
    fn read_mesh_3d_works() -> Result<(), StrError> {
        let mesh = read_mesh("./data/meshes/ok2.msh")?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 3\n\
             npoint = 12\n\
             ncell = 2\n\
             n_boundary_point = 0\n\
             n_boundary_edge = 0\n\
             n_boundary_face = 0\n\
             \n\
             points\n\
             i:0 x:[0.0, 0.0, 0.0] e:[] f:[]\n\
             i:1 x:[1.0, 0.0, 0.0] e:[] f:[]\n\
             i:2 x:[1.0, 1.0, 0.0] e:[] f:[]\n\
             i:3 x:[0.0, 1.0, 0.0] e:[] f:[]\n\
             i:4 x:[0.0, 0.0, 1.0] e:[] f:[]\n\
             i:5 x:[1.0, 0.0, 1.0] e:[] f:[]\n\
             i:6 x:[1.0, 1.0, 1.0] e:[] f:[]\n\
             i:7 x:[0.0, 1.0, 1.0] e:[] f:[]\n\
             i:8 x:[0.0, 0.0, 2.0] e:[] f:[]\n\
             i:9 x:[1.0, 0.0, 2.0] e:[] f:[]\n\
             i:10 x:[1.0, 1.0, 2.0] e:[] f:[]\n\
             i:11 x:[0.0, 1.0, 2.0] e:[] f:[]\n\
             \n\
             cells\n\
             i:0 a:1 n:3 p:[0, 1, 2, 3, 4, 5, 6, 7]\n\
             i:1 a:0 n:3 p:[4, 5, 6, 7, 8, 9, 10, 11]\n\
             \n\
             boundary_points\n\
             \n\
             boundary_edges\n\
             \n\
             boundary_faces\n"
        );
        Ok(())
    }

    #[test]
    fn parse_mesh_handle_wrong_data() -> Result<(), StrError> {
        assert_eq!(
            parse_mesh(
                "# header\n\
                 # ndim npoint ncell\n"
            )
            .err(),
            Some("text string is empty or header is missing")
        );

        assert_eq!(
            parse_mesh(
                "# header\n\
                 # ndim npoint ncell\n\
                 2 4 1\n\
                 \n\
                 # points\n\
                 # id x y\n\
                 0 0.0 0.0\n"
            )
            .err(),
            Some("not all points have been found")
        );

        assert_eq!(
            parse_mesh(
                "# header\n\
                 # ndim npoint ncell\n\
                 2 6 2\n\
                 \n\
                 # points\n\
                 # id x y\n\
                 0  0.0 0.0\n\
                 1  1.0 0.0\n\
                 2  1.0 1.0\n\
                 3  0.0 1.0\n\
                 4  2.0 0.0\n\
                 5  2.0 1.0\n\
                 \n\
                 # cells\n\
                 # idx attribute point_ids...\n\
                 0 1  2 4  0 1 2 3\n"
            )
            .err(),
            Some("not all cells have been found")
        );
        Ok(())
    }

    #[test]
    fn parse_mesh_2d_works() -> Result<(), StrError> {
        let mesh = parse_mesh(
            r"# header
            # ndim npoint ncell
            2 6 2
            
            # points
            # id x y
            0  0.0 0.0
            1  1.0 0.0
            2  1.0 1.0
            3  0.0 1.0
            4  2.0 0.0
            5  2.0 1.0
            
            # cells
            # id attribute ndim npoint point_ids...
            0 1  2 4  0 1 2 3
            1 0  2 4  1 4 5 2",
        )?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 2\n\
             npoint = 6\n\
             ncell = 2\n\
             n_boundary_point = 0\n\
             n_boundary_edge = 0\n\
             n_boundary_face = 0\n\
             \n\
             points\n\
             i:0 x:[0.0, 0.0] e:[] f:[]\n\
             i:1 x:[1.0, 0.0] e:[] f:[]\n\
             i:2 x:[1.0, 1.0] e:[] f:[]\n\
             i:3 x:[0.0, 1.0] e:[] f:[]\n\
             i:4 x:[2.0, 0.0] e:[] f:[]\n\
             i:5 x:[2.0, 1.0] e:[] f:[]\n\
             \n\
             cells\n\
             i:0 a:1 n:2 p:[0, 1, 2, 3]\n\
             i:1 a:0 n:2 p:[1, 4, 5, 2]\n\
             \n\
             boundary_points\n\
             \n\
             boundary_edges\n\
             \n\
             boundary_faces\n"
        );
        Ok(())
    }

    #[test]
    fn parse_mesh_3d_works() -> Result<(), StrError> {
        let mesh = parse_mesh(
            r"# header
            # ndim npoint ncell
            3 12 2
            
            # points
            # id x y z
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
            # id attribute ndim npoint point_ids...
            0 1  3 8  0 1 2 3 4 5  6  7
            1 0  3 8  4 5 6 7 8 9 10 11",
        )?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 3\n\
             npoint = 12\n\
             ncell = 2\n\
             n_boundary_point = 0\n\
             n_boundary_edge = 0\n\
             n_boundary_face = 0\n\
             \n\
             points\n\
             i:0 x:[0.0, 0.0, 0.0] e:[] f:[]\n\
             i:1 x:[1.0, 0.0, 0.0] e:[] f:[]\n\
             i:2 x:[1.0, 1.0, 0.0] e:[] f:[]\n\
             i:3 x:[0.0, 1.0, 0.0] e:[] f:[]\n\
             i:4 x:[0.0, 0.0, 1.0] e:[] f:[]\n\
             i:5 x:[1.0, 0.0, 1.0] e:[] f:[]\n\
             i:6 x:[1.0, 1.0, 1.0] e:[] f:[]\n\
             i:7 x:[0.0, 1.0, 1.0] e:[] f:[]\n\
             i:8 x:[0.0, 0.0, 2.0] e:[] f:[]\n\
             i:9 x:[1.0, 0.0, 2.0] e:[] f:[]\n\
             i:10 x:[1.0, 1.0, 2.0] e:[] f:[]\n\
             i:11 x:[0.0, 1.0, 2.0] e:[] f:[]\n\
             \n\
             cells\n\
             i:0 a:1 n:3 p:[0, 1, 2, 3, 4, 5, 6, 7]\n\
             i:1 a:0 n:3 p:[4, 5, 6, 7, 8, 9, 10, 11]\n\
             \n\
             boundary_points\n\
             \n\
             boundary_edges\n\
             \n\
             boundary_faces\n"
        );
        Ok(())
    }
}
