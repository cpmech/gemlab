use crate::Mesh;
use std::fs::File;
use std::io::{BufRead, BufReader};

struct ReadMeshData {
    ndim: usize,
    npoint: usize,
    ncell: usize,
    current_npoint: usize,
    current_ncell: usize,
}

impl ReadMeshData {
    fn new() -> Self {
        ReadMeshData {
            ndim: 0,
            npoint: 0,
            ncell: 0,
            current_npoint: 0,
            current_ncell: 0,
        }
    }

    fn parse_sizes(&mut self, line: &String) -> Result<bool, &'static str> {
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

    fn parse_point(&mut self, mesh: &mut Mesh, line: &String) -> Result<bool, &'static str> {
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
            Some(v) => mesh.points[i].group = v.parse().map_err(|_| "cannot parse point group")?,
            None => return Err("cannot read point group"),
        };

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

    fn parse_cell(&mut self, mesh: &mut Mesh, line: &String) -> Result<bool, &'static str> {
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
            Some(v) => mesh.cells[i].group = v.parse().map_err(|_| "cannot parse cell group")?,
            None => return Err("cannot read cell group"),
        };

        loop {
            match data.next() {
                Some(v) => {
                    let point_id: usize = v.parse().map_err(|_| "cannot parse point id of cell (connectivity)")?;
                    mesh.cells[i].point_ids.push(point_id);
                }
                None => break,
            }
        }

        self.current_ncell += 1; // next cell

        Ok(true) // returns true == parsed
    }
}

pub fn read_mesh(filepath: &String) -> Result<Mesh, &'static str> {
    let input = File::open(filepath).map_err(|_| "cannot open file")?;
    let buffered = BufReader::new(input);
    let mut lines_iter = buffered.lines();

    // auxiliary data structure
    let mut data = ReadMeshData::new();

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
    let mut mesh = Mesh::new_zeroed(data.ndim, data.npoint, data.ncell)?;

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{read_mesh, Mesh, ReadMeshData};

    #[test]
    fn parse_sizes_captures_errors() -> Result<(), &'static str> {
        let mut data = ReadMeshData::new();

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
    fn parse_point_captures_errors() -> Result<(), &'static str> {
        let mut data = ReadMeshData::new();
        data.ndim = 3;
        data.npoint = 2;
        data.ncell = 1;

        let mut mesh = Mesh::new_zeroed(data.ndim, data.npoint, data.ncell)?;

        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" wrong \n")).err(),
            Some("cannot parse point id")
        );

        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 2 1 0.0 0.0 0.0 \n")).err(),
            Some("the id and index of points must equal each other")
        );

        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0 \n")).err(),
            Some("cannot read point group")
        );
        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0 wrong")).err(),
            Some("cannot parse point group")
        );

        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0 1   \n")).err(),
            Some("cannot read point x coordinate")
        );
        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0 1  wrong")).err(),
            Some("cannot parse point x coordinate")
        );

        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0 1 0.0  \n")).err(),
            Some("cannot read point y coordinate")
        );
        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0 1 0.0 wrong")).err(),
            Some("cannot parse point y coordinate")
        );

        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0 1 0.0 0.0 \n")).err(),
            Some("cannot read point z coordinate")
        );
        assert_eq!(
            data.parse_point(&mut mesh, &String::from(" 0 1 0.0 0.0 wrong")).err(),
            Some("cannot parse point z coordinate")
        );
        Ok(())
    }

    #[test]
    fn read_mesh_handle_wrong_files() -> Result<(), &'static str> {
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
        /*
        assert_eq!(
            read_mesh(&String::from("./data/meshes/bad_many_points.msh")).err(),
            Some("there are more points than specified")
        );
        */
        Ok(())
    }

    #[test]
    fn read_mesh_2d_works() -> Result<(), &'static str> {
        let filepath = "./data/meshes/ok1.msh".to_string();
        let mesh = read_mesh(&filepath)?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 2\n\
             npoint = 4\n\
             ncell = 1\n\
             n_boundary_point = 0\n\
             n_boundary_edge = 0\n\
             n_boundary_face = 0\n\
             \n\
             points\n\
             i:0 g:1 x:[0.0, 0.0] c:[]\n\
             i:1 g:1 x:[1.0, 0.0] c:[]\n\
             i:2 g:1 x:[1.0, 1.0] c:[]\n\
             i:3 g:1 x:[0.0, 1.0] c:[]\n\
             \n\
             cells\n\
             i:0 g:1 p:[0, 1, 2, 3] e:[] f:[]\n\
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
    fn read_mesh_3d_works() -> Result<(), &'static str> {
        let filepath = "./data/meshes/ok2.msh".to_string();
        let mesh = read_mesh(&filepath)?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 3\n\
             npoint = 8\n\
             ncell = 1\n\
             n_boundary_point = 0\n\
             n_boundary_edge = 0\n\
             n_boundary_face = 0\n\
             \n\
             points\n\
             i:0 g:1 x:[0.0, 0.0, 0.0] c:[]\n\
             i:1 g:1 x:[1.0, 0.0, 0.0] c:[]\n\
             i:2 g:1 x:[1.0, 1.0, 0.0] c:[]\n\
             i:3 g:1 x:[0.0, 1.0, 0.0] c:[]\n\
             i:4 g:1 x:[0.0, 0.0, 1.0] c:[]\n\
             i:5 g:1 x:[1.0, 0.0, 1.0] c:[]\n\
             i:6 g:1 x:[1.0, 1.0, 1.0] c:[]\n\
             i:7 g:1 x:[0.0, 1.0, 1.0] c:[]\n\
             \n\
             cells\n\
             i:0 g:1 p:[0, 1, 2, 3, 4, 5, 6, 7] e:[] f:[]\n\
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
