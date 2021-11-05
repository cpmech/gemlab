use crate::Mesh;
use std::fs::File;
use std::io::{BufRead, BufReader};

struct ReadMeshData {
    ndim: usize,
    npoint: usize,
    ncell: usize,
    current_npoint: usize,
    current_ncell: usize,
    point_id: usize,
    point_group: usize,
    point_x: f64,
    point_y: f64,
    point_z: f64,
}

impl ReadMeshData {
    fn new() -> Self {
        ReadMeshData {
            ndim: 0,
            npoint: 0,
            ncell: 0,
            current_npoint: 0,
            current_ncell: 0,
            point_id: 0,
            point_group: 0,
            point_x: 0.0,
            point_y: 0.0,
            point_z: 0.0,
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

    fn parse_point(&mut self, line: &String) -> Result<bool, &'static str> {
        let maybe_data = line.trim_start().trim_end_matches("\n");
        if maybe_data.starts_with("#") || maybe_data == "" {
            return Ok(false); // ignore comments or empty lines
        }

        if self.current_npoint == self.npoint {
            return Err("there are more points than specified");
        }

        let mut data = maybe_data.split_whitespace();

        self.point_id = data
            .next()
            .unwrap() // must panic because no error expected here
            .parse()
            .map_err(|_| "cannot parse point id")?;

        if self.point_id != self.current_npoint {
            return Err("the id and index of points must equal each other");
        }

        match data.next() {
            Some(v) => self.point_group = v.parse().map_err(|_| "cannot parse point group")?,
            None => return Err("cannot read point group"),
        };

        match data.next() {
            Some(v) => self.point_x = v.parse().map_err(|_| "cannot parse point x coordinate")?,
            None => return Err("cannot read point x coordinate"),
        };

        match data.next() {
            Some(v) => self.point_y = v.parse().map_err(|_| "cannot parse point y coordinate")?,
            None => return Err("cannot read point y coordinate"),
        };

        if self.ndim == 3 {
            match data.next() {
                Some(v) => self.point_z = v.parse().map_err(|_| "cannot parse point z coordinate")?,
                None => return Err("cannot read point z coordinate"),
            };
        }

        self.current_npoint += 1; // next point

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
                if data.parse_point(&line)? {
                    let index = data.point_id;
                    mesh.points[index].id = data.point_id;
                    mesh.points[index].group = data.point_group;
                    mesh.points[index].coords[0] = data.point_x;
                    mesh.points[index].coords[1] = data.point_y;
                    if data.ndim == 3 {
                        mesh.points[data.point_id].coords[2] = data.point_z;
                    }
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

    // done
    Ok(mesh)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{read_mesh, ReadMeshData};

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

        data.current_npoint = 2;

        assert_eq!(
            data.parse_point(&String::from("0 1 0.0 0.0 0.0\n")).err(),
            Some("there are more points than specified")
        );

        data.current_npoint = 0;

        assert_eq!(
            data.parse_point(&String::from(" wrong \n")).err(),
            Some("cannot parse point id")
        );

        assert_eq!(
            data.parse_point(&String::from(" 0 \n")).err(),
            Some("cannot read point group")
        );
        assert_eq!(
            data.parse_point(&String::from(" 0 wrong")).err(),
            Some("cannot parse point group")
        );

        assert_eq!(
            data.parse_point(&String::from(" 0 1   \n")).err(),
            Some("cannot read point x coordinate")
        );
        assert_eq!(
            data.parse_point(&String::from(" 0 1  wrong")).err(),
            Some("cannot parse point x coordinate")
        );

        assert_eq!(
            data.parse_point(&String::from(" 0 1 0.0  \n")).err(),
            Some("cannot read point y coordinate")
        );
        assert_eq!(
            data.parse_point(&String::from(" 0 1 0.0 wrong")).err(),
            Some("cannot parse point y coordinate")
        );

        assert_eq!(
            data.parse_point(&String::from(" 0 1 0.0 0.0 \n")).err(),
            Some("cannot read point z coordinate")
        );
        assert_eq!(
            data.parse_point(&String::from(" 0 1 0.0 0.0 wrong")).err(),
            Some("cannot parse point z coordinate")
        );
        Ok(())
    }

    #[test]
    fn read_mesh_handle_wrong_files() -> Result<(), &'static str> {
        assert_eq!(read_mesh(&String::from("__wrong__")).err(), Some("cannot open file"));
        assert_eq!(
            read_mesh(&String::from("./data/meshes/bad_empty.msh")).err(),
            Some("file is empty")
        );
        assert_eq!(
            read_mesh(&String::from("./data/meshes/bad_missing_points.msh")).err(),
            Some("not all points have been found")
        );
        assert_eq!(
            read_mesh(&String::from("./data/meshes/bad_many_points.msh")).err(),
            Some("there are more points than specified")
        );
        Ok(())
    }

    #[test]
    fn read_mesh_works() -> Result<(), &'static str> {
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
             i:0 g:0 p:[] e:[] f:[]\n\
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
