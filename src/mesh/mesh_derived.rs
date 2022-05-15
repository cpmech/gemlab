/// Aliases (usize,usize) as the key of Edge
///
/// # Note
///
/// Since the local numbering scheme runs over "corners" first, we can compare
/// edges using only two points; i.e., the middle points don't matter.
pub type EdgeKey = (usize, usize);

/// Aliases (usize,usize,usize,usize) as the key of Face
///
/// # Note
///
/// If all faces have at most 3 points, the fourth entry in the key will be equal to the total number of points.
/// In this way, we can compare 4-node (or more nodes) faces with each other, since that the local numbering
/// scheme runs over the "corners" first; i.e., the middle points don't matter.
pub type FaceKey = (usize, usize, usize, usize);

/// Holds point data
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct PointDerived {
    /// Identification number which equals the index of the point in the mesh
    pub id: PointId,

    /// Set of boundary edges sharing this point
    pub shared_by_boundary_edges: HashSet<EdgeKey>,

    /// Set of boundary faces sharing this point
    pub shared_by_boundary_faces: HashSet<FaceKey>,
}

/// Holds edge data (derived data structure)
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct EdgeDerived {
    /// List of points defining this edge; in the right order (unsorted)
    pub points: Vec<PointId>,

    /// Set of 2D cells sharing this edge (to find the boundary)
    ///
    /// **2D mesh only**
    pub shared_by_2d_cells: HashSet<CellId>,

    /// Set of boundary faces sharing this edge
    pub shared_by_boundary_faces: HashSet<FaceKey>,
}

/// Holds face data (derived data structure)
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct FaceDerived {
    /// List of points defining this face; in the right order (unsorted)
    pub points: Vec<PointId>,

    /// Set of cells sharing this face
    pub shared_by_cells: HashSet<CellId>,
}

/// Holds mesh data (derived)
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct MeshDerived {
    /// Space dimension of the mesh
    ///
    /// The mesh's ndim may be different that an cell's ndim.
    /// For example, a 3D mesh may contain 1D lines or 2D triangles.
    ///
    /// **raw data**
    pub space_ndim: usize,

    /// All points in the mesh
    ///
    /// **raw data**
    pub points: Vec<Point>,

    /// All cells (aka geometric shape, polygon, polyhedra) in the mesh
    ///
    /// **raw data**
    pub cells: Vec<Cell>,

    /// Set of points on the boundaries
    ///
    /// Note: a boundary point belongs to a boundary edge or a boundary face
    ///
    /// (derived property)
    pub boundary_points: HashSet<PointId>,

    /// Set of edges on the boundaries
    ///
    /// Note: In 2D, a boundary edge is such that it is shared by one 2D cell only (1D cells are ignored)
    ///
    /// Note: In 3D, a boundary edge belongs to a boundary face
    ///
    /// (derived property)
    pub boundary_edges: HashMap<EdgeKey, Edge>,

    /// Set of faces on the boundaries
    ///
    /// Note: A boundary face is such that it is shared by one 3D cell only
    ///
    /// (derived property)
    pub boundary_faces: HashMap<FaceKey, Face>,

    /// Min coordinates (space_ndim)
    ///
    /// (derived property)
    pub coords_min: Vec<f64>,

    /// Max coordinates (space_ndim)
    ///
    /// (derived property)
    pub coords_max: Vec<f64>,

    /// Difference max minus min coordinates (space_ndim), i.e., bounding box of Mesh
    ///
    /// (derived property)
    pub coords_delta: Vec<f64>,

    /// Allows searching boundary points using their coordinates
    grid_boundary_points: GridSearch,

    /// Indicates whether the derived variables and maps have been computed or not
    derived_props_computed: bool,
}

impl MeshDerived {
    /// Returns a new empty mesh
    pub(super) fn new(space_ndim: usize) -> Result<Self, StrError> {
        if space_ndim < 2 || space_ndim > 3 {
            return Err("space_ndim must be 2 or 3");
        }
        Ok(Mesh {
            space_ndim,
            points: Vec::new(),
            cells: Vec::new(),
            boundary_points: HashSet::new(),
            boundary_edges: HashMap::new(),
            boundary_faces: HashMap::new(),
            coords_min: Vec::new(),
            coords_max: Vec::new(),
            coords_delta: Vec::new(),
            grid_boundary_points: GridSearch::new(space_ndim)?,
            derived_props_computed: false,
        })
    }

    /// Computes derived properties such as boundaries and limits
    ///
    /// # Note
    ///
    /// This function is automatically called by `from_text` and `from_text_file`.
    pub fn compute_derived_props(&mut self) -> Result<(), StrError> {
        // derived props
        self.derived_props_computed = false;
        self.boundary_points.clear();
        self.boundary_edges.clear();
        self.boundary_faces.clear();
        if self.space_ndim == 2 {
            self.compute_derived_props_2d()?;
        } else {
            self.compute_derived_props_3d()?;
        }

        // limits
        self.compute_limits()?;

        // compute number of divisions for GridSearch
        let mut index_long = 0;
        let mut delta_long = self.coords_delta[index_long];
        for i in 1..self.space_ndim {
            if self.coords_delta[i] > delta_long {
                index_long = i;
                delta_long = self.coords_delta[i];
            }
        }
        let mut ndiv = vec![0; self.space_ndim];
        for i in 0..self.space_ndim {
            if i == index_long {
                ndiv[i] = GRID_SEARCH_NDIV_LONG;
            } else {
                ndiv[i] = ((self.coords_delta[i] / delta_long) * (GRID_SEARCH_NDIV_LONG as f64)) as usize;
            }
            ndiv[i] = usize::min(GRID_SEARCH_NDIV_MAX, usize::max(GRID_SEARCH_NDIV_MIN, ndiv[i]));
        }

        // initialize GridSearch and add all boundary points to it
        self.grid_boundary_points
            .initialize(&ndiv, &self.coords_min, &self.coords_max)?;
        for point_id in &self.boundary_points {
            self.grid_boundary_points
                .insert(*point_id, &self.points[*point_id].coords)?;
        }

        // done
        self.derived_props_computed = true;
        Ok(())
    }

    /// Computes derived properties of 2D mesh
    fn compute_derived_props_2d(&mut self) -> Result<(), StrError> {
        // maps all edge keys to (cell_id, e) where e is the cell's local edge index
        let mut all_edges: HashMap<EdgeKey, Vec<(CellId, usize)>> = HashMap::new(); // (edge_key) => [(cell_id,e)]

        // maps all cell shapes to a Shape instance
        let mut all_shapes: HashMap<(usize, usize), Shape> = HashMap::new(); // (geo_ndim,nnode) => Shape

        // loop over all cells
        for cell in &mut self.cells {
            let nnode = cell.points.len();

            // handle 1D shapes
            if cell.geo_ndim != 2 {
                for m in 0..nnode {
                    self.boundary_points.insert(cell.points[m]);
                }
                continue; // skip 1D line in 2D because it's not a boundary edge
            }

            // get or allocate Shape
            let cell_shape =
                all_shapes
                    .entry((cell.geo_ndim, nnode))
                    .or_insert(Shape::new(self.space_ndim, cell.geo_ndim, nnode)?);

            // set the cell node coordinates in the shape object
            for m in 0..nnode {
                let point_id = cell.points[m];
                for j in 0..self.space_ndim {
                    cell_shape.set_node(point_id, m, j, self.points[point_id].coords[j])?;
                }
            }

            // check if the determinant of Jacobian is positive => counterclockwise nodes
            let det_jac = cell_shape.calc_jacobian(&[0.0, 0.0, 0.0])?;
            if det_jac < 0.0 {
                return Err("a cell has incorrect ordering of nodes");
            }

            // set information about all edges
            for e in 0..cell_shape.nedge {
                // define edge key (sorted point ids)
                let mut edge_key: EdgeKey = (
                    cell.points[cell_shape.get_edge_node_id(e, 0)],
                    cell.points[cell_shape.get_edge_node_id(e, 1)],
                );
                sort2(&mut edge_key);

                // configure edge_key => (cell_id, e)
                let edge_data = all_edges.entry(edge_key).or_insert(Vec::new());
                edge_data.push((cell.id, e));
            }
        }

        // loop over all edges
        for (edge_key, cell_ids_and_es) in &all_edges {
            // skip inner edges (those shared by multiple cells)
            if cell_ids_and_es.len() != 1 {
                continue;
            }

            // edge data
            let (cell_id, e) = cell_ids_and_es[0];
            let cell = &self.cells[cell_id];
            let cell_shape = all_shapes.get(&(cell.geo_ndim, cell.points.len())).unwrap(); // must exist due to previous loop
            let edge_nnode = cell_shape.edge_nnode;
            let mut edge_points: Vec<PointId> = vec![0; edge_nnode];

            // loop over all edge nodes
            for i in 0..edge_nnode {
                // configure edge nodes
                edge_points[i] = cell.points[cell_shape.get_edge_node_id(e, i)];

                // set boundary points
                self.points[edge_points[i]].shared_by_boundary_edges.insert(*edge_key);
                self.boundary_points.insert(edge_points[i]);
            }

            // append new boundary edge
            self.boundary_edges.insert(
                *edge_key,
                Edge {
                    points: edge_points,
                    shared_by_2d_cells: HashSet::from([cell_id]),
                    shared_by_boundary_faces: HashSet::new(),
                },
            );
        }
        Ok(())
    }

    /// Computes derived properties of 3D mesh
    fn compute_derived_props_3d(&mut self) -> Result<(), StrError> {
        // maps all face keys to (cell_id, f) where f is the cell's local face index
        let mut all_faces: HashMap<FaceKey, Vec<(CellId, usize)>> = HashMap::new(); // (face_key) => [(cell_id,f)]

        // maps all cell shapes to a Shape instance
        let mut all_shapes: HashMap<(usize, usize), Shape> = HashMap::new(); // (geo_ndim,nnode) => Shape

        // loop over all cells
        for cell in &mut self.cells {
            let nnode = cell.points.len();

            // handle 1D and 2D shapes
            if cell.geo_ndim != 3 {
                for m in 0..nnode {
                    self.boundary_points.insert(cell.points[m]);
                }
                continue; // skip 1D line or 2D shape in 3D because they don't have faces
            }

            // get or allocate Shape
            let cell_shape =
                all_shapes
                    .entry((cell.geo_ndim, nnode))
                    .or_insert(Shape::new(self.space_ndim, cell.geo_ndim, nnode)?);

            // set the cell node coordinates in the shape object
            for m in 0..nnode {
                let point_id = cell.points[m];
                for j in 0..self.space_ndim {
                    cell_shape.set_node(point_id, m, j, self.points[point_id].coords[j])?;
                }
            }

            // check if the determinant of Jacobian is positive => counterclockwise nodes
            let det_jac = cell_shape.calc_jacobian(&[0.0, 0.0, 0.0])?;
            if det_jac < 0.0 {
                return Err("a cell has incorrect ordering of nodes");
            }

            // set information about all faces
            for f in 0..cell_shape.nface {
                // define face key (sorted ids)
                let mut face_key: FaceKey = if cell_shape.face_nnode > 3 {
                    (
                        cell.points[cell_shape.get_face_node_id(f, 0)],
                        cell.points[cell_shape.get_face_node_id(f, 1)],
                        cell.points[cell_shape.get_face_node_id(f, 2)],
                        cell.points[cell_shape.get_face_node_id(f, 3)],
                    )
                } else {
                    (
                        cell.points[cell_shape.get_face_node_id(f, 0)],
                        cell.points[cell_shape.get_face_node_id(f, 1)],
                        cell.points[cell_shape.get_face_node_id(f, 2)],
                        self.points.len(),
                    )
                };
                sort4(&mut face_key);

                // configure face_key => (cell_id, f)
                let face_data = all_faces.entry(face_key).or_insert(Vec::new());
                face_data.push((cell.id, f));
            }
        }

        // sort face keys just so the next loop is deterministic
        let mut face_keys: Vec<_> = all_faces.keys().collect();
        face_keys.sort();

        // loop over all faces
        for face_key in face_keys {
            let cell_ids_and_fs = all_faces.get(face_key).unwrap();
            // skip inner faces (those shared by multiple cells)
            if cell_ids_and_fs.len() != 1 {
                continue;
            }

            // face data
            let (cell_id, f) = cell_ids_and_fs[0];
            let cell = &self.cells[cell_id];
            let cell_shape = all_shapes.get(&(cell.geo_ndim, cell.points.len())).unwrap(); // must exist due to previous loop
            let face_nnode = cell_shape.face_nnode;
            let mut face_points: Vec<PointId> = vec![0; face_nnode];
            let face_shape = Shape::new(self.space_ndim, 2, face_nnode)?;

            // loop over all face nodes
            for i in 0..face_nnode {
                // configure face nodes
                face_points[i] = cell.points[cell_shape.get_face_node_id(f, i)];

                // set boundary points
                self.points[face_points[i]].shared_by_boundary_faces.insert(*face_key);
                self.boundary_points.insert(face_points[i]);
            }

            // loop over all face edges
            for e in 0..face_shape.nedge {
                // define edge key (sorted point ids)
                let mut edge_key: EdgeKey = (
                    face_points[face_shape.get_edge_node_id(e, 0)],
                    face_points[face_shape.get_edge_node_id(e, 1)],
                );
                sort2(&mut edge_key);

                // handle boundary edge
                match self.boundary_edges.get_mut(&edge_key) {
                    Some(edge) => {
                        // set boundary edge information
                        edge.shared_by_boundary_faces.insert(*face_key);
                    }
                    None => {
                        // edge data
                        let edge_nnode = face_shape.edge_nnode;
                        let mut edge_points: Vec<PointId> = vec![0; edge_nnode];

                        // loop over all edge nodes
                        for i in 0..edge_nnode {
                            // configure edge nodes
                            edge_points[i] = face_points[face_shape.get_edge_node_id(e, i)];

                            // set boundary points
                            self.points[edge_points[i]].shared_by_boundary_edges.insert(edge_key);
                        }

                        // append new boundary edge
                        self.boundary_edges.insert(
                            edge_key,
                            Edge {
                                points: edge_points,
                                shared_by_2d_cells: HashSet::new(),
                                shared_by_boundary_faces: HashSet::from([*face_key]),
                            },
                        );
                    }
                }
            }

            // append new boundary face
            self.boundary_faces.insert(
                *face_key,
                Face {
                    points: face_points,
                    shared_by_cells: HashSet::from([cell_id]),
                },
            );
        }
        Ok(())
    }

    /// Computes the range of coordinates
    fn compute_limits(&mut self) -> Result<(), StrError> {
        self.coords_min = vec![f64::MAX; self.space_ndim];
        self.coords_max = vec![f64::MIN; self.space_ndim];
        self.coords_delta = vec![0.0; self.space_ndim];
        for point in &self.points {
            for i in 0..self.space_ndim {
                if point.coords[i] < self.coords_min[i] {
                    self.coords_min[i] = point.coords[i];
                }
                if point.coords[i] > self.coords_max[i] {
                    self.coords_max[i] = point.coords[i];
                }
            }
        }
        for i in 0..self.space_ndim {
            if self.coords_min[i] >= self.coords_max[i] {
                return Err("mesh limits are invalid");
            }
            self.coords_delta[i] = self.coords_max[i] - self.coords_min[i];
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {

    #[test]
    fn new_works() -> Result<(), StrError> {
        let mesh = Mesh::new(2)?;
        assert_eq!(mesh.space_ndim, 2);
        assert_eq!(mesh.points.len(), 0);
        assert_eq!(mesh.cells.len(), 0);
        assert_eq!(mesh.boundary_points.len(), 0);
        assert_eq!(mesh.boundary_edges.len(), 0);
        assert_eq!(mesh.boundary_faces.len(), 0);
        assert_eq!(mesh.coords_min.len(), 0);
        assert_eq!(mesh.coords_max.len(), 0);
        assert_eq!(
            format!("{}", mesh.grid_boundary_points),
            "ids = []\n\
             nitem = 0\n\
             ncontainer = 0\n\
             ndiv = [0, 0]\n"
        );
        assert_eq!(mesh.derived_props_computed, false);
        Ok(())
    }

    #[test]
    fn compute_limits_fails_on_wrong_data() -> Result<(), StrError> {
        let mut mesh = Mesh::new(2)?;
        mesh.points.push(Point {
            id: 0,
            coords: vec![0.0, 0.0],
            shared_by_boundary_edges: HashSet::new(),
            shared_by_boundary_faces: HashSet::new(),
        });
        assert_eq!(mesh.compute_limits().err(), Some("mesh limits are invalid"));
        Ok(())
    }

    #[test]
    fn normals_are_correct_2d() -> Result<(), StrError> {
        //
        //  3--------2--------5
        //  |        |        |
        //  |        |        |
        //  |        |        |
        //  0--------1--------4
        //
        let mesh = Mesh::from_text_file("./data/meshes/ok1.msh")?;

        // the norm of the normal vector should be equal to 0.5 = edge_length / 2.0
        // where 2.0 corresponds to the edge_length in the reference system
        let l = 0.5; // norm of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 4)], [0.0, -l]),
            // right
            (vec![(4, 5)], [l, 0.0]),
            // top
            (vec![(2, 3), (2, 5)], [0.0, l]),
            // left
            (vec![(0, 3)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let mut normal = Vector::new(mesh.space_ndim);
        let ksi = &[0.0, 0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let mut edge_shape = mesh.alloc_shape_boundary_edge(edge_key)?;
                edge_shape.calc_boundary_normal(&mut normal, ksi)?;
                assert_vec_approx_eq!(normal.as_data(), solution, 1e-15);
            }
        }
        Ok(())
    }

    #[test]
    fn normals_are_correct_3d() -> Result<(), StrError> {
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

        // the norm of the normal vector should be equal to 0.25 = face_area / 4.0
        // where 4.0 corresponds to the face_area in the reference system
        let l = 0.25; // norm of normal vector

        // face keys and correct normal vectors (solutions)
        let face_keys_and_solutions = [
            // behind
            (vec![(0, 3, 4, 7), (4, 7, 8, 11)], [-l, 0.0, 0.0]),
            // front
            (vec![(1, 2, 5, 6), (5, 6, 9, 10)], [l, 0.0, 0.0]),
            // left
            (vec![(0, 1, 4, 5), (4, 5, 8, 9)], [0.0, -l, 0.0]),
            // right
            (vec![(2, 3, 6, 7), (6, 7, 10, 11)], [0.0, l, 0.0]),
            // bottom
            (vec![(0, 1, 2, 3)], [0.0, 0.0, -l]),
            // top
            (vec![(8, 9, 10, 11)], [0.0, 0.0, l]),
        ];

        // check if the normal vectors at boundary are outward
        let mut normal = Vector::new(mesh.space_ndim);
        let ksi = &[0.0, 0.0, 0.0];
        for (face_keys, solution) in &face_keys_and_solutions {
            for face_key in face_keys {
                let mut face_shape = mesh.alloc_shape_boundary_face(face_key)?;
                face_shape.calc_boundary_normal(&mut normal, ksi)?;
                assert_vec_approx_eq!(normal.as_data(), solution, 1e-15);
            }
        }
        Ok(())
    }
}
