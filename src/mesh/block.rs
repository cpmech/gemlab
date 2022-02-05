use super::{Cell, Constraint, Mesh, Point};
use crate::shapes::Shape;
use crate::util::{AsArray2D, GridSearch};
use crate::StrError;
use russell_lab::Vector;
use std::collections::HashSet;

/// Defines a polygon on polyhedron that can be split into smaller shapes
///
/// # Block shape
///
/// ## 2D: Local indices of points and edges
///
/// ```text
///             Points (4 or 8)                        Edges (4)
///   y
///   |        3      6      2                             2
///   +--x      @-----@-----@                        +-----------+
///             |           |                        |           |
///             |           |                        |           |
///           7 @           @ 5                     3|           |1
///             |           |                        |           |
///             |           |                        |           |
///             @-----@-----@                        +-----------+
///            0      4      1                             0
/// ```
///
/// ## 3D: Local indices of points, edges, and faces
///
/// ```text
///                   Points (8 or 20)                   Faces (6)
///    z
///    |           4        15        7
///   ,+--y         @-------@--------@                 +----------------+
/// x'            ,'|              ,'|               ,'|              ,'|
///          12 @'  |         14 ,'  |             ,'  |  ___       ,'  |
///           ,'    |16        ,@    |19         ,'    |,'5,'  [0],'    |
///     5   ,'      @      6 ,'      @         ,'      |~~~     ,'      |
///       @'=======@=======@'        |       +'===============+'  ,'|   |
///       |      13 |      |         |       |   ,'|   |      |   |3|   |
///       |         |      |  11     |       |   |2|   |      |   |,'   |
///    17 |       0 @- - - | @- - - -@       |   |,'   +- - - | +- - - -+
///       @       ,'       @       ,' 3      |       ,'       |       ,'
///       |   8 @'      18 |     ,'          |     ,' [1]  ___|     ,'
///       |   ,'           |   ,@ 10         |   ,'      ,'4,'|   ,'
///       | ,'             | ,'              | ,'        ~~~  | ,'
///       @-------@--------@'                +----------------+'
///     1         9         2
///
///      
///                  Edges (12)
///      
///                         7
///                 +----------------+
///               ,'|              ,'|
///           4 ,'  |8          6,'  |
///           ,'    |          ,'    |   
///         ,'      |   5    ,'      |11
///       +'===============+'        |
///       |         |      |         |
///       |         |      |  3      |
///       |         .- - - | -  - - -+
///      9|       ,'       |       ,'
///       |    0,'         |10   ,'
///       |   ,'           |   ,' 2
///       | ,'             | ,'
///       +----------------+'
///               1
///  
pub struct Block {
    attribute_id: usize, // attribute ID of all elements in this block
    space_ndim: usize,   // space dimension
    // coords: Matrix,           // coordinates of points (npoint, ndim)
    ndiv: Vec<usize>,         // number of divisions along each dim (ndim)
    delta_ksi: Vec<Vec<f64>>, // delta ksi along each dim (ndim, {ndiv[0],ndiv[1],ndiv[2]})

    edge_constraints: Vec<Option<Constraint>>, // constraints (nedge)
    face_constraints: Vec<Option<Constraint>>, // constraints (nface)

    // shape and interpolation functions
    shape: Shape,

    // grid to search reference coordinates
    grid_ksi: GridSearch,
}

impl Block {
    // constants
    const NAT_LENGTH: f64 = 2.0; // length of shape along each direction in reference coords space
    const NAT_TOLERANCE: f64 = 1e-4; // tolerance to compare coordinates in the reference space

    // valid output npoint
    const VALID_OUTPUT_NPOINT_2D: [usize; 6] = [4, 8, 9, 12, 16, 17];
    const VALID_OUTPUT_NPOINT_3D: [usize; 2] = [8, 20];

    /// Creates a new Block with default options
    pub fn new(space_ndim: usize) -> Result<Self, StrError> {
        // check
        if space_ndim < 2 || space_ndim > 3 {
            return Err("space_ndim must be 2 or 3");
        }

        // shape
        let geo_ndim = space_ndim;
        let nnode = if geo_ndim == 2 { 8 } else { 20 };
        let shape = Shape::new(space_ndim, geo_ndim, nnode)?;

        // constants
        const NDIV: usize = 2;
        const GRID_NDIV: usize = 20;
        const GRID_MIN: f64 = -1.0;
        const GRID_MAX: f64 = 1.0;

        // grid
        let mut grid_ksi = GridSearch::new(space_ndim)?;
        grid_ksi.set_tolerances(&vec![Block::NAT_TOLERANCE; space_ndim])?;
        grid_ksi.initialize(
            &vec![GRID_NDIV; space_ndim],
            &vec![GRID_MIN; space_ndim],
            &vec![GRID_MAX; space_ndim],
        )?;

        // done
        Ok(Block {
            attribute_id: 1,
            space_ndim,
            ndiv: vec![NDIV; space_ndim],
            delta_ksi: vec![vec![1.0; NDIV]; space_ndim],
            edge_constraints: vec![None; shape.nedge],
            face_constraints: vec![None; shape.nface],
            shape,
            grid_ksi,
        })
    }

    /// Sets group
    pub fn set_attribute_id(&mut self, attribute_id: usize) -> &mut Self {
        self.attribute_id = attribute_id;
        self
    }

    /// Sets the coordinates of all 2D or 3D nodes
    ///
    /// # Input (2D)
    ///
    /// * `coords` -- (4 or 8, 2) matrix with all coordinates
    ///
    /// # Input (3D)
    ///
    /// * `coords` -- (8 or 20, 3) matrix with all coordinates
    #[rustfmt::skip]
    pub fn set_coords<'a, T, U>(&mut self, coords: &'a T) -> Result<(), StrError>
    where
        T: AsArray2D<'a, U>,
        U: 'a + Into<f64>,
    {
        // check
        let (nrow, ncol) = coords.size();
        if ncol != self.space_ndim {
            return Err("number of columns must be equal to space_ndim");
        }
        if self.space_ndim == 2 {
            if nrow != 4 && nrow != 8 {
                return Err("in 2D, the number of rows must be either 4 or 8");
            }
        } else {
            if nrow != 8 && nrow != 20 {
                return Err("in 3D, the number of rows must be either 8 or 20");
            }
        }

        // set corner nodes
        for m in 0..nrow {
            for j in 0..self.space_ndim {
                self.shape.set_node(m, j, coords.at(m, j).into())?
            }
        }

        // generate mid nodes
        if self.space_ndim == 2 && nrow == 4 {
            for j in 0..self.space_ndim {
                self.shape.set_node(4, j, (coords.at(0, j).into() + coords.at(1, j).into()) / 2.0)?;
                self.shape.set_node(5, j, (coords.at(1, j).into() + coords.at(2, j).into()) / 2.0)?;
                self.shape.set_node(6, j, (coords.at(2, j).into() + coords.at(3, j).into()) / 2.0)?;
                self.shape.set_node(7, j, (coords.at(3, j).into() + coords.at(0, j).into()) / 2.0)?;
            }
        }
        if self.space_ndim == 3 && nrow == 8 {
            for j in 0..self.space_ndim {
                self.shape.set_node( 8, j, (coords.at(0, j).into() + coords.at(1, j).into()) / 2.0)?;
                self.shape.set_node( 9, j, (coords.at(1, j).into() + coords.at(2, j).into()) / 2.0)?;
                self.shape.set_node(10, j, (coords.at(2, j).into() + coords.at(3, j).into()) / 2.0)?;
                self.shape.set_node(11, j, (coords.at(3, j).into() + coords.at(0, j).into()) / 2.0)?;

                self.shape.set_node(12, j, (coords.at(4, j).into() + coords.at(5, j).into()) / 2.0)?;
                self.shape.set_node(13, j, (coords.at(5, j).into() + coords.at(6, j).into()) / 2.0)?;
                self.shape.set_node(14, j, (coords.at(6, j).into() + coords.at(7, j).into()) / 2.0)?;
                self.shape.set_node(15, j, (coords.at(7, j).into() + coords.at(4, j).into()) / 2.0)?;

                self.shape.set_node(16, j, (coords.at(0, j).into() + coords.at(4, j).into()) / 2.0)?;
                self.shape.set_node(17, j, (coords.at(1, j).into() + coords.at(5, j).into()) / 2.0)?;
                self.shape.set_node(18, j, (coords.at(2, j).into() + coords.at(6, j).into()) / 2.0)?;
                self.shape.set_node(19, j, (coords.at(3, j).into() + coords.at(7, j).into()) / 2.0)?;
            }
        }
        Ok(())
    }

    /// Sets the number of equal divisions
    ///
    /// For each direction:
    ///
    /// ```text
    /// Δξᵐ = wᵐ ⋅ L / Σ_m wᵐ
    /// ```
    ///
    /// where `L=2` is the edge-length (in reference coordinates) and `wᵐ` are
    /// the weights for each division `m`.
    pub fn set_ndiv(&mut self, ndiv: &[usize]) -> &mut Self {
        assert_eq!(ndiv.len(), self.space_ndim);
        for i in 0..self.space_ndim {
            assert!(ndiv[i] > 0);
            self.ndiv[i] = ndiv[i];
            let w = 1.0;
            let sum_w = ndiv[i] as f64;
            self.delta_ksi[i] = vec![w * Block::NAT_LENGTH / sum_w; ndiv[i]];
        }
        self
    }

    /// Sets constraint to an edge of this block
    ///
    /// # Input
    ///
    /// * `e` -- index of edge in [0, nedge-1]
    /// * `constraint` -- the constraint
    pub fn set_edge_constraint(&mut self, e: usize, constraint: Constraint) -> &mut Self {
        self.edge_constraints[e] = Some(constraint);
        self
    }

    /// Sets constraint to a face of this block
    ///
    /// # Input
    ///
    /// * `f` -- index of face in [0, nface-1]
    /// * `constraint` -- the constraint
    pub fn set_face_constraint(&mut self, f: usize, constraint: Constraint) -> &mut Self {
        self.face_constraints[f] = Some(constraint);
        self
    }

    /// Subdivide block into vertices and cells (mesh)
    ///
    /// # Input
    ///
    /// Valid output_npoint:
    ///
    /// * 2D: [4, 8, 9, 12, 16, 17]
    /// * 3D: [8, 20]
    pub fn subdivide(&mut self, output_npoint: usize) -> Result<Mesh, StrError> {
        // check
        let space_ndim = self.space_ndim;
        if space_ndim == 2 {
            if !Block::VALID_OUTPUT_NPOINT_2D.contains(&output_npoint) {
                return Err("output_npoint is invalid");
            }
        } else {
            if !Block::VALID_OUTPUT_NPOINT_3D.contains(&output_npoint) {
                return Err("output_npoint is invalid");
            }
        }

        // resulting mesh
        let mut mesh = Mesh::new(space_ndim)?;

        // shape object corresponding to the new (output) cells
        let shape_out = Shape::new(space_ndim, space_ndim, output_npoint)?;

        //          BLOCK                      OUTPUT CELL
        // (-1,+1)-----------(+1,+1)    (-1,+1)-----------(+1,+1)
        //    |                 |   _.~'   |                 |
        //    |                 _.~'       |                 |
        //    |      o------o ~'|   scale  |                 |
        //    |      |      |   | <<  &    |                 |
        //    |      |      |   | translate|                 |
        //    |      o------o .___         |                 |
        //    |                 | '~~~._   |                 |
        // (-1,-1)-----------(+1,-1)    (-1,-1)-----------(+1,-1)

        // center of new cell in reference space (natural coordinates)
        let mut center = vec![0.0; space_ndim];

        // size ratios between new cell and block in natural coordinates
        let mut scale = vec![0.0; space_ndim];

        // natural coordinates of new points
        let mut ksi = vec![0.0; space_ndim];

        // real coordinates of new points
        let mut x = Vector::new(space_ndim);

        // number of divisions along each direction
        let (nx, ny, nz) = (
            self.ndiv[0],
            self.ndiv[1],
            if space_ndim == 2 { 1 } else { self.ndiv[2] },
        );

        // for each z-division
        if space_ndim == 3 {
            center[2] = -1.0 + self.delta_ksi[2][0] / 2.0;
        }
        for k in 0..nz {
            if space_ndim == 3 {
                scale[2] = self.delta_ksi[2][k] / Block::NAT_LENGTH;
            }

            // for each y-division
            center[1] = -1.0 + self.delta_ksi[1][0] / 2.0;
            for j in 0..ny {
                scale[1] = self.delta_ksi[1][j] / Block::NAT_LENGTH;

                // for each x-division
                center[0] = -1.0 + self.delta_ksi[0][0] / 2.0;
                for i in 0..nx {
                    scale[0] = self.delta_ksi[0][i] / Block::NAT_LENGTH;

                    // new cell id
                    let cell_id = mesh.cells.len();

                    // for each node of the new (output) cell
                    // (there may be more nodes than the block; e.g., internal nodes)
                    let mut points = vec![0; shape_out.nnode];
                    for m in 0..shape_out.nnode {
                        // reference natural coordinates of the new cell nodes
                        let ksi_ref = shape_out.get_reference_coords(m);

                        // scale and translate the reference coordinates
                        for a in 0..space_ndim {
                            ksi[a] = center[a] + scale[a] * ksi_ref[a];
                        }

                        // get existent point id or create new point
                        let point_id = match self.grid_ksi.find(&ksi)? {
                            Some(point_id) => point_id,
                            None => {
                                // insert point id in grid-search
                                let point_id = mesh.points.len();
                                self.grid_ksi.insert(point_id, &ksi)?;

                                // compute real coordinates of point
                                self.shape.calc_coords(&mut x, &ksi)?;

                                // add new point to mesh
                                let mut shared_by_cells = HashSet::new();
                                shared_by_cells.insert(cell_id);
                                mesh.points.push(Point {
                                    id: point_id,
                                    coords: x.as_data().clone(),
                                    shared_by_boundary_edges: HashSet::new(),
                                    shared_by_boundary_faces: HashSet::new(),
                                });
                                point_id
                            }
                        };
                        points[m] = point_id;
                    }

                    // new cell
                    let cell = Cell {
                        id: cell_id,
                        attribute_id: self.attribute_id,
                        geo_ndim: shape_out.geo_ndim,
                        points,
                        shape: Shape::new(space_ndim, space_ndim, output_npoint)?,
                    };
                    mesh.cells.push(cell);

                    // next x-center
                    center[0] += self.delta_ksi[0][i];
                }

                // next y-center
                center[1] += self.delta_ksi[1][j];
            }

            // next z-center
            if space_ndim == 3 {
                center[2] += self.delta_ksi[2][k];
            }
        }

        // done
        mesh.compute_derived_props()?;
        Ok(mesh)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Block, StrError};
    use crate::geometry::Circle;
    use crate::mesh::Constraint;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;

    #[test]
    fn new_fails_on_wrong_input() {
        let b2d = Block::new(1);
        assert_eq!(b2d.err(), Some("space_ndim must be 2 or 3"));
    }

    #[test]
    fn new_works() -> Result<(), StrError> {
        let b2d = Block::new(2)?;
        assert_eq!(b2d.attribute_id, 1);
        assert_eq!(b2d.space_ndim, 2);
        assert_eq!(
            format!("{}", b2d.shape.coords_transp),
            "┌                 ┐\n\
             │ 0 0 0 0 0 0 0 0 │\n\
             │ 0 0 0 0 0 0 0 0 │\n\
             └                 ┘"
        );
        assert_eq!(b2d.ndiv, &[2, 2]);
        assert_eq!(format!("{:?}", b2d.delta_ksi), "[[1.0, 1.0], [1.0, 1.0]]");
        assert_eq!(b2d.shape.nnode, 8);

        let b3d = Block::new(3)?;
        assert_eq!(b3d.attribute_id, 1);
        assert_eq!(b3d.space_ndim, 3);
        assert_eq!(
            format!("{}", b3d.shape.coords_transp),
            "┌                                         ┐\n\
             │ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 │\n\
             │ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 │\n\
             │ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 │\n\
             └                                         ┘"
        );
        assert_eq!(b3d.ndiv, &[2, 2, 2]);
        assert_eq!(format!("{:?}", b3d.delta_ksi), "[[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]");
        assert_eq!(b3d.shape.nnode, 20);
        Ok(())
    }

    #[test]
    fn set_attribute_id_works() -> Result<(), StrError> {
        let mut block = Block::new(2)?;
        block.set_attribute_id(2);
        assert_eq!(block.attribute_id, 2);
        Ok(())
    }

    #[test]
    fn set_coords_fail_on_wrong_input() -> Result<(), StrError> {
        let mut b2d = Block::new(2)?;
        let wrong: [[f64; 1]; 1] = [[0.0]];
        assert_eq!(
            b2d.set_coords(&wrong).err(),
            Some("number of columns must be equal to space_ndim")
        );
        let wrong: [[f64; 2]; 1] = [[0.0, 0.0]];
        assert_eq!(
            b2d.set_coords(&wrong).err(),
            Some("in 2D, the number of rows must be either 4 or 8")
        );
        let mut b3d = Block::new(3)?;
        let wrong: [[f64; 3]; 1] = [[0.0, 0.0, 0.0]];
        assert_eq!(
            b3d.set_coords(&wrong).err(),
            Some("in 3D, the number of rows must be either 8 or 20")
        );
        Ok(())
    }

    #[test]
    fn set_coords_works() -> Result<(), StrError> {
        let mut b2d = Block::new(2)?;
        #[rustfmt::skip]
        b2d.set_coords(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ])?;
        assert_eq!(
            format!("{}", b2d.shape.coords_transp),
            "┌                 ┐\n\
             │ 0 2 2 0 1 2 1 0 │\n\
             │ 0 0 2 2 0 1 2 1 │\n\
             └                 ┘"
        );

        let mut b3d = Block::new(3)?;
        #[rustfmt::skip]
        b3d.set_coords(&[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [2.0, 0.0, 2.0],
            [2.0, 2.0, 2.0],
            [0.0, 2.0, 2.0],
        ])?;
        assert_eq!(
            format!("{}", b3d.shape.coords_transp),
            "┌                                         ┐\n\
             │ 0 2 2 0 0 2 2 0 1 2 1 0 1 2 1 0 0 2 2 0 │\n\
             │ 0 0 2 2 0 0 2 2 0 1 2 1 0 1 2 1 0 0 2 2 │\n\
             │ 0 0 0 0 2 2 2 2 0 0 0 0 2 2 2 2 1 1 1 1 │\n\
             └                                         ┘"
        );
        Ok(())
    }

    #[test]
    fn set_ndiv_works() -> Result<(), StrError> {
        let mut block = Block::new(2)?;
        block.set_ndiv(&[2, 4]);
        assert_eq!(block.ndiv, &[2, 4]);
        assert_eq!(format!("{:?}", block.delta_ksi), "[[1.0, 1.0], [0.5, 0.5, 0.5, 0.5]]");
        Ok(())
    }

    #[test]
    fn set_edge_constraint_works() -> Result<(), StrError> {
        let mut block = Block::new(2)?;
        let constraint = Constraint::Arc(Circle {
            center: [-1.0, -1.0],
            radius: 2.0,
            tolerance: 1e-3,
        });
        block.set_edge_constraint(0, constraint);
        assert_eq!(block.edge_constraints[0], Some(constraint));
        Ok(())
    }

    #[test]
    fn set_face_constraint_works() -> Result<(), StrError> {
        let mut block = Block::new(3)?;
        let constraint = Constraint::ArcX(Circle {
            center: [-1.0, -1.0],
            radius: 2.0,
            tolerance: 1e-3,
        });
        block.set_face_constraint(0, constraint);
        assert_eq!(block.face_constraints[0], Some(constraint));
        Ok(())
    }

    #[test]
    fn subdivide_fails_on_wrong_input() -> Result<(), StrError> {
        let mut b2d = Block::new(2)?;
        assert_eq!(b2d.subdivide(1).err(), Some("output_npoint is invalid"));
        let mut b3d = Block::new(3)?;
        assert_eq!(b3d.subdivide(1).err(), Some("output_npoint is invalid"));
        Ok(())
    }

    #[test]
    fn subdivide_2d_works() -> Result<(), StrError> {
        let mut block = Block::new(2)?;
        #[rustfmt::skip]
        block.set_coords(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ])?;
        let mut mesh = block.subdivide(4)?;
        //
        //   7---------------6---------------8
        //   |               |               |
        //   |               |               |
        //   |      [2]      |      [3]      |
        //   |               |               |
        //   |               |               |
        //   3---------------2---------------5
        //   |               |               |
        //   |               |               |
        //   |      [0]      |      [1]      |
        //   |               |               |
        //   |               |               |
        //   0---------------1---------------4
        //
        assert_eq!(
            format!("{}", mesh),
            "SUMMARY\n\
             =======\n\
             space_ndim = 2\n\
             npoint = 9\n\
             ncell = 4\n\
             n_boundary_point = 8\n\
             n_boundary_edge = 8\n\
             n_boundary_face = 0\n\
             \n\
             POINTS\n\
             ======\n\
             i:0 x:[0.0, 0.0] e:[(0, 1), (0, 3)] f:[]\n\
             i:1 x:[1.0, 0.0] e:[(0, 1), (1, 4)] f:[]\n\
             i:2 x:[1.0, 1.0] e:[] f:[]\n\
             i:3 x:[0.0, 1.0] e:[(0, 3), (3, 7)] f:[]\n\
             i:4 x:[2.0, 0.0] e:[(1, 4), (4, 5)] f:[]\n\
             i:5 x:[2.0, 1.0] e:[(4, 5), (5, 8)] f:[]\n\
             i:6 x:[1.0, 2.0] e:[(6, 7), (6, 8)] f:[]\n\
             i:7 x:[0.0, 2.0] e:[(3, 7), (6, 7)] f:[]\n\
             i:8 x:[2.0, 2.0] e:[(5, 8), (6, 8)] f:[]\n\
             \n\
             CELLS\n\
             =====\n\
             i:0 a:1 g:2 p:[0, 1, 2, 3]\n\
             i:1 a:1 g:2 p:[1, 4, 5, 2]\n\
             i:2 a:1 g:2 p:[3, 2, 6, 7]\n\
             i:3 a:1 g:2 p:[2, 5, 8, 6]\n\
             \n\
             BOUNDARY POINTS\n\
             ===============\n\
             [0, 1, 3, 4, 5, 6, 7, 8]\n\
             \n\
             BOUNDARY EDGES\n\
             ==============\n\
             k:(0,1) p:[1, 0] c:[0] f:[]\n\
             k:(0,3) p:[0, 3] c:[0] f:[]\n\
             k:(1,4) p:[4, 1] c:[1] f:[]\n\
             k:(3,7) p:[3, 7] c:[2] f:[]\n\
             k:(4,5) p:[5, 4] c:[1] f:[]\n\
             k:(5,8) p:[8, 5] c:[3] f:[]\n\
             k:(6,7) p:[7, 6] c:[2] f:[]\n\
             k:(6,8) p:[6, 8] c:[3] f:[]\n\
             \n\
             BOUNDARY FACES\n\
             ==============\n"
        );

        // the magnitude of the normal vector should be equal to edge_length / 2.0 = 1.0 / 2.0
        // where 2.0 corresponds to the edge_length in the reference system
        let l = 0.5; // magnitude of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = vec![
            // bottom
            (vec![(0, 1), (1, 4)], [0.0, -l]),
            // right
            (vec![(4, 5), (5, 8)], [l, 0.0]),
            // top
            (vec![(6, 7), (6, 8)], [0.0, l]),
            // left
            (vec![(0, 3), (3, 7)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let mut normal = Vector::new(mesh.space_ndim);
        let ksi = &[0.0, 0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let edge = mesh.boundary_edges.get_mut(edge_key).unwrap();
                assert_eq!(edge.points.len(), 2);
                edge.shape.calc_boundary_normal(&mut normal, ksi)?;
                assert_vec_approx_eq!(normal.as_data(), solution, 1e-15);
            }
        }
        Ok(())
    }

    #[test]
    fn subdivide_2d_o2_works() -> Result<(), StrError> {
        let mut block = Block::new(2)?;
        #[rustfmt::skip]
        block.set_coords(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ])?;
        let mut mesh = block.subdivide(8)?;
        //
        //  14------16------13------20------18
        //   |               |               |
        //   |               |               |
        //  17      [2]     15      [3]     19
        //   |               |               |
        //   |               |               |
        //   3-------6-------2------12-------9
        //   |               |               |
        //   |               |               |
        //   7      [0]      5      [1]     11
        //   |               |               |
        //   |               |               |
        //   0-------4-------1------10-------8
        //
        assert_eq!(
            format!("{}", mesh),
            "SUMMARY\n\
             =======\n\
             space_ndim = 2\n\
             npoint = 21\n\
             ncell = 4\n\
             n_boundary_point = 16\n\
             n_boundary_edge = 8\n\
             n_boundary_face = 0\n\
             \n\
             POINTS\n\
             ======\n\
             i:0 x:[0.0, 0.0] e:[(0, 1), (0, 3)] f:[]\n\
             i:1 x:[1.0, 0.0] e:[(0, 1), (1, 8)] f:[]\n\
             i:2 x:[1.0, 1.0] e:[] f:[]\n\
             i:3 x:[0.0, 1.0] e:[(0, 3), (3, 14)] f:[]\n\
             i:4 x:[0.5, 0.0] e:[(0, 1)] f:[]\n\
             i:5 x:[1.0, 0.5] e:[] f:[]\n\
             i:6 x:[0.5, 1.0] e:[] f:[]\n\
             i:7 x:[0.0, 0.5] e:[(0, 3)] f:[]\n\
             i:8 x:[2.0, 0.0] e:[(1, 8), (8, 9)] f:[]\n\
             i:9 x:[2.0, 1.0] e:[(8, 9), (9, 18)] f:[]\n\
             i:10 x:[1.5, 0.0] e:[(1, 8)] f:[]\n\
             i:11 x:[2.0, 0.5] e:[(8, 9)] f:[]\n\
             i:12 x:[1.5, 1.0] e:[] f:[]\n\
             i:13 x:[1.0, 2.0] e:[(13, 14), (13, 18)] f:[]\n\
             i:14 x:[0.0, 2.0] e:[(3, 14), (13, 14)] f:[]\n\
             i:15 x:[1.0, 1.5] e:[] f:[]\n\
             i:16 x:[0.5, 2.0] e:[(13, 14)] f:[]\n\
             i:17 x:[0.0, 1.5] e:[(3, 14)] f:[]\n\
             i:18 x:[2.0, 2.0] e:[(9, 18), (13, 18)] f:[]\n\
             i:19 x:[2.0, 1.5] e:[(9, 18)] f:[]\n\
             i:20 x:[1.5, 2.0] e:[(13, 18)] f:[]\n\
             \n\
             CELLS\n\
             =====\n\
             i:0 a:1 g:2 p:[0, 1, 2, 3, 4, 5, 6, 7]\n\
             i:1 a:1 g:2 p:[1, 8, 9, 2, 10, 11, 12, 5]\n\
             i:2 a:1 g:2 p:[3, 2, 13, 14, 6, 15, 16, 17]\n\
             i:3 a:1 g:2 p:[2, 9, 18, 13, 12, 19, 20, 15]\n\
             \n\
             BOUNDARY POINTS\n\
             ===============\n\
             [0, 1, 3, 4, 7, 8, 9, 10, 11, 13, 14, 16, 17, 18, 19, 20]\n\
             \n\
             BOUNDARY EDGES\n\
             ==============\n\
             k:(0,1) p:[1, 0, 4] c:[0] f:[]\n\
             k:(0,3) p:[0, 3, 7] c:[0] f:[]\n\
             k:(1,8) p:[8, 1, 10] c:[1] f:[]\n\
             k:(3,14) p:[3, 14, 17] c:[2] f:[]\n\
             k:(8,9) p:[9, 8, 11] c:[1] f:[]\n\
             k:(9,18) p:[18, 9, 19] c:[3] f:[]\n\
             k:(13,14) p:[14, 13, 16] c:[2] f:[]\n\
             k:(13,18) p:[13, 18, 20] c:[3] f:[]\n\
             \n\
             BOUNDARY FACES\n\
             ==============\n"
        );

        // the magnitude of the normal vector should be equal to edge_length / 2.0 = 1.0 / 2.0
        // where 2.0 corresponds to the edge_length in the reference system
        let l = 0.5; // magnitude of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 8)], [0.0, -l]),
            // right
            (vec![(8, 9), (9, 18)], [l, 0.0]),
            // top
            (vec![(13, 14), (13, 18)], [0.0, l]),
            // left
            (vec![(0, 3), (3, 14)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let mut normal = Vector::new(mesh.space_ndim);
        let ksi = &[0.0, 0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let edge = mesh.boundary_edges.get_mut(edge_key).unwrap();
                assert_eq!(edge.points.len(), 3);
                edge.shape.calc_boundary_normal(&mut normal, ksi)?;
                assert_vec_approx_eq!(normal.as_data(), solution, 1e-15);
            }
        }
        Ok(())
    }

    #[test]
    fn subdivide_2d_qua9_works() -> Result<(), StrError> {
        //
        //  16------18------15------23------21
        //   |               |               |
        //   |               |               |
        //  19      20      17      24      22
        //   |               |               |
        //   |               |               |
        //   3-------6-------2------13------10
        //   |               |               |
        //   |               |               |
        //   7       8       5      14      12
        //   |               |               |
        //   |               |               |
        //   0-------4-------1------11-------9
        //
        let mut block = Block::new(2)?;
        #[rustfmt::skip]
        block.set_coords(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ])?;
        let mut mesh = block.subdivide(9)?;
        // check
        assert_eq!(mesh.points.len(), 25);
        assert_eq!(mesh.points[24].coords, &[1.5, 1.5]);
        assert_eq!(mesh.cells[3].points, &[2, 10, 21, 15, 13, 22, 23, 17, 24]);

        // the magnitude of the normal vector should be equal to edge_length / 2.0 = 1.0 / 2.0
        // where 2.0 corresponds to the edge_length in the reference system
        let l = 0.5; // magnitude of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 9)], [0.0, -l]),
            // right
            (vec![(9, 10), (10, 21)], [l, 0.0]),
            // top
            (vec![(15, 16), (15, 21)], [0.0, l]),
            // left
            (vec![(0, 3), (3, 16)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let mut normal = Vector::new(mesh.space_ndim);
        let ksi = &[0.0, 0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let edge = mesh.boundary_edges.get_mut(edge_key).unwrap();
                assert_eq!(edge.points.len(), 3);
                edge.shape.calc_boundary_normal(&mut normal, ksi)?;
                assert_vec_approx_eq!(normal.as_data(), solution, 1e-14);
            }
        }
        Ok(())
    }

    #[test]
    fn subdivide_2d_qua12_works() -> Result<(), StrError> {
        //
        //  21---26---23----20---32---30----28
        //   |               |               |
        //  24              25              31
        //   |               |               |
        //  27              22              29
        //   |               |               |
        //   3---10-----6----2---19---16----13
        //   |               |               |
        //   7               9              18
        //   |               |               |
        //  11               5              15
        //   |               |               |
        //   0----4-----8----1---14---17----12
        //
        let mut block = Block::new(2)?;
        #[rustfmt::skip]
        block.set_coords(&[
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [0.0, 3.0],
        ])?;
        let mut mesh = block.subdivide(12)?;

        // check
        assert_eq!(mesh.points.len(), 33);
        assert_vec_approx_eq!(mesh.points[9].coords, &[1.5, 1.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[15].coords, &[3.0, 0.5], 1e-15);
        assert_vec_approx_eq!(mesh.points[24].coords, &[0.0, 2.5], 1e-15);
        assert_vec_approx_eq!(mesh.points[32].coords, &[2.0, 3.0], 1e-15);
        assert_eq!(mesh.cells[0].points, &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
        assert_eq!(mesh.cells[1].points, &[1, 12, 13, 2, 14, 15, 16, 9, 17, 18, 19, 5]);
        assert_eq!(mesh.cells[2].points, &[3, 2, 20, 21, 10, 22, 23, 24, 6, 25, 26, 27]);
        assert_eq!(mesh.cells[3].points, &[2, 13, 28, 20, 19, 29, 30, 25, 16, 31, 32, 22]);

        // the magnitude of the normal vector should be equal to edge_length / 2.0 = 1.5 / 2.0
        // where 2.0 corresponds to the edge_length in the reference system
        let l = 0.75; // magnitude of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 12)], [0.0, -l]),
            // right
            (vec![(12, 13), (13, 28)], [l, 0.0]),
            // top
            (vec![(20, 21), (20, 28)], [0.0, l]),
            // left
            (vec![(0, 3), (3, 21)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let mut normal = Vector::new(mesh.space_ndim);
        let ksi = &[0.0, 0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let edge = mesh.boundary_edges.get_mut(edge_key).unwrap();
                assert_eq!(edge.points.len(), 4);
                edge.shape.calc_boundary_normal(&mut normal, ksi)?;
                assert_vec_approx_eq!(normal.as_data(), solution, 1e-14);
            }
        }
        Ok(())
    }

    #[test]
    fn subdivide_2d_qua16_works() -> Result<(), StrError> {
        //
        //  29---34----31---28---44---42----40
        //   |               |               |
        //  32   39    38   33   48   47    43
        //   |               |               |
        //  35   36    37   30   45   46    41
        //   |               |               |
        //   3---10-----6----2---23---20----17
        //   |               |               |
        //   7   15    14    9   27   26    22
        //   |               |               |
        //  11   12    13    5   24   25    19
        //   |               |               |
        //   0----4-----8----1---18---21----16
        //
        let mut block = Block::new(2)?;
        #[rustfmt::skip]
        block.set_coords(&[
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [0.0, 3.0],
        ])?;
        let mut mesh = block.subdivide(16)?;

        // check
        assert_eq!(mesh.points.len(), 49);
        assert_vec_approx_eq!(mesh.points[4].coords, &[0.5, 0.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[8].coords, &[1.0, 0.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[10].coords, &[0.5, 1.5], 1e-15);
        assert_vec_approx_eq!(mesh.points[13].coords, &[1.0, 0.5], 1e-15);
        assert_vec_approx_eq!(mesh.points[17].coords, &[3.0, 1.5], 1e-15);
        assert_vec_approx_eq!(mesh.points[24].coords, &[2.0, 0.5], 1e-15);
        assert_vec_approx_eq!(mesh.points[26].coords, &[2.5, 1.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[37].coords, &[1.0, 2.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[40].coords, &[3.0, 3.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[48].coords, &[2.0, 2.5], 1e-15);
        assert_vec_approx_eq!(mesh.points[43].coords, &[3.0, 2.5], 1e-15);
        assert_eq!(
            mesh.cells[0].points,
            &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        );
        assert_eq!(
            mesh.cells[1].points,
            &[1, 16, 17, 2, 18, 19, 20, 9, 21, 22, 23, 5, 24, 25, 26, 27]
        );
        assert_eq!(
            mesh.cells[2].points,
            &[3, 2, 28, 29, 10, 30, 31, 32, 6, 33, 34, 35, 36, 37, 38, 39]
        );
        assert_eq!(
            mesh.cells[3].points,
            &[2, 17, 40, 28, 23, 41, 42, 33, 20, 43, 44, 30, 45, 46, 47, 48]
        );

        // the magnitude of the normal vector should be equal to edge_length / 2.0 = 1.5 / 2.0
        // where 2.0 corresponds to the edge_length in the reference system
        let l = 0.75; // magnitude of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 16)], [0.0, -l]),
            // right
            (vec![(16, 17), (17, 40)], [l, 0.0]),
            // top
            (vec![(28, 29), (28, 40)], [0.0, l]),
            // left
            (vec![(0, 3), (3, 29)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let mut normal = Vector::new(mesh.space_ndim);
        let ksi = &[0.0, 0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let edge = mesh.boundary_edges.get_mut(edge_key).unwrap();
                assert_eq!(edge.points.len(), 4);
                edge.shape.calc_boundary_normal(&mut normal, ksi)?;
                assert_vec_approx_eq!(normal.as_data(), solution, 1e-14);
            }
        }
        Ok(())
    }

    #[test]
    fn subdivide_2d_qua17_works() -> Result<(), StrError> {
        //
        //  30---38---35---32---29---47---45---43---41
        //   |                   |                   |
        //  33                  37                  46
        //   |                   |                   |
        //  36        40        34        48        44
        //   |                   |                   |
        //  39                  31                  42
        //   |                   |                   |
        //   3---14---10----6----2---27---24---21---18
        //   |                   |                   |
        //   7                  13                  26
        //   |                   |                   |
        //  11        16         9        28        23
        //   |                   |                   |
        //  15                   5                  20
        //   |                   |                   |
        //   0----4----8---12----1---19---22---25---17
        //
        let mut block = Block::new(2)?;
        #[rustfmt::skip]
        block.set_coords(&[
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
        ])?;
        let mut mesh = block.subdivide(17)?;

        // check
        assert_eq!(mesh.points.len(), 49);
        assert_vec_approx_eq!(mesh.points[4].coords, &[0.5, 0.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[5].coords, &[2.0, 0.5], 1e-15);
        assert_vec_approx_eq!(mesh.points[6].coords, &[1.5, 2.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[11].coords, &[0.0, 1.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[12].coords, &[1.5, 0.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[16].coords, &[1.0, 1.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[19].coords, &[2.5, 0.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[23].coords, &[4.0, 1.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[27].coords, &[2.5, 2.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[28].coords, &[3.0, 1.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[32].coords, &[1.5, 4.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[40].coords, &[1.0, 3.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[41].coords, &[4.0, 4.0], 1e-15);
        assert_vec_approx_eq!(mesh.points[46].coords, &[4.0, 3.5], 1e-15);
        assert_vec_approx_eq!(mesh.points[48].coords, &[3.0, 3.0], 1e-15);
        assert_eq!(
            mesh.cells[0].points,
            &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        );
        assert_eq!(
            mesh.cells[1].points,
            &[1, 17, 18, 2, 19, 20, 21, 13, 22, 23, 24, 9, 25, 26, 27, 5, 28]
        );
        assert_eq!(
            mesh.cells[2].points,
            &[3, 2, 29, 30, 14, 31, 32, 33, 10, 34, 35, 36, 6, 37, 38, 39, 40]
        );
        assert_eq!(
            mesh.cells[3].points,
            &[2, 18, 41, 29, 27, 42, 43, 37, 24, 44, 45, 34, 21, 46, 47, 31, 48]
        );

        // the magnitude of the normal vector should be equal to edge_length / 2.0 = 2.0 / 2.0
        // where 2.0 corresponds to the edge_length in the reference system
        let l = 1.0; // magnitude of normal vector

        // edge keys and correct normal vectors (solutions)
        let edge_keys_and_solutions = [
            // bottom
            (vec![(0, 1), (1, 17)], [0.0, -l]),
            // right
            (vec![(17, 18), (18, 41)], [l, 0.0]),
            // top
            (vec![(29, 30), (29, 41)], [0.0, l]),
            // left
            (vec![(0, 3), (3, 30)], [-l, 0.0]),
        ];

        // check if the normal vectors at boundary are outward
        let mut normal = Vector::new(mesh.space_ndim);
        let ksi = &[0.0, 0.0, 0.0];
        for (edge_keys, solution) in &edge_keys_and_solutions {
            for edge_key in edge_keys {
                let edge = mesh.boundary_edges.get_mut(edge_key).unwrap();
                assert_eq!(edge.points.len(), 5);
                edge.shape.calc_boundary_normal(&mut normal, ksi)?;
                assert_vec_approx_eq!(normal.as_data(), solution, 1e-14);
            }
        }
        Ok(())
    }

    #[test]
    fn subdivide_3d_works() -> Result<(), StrError> {
        let mut block = Block::new(3)?;
        #[rustfmt::skip]
        block.set_coords(&[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 4.0],
            [2.0, 0.0, 4.0],
            [2.0, 2.0, 4.0],
            [0.0, 2.0, 4.0],
        ])?;
        let mesh = block.subdivide(8)?;
        //
        //              18------------------21------------------25
        //              /.                  /.                  /|
        //             / .                 / .                 / |
        //            /  .                /  .                /  |
        //           /   .               /   .               /   |
        //          /    .              /    .              /    |
        //        19------------------20------------------24     |
        //        /.     .            /.     .            /|     |
        //       / .     .           / .     .           / |     |
        //      /  .     .          /  .     .          /  |     |
        //     /   .     .         /   .     .         /   |     |
        //    /    .     .        /    .     .        /    |     |
        //  22==================23==================26     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |     .     4 - - - | - - . - - 7 - - - | - - | - -15
        //   |     .    /.       |     .    /.       |     |    /|
        //   |     .   / .       |     .   / .       |     |   / |
        //   |     .  /  .       |     .  /  .       |     |  /  |
        //   |     . /   .       |     . /   .       |     | /   |
        //   |     ./    .       |     ./    .       |     |/    |
        //   |     5 - - - - - - | - - 6 - - - - - - | - -14     |
        //   |    /.     .       |    /.     .       |    /|     |
        //   |   / .     .       |   / .     .       |   / |     |
        //   |  /  .     .       |  /  .     .       |  /  |     |
        //   | /   .     .       | /   .     .       | /   |     |
        //   |/    .     .       |/    .     .       |/    |     |
        //  10==================11==================17     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |     .     0 - - - | - - . - - 3 - - - | - - | - -13
        //   |     .    /        |     .    /        |     |    /
        //   |     .   /         |     .   /         |     |   /
        //   |     .  /          |     .  /          |     |  /
        //   |     . /           |     . /           |     | /
        //   |     ./            |     ./            |     |/
        //   |     1 - - - - - - | - - 2 - - - - - - | - -12
        //   |    /              |    /              |    /
        //   |   /               |   /               |   /
        //   |  /                |  /                |  /
        //   | /                 | /                 | /
        //   |/                  |/                  |/
        //   8===================9==================16
        //
        assert_eq!(
            format!("{}", mesh),
            "SUMMARY\n\
             =======\n\
             space_ndim = 3\n\
             npoint = 27\n\
             ncell = 8\n\
             n_boundary_point = 26\n\
             n_boundary_edge = 48\n\
             n_boundary_face = 24\n\
             \n\
             POINTS\n\
             ======\n\
             i:0 x:[0.0, 0.0, 0.0] e:[(0, 1), (0, 3), (0, 4)] f:[(0, 1, 2, 3), (0, 1, 4, 5), (0, 3, 4, 7)]\n\
             i:1 x:[1.0, 0.0, 0.0] e:[(0, 1), (1, 2), (1, 5), (1, 8)] f:[(0, 1, 2, 3), (0, 1, 4, 5), (1, 2, 8, 9), (1, 5, 8, 10)]\n\
             i:2 x:[1.0, 1.0, 0.0] e:[(1, 2), (2, 3), (2, 9), (2, 12)] f:[(0, 1, 2, 3), (1, 2, 8, 9), (2, 3, 12, 13), (2, 9, 12, 16)]\n\
             i:3 x:[0.0, 1.0, 0.0] e:[(0, 3), (2, 3), (3, 7), (3, 13)] f:[(0, 1, 2, 3), (0, 3, 4, 7), (2, 3, 12, 13), (3, 7, 13, 15)]\n\
             i:4 x:[0.0, 0.0, 2.0] e:[(0, 4), (4, 5), (4, 7), (4, 18)] f:[(0, 1, 4, 5), (0, 3, 4, 7), (4, 5, 18, 19), (4, 7, 18, 21)]\n\
             i:5 x:[1.0, 0.0, 2.0] e:[(1, 5), (4, 5), (5, 10), (5, 19)] f:[(0, 1, 4, 5), (1, 5, 8, 10), (4, 5, 18, 19), (5, 10, 19, 22)]\n\
             i:6 x:[1.0, 1.0, 2.0] e:[] f:[]\n\
             i:7 x:[0.0, 1.0, 2.0] e:[(3, 7), (4, 7), (7, 15), (7, 21)] f:[(0, 3, 4, 7), (3, 7, 13, 15), (4, 7, 18, 21), (7, 15, 21, 25)]\n\
             i:8 x:[2.0, 0.0, 0.0] e:[(1, 8), (8, 9), (8, 10)] f:[(1, 2, 8, 9), (1, 5, 8, 10), (8, 9, 10, 11)]\n\
             i:9 x:[2.0, 1.0, 0.0] e:[(2, 9), (8, 9), (9, 11), (9, 16)] f:[(1, 2, 8, 9), (2, 9, 12, 16), (8, 9, 10, 11), (9, 11, 16, 17)]\n\
             i:10 x:[2.0, 0.0, 2.0] e:[(5, 10), (8, 10), (10, 11), (10, 22)] f:[(1, 5, 8, 10), (5, 10, 19, 22), (8, 9, 10, 11), (10, 11, 22, 23)]\n\
             i:11 x:[2.0, 1.0, 2.0] e:[(9, 11), (10, 11), (11, 17), (11, 23)] f:[(8, 9, 10, 11), (9, 11, 16, 17), (10, 11, 22, 23), (11, 17, 23, 26)]\n\
             i:12 x:[1.0, 2.0, 0.0] e:[(2, 12), (12, 13), (12, 14), (12, 16)] f:[(2, 3, 12, 13), (2, 9, 12, 16), (12, 13, 14, 15), (12, 14, 16, 17)]\n\
             i:13 x:[0.0, 2.0, 0.0] e:[(3, 13), (12, 13), (13, 15)] f:[(2, 3, 12, 13), (3, 7, 13, 15), (12, 13, 14, 15)]\n\
             i:14 x:[1.0, 2.0, 2.0] e:[(12, 14), (14, 15), (14, 17), (14, 24)] f:[(12, 13, 14, 15), (12, 14, 16, 17), (14, 15, 24, 25), (14, 17, 24, 26)]\n\
             i:15 x:[0.0, 2.0, 2.0] e:[(7, 15), (13, 15), (14, 15), (15, 25)] f:[(3, 7, 13, 15), (7, 15, 21, 25), (12, 13, 14, 15), (14, 15, 24, 25)]\n\
             i:16 x:[2.0, 2.0, 0.0] e:[(9, 16), (12, 16), (16, 17)] f:[(2, 9, 12, 16), (9, 11, 16, 17), (12, 14, 16, 17)]\n\
             i:17 x:[2.0, 2.0, 2.0] e:[(11, 17), (14, 17), (16, 17), (17, 26)] f:[(9, 11, 16, 17), (11, 17, 23, 26), (12, 14, 16, 17), (14, 17, 24, 26)]\n\
             i:18 x:[0.0, 0.0, 4.0] e:[(4, 18), (18, 19), (18, 21)] f:[(4, 5, 18, 19), (4, 7, 18, 21), (18, 19, 20, 21)]\n\
             i:19 x:[1.0, 0.0, 4.0] e:[(5, 19), (18, 19), (19, 20), (19, 22)] f:[(4, 5, 18, 19), (5, 10, 19, 22), (18, 19, 20, 21), (19, 20, 22, 23)]\n\
             i:20 x:[1.0, 1.0, 4.0] e:[(19, 20), (20, 21), (20, 23), (20, 24)] f:[(18, 19, 20, 21), (19, 20, 22, 23), (20, 21, 24, 25), (20, 23, 24, 26)]\n\
             i:21 x:[0.0, 1.0, 4.0] e:[(7, 21), (18, 21), (20, 21), (21, 25)] f:[(4, 7, 18, 21), (7, 15, 21, 25), (18, 19, 20, 21), (20, 21, 24, 25)]\n\
             i:22 x:[2.0, 0.0, 4.0] e:[(10, 22), (19, 22), (22, 23)] f:[(5, 10, 19, 22), (10, 11, 22, 23), (19, 20, 22, 23)]\n\
             i:23 x:[2.0, 1.0, 4.0] e:[(11, 23), (20, 23), (22, 23), (23, 26)] f:[(10, 11, 22, 23), (11, 17, 23, 26), (19, 20, 22, 23), (20, 23, 24, 26)]\n\
             i:24 x:[1.0, 2.0, 4.0] e:[(14, 24), (20, 24), (24, 25), (24, 26)] f:[(14, 15, 24, 25), (14, 17, 24, 26), (20, 21, 24, 25), (20, 23, 24, 26)]\n\
             i:25 x:[0.0, 2.0, 4.0] e:[(15, 25), (21, 25), (24, 25)] f:[(7, 15, 21, 25), (14, 15, 24, 25), (20, 21, 24, 25)]\n\
             i:26 x:[2.0, 2.0, 4.0] e:[(17, 26), (23, 26), (24, 26)] f:[(11, 17, 23, 26), (14, 17, 24, 26), (20, 23, 24, 26)]\n\
             \n\
             CELLS\n\
             =====\n\
             i:0 a:1 g:3 p:[0, 1, 2, 3, 4, 5, 6, 7]\n\
             i:1 a:1 g:3 p:[1, 8, 9, 2, 5, 10, 11, 6]\n\
             i:2 a:1 g:3 p:[3, 2, 12, 13, 7, 6, 14, 15]\n\
             i:3 a:1 g:3 p:[2, 9, 16, 12, 6, 11, 17, 14]\n\
             i:4 a:1 g:3 p:[4, 5, 6, 7, 18, 19, 20, 21]\n\
             i:5 a:1 g:3 p:[5, 10, 11, 6, 19, 22, 23, 20]\n\
             i:6 a:1 g:3 p:[7, 6, 14, 15, 21, 20, 24, 25]\n\
             i:7 a:1 g:3 p:[6, 11, 17, 14, 20, 23, 26, 24]\n\
             \n\
             BOUNDARY POINTS\n\
             ===============\n\
             [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]\n\
             \n\
             BOUNDARY EDGES\n\
             ==============\n\
             k:(0,1) p:[0, 1] c:[] f:[(0, 1, 2, 3), (0, 1, 4, 5)]\n\
             k:(0,3) p:[3, 0] c:[] f:[(0, 1, 2, 3), (0, 3, 4, 7)]\n\
             k:(0,4) p:[0, 4] c:[] f:[(0, 1, 4, 5), (0, 3, 4, 7)]\n\
             k:(1,2) p:[1, 2] c:[] f:[(0, 1, 2, 3), (1, 2, 8, 9)]\n\
             k:(1,5) p:[5, 1] c:[] f:[(0, 1, 4, 5), (1, 5, 8, 10)]\n\
             k:(1,8) p:[1, 8] c:[] f:[(1, 2, 8, 9), (1, 5, 8, 10)]\n\
             k:(2,3) p:[2, 3] c:[] f:[(0, 1, 2, 3), (2, 3, 12, 13)]\n\
             k:(2,9) p:[9, 2] c:[] f:[(1, 2, 8, 9), (2, 9, 12, 16)]\n\
             k:(2,12) p:[2, 12] c:[] f:[(2, 3, 12, 13), (2, 9, 12, 16)]\n\
             k:(3,7) p:[3, 7] c:[] f:[(0, 3, 4, 7), (3, 7, 13, 15)]\n\
             k:(3,13) p:[13, 3] c:[] f:[(2, 3, 12, 13), (3, 7, 13, 15)]\n\
             k:(4,5) p:[4, 5] c:[] f:[(0, 1, 4, 5), (4, 5, 18, 19)]\n\
             k:(4,7) p:[7, 4] c:[] f:[(0, 3, 4, 7), (4, 7, 18, 21)]\n\
             k:(4,18) p:[4, 18] c:[] f:[(4, 5, 18, 19), (4, 7, 18, 21)]\n\
             k:(5,10) p:[5, 10] c:[] f:[(1, 5, 8, 10), (5, 10, 19, 22)]\n\
             k:(5,19) p:[19, 5] c:[] f:[(4, 5, 18, 19), (5, 10, 19, 22)]\n\
             k:(7,15) p:[15, 7] c:[] f:[(3, 7, 13, 15), (7, 15, 21, 25)]\n\
             k:(7,21) p:[7, 21] c:[] f:[(4, 7, 18, 21), (7, 15, 21, 25)]\n\
             k:(8,9) p:[8, 9] c:[] f:[(1, 2, 8, 9), (8, 9, 10, 11)]\n\
             k:(8,10) p:[10, 8] c:[] f:[(1, 5, 8, 10), (8, 9, 10, 11)]\n\
             k:(9,11) p:[11, 9] c:[] f:[(8, 9, 10, 11), (9, 11, 16, 17)]\n\
             k:(9,16) p:[9, 16] c:[] f:[(2, 9, 12, 16), (9, 11, 16, 17)]\n\
             k:(10,11) p:[10, 11] c:[] f:[(8, 9, 10, 11), (10, 11, 22, 23)]\n\
             k:(10,22) p:[22, 10] c:[] f:[(5, 10, 19, 22), (10, 11, 22, 23)]\n\
             k:(11,17) p:[11, 17] c:[] f:[(9, 11, 16, 17), (11, 17, 23, 26)]\n\
             k:(11,23) p:[23, 11] c:[] f:[(10, 11, 22, 23), (11, 17, 23, 26)]\n\
             k:(12,13) p:[12, 13] c:[] f:[(2, 3, 12, 13), (12, 13, 14, 15)]\n\
             k:(12,14) p:[12, 14] c:[] f:[(12, 13, 14, 15), (12, 14, 16, 17)]\n\
             k:(12,16) p:[16, 12] c:[] f:[(2, 9, 12, 16), (12, 14, 16, 17)]\n\
             k:(13,15) p:[13, 15] c:[] f:[(3, 7, 13, 15), (12, 13, 14, 15)]\n\
             k:(14,15) p:[14, 15] c:[] f:[(12, 13, 14, 15), (14, 15, 24, 25)]\n\
             k:(14,17) p:[17, 14] c:[] f:[(12, 14, 16, 17), (14, 17, 24, 26)]\n\
             k:(14,24) p:[14, 24] c:[] f:[(14, 15, 24, 25), (14, 17, 24, 26)]\n\
             k:(15,25) p:[15, 25] c:[] f:[(7, 15, 21, 25), (14, 15, 24, 25)]\n\
             k:(16,17) p:[17, 16] c:[] f:[(9, 11, 16, 17), (12, 14, 16, 17)]\n\
             k:(17,26) p:[26, 17] c:[] f:[(11, 17, 23, 26), (14, 17, 24, 26)]\n\
             k:(18,19) p:[18, 19] c:[] f:[(4, 5, 18, 19), (18, 19, 20, 21)]\n\
             k:(18,21) p:[21, 18] c:[] f:[(4, 7, 18, 21), (18, 19, 20, 21)]\n\
             k:(19,20) p:[20, 19] c:[] f:[(18, 19, 20, 21), (19, 20, 22, 23)]\n\
             k:(19,22) p:[19, 22] c:[] f:[(5, 10, 19, 22), (19, 20, 22, 23)]\n\
             k:(20,21) p:[21, 20] c:[] f:[(18, 19, 20, 21), (20, 21, 24, 25)]\n\
             k:(20,23) p:[20, 23] c:[] f:[(19, 20, 22, 23), (20, 23, 24, 26)]\n\
             k:(20,24) p:[24, 20] c:[] f:[(20, 21, 24, 25), (20, 23, 24, 26)]\n\
             k:(21,25) p:[25, 21] c:[] f:[(7, 15, 21, 25), (20, 21, 24, 25)]\n\
             k:(22,23) p:[22, 23] c:[] f:[(10, 11, 22, 23), (19, 20, 22, 23)]\n\
             k:(23,26) p:[23, 26] c:[] f:[(11, 17, 23, 26), (20, 23, 24, 26)]\n\
             k:(24,25) p:[24, 25] c:[] f:[(14, 15, 24, 25), (20, 21, 24, 25)]\n\
             k:(24,26) p:[26, 24] c:[] f:[(14, 17, 24, 26), (20, 23, 24, 26)]\n\
             \n\
             BOUNDARY FACES\n\
             ==============\n\
             k:(0,1,2,3) p:[0, 3, 2, 1] c:[0]\n\
             k:(0,1,4,5) p:[0, 1, 5, 4] c:[0]\n\
             k:(0,3,4,7) p:[0, 4, 7, 3] c:[0]\n\
             k:(1,2,8,9) p:[1, 2, 9, 8] c:[1]\n\
             k:(1,5,8,10) p:[1, 8, 10, 5] c:[1]\n\
             k:(2,3,12,13) p:[3, 13, 12, 2] c:[2]\n\
             k:(2,9,12,16) p:[2, 12, 16, 9] c:[3]\n\
             k:(3,7,13,15) p:[3, 7, 15, 13] c:[2]\n\
             k:(4,5,18,19) p:[4, 5, 19, 18] c:[4]\n\
             k:(4,7,18,21) p:[4, 18, 21, 7] c:[4]\n\
             k:(5,10,19,22) p:[5, 10, 22, 19] c:[5]\n\
             k:(7,15,21,25) p:[7, 21, 25, 15] c:[6]\n\
             k:(8,9,10,11) p:[8, 9, 11, 10] c:[1]\n\
             k:(9,11,16,17) p:[9, 16, 17, 11] c:[3]\n\
             k:(10,11,22,23) p:[10, 11, 23, 22] c:[5]\n\
             k:(11,17,23,26) p:[11, 17, 26, 23] c:[7]\n\
             k:(12,13,14,15) p:[12, 13, 15, 14] c:[2]\n\
             k:(12,14,16,17) p:[16, 12, 14, 17] c:[3]\n\
             k:(14,15,24,25) p:[14, 15, 25, 24] c:[6]\n\
             k:(14,17,24,26) p:[17, 14, 24, 26] c:[7]\n\
             k:(18,19,20,21) p:[18, 19, 20, 21] c:[4]\n\
             k:(19,20,22,23) p:[19, 22, 23, 20] c:[5]\n\
             k:(20,21,24,25) p:[21, 20, 24, 25] c:[6]\n\
             k:(20,23,24,26) p:[20, 23, 26, 24] c:[7]\n"
        );
        Ok(())
    }

    #[test]
    fn subdivide_3d_o2_works() -> Result<(), StrError> {
        let mut block = Block::new(3)?;
        #[rustfmt::skip]
        block.set_coords(&[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 4.0],
            [2.0, 0.0, 4.0],
            [2.0, 2.0, 4.0],
            [0.0, 2.0, 4.0],
        ])?;
        let mesh = block.subdivide(20)?;
        //
        //              51--------58--------54--------74--------71
        //              /.                  /.                  /|
        //             / .                 / .                 / |
        //           55  .               57  .               73  |
        //           /   .               /   .               /   |
        //          /    .              /    .              /    |
        //        52--------56--------53--------72--------70     |
        //        /.     .            /.     .            /|     |
        //       / .    59           / .    62           / |    76
        //     65  .     .         67  .     .         79  |     |
        //     /   .     .         /   .     .         /   |     |
        //    /    .     .        /    .     .        /    |     |
        //  63========66========64========78========77     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |    60     .       |    61     .       |    75     |
        //   |     .     4 - - - |15 - . - - 7 - - - |41 - | - -35
        //   |     .    /.       |     .    /.       |     |    /|
        //   |     .   / .       |     .   / .       |     |   / |
        //   |     . 12  .       |     . 14  .       |     | 40  |
        //   |     . /   .       |     . /   .       |     | /   |
        //  68     ./    .      69     ./    .      80     |/    |
        //   |     5 - - - -13 - | - - 6 - - - -39 - | - -34     |
        //   |    /.     .       |    /.     .       |    /|     |
        //   |   / .    16       |   / .    19       |   / |    43
        //   | 27  .     .       | 29  .     .       | 49  |     |
        //   | /   .     .       | /   .     .       | /   |     |
        //   |/    .     .       |/    .     .       |/    |     |
        //  22========28========23========48========45     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |    17     .       |    18     .       |    42     |
        //   |     .     0 - - - |11 - . - - 3 - - - |38 - | - -33
        //   |     .    /        |     .    /        |     |    /
        //   |     .   /         |     .   /         |     |   /
        //   |     .  8          |     . 10          |     | 37
        //   |     . /           |     . /           |     | /
        //  30     ./           31     ./           50     |/
        //   |     1 - - - - 9 - | - - 2 - - - -36 - | - -32
        //   |    /              |    /              |    /
        //   |   /               |   /               |   /
        //   | 24                | 26                | 47
        //   | /                 | /                 | /
        //   |/                  |/                  |/
        //  20========25========21========46========44
        //
        assert_eq!(
            format!("{}", mesh),
            "SUMMARY\n\
             =======\n\
             space_ndim = 3\n\
             npoint = 81\n\
             ncell = 8\n\
             n_boundary_point = 74\n\
             n_boundary_edge = 48\n\
             n_boundary_face = 24\n\
             \n\
             POINTS\n\
             ======\n\
             i:0 x:[0.0, 0.0, 0.0] e:[(0, 1), (0, 3), (0, 4)] f:[(0, 1, 2, 3), (0, 1, 4, 5), (0, 3, 4, 7)]\n\
             i:1 x:[1.0, 0.0, 0.0] e:[(0, 1), (1, 2), (1, 5), (1, 20)] f:[(0, 1, 2, 3), (0, 1, 4, 5), (1, 2, 20, 21), (1, 5, 20, 22)]\n\
             i:2 x:[1.0, 1.0, 0.0] e:[(1, 2), (2, 3), (2, 21), (2, 32)] f:[(0, 1, 2, 3), (1, 2, 20, 21), (2, 3, 32, 33), (2, 21, 32, 44)]\n\
             i:3 x:[0.0, 1.0, 0.0] e:[(0, 3), (2, 3), (3, 7), (3, 33)] f:[(0, 1, 2, 3), (0, 3, 4, 7), (2, 3, 32, 33), (3, 7, 33, 35)]\n\
             i:4 x:[0.0, 0.0, 2.0] e:[(0, 4), (4, 5), (4, 7), (4, 51)] f:[(0, 1, 4, 5), (0, 3, 4, 7), (4, 5, 51, 52), (4, 7, 51, 54)]\n\
             i:5 x:[1.0, 0.0, 2.0] e:[(1, 5), (4, 5), (5, 22), (5, 52)] f:[(0, 1, 4, 5), (1, 5, 20, 22), (4, 5, 51, 52), (5, 22, 52, 63)]\n\
             i:6 x:[1.0, 1.0, 2.0] e:[] f:[]\n\
             i:7 x:[0.0, 1.0, 2.0] e:[(3, 7), (4, 7), (7, 35), (7, 54)] f:[(0, 3, 4, 7), (3, 7, 33, 35), (4, 7, 51, 54), (7, 35, 54, 71)]\n\
             i:8 x:[0.5, 0.0, 0.0] e:[(0, 1)] f:[(0, 1, 2, 3), (0, 1, 4, 5)]\n\
             i:9 x:[1.0, 0.5, 0.0] e:[(1, 2)] f:[(0, 1, 2, 3), (1, 2, 20, 21)]\n\
             i:10 x:[0.5, 1.0, 0.0] e:[(2, 3)] f:[(0, 1, 2, 3), (2, 3, 32, 33)]\n\
             i:11 x:[0.0, 0.5, 0.0] e:[(0, 3)] f:[(0, 1, 2, 3), (0, 3, 4, 7)]\n\
             i:12 x:[0.5, 0.0, 2.0] e:[(4, 5)] f:[(0, 1, 4, 5), (4, 5, 51, 52)]\n\
             i:13 x:[1.0, 0.5, 2.0] e:[] f:[]\n\
             i:14 x:[0.5, 1.0, 2.0] e:[] f:[]\n\
             i:15 x:[0.0, 0.5, 2.0] e:[(4, 7)] f:[(0, 3, 4, 7), (4, 7, 51, 54)]\n\
             i:16 x:[0.0, 0.0, 1.0] e:[(0, 4)] f:[(0, 1, 4, 5), (0, 3, 4, 7)]\n\
             i:17 x:[1.0, 0.0, 1.0] e:[(1, 5)] f:[(0, 1, 4, 5), (1, 5, 20, 22)]\n\
             i:18 x:[1.0, 1.0, 1.0] e:[] f:[]\n\
             i:19 x:[0.0, 1.0, 1.0] e:[(3, 7)] f:[(0, 3, 4, 7), (3, 7, 33, 35)]\n\
             i:20 x:[2.0, 0.0, 0.0] e:[(1, 20), (20, 21), (20, 22)] f:[(1, 2, 20, 21), (1, 5, 20, 22), (20, 21, 22, 23)]\n\
             i:21 x:[2.0, 1.0, 0.0] e:[(2, 21), (20, 21), (21, 23), (21, 44)] f:[(1, 2, 20, 21), (2, 21, 32, 44), (20, 21, 22, 23), (21, 23, 44, 45)]\n\
             i:22 x:[2.0, 0.0, 2.0] e:[(5, 22), (20, 22), (22, 23), (22, 63)] f:[(1, 5, 20, 22), (5, 22, 52, 63), (20, 21, 22, 23), (22, 23, 63, 64)]\n\
             i:23 x:[2.0, 1.0, 2.0] e:[(21, 23), (22, 23), (23, 45), (23, 64)] f:[(20, 21, 22, 23), (21, 23, 44, 45), (22, 23, 63, 64), (23, 45, 64, 77)]\n\
             i:24 x:[1.5, 0.0, 0.0] e:[(1, 20)] f:[(1, 2, 20, 21), (1, 5, 20, 22)]\n\
             i:25 x:[2.0, 0.5, 0.0] e:[(20, 21)] f:[(1, 2, 20, 21), (20, 21, 22, 23)]\n\
             i:26 x:[1.5, 1.0, 0.0] e:[(2, 21)] f:[(1, 2, 20, 21), (2, 21, 32, 44)]\n\
             i:27 x:[1.5, 0.0, 2.0] e:[(5, 22)] f:[(1, 5, 20, 22), (5, 22, 52, 63)]\n\
             i:28 x:[2.0, 0.5, 2.0] e:[(22, 23)] f:[(20, 21, 22, 23), (22, 23, 63, 64)]\n\
             i:29 x:[1.5, 1.0, 2.0] e:[] f:[]\n\
             i:30 x:[2.0, 0.0, 1.0] e:[(20, 22)] f:[(1, 5, 20, 22), (20, 21, 22, 23)]\n\
             i:31 x:[2.0, 1.0, 1.0] e:[(21, 23)] f:[(20, 21, 22, 23), (21, 23, 44, 45)]\n\
             i:32 x:[1.0, 2.0, 0.0] e:[(2, 32), (32, 33), (32, 34), (32, 44)] f:[(2, 3, 32, 33), (2, 21, 32, 44), (32, 33, 34, 35), (32, 34, 44, 45)]\n\
             i:33 x:[0.0, 2.0, 0.0] e:[(3, 33), (32, 33), (33, 35)] f:[(2, 3, 32, 33), (3, 7, 33, 35), (32, 33, 34, 35)]\n\
             i:34 x:[1.0, 2.0, 2.0] e:[(32, 34), (34, 35), (34, 45), (34, 70)] f:[(32, 33, 34, 35), (32, 34, 44, 45), (34, 35, 70, 71), (34, 45, 70, 77)]\n\
             i:35 x:[0.0, 2.0, 2.0] e:[(7, 35), (33, 35), (34, 35), (35, 71)] f:[(3, 7, 33, 35), (7, 35, 54, 71), (32, 33, 34, 35), (34, 35, 70, 71)]\n\
             i:36 x:[1.0, 1.5, 0.0] e:[(2, 32)] f:[(2, 3, 32, 33), (2, 21, 32, 44)]\n\
             i:37 x:[0.5, 2.0, 0.0] e:[(32, 33)] f:[(2, 3, 32, 33), (32, 33, 34, 35)]\n\
             i:38 x:[0.0, 1.5, 0.0] e:[(3, 33)] f:[(2, 3, 32, 33), (3, 7, 33, 35)]\n\
             i:39 x:[1.0, 1.5, 2.0] e:[] f:[]\n\
             i:40 x:[0.5, 2.0, 2.0] e:[(34, 35)] f:[(32, 33, 34, 35), (34, 35, 70, 71)]\n\
             i:41 x:[0.0, 1.5, 2.0] e:[(7, 35)] f:[(3, 7, 33, 35), (7, 35, 54, 71)]\n\
             i:42 x:[1.0, 2.0, 1.0] e:[(32, 34)] f:[(32, 33, 34, 35), (32, 34, 44, 45)]\n\
             i:43 x:[0.0, 2.0, 1.0] e:[(33, 35)] f:[(3, 7, 33, 35), (32, 33, 34, 35)]\n\
             i:44 x:[2.0, 2.0, 0.0] e:[(21, 44), (32, 44), (44, 45)] f:[(2, 21, 32, 44), (21, 23, 44, 45), (32, 34, 44, 45)]\n\
             i:45 x:[2.0, 2.0, 2.0] e:[(23, 45), (34, 45), (44, 45), (45, 77)] f:[(21, 23, 44, 45), (23, 45, 64, 77), (32, 34, 44, 45), (34, 45, 70, 77)]\n\
             i:46 x:[2.0, 1.5, 0.0] e:[(21, 44)] f:[(2, 21, 32, 44), (21, 23, 44, 45)]\n\
             i:47 x:[1.5, 2.0, 0.0] e:[(32, 44)] f:[(2, 21, 32, 44), (32, 34, 44, 45)]\n\
             i:48 x:[2.0, 1.5, 2.0] e:[(23, 45)] f:[(21, 23, 44, 45), (23, 45, 64, 77)]\n\
             i:49 x:[1.5, 2.0, 2.0] e:[(34, 45)] f:[(32, 34, 44, 45), (34, 45, 70, 77)]\n\
             i:50 x:[2.0, 2.0, 1.0] e:[(44, 45)] f:[(21, 23, 44, 45), (32, 34, 44, 45)]\n\
             i:51 x:[0.0, 0.0, 4.0] e:[(4, 51), (51, 52), (51, 54)] f:[(4, 5, 51, 52), (4, 7, 51, 54), (51, 52, 53, 54)]\n\
             i:52 x:[1.0, 0.0, 4.0] e:[(5, 52), (51, 52), (52, 53), (52, 63)] f:[(4, 5, 51, 52), (5, 22, 52, 63), (51, 52, 53, 54), (52, 53, 63, 64)]\n\
             i:53 x:[1.0, 1.0, 4.0] e:[(52, 53), (53, 54), (53, 64), (53, 70)] f:[(51, 52, 53, 54), (52, 53, 63, 64), (53, 54, 70, 71), (53, 64, 70, 77)]\n\
             i:54 x:[0.0, 1.0, 4.0] e:[(7, 54), (51, 54), (53, 54), (54, 71)] f:[(4, 7, 51, 54), (7, 35, 54, 71), (51, 52, 53, 54), (53, 54, 70, 71)]\n\
             i:55 x:[0.5, 0.0, 4.0] e:[(51, 52)] f:[(4, 5, 51, 52), (51, 52, 53, 54)]\n\
             i:56 x:[1.0, 0.5, 4.0] e:[(52, 53)] f:[(51, 52, 53, 54), (52, 53, 63, 64)]\n\
             i:57 x:[0.5, 1.0, 4.0] e:[(53, 54)] f:[(51, 52, 53, 54), (53, 54, 70, 71)]\n\
             i:58 x:[0.0, 0.5, 4.0] e:[(51, 54)] f:[(4, 7, 51, 54), (51, 52, 53, 54)]\n\
             i:59 x:[0.0, 0.0, 3.0] e:[(4, 51)] f:[(4, 5, 51, 52), (4, 7, 51, 54)]\n\
             i:60 x:[1.0, 0.0, 3.0] e:[(5, 52)] f:[(4, 5, 51, 52), (5, 22, 52, 63)]\n\
             i:61 x:[1.0, 1.0, 3.0] e:[] f:[]\n\
             i:62 x:[0.0, 1.0, 3.0] e:[(7, 54)] f:[(4, 7, 51, 54), (7, 35, 54, 71)]\n\
             i:63 x:[2.0, 0.0, 4.0] e:[(22, 63), (52, 63), (63, 64)] f:[(5, 22, 52, 63), (22, 23, 63, 64), (52, 53, 63, 64)]\n\
             i:64 x:[2.0, 1.0, 4.0] e:[(23, 64), (53, 64), (63, 64), (64, 77)] f:[(22, 23, 63, 64), (23, 45, 64, 77), (52, 53, 63, 64), (53, 64, 70, 77)]\n\
             i:65 x:[1.5, 0.0, 4.0] e:[(52, 63)] f:[(5, 22, 52, 63), (52, 53, 63, 64)]\n\
             i:66 x:[2.0, 0.5, 4.0] e:[(63, 64)] f:[(22, 23, 63, 64), (52, 53, 63, 64)]\n\
             i:67 x:[1.5, 1.0, 4.0] e:[(53, 64)] f:[(52, 53, 63, 64), (53, 64, 70, 77)]\n\
             i:68 x:[2.0, 0.0, 3.0] e:[(22, 63)] f:[(5, 22, 52, 63), (22, 23, 63, 64)]\n\
             i:69 x:[2.0, 1.0, 3.0] e:[(23, 64)] f:[(22, 23, 63, 64), (23, 45, 64, 77)]\n\
             i:70 x:[1.0, 2.0, 4.0] e:[(34, 70), (53, 70), (70, 71), (70, 77)] f:[(34, 35, 70, 71), (34, 45, 70, 77), (53, 54, 70, 71), (53, 64, 70, 77)]\n\
             i:71 x:[0.0, 2.0, 4.0] e:[(35, 71), (54, 71), (70, 71)] f:[(7, 35, 54, 71), (34, 35, 70, 71), (53, 54, 70, 71)]\n\
             i:72 x:[1.0, 1.5, 4.0] e:[(53, 70)] f:[(53, 54, 70, 71), (53, 64, 70, 77)]\n\
             i:73 x:[0.5, 2.0, 4.0] e:[(70, 71)] f:[(34, 35, 70, 71), (53, 54, 70, 71)]\n\
             i:74 x:[0.0, 1.5, 4.0] e:[(54, 71)] f:[(7, 35, 54, 71), (53, 54, 70, 71)]\n\
             i:75 x:[1.0, 2.0, 3.0] e:[(34, 70)] f:[(34, 35, 70, 71), (34, 45, 70, 77)]\n\
             i:76 x:[0.0, 2.0, 3.0] e:[(35, 71)] f:[(7, 35, 54, 71), (34, 35, 70, 71)]\n\
             i:77 x:[2.0, 2.0, 4.0] e:[(45, 77), (64, 77), (70, 77)] f:[(23, 45, 64, 77), (34, 45, 70, 77), (53, 64, 70, 77)]\n\
             i:78 x:[2.0, 1.5, 4.0] e:[(64, 77)] f:[(23, 45, 64, 77), (53, 64, 70, 77)]\n\
             i:79 x:[1.5, 2.0, 4.0] e:[(70, 77)] f:[(34, 45, 70, 77), (53, 64, 70, 77)]\n\
             i:80 x:[2.0, 2.0, 3.0] e:[(45, 77)] f:[(23, 45, 64, 77), (34, 45, 70, 77)]\n\
             \n\
             CELLS\n\
             =====\n\
             i:0 a:1 g:3 p:[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]\n\
             i:1 a:1 g:3 p:[1, 20, 21, 2, 5, 22, 23, 6, 24, 25, 26, 9, 27, 28, 29, 13, 17, 30, 31, 18]\n\
             i:2 a:1 g:3 p:[3, 2, 32, 33, 7, 6, 34, 35, 10, 36, 37, 38, 14, 39, 40, 41, 19, 18, 42, 43]\n\
             i:3 a:1 g:3 p:[2, 21, 44, 32, 6, 23, 45, 34, 26, 46, 47, 36, 29, 48, 49, 39, 18, 31, 50, 42]\n\
             i:4 a:1 g:3 p:[4, 5, 6, 7, 51, 52, 53, 54, 12, 13, 14, 15, 55, 56, 57, 58, 59, 60, 61, 62]\n\
             i:5 a:1 g:3 p:[5, 22, 23, 6, 52, 63, 64, 53, 27, 28, 29, 13, 65, 66, 67, 56, 60, 68, 69, 61]\n\
             i:6 a:1 g:3 p:[7, 6, 34, 35, 54, 53, 70, 71, 14, 39, 40, 41, 57, 72, 73, 74, 62, 61, 75, 76]\n\
             i:7 a:1 g:3 p:[6, 23, 45, 34, 53, 64, 77, 70, 29, 48, 49, 39, 67, 78, 79, 72, 61, 69, 80, 75]\n\
             \n\
             BOUNDARY POINTS\n\
             ===============\n\
             [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 35, 36, 37, 38, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80]\n\
             \n\
             BOUNDARY EDGES\n\
             ==============\n\
             k:(0,1) p:[0, 1, 8] c:[] f:[(0, 1, 2, 3), (0, 1, 4, 5)]\n\
             k:(0,3) p:[3, 0, 11] c:[] f:[(0, 1, 2, 3), (0, 3, 4, 7)]\n\
             k:(0,4) p:[0, 4, 16] c:[] f:[(0, 1, 4, 5), (0, 3, 4, 7)]\n\
             k:(1,2) p:[1, 2, 9] c:[] f:[(0, 1, 2, 3), (1, 2, 20, 21)]\n\
             k:(1,5) p:[5, 1, 17] c:[] f:[(0, 1, 4, 5), (1, 5, 20, 22)]\n\
             k:(1,20) p:[1, 20, 24] c:[] f:[(1, 2, 20, 21), (1, 5, 20, 22)]\n\
             k:(2,3) p:[2, 3, 10] c:[] f:[(0, 1, 2, 3), (2, 3, 32, 33)]\n\
             k:(2,21) p:[21, 2, 26] c:[] f:[(1, 2, 20, 21), (2, 21, 32, 44)]\n\
             k:(2,32) p:[2, 32, 36] c:[] f:[(2, 3, 32, 33), (2, 21, 32, 44)]\n\
             k:(3,7) p:[3, 7, 19] c:[] f:[(0, 3, 4, 7), (3, 7, 33, 35)]\n\
             k:(3,33) p:[33, 3, 38] c:[] f:[(2, 3, 32, 33), (3, 7, 33, 35)]\n\
             k:(4,5) p:[4, 5, 12] c:[] f:[(0, 1, 4, 5), (4, 5, 51, 52)]\n\
             k:(4,7) p:[7, 4, 15] c:[] f:[(0, 3, 4, 7), (4, 7, 51, 54)]\n\
             k:(4,51) p:[4, 51, 59] c:[] f:[(4, 5, 51, 52), (4, 7, 51, 54)]\n\
             k:(5,22) p:[5, 22, 27] c:[] f:[(1, 5, 20, 22), (5, 22, 52, 63)]\n\
             k:(5,52) p:[52, 5, 60] c:[] f:[(4, 5, 51, 52), (5, 22, 52, 63)]\n\
             k:(7,35) p:[35, 7, 41] c:[] f:[(3, 7, 33, 35), (7, 35, 54, 71)]\n\
             k:(7,54) p:[7, 54, 62] c:[] f:[(4, 7, 51, 54), (7, 35, 54, 71)]\n\
             k:(20,21) p:[20, 21, 25] c:[] f:[(1, 2, 20, 21), (20, 21, 22, 23)]\n\
             k:(20,22) p:[22, 20, 30] c:[] f:[(1, 5, 20, 22), (20, 21, 22, 23)]\n\
             k:(21,23) p:[23, 21, 31] c:[] f:[(20, 21, 22, 23), (21, 23, 44, 45)]\n\
             k:(21,44) p:[21, 44, 46] c:[] f:[(2, 21, 32, 44), (21, 23, 44, 45)]\n\
             k:(22,23) p:[22, 23, 28] c:[] f:[(20, 21, 22, 23), (22, 23, 63, 64)]\n\
             k:(22,63) p:[63, 22, 68] c:[] f:[(5, 22, 52, 63), (22, 23, 63, 64)]\n\
             k:(23,45) p:[23, 45, 48] c:[] f:[(21, 23, 44, 45), (23, 45, 64, 77)]\n\
             k:(23,64) p:[64, 23, 69] c:[] f:[(22, 23, 63, 64), (23, 45, 64, 77)]\n\
             k:(32,33) p:[32, 33, 37] c:[] f:[(2, 3, 32, 33), (32, 33, 34, 35)]\n\
             k:(32,34) p:[32, 34, 42] c:[] f:[(32, 33, 34, 35), (32, 34, 44, 45)]\n\
             k:(32,44) p:[44, 32, 47] c:[] f:[(2, 21, 32, 44), (32, 34, 44, 45)]\n\
             k:(33,35) p:[33, 35, 43] c:[] f:[(3, 7, 33, 35), (32, 33, 34, 35)]\n\
             k:(34,35) p:[34, 35, 40] c:[] f:[(32, 33, 34, 35), (34, 35, 70, 71)]\n\
             k:(34,45) p:[45, 34, 49] c:[] f:[(32, 34, 44, 45), (34, 45, 70, 77)]\n\
             k:(34,70) p:[34, 70, 75] c:[] f:[(34, 35, 70, 71), (34, 45, 70, 77)]\n\
             k:(35,71) p:[35, 71, 76] c:[] f:[(7, 35, 54, 71), (34, 35, 70, 71)]\n\
             k:(44,45) p:[45, 44, 50] c:[] f:[(21, 23, 44, 45), (32, 34, 44, 45)]\n\
             k:(45,77) p:[77, 45, 80] c:[] f:[(23, 45, 64, 77), (34, 45, 70, 77)]\n\
             k:(51,52) p:[51, 52, 55] c:[] f:[(4, 5, 51, 52), (51, 52, 53, 54)]\n\
             k:(51,54) p:[54, 51, 58] c:[] f:[(4, 7, 51, 54), (51, 52, 53, 54)]\n\
             k:(52,53) p:[53, 52, 56] c:[] f:[(51, 52, 53, 54), (52, 53, 63, 64)]\n\
             k:(52,63) p:[52, 63, 65] c:[] f:[(5, 22, 52, 63), (52, 53, 63, 64)]\n\
             k:(53,54) p:[54, 53, 57] c:[] f:[(51, 52, 53, 54), (53, 54, 70, 71)]\n\
             k:(53,64) p:[53, 64, 67] c:[] f:[(52, 53, 63, 64), (53, 64, 70, 77)]\n\
             k:(53,70) p:[70, 53, 72] c:[] f:[(53, 54, 70, 71), (53, 64, 70, 77)]\n\
             k:(54,71) p:[71, 54, 74] c:[] f:[(7, 35, 54, 71), (53, 54, 70, 71)]\n\
             k:(63,64) p:[63, 64, 66] c:[] f:[(22, 23, 63, 64), (52, 53, 63, 64)]\n\
             k:(64,77) p:[64, 77, 78] c:[] f:[(23, 45, 64, 77), (53, 64, 70, 77)]\n\
             k:(70,71) p:[70, 71, 73] c:[] f:[(34, 35, 70, 71), (53, 54, 70, 71)]\n\
             k:(70,77) p:[77, 70, 79] c:[] f:[(34, 45, 70, 77), (53, 64, 70, 77)]\n\
             \n\
             BOUNDARY FACES\n\
             ==============\n\
             k:(0,1,2,3) p:[0, 3, 2, 1, 11, 10, 9, 8] c:[0]\n\
             k:(0,1,4,5) p:[0, 1, 5, 4, 8, 17, 12, 16] c:[0]\n\
             k:(0,3,4,7) p:[0, 4, 7, 3, 16, 15, 19, 11] c:[0]\n\
             k:(1,2,20,21) p:[1, 2, 21, 20, 9, 26, 25, 24] c:[1]\n\
             k:(1,5,20,22) p:[1, 20, 22, 5, 24, 30, 27, 17] c:[1]\n\
             k:(2,3,32,33) p:[3, 33, 32, 2, 38, 37, 36, 10] c:[2]\n\
             k:(2,21,32,44) p:[2, 32, 44, 21, 36, 47, 46, 26] c:[3]\n\
             k:(3,7,33,35) p:[3, 7, 35, 33, 19, 41, 43, 38] c:[2]\n\
             k:(4,5,51,52) p:[4, 5, 52, 51, 12, 60, 55, 59] c:[4]\n\
             k:(4,7,51,54) p:[4, 51, 54, 7, 59, 58, 62, 15] c:[4]\n\
             k:(5,22,52,63) p:[5, 22, 63, 52, 27, 68, 65, 60] c:[5]\n\
             k:(7,35,54,71) p:[7, 54, 71, 35, 62, 74, 76, 41] c:[6]\n\
             k:(20,21,22,23) p:[20, 21, 23, 22, 25, 31, 28, 30] c:[1]\n\
             k:(21,23,44,45) p:[21, 44, 45, 23, 46, 50, 48, 31] c:[3]\n\
             k:(22,23,63,64) p:[22, 23, 64, 63, 28, 69, 66, 68] c:[5]\n\
             k:(23,45,64,77) p:[23, 45, 77, 64, 48, 80, 78, 69] c:[7]\n\
             k:(32,33,34,35) p:[32, 33, 35, 34, 37, 43, 40, 42] c:[2]\n\
             k:(32,34,44,45) p:[44, 32, 34, 45, 47, 42, 49, 50] c:[3]\n\
             k:(34,35,70,71) p:[34, 35, 71, 70, 40, 76, 73, 75] c:[6]\n\
             k:(34,45,70,77) p:[45, 34, 70, 77, 49, 75, 79, 80] c:[7]\n\
             k:(51,52,53,54) p:[51, 52, 53, 54, 55, 56, 57, 58] c:[4]\n\
             k:(52,53,63,64) p:[52, 63, 64, 53, 65, 66, 67, 56] c:[5]\n\
             k:(53,54,70,71) p:[54, 53, 70, 71, 57, 72, 73, 74] c:[6]\n\
             k:(53,64,70,77) p:[53, 64, 77, 70, 67, 78, 79, 72] c:[7]\n"
        );
        Ok(())
    }
}
