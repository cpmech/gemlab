use super::{Cell, Constraint, Mesh, Point};
use crate::shapes::{Shape, StateOfShape};
use crate::util::{AsArray2D, GridSearch, GsNdiv, GsTol};
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
    attribute_id: usize,      // attribute ID of all elements in this block
    space_ndim: usize,        // space dimension
    ndiv: Vec<usize>,         // number of divisions along each dim (ndim)
    delta_ksi: Vec<Vec<f64>>, // delta ksi along each dim (ndim, {ndiv[0],ndiv[1],ndiv[2]})

    edge_constraints: Vec<Option<Constraint>>, // constraints (nedge)
    face_constraints: Vec<Option<Constraint>>, // constraints (nface)

    // shape and state
    shape: Shape,
    state: StateOfShape,

    // grid to search reference coordinates
    grid_ksi: GridSearch,
}

impl Block {
    // constants
    const NAT_LENGTH: f64 = 2.0; // length of shape along each direction in reference coords space

    // valid output npoint
    const VALID_OUTPUT_NPOINT_2D: [usize; 6] = [4, 8, 9, 12, 16, 17];
    const VALID_OUTPUT_NPOINT_3D: [usize; 2] = [8, 20];

    /// Allocate a new instance
    ///
    /// # Input (2D)
    ///
    /// * `coords` -- (nrow = 4 or 8, ncol = space_ndim = 2) matrix with all coordinates
    ///
    /// # Input (3D)
    ///
    /// * `coords` -- (nrow = 8 or 20, ncol = space_ndim = 3) matrix with all coordinates
    pub fn new<'a, T>(coords: &'a T) -> Result<Self, StrError>
    where
        T: AsArray2D<'a, f64>,
    {
        // check
        let (nrow, space_ndim) = coords.size();
        if space_ndim < 2 || space_ndim > 3 {
            return Err("ncol = space_ndim must be 2 or 3");
        }
        if space_ndim == 2 {
            if nrow != 4 && nrow != 8 {
                return Err("in 2D, nrow must be 4 or 8");
            }
        } else {
            if nrow != 8 && nrow != 20 {
                return Err("in 3D, nrow must be 8 or 20");
            }
        }

        // shape
        let geo_ndim = space_ndim;
        let nnode = if geo_ndim == 2 { 8 } else { 20 }; // Qua8 or Hex20
        let shape = Shape::new(space_ndim, geo_ndim, nnode).unwrap(); // should not fail here

        // state
        let state = if space_ndim == 2 {
            if nrow == 8 {
                // all vertices given, ok
                StateOfShape::new(geo_ndim, coords).unwrap() // should not fail here
            } else {
                // copy "corner" vertices
                let mut xx = vec![vec![0.0; space_ndim]; nnode];
                for m in 0..nrow {
                    for j in 0..space_ndim {
                        xx[m][j] = coords.at(m, j);
                    }
                }
                // generate mid vertices
                for j in 0..space_ndim {
                    xx[4][j] = (coords.at(0, j) + coords.at(1, j)) / 2.0;
                    xx[5][j] = (coords.at(1, j) + coords.at(2, j)) / 2.0;
                    xx[6][j] = (coords.at(2, j) + coords.at(3, j)) / 2.0;
                    xx[7][j] = (coords.at(3, j) + coords.at(0, j)) / 2.0;
                }
                StateOfShape::new(geo_ndim, &xx)?
            }
        } else {
            if nrow == 20 {
                // all vertices given, ok
                StateOfShape::new(geo_ndim, coords).unwrap() // should not fail here
            } else {
                // copy "corner" vertices
                let mut xx = vec![vec![0.0; space_ndim]; nnode];
                for m in 0..nrow {
                    for j in 0..space_ndim {
                        xx[m][j] = coords.at(m, j);
                    }
                }
                // generate mid vertices
                for j in 0..space_ndim {
                    xx[8][j] = (coords.at(0, j) + coords.at(1, j)) / 2.0;
                    xx[9][j] = (coords.at(1, j) + coords.at(2, j)) / 2.0;
                    xx[10][j] = (coords.at(2, j) + coords.at(3, j)) / 2.0;
                    xx[11][j] = (coords.at(3, j) + coords.at(0, j)) / 2.0;

                    xx[12][j] = (coords.at(4, j) + coords.at(5, j)) / 2.0;
                    xx[13][j] = (coords.at(5, j) + coords.at(6, j)) / 2.0;
                    xx[14][j] = (coords.at(6, j) + coords.at(7, j)) / 2.0;
                    xx[15][j] = (coords.at(7, j) + coords.at(4, j)) / 2.0;

                    xx[16][j] = (coords.at(0, j) + coords.at(4, j)) / 2.0;
                    xx[17][j] = (coords.at(1, j) + coords.at(5, j)) / 2.0;
                    xx[18][j] = (coords.at(2, j) + coords.at(6, j)) / 2.0;
                    xx[19][j] = (coords.at(3, j) + coords.at(7, j)) / 2.0;
                }
                StateOfShape::new(geo_ndim, &xx)?
            }
        };

        // grid search
        let grid_ksi = GridSearch::new(
            &vec![-1.0; space_ndim],
            &vec![1.0; space_ndim],
            GsNdiv::Default,
            GsTol::Default,
        )?;

        // done
        const NDIV: usize = 2;
        Ok(Block {
            attribute_id: 1,
            space_ndim,
            ndiv: vec![NDIV; space_ndim],
            delta_ksi: vec![vec![1.0; NDIV]; space_ndim],
            edge_constraints: vec![None; shape.nedge],
            face_constraints: vec![None; shape.nface],
            shape,
            state,
            grid_ksi,
        })
    }

    /// Sets group
    pub fn set_attribute_id(&mut self, attribute_id: usize) -> &mut Self {
        self.attribute_id = attribute_id;
        self
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
        let mut mesh = Mesh {
            space_ndim,
            points: Vec::new(),
            cells: Vec::new(),
        };

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
                        let ksi_ref = shape_out.reference_coords(m);

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
                                self.shape.calc_coords(&mut x, &mut self.state, &ksi)?;

                                // add new point to mesh
                                let mut shared_by_cells = HashSet::new();
                                shared_by_cells.insert(cell_id);
                                mesh.points.push(Point {
                                    id: point_id,
                                    coords: x.as_data().clone(),
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
        Ok(mesh)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Block, StrError};
    use crate::geometry::Circle;
    use crate::mesh::{Constraint, Samples};
    use russell_chk::assert_vec_approx_eq;

    #[test]
    fn new_fails_on_wrong_input() {
        assert_eq!(Block::new(&[[]]).err(), Some("ncol = space_ndim must be 2 or 3"));
        assert_eq!(Block::new(&[[0.0, 0.0]]).err(), Some("in 2D, nrow must be 4 or 8"));
        assert_eq!(
            Block::new(&[[0.0, 0.0, 0.0]]).err(),
            Some("in 3D, nrow must be 8 or 20")
        );
    }

    #[test]
    fn new_works() -> Result<(), StrError> {
        let b2d = Block::new(&[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]])?;
        assert_eq!(b2d.attribute_id, 1);
        assert_eq!(b2d.space_ndim, 2);
        assert_eq!(
            format!("{}", b2d.state.coords_transp),
            "┌                 ┐\n\
             │ 0 2 2 0 1 2 1 0 │\n\
             │ 0 0 2 2 0 1 2 1 │\n\
             └                 ┘"
        );
        assert_eq!(b2d.ndiv, &[2, 2]);
        assert_eq!(format!("{:?}", b2d.delta_ksi), "[[1.0, 1.0], [1.0, 1.0]]");
        assert_eq!(b2d.shape.nnode, 8);

        let b3d = Block::new(&[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [2.0, 0.0, 2.0],
            [2.0, 2.0, 2.0],
            [0.0, 2.0, 2.0],
        ])?;
        assert_eq!(b3d.attribute_id, 1);
        assert_eq!(b3d.space_ndim, 3);
        assert_eq!(
            format!("{}", b3d.state.coords_transp),
            "┌                                         ┐\n\
             │ 0 2 2 0 0 2 2 0 1 2 1 0 1 2 1 0 0 2 2 0 │\n\
             │ 0 0 2 2 0 0 2 2 0 1 2 1 0 1 2 1 0 0 2 2 │\n\
             │ 0 0 0 0 2 2 2 2 0 0 0 0 2 2 2 2 1 1 1 1 │\n\
             └                                         ┘"
        );
        assert_eq!(b3d.ndiv, &[2, 2, 2]);
        assert_eq!(format!("{:?}", b3d.delta_ksi), "[[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]");
        assert_eq!(b3d.shape.nnode, 20);
        Ok(())
    }

    #[test]
    fn set_attribute_id_works() -> Result<(), StrError> {
        let mut block = Block::new(&[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]])?;
        block.set_attribute_id(2);
        assert_eq!(block.attribute_id, 2);
        Ok(())
    }

    #[test]
    fn set_ndiv_works() -> Result<(), StrError> {
        let mut block = Block::new(&[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]])?;
        block.set_ndiv(&[2, 4]);
        assert_eq!(block.ndiv, &[2, 4]);
        assert_eq!(format!("{:?}", block.delta_ksi), "[[1.0, 1.0], [0.5, 0.5, 0.5, 0.5]]");
        Ok(())
    }

    #[test]
    fn set_edge_constraint_works() -> Result<(), StrError> {
        let mut block = Block::new(&[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]])?;
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
        let mut block = Block::new(&[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [2.0, 0.0, 2.0],
            [2.0, 2.0, 2.0],
            [0.0, 2.0, 2.0],
        ])?;
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
        let mut b2d = Block::new(&[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]])?;
        assert_eq!(b2d.subdivide(1).err(), Some("output_npoint is invalid"));
        let mut b3d = Block::new(&[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [2.0, 0.0, 2.0],
            [2.0, 2.0, 2.0],
            [0.0, 2.0, 2.0],
        ])?;
        assert_eq!(b3d.subdivide(1).err(), Some("output_npoint is invalid"));
        Ok(())
    }

    #[test]
    fn subdivide_2d_qua4_works() -> Result<(), StrError> {
        // 7---------------6---------------8
        // |               |               |
        // |               |               |
        // |      [2]      |      [3]      |
        // |               |               |
        // |               |               |
        // 3---------------2---------------5
        // |               |               |
        // |               |               |
        // |      [0]      |      [1]      |
        // |               |               |
        // |               |               |
        // 0---------------1---------------4
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ])?;
        let mesh = block.subdivide(4)?;
        //
        let correct = Samples::block_2d_four_qua4();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        Ok(())
    }

    #[test]
    fn subdivide_2d_qua8_works() -> Result<(), StrError> {
        // 14------16------13------20------18
        //  |               |               |
        //  |               |               |
        // 17      [2]     15      [3]     19
        //  |               |               |
        //  |               |               |
        //  3-------6-------2------12-------9
        //  |               |               |
        //  |               |               |
        //  7      [0]      5      [1]     11
        //  |               |               |
        //  |               |               |
        //  0-------4-------1------10-------8
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ])?;
        let mesh = block.subdivide(8)?;
        let correct = Samples::block_2d_four_qua8();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        Ok(())
    }

    #[test]
    fn subdivide_2d_qua9_works() -> Result<(), StrError> {
        // 16------18------15------23------21
        //  |               |               |
        //  |               |               |
        // 19      20      17      24      22
        //  |               |               |
        //  |               |               |
        //  3-------6-------2------13------10
        //  |               |               |
        //  |               |               |
        //  7       8       5      14      12
        //  |               |               |
        //  |               |               |
        //  0-------4-------1------11-------9
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ])?;
        let mesh = block.subdivide(9)?;
        let correct = Samples::block_2d_four_qua9();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        Ok(())
    }

    #[test]
    fn subdivide_2d_qua12_works() -> Result<(), StrError> {
        // 21---26---23----20---32---30----28
        //  |               |               |
        // 24              25              31
        //  |               |               |
        // 27              22              29
        //  |               |               |
        //  3---10-----6----2---19---16----13
        //  |               |               |
        //  7               9              18
        //  |               |               |
        // 11               5              15
        //  |               |               |
        //  0----4-----8----1---14---17----12
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [0.0, 3.0],
        ])?;
        let mesh = block.subdivide(12)?;
        let correct = Samples::block_2d_four_qua12();
        assert_eq!(mesh.points.len(), correct.points.len());
        assert_eq!(mesh.cells.len(), correct.cells.len());
        for point in &correct.points {
            assert_vec_approx_eq!(point.coords, correct.points[point.id].coords, 1e-15);
        }
        for cell in &correct.cells {
            assert_eq!(cell.points, correct.cells[cell.id].points);
        }
        Ok(())
    }

    #[test]
    fn subdivide_2d_qua16_works() -> Result<(), StrError> {
        // 29---34----31---28---44---42----40
        //  |               |               |
        // 32   39    38   33   48   47    43
        //  |               |               |
        // 35   36    37   30   45   46    41
        //  |               |               |
        //  3---10-----6----2---23---20----17
        //  |               |               |
        //  7   15    14    9   27   26    22
        //  |               |               |
        // 11   12    13    5   24   25    19
        //  |               |               |
        //  0----4-----8----1---18---21----16
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [0.0, 3.0],
        ])?;
        let mesh = block.subdivide(16)?;
        let correct = Samples::block_2d_four_qua16();
        assert_eq!(mesh.points.len(), correct.points.len());
        assert_eq!(mesh.cells.len(), correct.cells.len());
        for point in &correct.points {
            assert_vec_approx_eq!(point.coords, correct.points[point.id].coords, 1e-15);
        }
        for cell in &correct.cells {
            assert_eq!(cell.points, correct.cells[cell.id].points);
        }
        Ok(())
    }

    #[test]
    fn subdivide_2d_qua17_works() -> Result<(), StrError> {
        // 30---38---35---32---29---47---45---43---41
        //  |                   |                   |
        // 33                  37                  46
        //  |                   |                   |
        // 36        40        34        48        44
        //  |                   |                   |
        // 39                  31                  42
        //  |                   |                   |
        //  3---14---10----6----2---27---24---21---18
        //  |                   |                   |
        //  7                  13                  26
        //  |                   |                   |
        // 11        16         9        28        23
        //  |                   |                   |
        // 15                   5                  20
        //  |                   |                   |
        //  0----4----8---12----1---19---22---25---17
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
        ])?;
        let mesh = block.subdivide(17)?;
        let correct = Samples::block_2d_four_qua17();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        Ok(())
    }

    #[test]
    fn subdivide_3d_works() -> Result<(), StrError> {
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
        #[rustfmt::skip]
        let mut block = Block::new(&[
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
        let correct = Samples::block_3d_eight_hex8();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        Ok(())
    }

    #[test]
    fn subdivide_3d_o2_works() -> Result<(), StrError> {
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
        #[rustfmt::skip]
        let mut block = Block::new(&[
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
        let correct = Samples::block_3d_eight_hex20();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        Ok(())
    }
}
