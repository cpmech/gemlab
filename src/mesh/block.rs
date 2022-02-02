use super::{Cell, Constraint, Mesh, Point};
use crate::shapes::Shape;
use crate::util::{AsArray2D, GridSearch};
use crate::StrError;
use russell_lab::{mat_vec_mul, Matrix, Vector};
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
    npoint: usize,            // number of points (8 or 20) (quadratic block internally)
    nedge: usize,             // number of edges (4 or 12)
    nface: usize,             // number of faces (4 or 6)
    coords: Matrix,           // coordinates of points (npoint, ndim)
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
    const NAT_MIN: f64 = -1.0; // min reference coordinate value
    const NAT_MAX: f64 = 1.0; // max reference coordinate value
    const NAT_LENGTH: f64 = 2.0; // length of shape along each direction in reference coords space
    const NAT_TOLERANCE: f64 = 1e-4; // tolerance to compare coordinates in the reference space

    // valid output npoint
    const VALID_OUTPUT_NPOINT_2D: [usize; 5] = [4, 8, 9, 12, 16];
    const VALID_OUTPUT_NPOINT_3D: [usize; 2] = [8, 20];

    /// Creates a new Block with default options
    pub fn new(space_ndim: usize) -> Result<Self, StrError> {
        // check
        if space_ndim < 2 || space_ndim > 3 {
            return Err("space_ndim must be 2 or 3");
        }

        // shape
        let geo_ndim = space_ndim;
        let npoint = if geo_ndim == 2 { 8 } else { 20 };
        let shape = Shape::new(space_ndim, geo_ndim, npoint)?;

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
            npoint,
            nedge: shape.nedge,
            nface: shape.nface,
            coords: Matrix::new(npoint, space_ndim),
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

    /// Sets the coordinates of all 2D or 3D points
    ///
    /// # Input (2D)
    ///
    /// * `coords` -- (4 or 8, 2) matrix with all coordinates
    ///
    /// # Input (3D)
    ///
    /// * `coords` -- (8 or 20, 3) matrix with all coordinates
    pub fn set_coords<'a, T, U>(&mut self, coords: &'a T) -> &mut Self
    where
        T: AsArray2D<'a, U>,
        U: 'a + Into<f64>,
    {
        // check
        let (nrow, ncol) = coords.size();
        assert_eq!(ncol, self.space_ndim);
        if self.space_ndim == 2 {
            assert!(nrow == 4 || nrow == 8);
        } else {
            assert!(nrow == 8 || nrow == 20);
        }
        // set corner points
        for i in 0..nrow {
            for j in 0..ncol {
                self.coords[i][j] = coords.at(i, j).into();
            }
        }
        // generate mid points
        if self.space_ndim == 2 && nrow == 4 {
            for i in 0..self.space_ndim {
                self.coords[4][i] = (self.coords[0][i] + self.coords[1][i]) / 2.0;
                self.coords[5][i] = (self.coords[1][i] + self.coords[2][i]) / 2.0;
                self.coords[6][i] = (self.coords[2][i] + self.coords[3][i]) / 2.0;
                self.coords[7][i] = (self.coords[3][i] + self.coords[0][i]) / 2.0;
            }
        }
        if self.space_ndim == 3 && nrow == 8 {
            for i in 0..self.space_ndim {
                self.coords[8][i] = (self.coords[0][i] + self.coords[1][i]) / 2.0;
                self.coords[9][i] = (self.coords[1][i] + self.coords[2][i]) / 2.0;
                self.coords[10][i] = (self.coords[2][i] + self.coords[3][i]) / 2.0;
                self.coords[11][i] = (self.coords[3][i] + self.coords[0][i]) / 2.0;

                self.coords[12][i] = (self.coords[4][i] + self.coords[5][i]) / 2.0;
                self.coords[13][i] = (self.coords[5][i] + self.coords[6][i]) / 2.0;
                self.coords[14][i] = (self.coords[6][i] + self.coords[7][i]) / 2.0;
                self.coords[15][i] = (self.coords[7][i] + self.coords[4][i]) / 2.0;

                self.coords[16][i] = (self.coords[0][i] + self.coords[4][i]) / 2.0;
                self.coords[17][i] = (self.coords[1][i] + self.coords[5][i]) / 2.0;
                self.coords[18][i] = (self.coords[2][i] + self.coords[6][i]) / 2.0;
                self.coords[19][i] = (self.coords[3][i] + self.coords[7][i]) / 2.0;
            }
        }
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
    /// * 2D: [4, 8, 9, 12, 16]
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

        // auxiliary variables
        let geo_ndim = 2;
        let shape_out = Shape::new(space_ndim, geo_ndim, output_npoint)?;

        // transformation matrix: scale and translate reference space
        //   _                                       _
        //  |  scale_x   0.0    0.0    translation_x  |
        //  |    0.0   scale_y  0.0    translation_y  |
        //  |    0.0     0.0   scale_z translation_z  |
        //  |_   0.0     0.0    0.0         1.0      _|
        let ndim = space_ndim;
        let mut transform = Matrix::identity(ndim + 1);

        // augmented reference coordinates [r,s,1] or [r,s,t,1]
        let mut ksi_aug = Vector::new(ndim + 1);
        ksi_aug[ndim] = 1.0;

        // augmented transformed nat-coordinates
        let mut ksi = Vector::new(ndim + 1);

        // center of shape in nat-coords
        let mut center = vec![0.0; ndim];

        // real point coordinates
        let mut x = Vector::new(ndim);

        // number of divisions along each direction
        let (nx, ny, nz) = (self.ndiv[0], self.ndiv[1], if ndim == 2 { 1 } else { self.ndiv[2] });

        // for each z-division
        if ndim == 3 {
            center[2] = -1.0 + self.delta_ksi[2][0] / 2.0;
        }
        for k in 0..nz {
            if ndim == 3 {
                transform[2][2] = self.delta_ksi[2][k] / Block::NAT_LENGTH; // scale
                transform[2][ndim] = center[2]; // translation
            }

            // for each y-division
            center[1] = -1.0 + self.delta_ksi[1][0] / 2.0;
            for j in 0..ny {
                transform[1][1] = self.delta_ksi[1][j] / Block::NAT_LENGTH; // scale
                transform[1][ndim] = center[1]; // translation

                // for each x-division
                center[0] = -1.0 + self.delta_ksi[0][0] / 2.0;
                for i in 0..nx {
                    transform[0][0] = self.delta_ksi[0][i] / Block::NAT_LENGTH; // scale
                    transform[0][ndim] = center[0]; // translation

                    // new cell id
                    let cell_id = mesh.cells.len();

                    // for each point
                    let mut points = vec![0; shape_out.nnode];
                    for m in 0..shape_out.nnode {
                        // transform reference coords: scale and translate
                        // shape_out.get_reference_coords(&mut ksi_aug, m); // TODO
                        mat_vec_mul(&mut ksi, 1.0, &transform, &ksi_aug)?;

                        // maybe append point to mesh
                        let index = self.maybe_append_new_point(&mut mesh, &mut x, &ksi, cell_id)?;
                        points[m] = index;
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
            if ndim == 3 {
                center[2] += self.delta_ksi[2][k];
            }
        }

        // done
        Ok(mesh)
    }

    /// Appends new point, if not existent yet
    ///
    /// # Output
    ///
    /// * `mesh` -- updated mesh with new point
    /// * `x` -- real coordinates of the new/existent point point
    ///
    /// # Input
    ///
    /// * `ksi_vec` -- (augmented) reference coordinates of the point (scaled and translated already)
    /// * `cell_id` -- ID of the cell adding this point (to set shared-by information)
    ///
    /// # Returns
    ///
    /// Returns the index of the new or existent point.
    fn maybe_append_new_point(
        &mut self,
        mesh: &mut Mesh,
        x: &mut Vector,
        ksi_vec: &Vector,
        cell_id: usize,
    ) -> Result<usize, StrError> {
        // skip if existent point
        let ndim = self.space_ndim;
        let ksi = &ksi_vec.as_data()[0..ndim];
        if let Some(index) = self.grid_ksi.find(ksi)? {
            return Ok(index);
        }

        // insert point in grid search object
        let index = mesh.points.len();
        self.grid_ksi.insert(index, ksi)?;

        // compute real coords
        // self.shape.calc_coords(x, ksi_vec)?; // TODO
        // self.shape.calc_ss_vec(ksi_vec);
        // self.shape.mul_interp_by_matrix(x, &self.coords)?;

        // add new point to mesh
        let mut shared_by_cells = HashSet::new();
        shared_by_cells.insert(cell_id);
        mesh.points.push(Point {
            id: index,
            coords: x.as_data().clone(),
            shared_by_boundary_edges: HashSet::new(),
            shared_by_boundary_faces: HashSet::new(),
        });

        // done
        Ok(index)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Block, StrError};
    use crate::{geometry::Circle, mesh::Constraint};

    #[test]
    fn new_fails_on_wrong_input() {
        let b2d = Block::new(1);
        assert_eq!(b2d.err(), Some("ndim must be 2 or 3"));
    }

    #[test]
    fn new_works() -> Result<(), StrError> {
        let b2d = Block::new(2)?;
        assert_eq!(b2d.attribute_id, 1);
        assert_eq!(b2d.space_ndim, 2);
        assert_eq!(b2d.npoint, 8);
        assert_eq!(b2d.nedge, 4);
        assert_eq!(b2d.nface, 0);
        assert_eq!(
            format!("{}", b2d.coords),
            "┌     ┐\n\
             │ 0 0 │\n\
             │ 0 0 │\n\
             │ 0 0 │\n\
             │ 0 0 │\n\
             │ 0 0 │\n\
             │ 0 0 │\n\
             │ 0 0 │\n\
             │ 0 0 │\n\
             └     ┘"
        );
        assert_eq!(b2d.ndiv, &[2, 2]);
        assert_eq!(format!("{:?}", b2d.delta_ksi), "[[1.0, 1.0], [1.0, 1.0]]");
        assert_eq!(b2d.shape.nnode, 8);

        let b3d = Block::new(3)?;
        assert_eq!(b3d.attribute_id, 1);
        assert_eq!(b3d.space_ndim, 3);
        assert_eq!(b3d.npoint, 20);
        assert_eq!(b3d.nedge, 12);
        assert_eq!(b3d.nface, 6);
        assert_eq!(
            format!("{}", b3d.coords),
            "┌       ┐\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             │ 0 0 0 │\n\
             └       ┘"
        );
        assert_eq!(b3d.ndiv, &[2, 2, 2]);
        assert_eq!(format!("{:?}", b3d.delta_ksi), "[[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]");
        assert_eq!(b3d.shape.nnode, 20);
        Ok(())
    }

    #[test]
    fn set_group_works() -> Result<(), StrError> {
        let mut block = Block::new(2)?;
        block.set_attribute_id(2);
        assert_eq!(block.attribute_id, 2);
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
        ]);
        assert_eq!(
            format!("{}", b2d.coords),
            "┌     ┐\n\
             │ 0 0 │\n\
             │ 2 0 │\n\
             │ 2 2 │\n\
             │ 0 2 │\n\
             │ 1 0 │\n\
             │ 2 1 │\n\
             │ 1 2 │\n\
             │ 0 1 │\n\
             └     ┘"
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
        ]);
        assert_eq!(
            format!("{}", b3d.coords),
            "┌       ┐\n\
             │ 0 0 0 │\n\
             │ 2 0 0 │\n\
             │ 2 2 0 │\n\
             │ 0 2 0 │\n\
             │ 0 0 2 │\n\
             │ 2 0 2 │\n\
             │ 2 2 2 │\n\
             │ 0 2 2 │\n\
             │ 1 0 0 │\n\
             │ 2 1 0 │\n\
             │ 1 2 0 │\n\
             │ 0 1 0 │\n\
             │ 1 0 2 │\n\
             │ 2 1 2 │\n\
             │ 1 2 2 │\n\
             │ 0 1 2 │\n\
             │ 0 0 1 │\n\
             │ 2 0 1 │\n\
             │ 2 2 1 │\n\
             │ 0 2 1 │\n\
             └       ┘"
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
    fn subdivide_2d_works() -> Result<(), StrError> {
        let mut block = Block::new(2)?;
        #[rustfmt::skip]
        block.set_coords(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ]);
        let mesh = block.subdivide(4)?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 2\n\
             npoint = 9\n\
             ncell = 4\n\
             n_boundary_point = 8\n\
             n_boundary_edge = 8\n\
             n_boundary_face = 0\n\
             \n\
             points\n\
             i:0 g:1 x:[0.0, 0.0] c:[0]\n\
             i:1 g:1 x:[1.0, 0.0] c:[0, 1]\n\
             i:2 g:1 x:[1.0, 1.0] c:[0, 1, 2, 3]\n\
             i:3 g:1 x:[0.0, 1.0] c:[0, 2]\n\
             i:4 g:1 x:[2.0, 0.0] c:[1]\n\
             i:5 g:1 x:[2.0, 1.0] c:[1, 3]\n\
             i:6 g:1 x:[1.0, 2.0] c:[2, 3]\n\
             i:7 g:1 x:[0.0, 2.0] c:[2]\n\
             i:8 g:1 x:[2.0, 2.0] c:[3]\n\
             \n\
             cells\n\
             i:0 g:1 n:2 p:[0, 1, 2, 3] e:[(0, 1), (0, 3)] f:[]\n\
             i:1 g:1 n:2 p:[1, 4, 5, 2] e:[(1, 4), (4, 5)] f:[]\n\
             i:2 g:1 n:2 p:[3, 2, 6, 7] e:[(3, 7), (6, 7)] f:[]\n\
             i:3 g:1 n:2 p:[2, 5, 8, 6] e:[(5, 8), (6, 8)] f:[]\n\
             \n\
             boundary_points\n\
             0 1 3 4 5 6 7 8 \n\
             \n\
             boundary_edges\n\
             g:1 k:(0,1) p:[0, 1] c:[0]\n\
             g:1 k:(0,3) p:[3, 0] c:[0]\n\
             g:1 k:(1,4) p:[1, 4] c:[1]\n\
             g:1 k:(3,7) p:[7, 3] c:[2]\n\
             g:1 k:(4,5) p:[4, 5] c:[1]\n\
             g:1 k:(5,8) p:[5, 8] c:[3]\n\
             g:1 k:(6,7) p:[6, 7] c:[2]\n\
             g:1 k:(6,8) p:[8, 6] c:[3]\n\
             \n\
             boundary_faces\n"
        );
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
        ]);
        let mesh = block.subdivide(8)?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 2\n\
             npoint = 21\n\
             ncell = 4\n\
             n_boundary_point = 16\n\
             n_boundary_edge = 8\n\
             n_boundary_face = 0\n\
             \n\
             points\n\
             i:0 g:1 x:[0.0, 0.0] c:[0]\n\
             i:1 g:1 x:[1.0, 0.0] c:[0, 1]\n\
             i:2 g:1 x:[1.0, 1.0] c:[0, 1, 2, 3]\n\
             i:3 g:1 x:[0.0, 1.0] c:[0, 2]\n\
             i:4 g:1 x:[0.5, 0.0] c:[0]\n\
             i:5 g:1 x:[1.0, 0.5] c:[0, 1]\n\
             i:6 g:1 x:[0.5, 1.0] c:[0, 2]\n\
             i:7 g:1 x:[0.0, 0.5] c:[0]\n\
             i:8 g:1 x:[2.0, 0.0] c:[1]\n\
             i:9 g:1 x:[2.0, 1.0] c:[1, 3]\n\
             i:10 g:1 x:[1.5, 0.0] c:[1]\n\
             i:11 g:1 x:[2.0, 0.5] c:[1]\n\
             i:12 g:1 x:[1.5, 1.0] c:[1, 3]\n\
             i:13 g:1 x:[1.0, 2.0] c:[2, 3]\n\
             i:14 g:1 x:[0.0, 2.0] c:[2]\n\
             i:15 g:1 x:[1.0, 1.5] c:[2, 3]\n\
             i:16 g:1 x:[0.5, 2.0] c:[2]\n\
             i:17 g:1 x:[0.0, 1.5] c:[2]\n\
             i:18 g:1 x:[2.0, 2.0] c:[3]\n\
             i:19 g:1 x:[2.0, 1.5] c:[3]\n\
             i:20 g:1 x:[1.5, 2.0] c:[3]\n\
             \n\
             cells\n\
             i:0 g:1 n:2 p:[0, 1, 2, 3, 4, 5, 6, 7] e:[(0, 1), (0, 3)] f:[]\n\
             i:1 g:1 n:2 p:[1, 8, 9, 2, 10, 11, 12, 5] e:[(1, 8), (8, 9)] f:[]\n\
             i:2 g:1 n:2 p:[3, 2, 13, 14, 6, 15, 16, 17] e:[(3, 14), (13, 14)] f:[]\n\
             i:3 g:1 n:2 p:[2, 9, 18, 13, 12, 19, 20, 15] e:[(9, 18), (13, 18)] f:[]\n\
             \n\
             boundary_points\n\
             0 1 3 4 7 8 9 10 11 13 14 16 17 18 19 20 \n\
             \n\
             boundary_edges\n\
             g:1 k:(0,1) p:[0, 1, 4] c:[0]\n\
             g:1 k:(0,3) p:[3, 0, 7] c:[0]\n\
             g:1 k:(1,8) p:[1, 8, 10] c:[1]\n\
             g:1 k:(3,14) p:[14, 3, 17] c:[2]\n\
             g:1 k:(8,9) p:[8, 9, 11] c:[1]\n\
             g:1 k:(9,18) p:[9, 18, 19] c:[3]\n\
             g:1 k:(13,14) p:[13, 14, 16] c:[2]\n\
             g:1 k:(13,18) p:[18, 13, 20] c:[3]\n\
             \n\
             boundary_faces\n"
        );
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
            [0.0, 0.0, 2.0],
            [2.0, 0.0, 2.0],
            [2.0, 2.0, 2.0],
            [0.0, 2.0, 2.0],
        ]);
        let mesh = block.subdivide(8)?;
        println!("{}", mesh);
        assert_eq!(
        format!("{}", mesh),
            "ndim = 3\n\
             npoint = 27\n\
             ncell = 8\n\
             n_boundary_point = 26\n\
             n_boundary_edge = 48\n\
             n_boundary_face = 24\n\
             \n\
             points\n\
             i:0 g:1 x:[0.0, 0.0, 0.0] c:[0]\n\
             i:1 g:1 x:[1.0, 0.0, 0.0] c:[0, 1]\n\
             i:2 g:1 x:[1.0, 1.0, 0.0] c:[0, 1, 2, 3]\n\
             i:3 g:1 x:[0.0, 1.0, 0.0] c:[0, 2]\n\
             i:4 g:1 x:[0.0, 0.0, 1.0] c:[0, 4]\n\
             i:5 g:1 x:[1.0, 0.0, 1.0] c:[0, 1, 4, 5]\n\
             i:6 g:1 x:[1.0, 1.0, 1.0] c:[0, 1, 2, 3, 4, 5, 6, 7]\n\
             i:7 g:1 x:[0.0, 1.0, 1.0] c:[0, 2, 4, 6]\n\
             i:8 g:1 x:[2.0, 0.0, 0.0] c:[1]\n\
             i:9 g:1 x:[2.0, 1.0, 0.0] c:[1, 3]\n\
             i:10 g:1 x:[2.0, 0.0, 1.0] c:[1, 5]\n\
             i:11 g:1 x:[2.0, 1.0, 1.0] c:[1, 3, 5, 7]\n\
             i:12 g:1 x:[1.0, 2.0, 0.0] c:[2, 3]\n\
             i:13 g:1 x:[0.0, 2.0, 0.0] c:[2]\n\
             i:14 g:1 x:[1.0, 2.0, 1.0] c:[2, 3, 6, 7]\n\
             i:15 g:1 x:[0.0, 2.0, 1.0] c:[2, 6]\n\
             i:16 g:1 x:[2.0, 2.0, 0.0] c:[3]\n\
             i:17 g:1 x:[2.0, 2.0, 1.0] c:[3, 7]\n\
             i:18 g:1 x:[0.0, 0.0, 2.0] c:[4]\n\
             i:19 g:1 x:[1.0, 0.0, 2.0] c:[4, 5]\n\
             i:20 g:1 x:[1.0, 1.0, 2.0] c:[4, 5, 6, 7]\n\
             i:21 g:1 x:[0.0, 1.0, 2.0] c:[4, 6]\n\
             i:22 g:1 x:[2.0, 0.0, 2.0] c:[5]\n\
             i:23 g:1 x:[2.0, 1.0, 2.0] c:[5, 7]\n\
             i:24 g:1 x:[1.0, 2.0, 2.0] c:[6, 7]\n\
             i:25 g:1 x:[0.0, 2.0, 2.0] c:[6]\n\
             i:26 g:1 x:[2.0, 2.0, 2.0] c:[7]\n\
             \n\
             cells\n\
             i:0 g:1 n:3 p:[0, 1, 2, 3, 4, 5, 6, 7] e:[(0, 1), (0, 3), (0, 4), (1, 2), (1, 5), (2, 3), (3, 7), (4, 5), (4, 7)] f:[(0, 1, 5), (0, 2, 3), (0, 4, 7)]\n\
             i:1 g:1 n:3 p:[1, 8, 9, 2, 5, 10, 11, 6] e:[(1, 2), (1, 5), (1, 8), (2, 9), (5, 10), (8, 9), (8, 10), (9, 11), (10, 11)] f:[(1, 2, 9), (1, 8, 10), (8, 9, 11)]\n\
             i:2 g:1 n:3 p:[3, 2, 12, 13, 7, 6, 14, 15] e:[(2, 3), (2, 12), (3, 7), (3, 13), (7, 15), (12, 13), (12, 14), (13, 15), (14, 15)] f:[(3, 7, 15), (3, 12, 13), (12, 13, 15)]\n\
             i:3 g:1 n:3 p:[2, 9, 16, 12, 6, 11, 17, 14] e:[(2, 9), (2, 12), (9, 11), (9, 16), (11, 17), (12, 14), (12, 16), (14, 17), (16, 17)] f:[(2, 12, 16), (9, 16, 17), (12, 14, 16)]\n\
             i:4 g:1 n:3 p:[4, 5, 6, 7, 18, 19, 20, 21] e:[(4, 5), (4, 7), (4, 18), (5, 19), (7, 21), (18, 19), (18, 21), (19, 20), (20, 21)] f:[(4, 5, 19), (4, 18, 21), (18, 19, 20)]\n\
             i:5 g:1 n:3 p:[5, 10, 11, 6, 19, 22, 23, 20] e:[(5, 10), (5, 19), (10, 11), (10, 22), (11, 23), (19, 20), (19, 22), (20, 23), (22, 23)] f:[(5, 10, 22), (10, 11, 23), (19, 22, 23)]\n\
             i:6 g:1 n:3 p:[7, 6, 14, 15, 21, 20, 24, 25] e:[(7, 15), (7, 21), (14, 15), (14, 24), (15, 25), (20, 21), (20, 24), (21, 25), (24, 25)] f:[(7, 21, 25), (14, 15, 25), (20, 21, 24)]\n\
             i:7 g:1 n:3 p:[6, 11, 17, 14, 20, 23, 26, 24] e:[(11, 17), (11, 23), (14, 17), (14, 24), (17, 26), (20, 23), (20, 24), (23, 26), (24, 26)] f:[(11, 17, 26), (14, 17, 24), (20, 23, 26)]\n\
             \n\
             boundary_points\n\
             0 1 2 3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 \n\
             \n\
             boundary_edges\n\
             g:1 k:(0,1) p:[0, 1] c:[0]\n\
             g:1 k:(0,3) p:[3, 0] c:[0]\n\
             g:1 k:(0,4) p:[0, 4] c:[0]\n\
             g:1 k:(1,2) p:[1, 2] c:[0, 1]\n\
             g:1 k:(1,5) p:[1, 5] c:[0, 1]\n\
             g:1 k:(1,8) p:[1, 8] c:[1]\n\
             g:1 k:(2,3) p:[2, 3] c:[0, 2]\n\
             g:1 k:(2,9) p:[9, 2] c:[1, 3]\n\
             g:1 k:(2,12) p:[2, 12] c:[2, 3]\n\
             g:1 k:(3,7) p:[3, 7] c:[0, 2]\n\
             g:1 k:(3,13) p:[13, 3] c:[2]\n\
             g:1 k:(4,5) p:[4, 5] c:[0, 4]\n\
             g:1 k:(4,7) p:[7, 4] c:[0, 4]\n\
             g:1 k:(4,18) p:[4, 18] c:[4]\n\
             g:1 k:(5,10) p:[5, 10] c:[1, 5]\n\
             g:1 k:(5,19) p:[5, 19] c:[4, 5]\n\
             g:1 k:(7,15) p:[15, 7] c:[2, 6]\n\
             g:1 k:(7,21) p:[7, 21] c:[4, 6]\n\
             g:1 k:(8,9) p:[8, 9] c:[1]\n\
             g:1 k:(8,10) p:[8, 10] c:[1]\n\
             g:1 k:(9,11) p:[9, 11] c:[1, 3]\n\
             g:1 k:(9,16) p:[9, 16] c:[3]\n\
             g:1 k:(10,11) p:[10, 11] c:[1, 5]\n\
             g:1 k:(10,22) p:[10, 22] c:[5]\n\
             g:1 k:(11,17) p:[11, 17] c:[3, 7]\n\
             g:1 k:(11,23) p:[11, 23] c:[5, 7]\n\
             g:1 k:(12,13) p:[12, 13] c:[2]\n\
             g:1 k:(12,14) p:[12, 14] c:[2, 3]\n\
             g:1 k:(12,16) p:[16, 12] c:[3]\n\
             g:1 k:(13,15) p:[13, 15] c:[2]\n\
             g:1 k:(14,15) p:[14, 15] c:[2, 6]\n\
             g:1 k:(14,17) p:[17, 14] c:[3, 7]\n\
             g:1 k:(14,24) p:[14, 24] c:[6, 7]\n\
             g:1 k:(15,25) p:[15, 25] c:[6]\n\
             g:1 k:(16,17) p:[16, 17] c:[3]\n\
             g:1 k:(17,26) p:[17, 26] c:[7]\n\
             g:1 k:(18,19) p:[18, 19] c:[4]\n\
             g:1 k:(18,21) p:[21, 18] c:[4]\n\
             g:1 k:(19,20) p:[19, 20] c:[4, 5]\n\
             g:1 k:(19,22) p:[19, 22] c:[5]\n\
             g:1 k:(20,21) p:[20, 21] c:[4, 6]\n\
             g:1 k:(20,23) p:[23, 20] c:[5, 7]\n\
             g:1 k:(20,24) p:[20, 24] c:[6, 7]\n\
             g:1 k:(21,25) p:[25, 21] c:[6]\n\
             g:1 k:(22,23) p:[22, 23] c:[5]\n\
             g:1 k:(23,26) p:[23, 26] c:[7]\n\
             g:1 k:(24,25) p:[24, 25] c:[6]\n\
             g:1 k:(24,26) p:[26, 24] c:[7]\n\
             \n\
             boundary_faces\n\
             g:1 k:(0,1,5) p:[0, 1, 5, 4] c:[0]\n\
             g:1 k:(0,2,3) p:[0, 3, 2, 1] c:[0]\n\
             g:1 k:(0,4,7) p:[0, 4, 7, 3] c:[0]\n\
             g:1 k:(1,2,9) p:[1, 2, 9, 8] c:[1]\n\
             g:1 k:(1,8,10) p:[1, 8, 10, 5] c:[1]\n\
             g:1 k:(2,12,16) p:[2, 12, 16, 9] c:[3]\n\
             g:1 k:(3,7,15) p:[3, 7, 15, 13] c:[2]\n\
             g:1 k:(3,12,13) p:[3, 13, 12, 2] c:[2]\n\
             g:1 k:(4,5,19) p:[4, 5, 19, 18] c:[4]\n\
             g:1 k:(4,18,21) p:[4, 18, 21, 7] c:[4]\n\
             g:1 k:(5,10,22) p:[5, 10, 22, 19] c:[5]\n\
             g:1 k:(7,21,25) p:[7, 21, 25, 15] c:[6]\n\
             g:1 k:(8,9,11) p:[8, 9, 11, 10] c:[1]\n\
             g:1 k:(9,16,17) p:[9, 16, 17, 11] c:[3]\n\
             g:1 k:(10,11,23) p:[10, 11, 23, 22] c:[5]\n\
             g:1 k:(11,17,26) p:[11, 17, 26, 23] c:[7]\n\
             g:1 k:(12,13,15) p:[12, 13, 15, 14] c:[2]\n\
             g:1 k:(12,14,16) p:[16, 12, 14, 17] c:[3]\n\
             g:1 k:(14,15,25) p:[14, 15, 25, 24] c:[6]\n\
             g:1 k:(14,17,24) p:[17, 14, 24, 26] c:[7]\n\
             g:1 k:(18,19,20) p:[18, 19, 20, 21] c:[4]\n\
             g:1 k:(19,22,23) p:[19, 22, 23, 20] c:[5]\n\
             g:1 k:(20,21,24) p:[21, 20, 24, 25] c:[6]\n\
             g:1 k:(20,23,26) p:[20, 23, 26, 24] c:[7]\n"
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
            [0.0, 0.0, 2.0],
            [2.0, 0.0, 2.0],
            [2.0, 2.0, 2.0],
            [0.0, 2.0, 2.0],
        ]);
        let mesh = block.subdivide(20)?;
        println!("{}", mesh);
        assert_eq!(
        format!("{}", mesh),
            "ndim = 3\n\
             npoint = 81\n\
             ncell = 8\n\
             n_boundary_point = 74\n\
             n_boundary_edge = 48\n\
             n_boundary_face = 24\n\
             \n\
             points\n\
             i:0 g:1 x:[0.0, 0.0, 0.0] c:[0]\n\
             i:1 g:1 x:[1.0, 0.0, 0.0] c:[0, 1]\n\
             i:2 g:1 x:[1.0, 1.0, 0.0] c:[0, 1, 2, 3]\n\
             i:3 g:1 x:[0.0, 1.0, 0.0] c:[0, 2]\n\
             i:4 g:1 x:[0.0, 0.0, 1.0] c:[0, 4]\n\
             i:5 g:1 x:[1.0, 0.0, 1.0] c:[0, 1, 4, 5]\n\
             i:6 g:1 x:[1.0, 1.0, 1.0] c:[0, 1, 2, 3, 4, 5, 6, 7]\n\
             i:7 g:1 x:[0.0, 1.0, 1.0] c:[0, 2, 4, 6]\n\
             i:8 g:1 x:[0.5, 0.0, 0.0] c:[0]\n\
             i:9 g:1 x:[1.0, 0.5, 0.0] c:[0, 1]\n\
             i:10 g:1 x:[0.5, 1.0, 0.0] c:[0, 2]\n\
             i:11 g:1 x:[0.0, 0.5, 0.0] c:[0]\n\
             i:12 g:1 x:[0.5, 0.0, 1.0] c:[0, 4]\n\
             i:13 g:1 x:[1.0, 0.5, 1.0] c:[0, 1, 4, 5]\n\
             i:14 g:1 x:[0.5, 1.0, 1.0] c:[0, 2, 4, 6]\n\
             i:15 g:1 x:[0.0, 0.5, 1.0] c:[0, 4]\n\
             i:16 g:1 x:[0.0, 0.0, 0.5] c:[0]\n\
             i:17 g:1 x:[1.0, 0.0, 0.5] c:[0, 1]\n\
             i:18 g:1 x:[1.0, 1.0, 0.5] c:[0, 1, 2, 3]\n\
             i:19 g:1 x:[0.0, 1.0, 0.5] c:[0, 2]\n\
             i:20 g:1 x:[2.0, 0.0, 0.0] c:[1]\n\
             i:21 g:1 x:[2.0, 1.0, 0.0] c:[1, 3]\n\
             i:22 g:1 x:[2.0, 0.0, 1.0] c:[1, 5]\n\
             i:23 g:1 x:[2.0, 1.0, 1.0] c:[1, 3, 5, 7]\n\
             i:24 g:1 x:[1.5, 0.0, 0.0] c:[1]\n\
             i:25 g:1 x:[2.0, 0.5, 0.0] c:[1]\n\
             i:26 g:1 x:[1.5, 1.0, 0.0] c:[1, 3]\n\
             i:27 g:1 x:[1.5, 0.0, 1.0] c:[1, 5]\n\
             i:28 g:1 x:[2.0, 0.5, 1.0] c:[1, 5]\n\
             i:29 g:1 x:[1.5, 1.0, 1.0] c:[1, 3, 5, 7]\n\
             i:30 g:1 x:[2.0, 0.0, 0.5] c:[1]\n\
             i:31 g:1 x:[2.0, 1.0, 0.5] c:[1, 3]\n\
             i:32 g:1 x:[1.0, 2.0, 0.0] c:[2, 3]\n\
             i:33 g:1 x:[0.0, 2.0, 0.0] c:[2]\n\
             i:34 g:1 x:[1.0, 2.0, 1.0] c:[2, 3, 6, 7]\n\
             i:35 g:1 x:[0.0, 2.0, 1.0] c:[2, 6]\n\
             i:36 g:1 x:[1.0, 1.5, 0.0] c:[2, 3]\n\
             i:37 g:1 x:[0.5, 2.0, 0.0] c:[2]\n\
             i:38 g:1 x:[0.0, 1.5, 0.0] c:[2]\n\
             i:39 g:1 x:[1.0, 1.5, 1.0] c:[2, 3, 6, 7]\n\
             i:40 g:1 x:[0.5, 2.0, 1.0] c:[2, 6]\n\
             i:41 g:1 x:[0.0, 1.5, 1.0] c:[2, 6]\n\
             i:42 g:1 x:[1.0, 2.0, 0.5] c:[2, 3]\n\
             i:43 g:1 x:[0.0, 2.0, 0.5] c:[2]\n\
             i:44 g:1 x:[2.0, 2.0, 0.0] c:[3]\n\
             i:45 g:1 x:[2.0, 2.0, 1.0] c:[3, 7]\n\
             i:46 g:1 x:[2.0, 1.5, 0.0] c:[3]\n\
             i:47 g:1 x:[1.5, 2.0, 0.0] c:[3]\n\
             i:48 g:1 x:[2.0, 1.5, 1.0] c:[3, 7]\n\
             i:49 g:1 x:[1.5, 2.0, 1.0] c:[3, 7]\n\
             i:50 g:1 x:[2.0, 2.0, 0.5] c:[3]\n\
             i:51 g:1 x:[0.0, 0.0, 2.0] c:[4]\n\
             i:52 g:1 x:[1.0, 0.0, 2.0] c:[4, 5]\n\
             i:53 g:1 x:[1.0, 1.0, 2.0] c:[4, 5, 6, 7]\n\
             i:54 g:1 x:[0.0, 1.0, 2.0] c:[4, 6]\n\
             i:55 g:1 x:[0.5, 0.0, 2.0] c:[4]\n\
             i:56 g:1 x:[1.0, 0.5, 2.0] c:[4, 5]\n\
             i:57 g:1 x:[0.5, 1.0, 2.0] c:[4, 6]\n\
             i:58 g:1 x:[0.0, 0.5, 2.0] c:[4]\n\
             i:59 g:1 x:[0.0, 0.0, 1.5] c:[4]\n\
             i:60 g:1 x:[1.0, 0.0, 1.5] c:[4, 5]\n\
             i:61 g:1 x:[1.0, 1.0, 1.5] c:[4, 5, 6, 7]\n\
             i:62 g:1 x:[0.0, 1.0, 1.5] c:[4, 6]\n\
             i:63 g:1 x:[2.0, 0.0, 2.0] c:[5]\n\
             i:64 g:1 x:[2.0, 1.0, 2.0] c:[5, 7]\n\
             i:65 g:1 x:[1.5, 0.0, 2.0] c:[5]\n\
             i:66 g:1 x:[2.0, 0.5, 2.0] c:[5]\n\
             i:67 g:1 x:[1.5, 1.0, 2.0] c:[5, 7]\n\
             i:68 g:1 x:[2.0, 0.0, 1.5] c:[5]\n\
             i:69 g:1 x:[2.0, 1.0, 1.5] c:[5, 7]\n\
             i:70 g:1 x:[1.0, 2.0, 2.0] c:[6, 7]\n\
             i:71 g:1 x:[0.0, 2.0, 2.0] c:[6]\n\
             i:72 g:1 x:[1.0, 1.5, 2.0] c:[6, 7]\n\
             i:73 g:1 x:[0.5, 2.0, 2.0] c:[6]\n\
             i:74 g:1 x:[0.0, 1.5, 2.0] c:[6]\n\
             i:75 g:1 x:[1.0, 2.0, 1.5] c:[6, 7]\n\
             i:76 g:1 x:[0.0, 2.0, 1.5] c:[6]\n\
             i:77 g:1 x:[2.0, 2.0, 2.0] c:[7]\n\
             i:78 g:1 x:[2.0, 1.5, 2.0] c:[7]\n\
             i:79 g:1 x:[1.5, 2.0, 2.0] c:[7]\n\
             i:80 g:1 x:[2.0, 2.0, 1.5] c:[7]\n\
             \n\
             cells\n\
             i:0 g:1 n:3 p:[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19] e:[(0, 1), (0, 3), (0, 4), (1, 2), (1, 5), (2, 3), (3, 7), (4, 5), (4, 7)] f:[(0, 1, 5), (0, 2, 3), (0, 4, 7)]\n\
             i:1 g:1 n:3 p:[1, 20, 21, 2, 5, 22, 23, 6, 24, 25, 26, 9, 27, 28, 29, 13, 17, 30, 31, 18] e:[(1, 2), (1, 5), (1, 20), (2, 21), (5, 22), (20, 21), (20, 22), (21, 23), (22, 23)] f:[(1, 2, 21), (1, 20, 22), (20, 21, 23)]\n\
             i:2 g:1 n:3 p:[3, 2, 32, 33, 7, 6, 34, 35, 10, 36, 37, 38, 14, 39, 40, 41, 19, 18, 42, 43] e:[(2, 3), (2, 32), (3, 7), (3, 33), (7, 35), (32, 33), (32, 34), (33, 35), (34, 35)] f:[(3, 7, 35), (3, 32, 33), (32, 33, 35)]\n\
             i:3 g:1 n:3 p:[2, 21, 44, 32, 6, 23, 45, 34, 26, 46, 47, 36, 29, 48, 49, 39, 18, 31, 50, 42] e:[(2, 21), (2, 32), (21, 23), (21, 44), (23, 45), (32, 34), (32, 44), (34, 45), (44, 45)] f:[(2, 32, 44), (21, 44, 45), (32, 34, 44)]\n\
             i:4 g:1 n:3 p:[4, 5, 6, 7, 51, 52, 53, 54, 12, 13, 14, 15, 55, 56, 57, 58, 59, 60, 61, 62] e:[(4, 5), (4, 7), (4, 51), (5, 52), (7, 54), (51, 52), (51, 54), (52, 53), (53, 54)] f:[(4, 5, 52), (4, 51, 54), (51, 52, 53)]\n\
             i:5 g:1 n:3 p:[5, 22, 23, 6, 52, 63, 64, 53, 27, 28, 29, 13, 65, 66, 67, 56, 60, 68, 69, 61] e:[(5, 22), (5, 52), (22, 23), (22, 63), (23, 64), (52, 53), (52, 63), (53, 64), (63, 64)] f:[(5, 22, 63), (22, 23, 64), (52, 63, 64)]\n\
             i:6 g:1 n:3 p:[7, 6, 34, 35, 54, 53, 70, 71, 14, 39, 40, 41, 57, 72, 73, 74, 62, 61, 75, 76] e:[(7, 35), (7, 54), (34, 35), (34, 70), (35, 71), (53, 54), (53, 70), (54, 71), (70, 71)] f:[(7, 54, 71), (34, 35, 71), (53, 54, 70)]\n\
             i:7 g:1 n:3 p:[6, 23, 45, 34, 53, 64, 77, 70, 29, 48, 49, 39, 67, 78, 79, 72, 61, 69, 80, 75] e:[(23, 45), (23, 64), (34, 45), (34, 70), (45, 77), (53, 64), (53, 70), (64, 77), (70, 77)] f:[(23, 45, 77), (34, 45, 70), (53, 64, 77)]\n\
             \n\
             boundary_points\n\
             0 1 2 3 4 5 7 8 9 10 11 12 15 16 17 19 20 21 22 23 24 25 26 27 28 30 31 32 33 34 35 36 37 38 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 \n\
             \n\
             boundary_edges\n\
             g:1 k:(0,1) p:[0, 1, 8] c:[0]\n\
             g:1 k:(0,3) p:[3, 0, 11] c:[0]\n\
             g:1 k:(0,4) p:[0, 4, 16] c:[0]\n\
             g:1 k:(1,2) p:[1, 2, 9] c:[0, 1]\n\
             g:1 k:(1,5) p:[1, 5, 17] c:[0, 1]\n\
             g:1 k:(1,20) p:[1, 20, 24] c:[1]\n\
             g:1 k:(2,3) p:[2, 3, 10] c:[0, 2]\n\
             g:1 k:(2,21) p:[21, 2, 26] c:[1, 3]\n\
             g:1 k:(2,32) p:[2, 32, 36] c:[2, 3]\n\
             g:1 k:(3,7) p:[3, 7, 19] c:[0, 2]\n\
             g:1 k:(3,33) p:[33, 3, 38] c:[2]\n\
             g:1 k:(4,5) p:[4, 5, 12] c:[0, 4]\n\
             g:1 k:(4,7) p:[7, 4, 15] c:[0, 4]\n\
             g:1 k:(4,51) p:[4, 51, 59] c:[4]\n\
             g:1 k:(5,22) p:[5, 22, 27] c:[1, 5]\n\
             g:1 k:(5,52) p:[5, 52, 60] c:[4, 5]\n\
             g:1 k:(7,35) p:[35, 7, 41] c:[2, 6]\n\
             g:1 k:(7,54) p:[7, 54, 62] c:[4, 6]\n\
             g:1 k:(20,21) p:[20, 21, 25] c:[1]\n\
             g:1 k:(20,22) p:[20, 22, 30] c:[1]\n\
             g:1 k:(21,23) p:[21, 23, 31] c:[1, 3]\n\
             g:1 k:(21,44) p:[21, 44, 46] c:[3]\n\
             g:1 k:(22,23) p:[22, 23, 28] c:[1, 5]\n\
             g:1 k:(22,63) p:[22, 63, 68] c:[5]\n\
             g:1 k:(23,45) p:[23, 45, 48] c:[3, 7]\n\
             g:1 k:(23,64) p:[23, 64, 69] c:[5, 7]\n\
             g:1 k:(32,33) p:[32, 33, 37] c:[2]\n\
             g:1 k:(32,34) p:[32, 34, 42] c:[2, 3]\n\
             g:1 k:(32,44) p:[44, 32, 47] c:[3]\n\
             g:1 k:(33,35) p:[33, 35, 43] c:[2]\n\
             g:1 k:(34,35) p:[34, 35, 40] c:[2, 6]\n\
             g:1 k:(34,45) p:[45, 34, 49] c:[3, 7]\n\
             g:1 k:(34,70) p:[34, 70, 75] c:[6, 7]\n\
             g:1 k:(35,71) p:[35, 71, 76] c:[6]\n\
             g:1 k:(44,45) p:[44, 45, 50] c:[3]\n\
             g:1 k:(45,77) p:[45, 77, 80] c:[7]\n\
             g:1 k:(51,52) p:[51, 52, 55] c:[4]\n\
             g:1 k:(51,54) p:[54, 51, 58] c:[4]\n\
             g:1 k:(52,53) p:[52, 53, 56] c:[4, 5]\n\
             g:1 k:(52,63) p:[52, 63, 65] c:[5]\n\
             g:1 k:(53,54) p:[53, 54, 57] c:[4, 6]\n\
             g:1 k:(53,64) p:[64, 53, 67] c:[5, 7]\n\
             g:1 k:(53,70) p:[53, 70, 72] c:[6, 7]\n\
             g:1 k:(54,71) p:[71, 54, 74] c:[6]\n\
             g:1 k:(63,64) p:[63, 64, 66] c:[5]\n\
             g:1 k:(64,77) p:[64, 77, 78] c:[7]\n\
             g:1 k:(70,71) p:[70, 71, 73] c:[6]\n\
             g:1 k:(70,77) p:[77, 70, 79] c:[7]\n\
             \n\
             boundary_faces\n\
             g:1 k:(0,1,5) p:[0, 1, 5, 4, 8, 17, 12, 16] c:[0]\n\
             g:1 k:(0,2,3) p:[0, 3, 2, 1, 11, 10, 9, 8] c:[0]\n\
             g:1 k:(0,4,7) p:[0, 4, 7, 3, 16, 15, 19, 11] c:[0]\n\
             g:1 k:(1,2,21) p:[1, 2, 21, 20, 9, 26, 25, 24] c:[1]\n\
             g:1 k:(1,20,22) p:[1, 20, 22, 5, 24, 30, 27, 17] c:[1]\n\
             g:1 k:(2,32,44) p:[2, 32, 44, 21, 36, 47, 46, 26] c:[3]\n\
             g:1 k:(3,7,35) p:[3, 7, 35, 33, 19, 41, 43, 38] c:[2]\n\
             g:1 k:(3,32,33) p:[3, 33, 32, 2, 38, 37, 36, 10] c:[2]\n\
             g:1 k:(4,5,52) p:[4, 5, 52, 51, 12, 60, 55, 59] c:[4]\n\
             g:1 k:(4,51,54) p:[4, 51, 54, 7, 59, 58, 62, 15] c:[4]\n\
             g:1 k:(5,22,63) p:[5, 22, 63, 52, 27, 68, 65, 60] c:[5]\n\
             g:1 k:(7,54,71) p:[7, 54, 71, 35, 62, 74, 76, 41] c:[6]\n\
             g:1 k:(20,21,23) p:[20, 21, 23, 22, 25, 31, 28, 30] c:[1]\n\
             g:1 k:(21,44,45) p:[21, 44, 45, 23, 46, 50, 48, 31] c:[3]\n\
             g:1 k:(22,23,64) p:[22, 23, 64, 63, 28, 69, 66, 68] c:[5]\n\
             g:1 k:(23,45,77) p:[23, 45, 77, 64, 48, 80, 78, 69] c:[7]\n\
             g:1 k:(32,33,35) p:[32, 33, 35, 34, 37, 43, 40, 42] c:[2]\n\
             g:1 k:(32,34,44) p:[44, 32, 34, 45, 47, 42, 49, 50] c:[3]\n\
             g:1 k:(34,35,71) p:[34, 35, 71, 70, 40, 76, 73, 75] c:[6]\n\
             g:1 k:(34,45,70) p:[45, 34, 70, 77, 49, 75, 79, 80] c:[7]\n\
             g:1 k:(51,52,53) p:[51, 52, 53, 54, 55, 56, 57, 58] c:[4]\n\
             g:1 k:(52,63,64) p:[52, 63, 64, 53, 65, 66, 67, 56] c:[5]\n\
             g:1 k:(53,54,70) p:[54, 53, 70, 71, 57, 72, 73, 74] c:[6]\n\
             g:1 k:(53,64,77) p:[53, 64, 77, 70, 67, 78, 79, 72] c:[7]\n"
        );
        Ok(())
    }
}
