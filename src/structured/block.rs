use crate::{
    new_shape_qua_or_hex, AsArray2D, Cell, Circle, Edge, Face, GridSearch, KeyEdge, KeyFace, KindQuaOrHex, Mesh, Point,
    Shape,
};
use russell_lab::{mat_vec_mul, Matrix, Vector};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Constraint {
    /// Arc
    Arc(Circle),

    /// Arc surface extruded along X
    ArcX(Circle),
}

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
    group: usize,             // group of all elements in this block
    ndim: usize,              // space dimension
    npoint: usize,            // number of points (8 or 20) (quadratic block internally)
    nedge: usize,             // number of edges (4 or 12)
    nface: usize,             // number of faces (4 or 6)
    coords: Matrix,           // coordinates of points (npoint, ndim)
    point_groups: Vec<usize>, // point groups (npoint)
    edge_groups: Vec<usize>,  // edge groups (nedge)
    face_groups: Vec<usize>,  // face groups (nface)
    ndiv: Vec<usize>,         // number of divisions along each dim (ndim)
    delta_ksi: Vec<Vec<f64>>, // delta ksi along each dim (ndim, {ndiv[0],ndiv[1],ndiv[2]})

    edge_constraints: Vec<Option<Constraint>>, // constraints (nedge)
    face_constraints: Vec<Option<Constraint>>, // constraints (nface)

    // shape and interpolation functions
    shape: Box<dyn Shape>,

    // grid to search natural coordinates
    grid_ksi: GridSearch,
}

impl Block {
    // constants
    const NAT_MIN: f64 = -1.0; // min natural coordinate value
    const NAT_MAX: f64 = 1.0; // max natural coordinate value
    const NAT_LENGTH: f64 = 2.0; // length of shape along each direction in natural coords space
    const NAT_TOLERANCE: f64 = 1e-4; // tolerance to compare coordinates in the natural space

    /// Creates a new Block with default options
    pub fn new(ndim: usize) -> Self {
        let shape = if ndim == 2 {
            new_shape_qua_or_hex(KindQuaOrHex::Qua8)
        } else {
            new_shape_qua_or_hex(KindQuaOrHex::Hex20)
        };
        let (npoint, nedge, nface) = (shape.get_npoint(), shape.get_nedge(), shape.get_nface());

        const NDIV: usize = 2;
        const GRID_NDIV: usize = 20;
        const GRID_MIN: f64 = -1.0;
        const GRID_MAX: f64 = 1.0;

        Block {
            group: 1,
            ndim,
            npoint,
            nedge,
            nface,
            coords: Matrix::new(npoint, ndim),
            point_groups: vec![0; npoint],
            edge_groups: vec![0; nedge],
            face_groups: vec![0; nface],
            ndiv: vec![NDIV; ndim],
            delta_ksi: vec![vec![1.0; NDIV]; ndim],
            edge_constraints: vec![None; nedge],
            face_constraints: vec![None; nface],
            shape,
            grid_ksi: GridSearch::new(
                &vec![GRID_NDIV; ndim],
                &vec![GRID_MIN; ndim],
                &vec![GRID_MAX; ndim],
                &vec![Block::NAT_TOLERANCE; ndim],
            )
            .unwrap(),
        }
    }

    /// Sets group
    pub fn set_group(&mut self, group: usize) -> &mut Self {
        self.group = group;
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
        assert_eq!(ncol, self.ndim);
        if self.ndim == 2 {
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
        if self.ndim == 2 && nrow == 4 {
            for i in 0..self.ndim {
                self.coords[4][i] = (self.coords[0][i] + self.coords[1][i]) / 2.0;
                self.coords[5][i] = (self.coords[1][i] + self.coords[2][i]) / 2.0;
                self.coords[6][i] = (self.coords[2][i] + self.coords[3][i]) / 2.0;
                self.coords[7][i] = (self.coords[3][i] + self.coords[0][i]) / 2.0;
            }
        }
        if self.ndim == 3 && nrow == 8 {
            for i in 0..self.ndim {
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

    /// Sets the group of a point
    ///
    /// # Input
    ///
    /// * `p` -- index of point in [0, npoint-1]
    pub fn set_point_group(&mut self, p: usize, group: usize) -> &mut Self {
        assert!(p < self.npoint);
        self.point_groups[p] = group;
        self
    }

    /// Sets the group of an edge
    ///
    /// # Input
    ///
    /// * `e` -- index of edge in [0, nedge-1]
    pub fn set_edge_group(&mut self, e: usize, group: usize) -> &mut Self {
        assert!(e < self.nedge);
        self.edge_groups[e] = group;
        self
    }

    /// Sets the group of a face
    ///
    /// # Input
    ///
    /// * `f` -- index of edge in [0, nface-1]
    pub fn set_face_group(&mut self, f: usize, group: usize) -> &mut Self {
        assert!(f < self.nface);
        self.face_groups[f] = group;
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
    /// where `L=2` is the edge-length (in natural coordinates) and `wᵐ` are
    /// the weights for each division `m`.
    pub fn set_ndiv(&mut self, ndiv: &[usize]) -> &mut Self {
        assert_eq!(ndiv.len(), self.ndim);
        for i in 0..self.ndim {
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
    pub fn subdivide(&mut self, output: KindQuaOrHex) -> Result<Mesh, &'static str> {
        // results
        let mut mesh = Mesh::new(self.ndim);

        // auxiliary variables
        let shape_out = new_shape_qua_or_hex(output);
        let npoint_out = shape_out.get_npoint();

        // transformation matrix: scale and translate natural space
        //   _                                       _
        //  |  scale_x   0.0    0.0    translation_x  |
        //  |    0.0   scale_y  0.0    translation_y  |
        //  |    0.0     0.0   scale_z translation_z  |
        //  |_   0.0     0.0    0.0         1.0      _|
        let mut transform = Matrix::identity(self.ndim + 1);

        // augmented natural coordinates [r,s,1] or [r,s,t,1]
        let mut ksi_aug = Vector::new(self.ndim + 1);
        ksi_aug[self.ndim] = 1.0;

        // augmented transformed nat-coordinates
        let mut ksi = Vector::new(self.ndim + 1);

        // center of shape in nat-coords
        let mut center = vec![0.0; self.ndim];

        // real point coordinates
        let mut x = Vector::new(self.ndim);

        // number of divisions along each direction
        let (nx, ny, nz) = (
            self.ndiv[0],
            self.ndiv[1],
            if self.ndim == 2 { 1 } else { self.ndiv[2] },
        );

        // for each z-division
        if self.ndim == 3 {
            center[2] = -1.0 + self.delta_ksi[2][0] / 2.0;
        }
        for k in 0..nz {
            if self.ndim == 3 {
                transform[2][2] = self.delta_ksi[2][k] / Block::NAT_LENGTH; // scale
                transform[2][self.ndim] = center[2]; // translation
            }

            // for each y-division
            center[1] = -1.0 + self.delta_ksi[1][0] / 2.0;
            for j in 0..ny {
                transform[1][1] = self.delta_ksi[1][j] / Block::NAT_LENGTH; // scale
                transform[1][self.ndim] = center[1]; // translation

                // for each x-division
                center[0] = -1.0 + self.delta_ksi[0][0] / 2.0;
                for i in 0..nx {
                    transform[0][0] = self.delta_ksi[0][i] / Block::NAT_LENGTH; // scale
                    transform[0][self.ndim] = center[0]; // translation

                    // new cell id
                    let cell_id = mesh.cells.len();

                    // for each point
                    let mut point_ids = vec![0; npoint_out];
                    for m in 0..npoint_out {
                        // transform natural coords: scale and translate
                        shape_out.get_ksi(&mut ksi_aug, m);
                        mat_vec_mul(&mut ksi, 1.0, &transform, &ksi_aug)?;

                        // maybe append point to mesh
                        let index = self.maybe_append_new_point(&mut mesh, &mut x, &ksi, cell_id)?;
                        point_ids[m] = index;

                        // mark boundary point
                        if self.is_boundary_point(&ksi) {
                            mesh.boundary_points.insert(index, true);
                        }
                    }

                    // new edge
                    let boundary_edge_ids = self.maybe_append_new_edge(&mut mesh, &shape_out, &point_ids, cell_id)?;

                    // new cell
                    let cell = Cell {
                        id: cell_id,
                        group: self.group,
                        point_ids,
                        boundary_edge_ids,
                        boundary_face_ids: Vec::new(),
                    };
                    mesh.cells.push(cell);

                    // next x-center
                    center[0] += self.delta_ksi[0][i];
                }

                // next y-center
                center[1] += self.delta_ksi[1][j];
            }

            // next z-center
            if self.ndim == 3 {
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
    /// * `ksi_vec` -- (augmented) natural coordinates of the point (scaled and translated already)
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
    ) -> Result<usize, &'static str> {
        // handle existent point
        let ksi = &ksi_vec.as_data()[0..self.ndim];
        if let Some(index) = self.grid_ksi.find(ksi)? {
            mesh.points[index].shared_by_cell_ids.push(cell_id);
            return Ok(index);
        }

        // insert point in grid search object
        let index = mesh.points.len();
        self.grid_ksi.insert(index, ksi)?;

        // compute real coords
        self.shape.calc_interp(ksi_vec);
        self.shape.mul_interp_by_matrix(x, &self.coords)?;

        // add new point to mesh
        mesh.points.push(Point {
            id: index,
            group: self.group,
            coords: x.as_data().clone(),
            shared_by_cell_ids: vec![cell_id],
        });

        // done
        Ok(index)
    }

    /// Appends new boundary edge, if not existent yet
    ///
    /// # Output
    ///
    /// * `mesh` -- updated mesh with new boundary edge
    ///
    /// # Input
    ///
    /// * `shape_out` -- shape object for the output cells
    /// * `point_ids` -- (just added) cell point ids
    /// * `cell_id` -- ID of the cell adding this point (to set shared-by information)
    ///
    /// # Returns
    ///
    /// Returns the IDs of the new or existent edge on boundary.
    fn maybe_append_new_edge(
        &mut self,
        mesh: &mut Mesh,
        shape_out: &Box<dyn Shape>,
        point_ids: &Vec<usize>,
        cell_id: usize,
    ) -> Result<Vec<usize>, &'static str> {
        // results
        let mut boundary_edge_ids = Vec::new();

        // loop over each edge
        for e in 0..shape_out.get_nedge() {
            // point ids
            let mut boundary_point_ids = Vec::new();
            for i in 0..shape_out.get_edge_npoint() {
                let local_point_id = self.shape.get_edge_local_point_id(e, i);
                let point_id = point_ids[local_point_id];
                if mesh.is_point_on_boundary(point_id) {
                    boundary_point_ids.push(point_id);
                }
                if boundary_point_ids.len() == 2 {
                    break;
                }
            }

            // skip if edge is not on boundary
            if boundary_point_ids.len() < 2 {
                continue;
            }

            // define key
            boundary_point_ids.sort(); // important
            let key = KeyEdge {
                a: boundary_point_ids[0],
                b: boundary_point_ids[1],
            };

            // existing item
            if mesh.boundary_edges.contains_key(&key) {
                let edge = mesh.boundary_edges.get_mut(&key).unwrap();
                edge.shared_by_cell_ids.push(cell_id);
                boundary_edge_ids.push(edge.id);

            // new item
            } else {
                let new_edge_id = mesh.boundary_edges.len();
                boundary_edge_ids.push(new_edge_id);
                mesh.boundary_edges.insert(
                    key,
                    Edge {
                        id: new_edge_id,
                        group: self.group,
                        point_ids: boundary_point_ids,
                        shared_by_cell_ids: vec![cell_id],
                    },
                );
            }
        }

        // done
        Ok(boundary_edge_ids)
    }

    fn is_boundary_point(&self, ksi: &Vector) -> bool {
        for i in 0..self.ndim {
            if f64::abs(ksi[i] - Block::NAT_MIN) <= Block::NAT_TOLERANCE {
                return true;
            }
            if f64::abs(ksi[i] - Block::NAT_MAX) <= Block::NAT_TOLERANCE {
                return true;
            }
        }
        false
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constraint_traits_work() {
        let constraint = Constraint::Arc(Circle {
            center: [2.0, 3.0],
            radius: 1.0,
            tolerance: 1e-2,
        });
        let clone = constraint.clone();
        assert_eq!(constraint, clone);
        assert_eq!(
            format!("{:?}", constraint),
            "Arc(Circle { center: [2.0, 3.0], radius: 1.0, tolerance: 0.01 })"
        );
    }

    #[test]
    fn new_works() {
        let b2d = Block::new(2);
        assert_eq!(b2d.group, 1);
        assert_eq!(b2d.ndim, 2);
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
        assert_eq!(b2d.point_groups, &[0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(b2d.edge_groups, &[0, 0, 0, 0]);
        assert_eq!(b2d.face_groups.len(), 0);
        assert_eq!(b2d.ndiv, &[2, 2]);
        assert_eq!(format!("{:?}", b2d.delta_ksi), "[[1.0, 1.0], [1.0, 1.0]]");
        assert_eq!(b2d.shape.get_npoint(), 8);

        let b3d = Block::new(3);
        assert_eq!(b3d.group, 1);
        assert_eq!(b3d.ndim, 3);
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
        assert_eq!(
            b3d.point_groups,
            &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        );
        assert_eq!(b3d.edge_groups, &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(b3d.face_groups, &[0, 0, 0, 0, 0, 0]);
        assert_eq!(b3d.ndiv, &[2, 2, 2]);
        assert_eq!(format!("{:?}", b3d.delta_ksi), "[[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]");
        assert_eq!(b3d.shape.get_npoint(), 20);
    }

    #[test]
    fn set_group_works() {
        let mut block = Block::new(2);
        block.set_group(2);
        assert_eq!(block.group, 2);
    }

    #[test]
    fn set_coords_works() {
        let mut b2d = Block::new(2);
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

        let mut b3d = Block::new(3);
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
    }

    #[test]
    fn set_point_group_works() {
        let mut block = Block::new(2);
        block.set_point_group(0, 111);
        assert_eq!(block.point_groups, &[111, 0, 0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn set_edge_group_works() {
        let mut block = Block::new(2);
        block.set_edge_group(0, 111);
        assert_eq!(block.edge_groups, &[111, 0, 0, 0]);
    }

    #[test]
    fn set_face_group_works() {
        let mut block = Block::new(3);
        block.set_face_group(0, 111);
        assert_eq!(block.face_groups, &[111, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn set_ndiv_works() {
        let mut block = Block::new(2);
        block.set_ndiv(&[2, 4]);
        assert_eq!(block.ndiv, &[2, 4]);
        assert_eq!(format!("{:?}", block.delta_ksi), "[[1.0, 1.0], [0.5, 0.5, 0.5, 0.5]]");
    }

    #[test]
    fn set_edge_constraint_works() {
        let mut block = Block::new(2);
        let constraint = Constraint::Arc(Circle {
            center: [-1.0, -1.0],
            radius: 2.0,
            tolerance: 1e-3,
        });
        block.set_edge_constraint(0, constraint);
        assert_eq!(block.edge_constraints[0], Some(constraint));
    }

    #[test]
    fn set_face_constraint_works() {
        let mut block = Block::new(3);
        let constraint = Constraint::ArcX(Circle {
            center: [-1.0, -1.0],
            radius: 2.0,
            tolerance: 1e-3,
        });
        block.set_face_constraint(0, constraint);
        assert_eq!(block.face_constraints[0], Some(constraint));
    }

    #[test]
    fn subdivide_2d_works() -> Result<(), &'static str> {
        let mut block = Block::new(2);
        #[rustfmt::skip]
        block.set_coords(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ]);
        let mesh = block.subdivide(KindQuaOrHex::Qua4)?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 2\n\
             npoint = 9\n\
             ncell = 4\n\
             n_boundary_point = 8\n\
             n_boundary_edge = 0\n\
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
             i:0 g:1 p:[0, 1, 2, 3] e:[] f:[]\n\
             i:1 g:1 p:[1, 4, 5, 2] e:[] f:[]\n\
             i:2 g:1 p:[3, 2, 6, 7] e:[] f:[]\n\
             i:3 g:1 p:[2, 5, 8, 6] e:[] f:[]\n\
             \n\
             boundary_points\n\
             0 1 3 4 5 6 7 8 \n\
             \n\
             boundary_edges\n\
             \n\
             boundary_faces\n"
        );
        Ok(())
    }

    #[test]
    fn subdivide_2d_o2_works() -> Result<(), &'static str> {
        let mut block = Block::new(2);
        #[rustfmt::skip]
        block.set_coords(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ]);
        let mesh = block.subdivide(KindQuaOrHex::Qua8)?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 2\n\
             npoint = 21\n\
             ncell = 4\n\
             n_boundary_point = 16\n\
             n_boundary_edge = 0\n\
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
             i:0 g:1 p:[0, 1, 2, 3, 4, 5, 6, 7] e:[] f:[]\n\
             i:1 g:1 p:[1, 8, 9, 2, 10, 11, 12, 5] e:[] f:[]\n\
             i:2 g:1 p:[3, 2, 13, 14, 6, 15, 16, 17] e:[] f:[]\n\
             i:3 g:1 p:[2, 9, 18, 13, 12, 19, 20, 15] e:[] f:[]\n\
             \n\
             boundary_points\n\
             0 1 3 4 7 8 9 10 11 13 14 16 17 18 19 20 \n\
             \n\
             boundary_edges\n\
             \n\
             boundary_faces\n"
        );
        Ok(())
    }

    #[test]
    fn subdivide_3d_works() -> Result<(), &'static str> {
        let mut block = Block::new(3);
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
        let mesh = block.subdivide(KindQuaOrHex::Hex8)?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 3\n\
             npoint = 27\n\
             ncell = 8\n\
             n_boundary_point = 26\n\
             n_boundary_edge = 0\n\
             n_boundary_face = 0\n\
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
             i:0 g:1 p:[0, 1, 2, 3, 4, 5, 6, 7] e:[] f:[]\n\
             i:1 g:1 p:[1, 8, 9, 2, 5, 10, 11, 6] e:[] f:[]\n\
             i:2 g:1 p:[3, 2, 12, 13, 7, 6, 14, 15] e:[] f:[]\n\
             i:3 g:1 p:[2, 9, 16, 12, 6, 11, 17, 14] e:[] f:[]\n\
             i:4 g:1 p:[4, 5, 6, 7, 18, 19, 20, 21] e:[] f:[]\n\
             i:5 g:1 p:[5, 10, 11, 6, 19, 22, 23, 20] e:[] f:[]\n\
             i:6 g:1 p:[7, 6, 14, 15, 21, 20, 24, 25] e:[] f:[]\n\
             i:7 g:1 p:[6, 11, 17, 14, 20, 23, 26, 24] e:[] f:[]\n\
             \n\
             boundary_points\n\
             0 1 2 3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 \n\
             \n\
             boundary_edges\n\
             \n\
             boundary_faces\n"
        );
        Ok(())
    }

    #[test]
    fn subdivide_3d_o2_works() -> Result<(), &'static str> {
        let mut block = Block::new(3);
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
        let mesh = block.subdivide(KindQuaOrHex::Hex20)?;
        println!("{}", mesh);
        assert_eq!(
            format!("{}", mesh),
            "ndim = 3\n\
             npoint = 81\n\
             ncell = 8\n\
             n_boundary_point = 74\n\
             n_boundary_edge = 0\n\
             n_boundary_face = 0\n\
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
             i:0 g:1 p:[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19] e:[] f:[]\n\
             i:1 g:1 p:[1, 20, 21, 2, 5, 22, 23, 6, 24, 25, 26, 9, 27, 28, 29, 13, 17, 30, 31, 18] e:[] f:[]\n\
             i:2 g:1 p:[3, 2, 32, 33, 7, 6, 34, 35, 10, 36, 37, 38, 14, 39, 40, 41, 19, 18, 42, 43] e:[] f:[]\n\
             i:3 g:1 p:[2, 21, 44, 32, 6, 23, 45, 34, 26, 46, 47, 36, 29, 48, 49, 39, 18, 31, 50, 42] e:[] f:[]\n\
             i:4 g:1 p:[4, 5, 6, 7, 51, 52, 53, 54, 12, 13, 14, 15, 55, 56, 57, 58, 59, 60, 61, 62] e:[] f:[]\n\
             i:5 g:1 p:[5, 22, 23, 6, 52, 63, 64, 53, 27, 28, 29, 13, 65, 66, 67, 56, 60, 68, 69, 61] e:[] f:[]\n\
             i:6 g:1 p:[7, 6, 34, 35, 54, 53, 70, 71, 14, 39, 40, 41, 57, 72, 73, 74, 62, 61, 75, 76] e:[] f:[]\n\
             i:7 g:1 p:[6, 23, 45, 34, 53, 64, 77, 70, 29, 48, 49, 39, 67, 78, 79, 72, 61, 69, 80, 75] e:[] f:[]\n\
             \n\
             boundary_points\n\
             0 1 2 3 4 5 7 8 9 10 11 12 15 16 17 19 20 21 22 23 24 25 26 27 28 30 31 32 33 34 35 36 37 38 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 \n\
             \n\
             boundary_edges\n\
             \n\
             boundary_faces\n"
        );
        Ok(())
    }
}
