use crate::AsArray2D;
use russell_lab::Matrix;
use std::collections::HashMap;

pub enum BlockKind {
    Qua4,
    Qua8,
    Hex8,
    Hex20,
}

pub enum Constraint {}

pub struct ConstraintData {
    kind: Constraint,
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
    npoint: usize,            // number of points (4, 8, or 20)
    nedge: usize,             // number of edges (4 or 12)
    nface: usize,             // number of faces (4 or 6)
    coords: Matrix,           // coordinates of points (npoint, ndim)
    point_groups: Vec<usize>, // point groups (npoint)
    edge_groups: Vec<usize>,  // edge groups (nedge)
    face_groups: Vec<usize>,  // face groups (nface)
    ndiv: Vec<usize>,         // number of divisions along each dim (ndim)
    weights: Matrix,          // weights along each dimension (ndim, ndiv)
    sum_weights: Vec<f64>,    // sum of weights along each dimension (ndim)

    // maps side to constraint (num_sides)
    constraints: HashMap<Constraint, ConstraintData>,
}

impl Block {
    /// Creates a new Block with default options
    pub fn new(kind: BlockKind) -> Self {
        let (ndim, npoint, nedge, nface) = match kind {
            BlockKind::Qua4 => (2, 4, 4, 0),
            BlockKind::Qua8 => (2, 8, 4, 0),
            BlockKind::Hex8 => (3, 8, 12, 6),
            BlockKind::Hex20 => (3, 20, 12, 6),
        };
        const NDIV: usize = 2;
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
            weights: Matrix::filled(ndim, NDIV, 1.0),
            sum_weights: vec![NDIV as f64; ndim],
            constraints: HashMap::new(),
        }
    }

    /// Sets group
    pub fn set_group(&mut self, group: usize) -> &mut Self {
        self.group = group;
        self
    }

    /// Sets 2D point
    ///
    /// # Input
    ///
    /// * `p` -- index of point in [0, npoint-1]
    /// * `group` -- group of point
    /// * `x`, `y` -- coordinates of point
    pub fn set_point_2d(&mut self, p: usize, group: usize, x: f64, y: f64) -> &mut Self {
        assert_eq!(self.ndim, 2);
        assert!(p < self.npoint);
        self.point_groups[p] = group;
        self.coords[p][0] = x;
        self.coords[p][1] = y;
        self
    }

    /// Sets 3D point
    ///
    /// # Input
    ///
    /// * `p` -- index of point in [0, npoint-1]
    /// * `group` -- group of point
    /// * `x`, `y`, `z` -- coordinates of point
    pub fn set_point_3d(&mut self, p: usize, group: usize, x: f64, y: f64, z: f64) -> &mut Self {
        assert_eq!(self.ndim, 3);
        assert!(p < self.npoint);
        self.point_groups[p] = group;
        self.coords[p][0] = x;
        self.coords[p][1] = y;
        self.coords[p][2] = z;
        self
    }

    /// Sets the coordinates of all 2D or 3D points
    ///
    /// # Input
    ///
    /// * `coords` -- (npoint, ndim) matrix with all coordinates
    pub fn set_coords<'a, T, U>(&mut self, coords: &'a T) -> &mut Self
    where
        T: AsArray2D<'a, U>,
        U: 'a + Into<f64>,
    {
        let (nrow, ncol) = coords.size();
        assert_eq!(nrow, self.npoint);
        assert_eq!(ncol, self.ndim);
        for i in 0..nrow {
            for j in 0..ncol {
                self.coords[i][j] = coords.at(i, j).into();
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
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_works() {
        let qua4 = Block::new(BlockKind::Qua4);
        assert_eq!(qua4.group, 1);
        assert_eq!(qua4.ndim, 2);
        assert_eq!(qua4.npoint, 4);
        assert_eq!(qua4.nedge, 4);
        assert_eq!(qua4.nface, 0);
        assert_eq!(
            format!("{}", qua4.coords),
            "┌     ┐\n\
             │ 0 0 │\n\
             │ 0 0 │\n\
             │ 0 0 │\n\
             │ 0 0 │\n\
             └     ┘"
        );
        assert_eq!(qua4.point_groups, &[0, 0, 0, 0]);
        assert_eq!(qua4.edge_groups, &[0, 0, 0, 0]);
        assert_eq!(qua4.face_groups, &[]);
        assert_eq!(qua4.ndiv, &[2, 2]);
        assert_eq!(
            format!("{}", qua4.weights),
            "┌     ┐\n\
             │ 1 1 │\n\
             │ 1 1 │\n\
             └     ┘"
        );
        assert_eq!(qua4.sum_weights, &[2.0, 2.0]);

        let qua8 = Block::new(BlockKind::Qua8);
        assert_eq!(qua8.group, 1);
        assert_eq!(qua8.ndim, 2);
        assert_eq!(qua8.npoint, 8);
        assert_eq!(qua8.nedge, 4);
        assert_eq!(qua8.nface, 0);
        assert_eq!(
            format!("{}", qua8.coords),
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
        assert_eq!(qua8.point_groups, &[0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(qua8.edge_groups, &[0, 0, 0, 0]);
        assert_eq!(qua8.face_groups, &[]);
        assert_eq!(qua8.ndiv, &[2, 2]);
        assert_eq!(
            format!("{}", qua8.weights),
            "┌     ┐\n\
             │ 1 1 │\n\
             │ 1 1 │\n\
             └     ┘"
        );
        assert_eq!(qua8.sum_weights, &[2.0, 2.0]);

        let hex8 = Block::new(BlockKind::Hex8);
        assert_eq!(hex8.group, 1);
        assert_eq!(hex8.ndim, 3);
        assert_eq!(hex8.npoint, 8);
        assert_eq!(hex8.nedge, 12);
        assert_eq!(hex8.nface, 6);
        assert_eq!(
            format!("{}", hex8.coords),
            "┌       ┐\n\
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
        assert_eq!(hex8.point_groups, &[0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(hex8.edge_groups, &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(hex8.face_groups, &[0, 0, 0, 0, 0, 0]);
        assert_eq!(hex8.ndiv, &[2, 2, 2]);
        assert_eq!(
            format!("{}", hex8.weights),
            "┌     ┐\n\
             │ 1 1 │\n\
             │ 1 1 │\n\
             │ 1 1 │\n\
             └     ┘"
        );
        assert_eq!(hex8.sum_weights, &[2.0, 2.0, 2.0]);

        let hex20 = Block::new(BlockKind::Hex20);
        assert_eq!(hex20.group, 1);
        assert_eq!(hex20.ndim, 3);
        assert_eq!(hex20.npoint, 20);
        assert_eq!(hex20.nedge, 12);
        assert_eq!(hex20.nface, 6);
        assert_eq!(
            format!("{}", hex20.coords),
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
            hex20.point_groups,
            &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        );
        assert_eq!(hex20.edge_groups, &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(hex20.face_groups, &[0, 0, 0, 0, 0, 0]);
        assert_eq!(hex20.ndiv, &[2, 2, 2]);
        assert_eq!(
            format!("{}", hex20.weights),
            "┌     ┐\n\
             │ 1 1 │\n\
             │ 1 1 │\n\
             │ 1 1 │\n\
             └     ┘"
        );
        assert_eq!(hex20.sum_weights, &[2.0, 2.0, 2.0]);
    }

    #[test]
    fn set_group_works() {
        let mut block = Block::new(BlockKind::Qua4);
        block.set_group(2);
        assert_eq!(block.group, 2);
    }

    #[test]
    fn set_point_2d_works() {
        let mut block = Block::new(BlockKind::Qua4);
        block
            .set_point_2d(0, 1, 0.0, 0.0)
            .set_point_2d(1, 2, 2.0, 0.0)
            .set_point_2d(2, 3, 2.0, 2.0)
            .set_point_2d(3, 4, 0.0, 2.0);
        assert_eq!(
            format!("{}", block.coords),
            "┌     ┐\n\
             │ 0 0 │\n\
             │ 2 0 │\n\
             │ 2 2 │\n\
             │ 0 2 │\n\
             └     ┘"
        );
    }

    #[test]
    fn set_point_3d_works() {
        let mut block = Block::new(BlockKind::Hex8);
        block
            .set_point_3d(0, 1, 0.0, 0.0, 0.0)
            .set_point_3d(1, 2, 2.0, 0.0, 0.0)
            .set_point_3d(2, 3, 2.0, 2.0, 0.0)
            .set_point_3d(3, 4, 0.0, 2.0, 0.0)
            .set_point_3d(4, 1, 0.0, 0.0, 2.0)
            .set_point_3d(5, 2, 2.0, 0.0, 2.0)
            .set_point_3d(6, 3, 2.0, 2.0, 2.0)
            .set_point_3d(7, 4, 0.0, 2.0, 2.0);
        assert_eq!(
            format!("{}", block.coords),
            "┌       ┐\n\
             │ 0 0 0 │\n\
             │ 2 0 0 │\n\
             │ 2 2 0 │\n\
             │ 0 2 0 │\n\
             │ 0 0 2 │\n\
             │ 2 0 2 │\n\
             │ 2 2 2 │\n\
             │ 0 2 2 │\n\
             └       ┘"
        );
    }

    #[test]
    fn set_coords_works() {
        let mut qua4 = Block::new(BlockKind::Qua4);
        #[rustfmt::skip]
        qua4.set_coords(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ]);
        assert_eq!(
            format!("{}", qua4.coords),
            "┌     ┐\n\
             │ 0 0 │\n\
             │ 2 0 │\n\
             │ 2 2 │\n\
             │ 0 2 │\n\
             └     ┘"
        );

        let mut hex8 = Block::new(BlockKind::Hex8);
        #[rustfmt::skip]
        hex8.set_coords(&[
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
            format!("{}", hex8.coords),
            "┌       ┐\n\
             │ 0 0 0 │\n\
             │ 2 0 0 │\n\
             │ 2 2 0 │\n\
             │ 0 2 0 │\n\
             │ 0 0 2 │\n\
             │ 2 0 2 │\n\
             │ 2 2 2 │\n\
             │ 0 2 2 │\n\
             └       ┘"
        );
    }

    #[test]
    fn set_point_group_works() {
        let mut block = Block::new(BlockKind::Qua4);
        block.set_point_group(0, 111);
        assert_eq!(block.point_groups, &[111, 0, 0, 0]);
    }

    #[test]
    fn set_edge_group_works() {
        let mut block = Block::new(BlockKind::Qua4);
        block.set_edge_group(0, 111);
        assert_eq!(block.edge_groups, &[111, 0, 0, 0]);
    }

    #[test]
    fn set_face_group_works() {
        let mut block = Block::new(BlockKind::Hex8);
        block.set_face_group(0, 111);
        assert_eq!(block.face_groups, &[111, 0, 0, 0, 0, 0]);
    }
}
