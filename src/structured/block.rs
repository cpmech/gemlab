use std::collections::HashMap;

enum BlockShape {}

enum Constraint {}

struct ConstraintData {
    kind: Constraint,
}

/// Defines a polygon on polyhedron that can be split into smaller shapes
///
/// # Geometry
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
    npoints: usize,           // number of points (4, 8, or 20)
    nedges: usize,            // number of edges (4 or 12)
    nfaces: usize,            // number of faces (4 or 6)
    coords: Vec<Vec<f64>>,    // coordinates of points (ndim, num_points)
    point_groups: Vec<usize>, // point groups (npoints)
    edge_groups: Vec<usize>,  // edge groups (nedges)
    face_groups: Vec<usize>,  // face groups (nfaces)
    num_div: Vec<usize>,      // number of divisions along each dim (ndim)
    weights: Vec<Vec<f64>>,   // weights along each dimension (ndim, num_weights)
    sum_weights: Vec<f64>,    // sum of weights along each dimension (ndim)

    // maps side to constraint (num_sides)
    constraints: HashMap<Constraint, ConstraintData>,
}

impl Block {
    /*
    pub fn new(group: usize, ndim: usize, npoints: usize) -> Self {
        Block {
            group,
            ndim,
            coords: vec![Vec::new(); ndim],
            point_groups: Vec::new(),
            edge_groups: Vec::new(),
            face_groups: Vec::new(),
            num_div: vec![0; ndim],
            weights: vec![Vec::new(); ndim],
            sum_weights: vec![0.0; ndim],
            constraints: HashMap::new(),
        }
    }
    */
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_works() {
        // let block = Block::new(1, 3);
        // assert_eq!(block.ndim, 3);
    }
}
