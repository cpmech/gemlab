use super::{Cell, Mesh, Point};

/// Holds samples of meshes
pub struct Samples;

impl Samples {
    /// Returns two quadrilaterals along the horizontal
    ///
    /// ```text
    /// 1.0  3-----------2-----------5
    ///      |           |           |
    ///      |    [0]    |    [1]    |  [*] indicates id
    ///      |    (1)    |    (2)    |  (*) indicates attribute_id
    ///      |           |           |
    /// 0.0  0-----------1-----------4
    ///     0.0         1.0         2.0
    /// ```
    #[rustfmt::skip]
    pub fn two_quads_horizontal() -> Mesh {
        Mesh {
            space_ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0] },
                Point { id: 3, coords: vec![0.0, 1.0] },
                Point { id: 4, coords: vec![2.0, 0.0] },
                Point { id: 5, coords: vec![2.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, geo_ndim: 2, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute_id: 2, geo_ndim: 2, points: vec![1, 4, 5, 2] },
            ],
        }
    }

    /// Returns two cubes along the vertical
    ///
    /// ```text
    ///           8-------------11
    ///          /.             /|
    ///         / .            / |
    ///        /  .           /  |
    ///       /   .          /   |       id = 1
    /// 2.0  9-------------10    |       attribute_id = 2
    ///      |    .         |    |
    ///      |    4---------|----7
    ///      |   /.         |   /|
    ///      |  / .         |  / |
    ///      | /  .         | /  |
    ///      |/   .         |/   |
    /// 1.0  5--------------6    |       id = 0
    ///      |    .         |    |       attribute_id = 1
    ///      |    0---------|----3  0.0
    ///      |   /          |   /
    ///      |  /           |  /
    ///      | /            | /
    ///      |/             |/
    /// 0.0  1--------------2   1.0
    ///     0.0            1.0
    /// ```
    #[rustfmt::skip]
    pub fn two_cubes_vertical() -> Mesh {
        Mesh {
            space_ndim: 3,
            points: vec![
                Point { id:  0, coords: vec![0.0, 0.0, 0.0] },
                Point { id:  1, coords: vec![1.0, 0.0, 0.0] },
                Point { id:  2, coords: vec![1.0, 1.0, 0.0] },
                Point { id:  3, coords: vec![0.0, 1.0, 0.0] },
                Point { id:  4, coords: vec![0.0, 0.0, 1.0] },
                Point { id:  5, coords: vec![1.0, 0.0, 1.0] },
                Point { id:  6, coords: vec![1.0, 1.0, 1.0] },
                Point { id:  7, coords: vec![0.0, 1.0, 1.0] },
                Point { id:  8, coords: vec![0.0, 0.0, 2.0] },
                Point { id:  9, coords: vec![1.0, 0.0, 2.0] },
                Point { id: 10, coords: vec![1.0, 1.0, 2.0] },
                Point { id: 11, coords: vec![0.0, 1.0, 2.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, geo_ndim: 3, points: vec![0,1,2,3, 4,5,6,7] },
                Cell { id: 1, attribute_id: 2, geo_ndim: 3, points: vec![4,5,6,7, 8,9,10,11] },
            ],
        }
    }

    /// Returns a mesh with mixed shapes in 2D
    ///
    /// ```text
    /// 1.0              4-----------3
    ///                  |           |
    ///                  |    [1]    |   [*] indicates id
    ///                  |    (2)    |   (*) indicates attribute_id
    ///                  |           |
    /// 0.0  0-----------1-----------2
    ///           [0]
    ///           (1)
    ///
    ///     0.0         1.0         2.0
    /// ```
    #[rustfmt::skip]
    pub fn mixed_shapes_2d() -> Mesh {
        Mesh {
            space_ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![2.0, 0.0] },
                Point { id: 3, coords: vec![2.0, 1.0] },
                Point { id: 4, coords: vec![1.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, geo_ndim: 1, points: vec![0, 1] },
                Cell { id: 1, attribute_id: 2, geo_ndim: 2, points: vec![1, 2, 3, 4] },
            ],
        }
    }

    /// Returns a mesh with mixed shapes in 3D
    ///
    /// ```text
    /// Note: the tetrahedron (2,8,3,6) is incompatible
    ///                       4------------7-----------10
    ///                      /.           /|            |
    ///                     / .          / |            |
    ///                    /  .         /  |            |
    ///                   /   .        /   |            |
    ///                  5------------6    |            |
    ///                  |    .       |`.  |            |
    ///                  |    0-------|--`.3------------9
    ///                  |   /        |   /`.          /
    ///                  |  /         |  /   `.       /
    ///                  | /          | /      `.    /
    ///                  |/           |/         `. /
    ///  12-----11-------1------------2------------8
    /// ```
    #[rustfmt::skip]
    pub fn mixed_shapes_3d() -> Mesh {
        Mesh {
            space_ndim: 3,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0, 0.0] },
                Point { id: 3, coords: vec![0.0, 1.0, 0.0] },
                Point { id: 4, coords: vec![0.0, 0.0, 1.0] },
                Point { id: 5, coords: vec![1.0, 0.0, 1.0] },
                Point { id: 6, coords: vec![1.0, 1.0, 1.0] },
                Point { id: 7, coords: vec![0.0, 1.0, 1.0] },
                Point { id: 8, coords: vec![1.0, 2.0, 0.0] },
                Point { id: 9, coords: vec![0.0, 2.0, 0.0] },
                Point { id:10, coords: vec![0.0, 2.0, 1.0] },
                Point { id:11, coords: vec![1.0,-0.5, 0.0] },
                Point { id:12, coords: vec![1.0,-1.0, 0.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, geo_ndim: 3, points: vec![0,1,2,3,4,5,6,7] },
                Cell { id: 1, attribute_id: 2, geo_ndim: 3, points: vec![2,8,3,6] },
                Cell { id: 2, attribute_id: 2, geo_ndim: 2, points: vec![3,9,10,7] },
                Cell { id: 3, attribute_id: 2, geo_ndim: 2, points: vec![8,9,3] },
                Cell { id: 4, attribute_id: 3, geo_ndim: 1, points: vec![1,11,12] },
            ],
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Samples;

    #[test]
    fn samples_work() {
        let mesh = Samples::two_quads_horizontal();
        assert_eq!(mesh.space_ndim, 2);
        assert_eq!(mesh.points.len(), 6);
        assert_eq!(mesh.cells.len(), 2);

        let mesh = Samples::two_cubes_vertical();
        assert_eq!(mesh.space_ndim, 3);
        assert_eq!(mesh.points.len(), 12);
        assert_eq!(mesh.cells.len(), 2);

        let mesh = Samples::mixed_shapes_2d();
        assert_eq!(mesh.space_ndim, 2);
        assert_eq!(mesh.points.len(), 5);
        assert_eq!(mesh.cells.len(), 2);

        let mesh = Samples::mixed_shapes_3d();
        assert_eq!(mesh.space_ndim, 3);
        assert_eq!(mesh.points.len(), 13);
        assert_eq!(mesh.cells.len(), 5);
    }
}
