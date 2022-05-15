use super::{Cell, Mesh, Point};

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
    }
}
