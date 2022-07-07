use super::{Cell, Mesh, Point};
use crate::shapes::GeoKind;

/// Holds samples of meshes
pub struct Samples;

impl Samples {
    /// Returns two quadrilaterals along the horizontal
    ///
    /// ```text
    ///          [#] indicates id
    ///      y   (#) indicates attribute_id
    ///      ↑
    /// 1.0  3-----------2-----------5
    ///      |           |           |
    ///      |    [0]    |    [1]    |
    ///      |    (1)    |    (2)    |
    ///      |           |           |
    /// 0.0  0-----------1-----------4  → x
    ///     0.0         1.0         2.0
    /// ```
    #[rustfmt::skip]
    pub fn two_quads_horizontal() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0] },
                Point { id: 3, coords: vec![0.0, 1.0] },
                Point { id: 4, coords: vec![2.0, 0.0] },
                Point { id: 5, coords: vec![2.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute_id: 2, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
            ],
        }
    }

    /// Returns two cubes along the vertical
    ///
    /// ```text
    ///           [#] indicates id
    ///           (#) indicates attribute_id
    ///
    ///       8-------------11  2.0
    ///      /.             /|
    ///     / .            / |
    ///    /  .           /  |
    ///   /   .          /   |
    ///  9-------------10    |
    ///  |    .         |    |
    ///  |    4---------|----7  1.0
    ///  |   /. [1]     |   /|
    ///  |  / . (2)     |  / |
    ///  | /  .         | /  |
    ///  |/   .         |/   |
    ///  5--------------6    |          z
    ///  |    .         |    |          ↑
    ///  |    0---------|----3  0.0     o → y
    ///  |   /  [0]     |   /          ↙
    ///  |  /   (1)     |  /          x
    ///  | /            | /
    ///  |/             |/
    ///  1--------------2   1.0
    /// 0.0            1.0
    /// ```
    #[rustfmt::skip]
    pub fn two_cubes_vertical() -> Mesh {
        Mesh {
            ndim: 3,
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
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Hex8, points: vec![0,1,2,3, 4,5,6,7] },
                Cell { id: 1, attribute_id: 2, kind: GeoKind::Hex8, points: vec![4,5,6,7, 8,9,10,11] },
            ],
        }
    }

    /// Returns four cubes making a vertical wall
    ///
    /// ```text
    ///           [#] indicates id
    ///           (#) indicates attribute_id
    ///
    ///       8-------------11-------------17  2.0
    ///      /.             /.             /|
    ///     / .            / .            / |
    ///    /  .           /  .           /  |
    ///   /   .          /   .          /   |
    ///  9=============10=============16    |
    ///  |    .         |    .         |    |
    ///  |    4---------|----7---------|---15  1.0
    ///  |   /. [1]     |   /. [3]     |   /|
    ///  |  / . (2)     |  / . (2)     |  / |
    ///  | /  .         | /  .         | /  |
    ///  |/   .         |/   .         |/   |
    ///  5--------------6-------------14    |          z
    ///  |    .         |    .         |    |          ↑
    ///  |    0---------|----3---------|---13  0.0     o → y
    ///  |   /  [0]     |   /  [2]     |   /          ↙
    ///  |  /   (1)     |  /   (1)     |  /          x
    ///  | /            | /            | /
    ///  |/             |/             |/
    ///  1--------------2-------------12   1.0
    /// 0.0            1.0            2.0
    /// ```
    #[rustfmt::skip]
    pub fn four_cubes_wall() -> Mesh {
        Mesh {
            ndim: 3,
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
                Point { id: 12, coords: vec![1.0, 2.0, 0.0] },
                Point { id: 13, coords: vec![0.0, 2.0, 0.0] },
                Point { id: 14, coords: vec![1.0, 2.0, 1.0] },
                Point { id: 15, coords: vec![0.0, 2.0, 1.0] },
                Point { id: 16, coords: vec![1.0, 2.0, 2.0] },
                Point { id: 17, coords: vec![0.0, 2.0, 2.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Hex8, points: vec![0,1, 2, 3,  4, 5, 6, 7] },
                Cell { id: 1, attribute_id: 2, kind: GeoKind::Hex8, points: vec![4,5, 6, 7,  8, 9,10,11] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Hex8, points: vec![3,2,12,13,  7, 6,14,15] },
                Cell { id: 3, attribute_id: 2, kind: GeoKind::Hex8, points: vec![7,6,14,15, 11,10,16,17] },
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
    /// 0.0  0-----------1-----------2-----------5
    ///           [0]                     [2]
    ///           (1)                     (1)
    ///
    ///     0.0         1.0         2.0         3.0
    /// ```
    #[rustfmt::skip]
    pub fn mixed_shapes_2d() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![2.0, 0.0] },
                Point { id: 3, coords: vec![2.0, 1.0] },
                Point { id: 4, coords: vec![1.0, 1.0] },
                Point { id: 5, coords: vec![3.0, 0.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
                Cell { id: 1, attribute_id: 2, kind: GeoKind::Qua4, points: vec![1, 2, 3, 4] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Lin2, points: vec![2, 5] },
            ],
        }
    }

    /// Returns a mesh with mixed shapes in 3D
    ///
    /// **Note:** There is a tetrahedron (Tet4) on with nodes (2,8,3,6).
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
            ndim: 3,
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
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Hex8, points: vec![0,1,2,3,4,5,6,7] },
                Cell { id: 1, attribute_id: 2, kind: GeoKind::Tet4, points: vec![2,8,3,6] },
                Cell { id: 2, attribute_id: 2, kind: GeoKind::Qua4, points: vec![3,9,10,7] },
                Cell { id: 3, attribute_id: 2, kind: GeoKind::Tri3, points: vec![8,9,3] },
                Cell { id: 4, attribute_id: 3, kind: GeoKind::Lin3, points: vec![1,12,11] },
            ],
        }
    }

    /// Returns a mesh with four qua4 cells
    ///
    /// ```text
    /// 2.0   7---------------6---------------8
    ///       |               |               |
    ///       |               |               |
    ///       |      [2]      |      [3]      |
    ///       |               |               |
    ///       |               |               |
    /// 1.0   3---------------2---------------5
    ///       |               |               |
    ///       |               |               |
    ///       |      [0]      |      [1]      |
    ///       |               |               |
    ///       |               |               |
    /// 0.0   0---------------1---------------4
    ///
    ///      0.0             1.0             2.0
    ///
    /// xmin = 0.0, xmax = 2.0
    /// ymin = 0.0, ymax = 2.0
    /// ```
    #[rustfmt::skip]
    pub fn block_2d_four_qua4() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, coords: vec![1.0, 0.0] },
                Point { id: 2, coords: vec![1.0, 1.0] },
                Point { id: 3, coords: vec![0.0, 1.0] },
                Point { id: 4, coords: vec![2.0, 0.0] },
                Point { id: 5, coords: vec![2.0, 1.0] },
                Point { id: 6, coords: vec![1.0, 2.0] },
                Point { id: 7, coords: vec![0.0, 2.0] },
                Point { id: 8, coords: vec![2.0, 2.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua4, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Qua4, points: vec![1, 4, 5, 2] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Qua4, points: vec![3, 2, 6, 7] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Qua4, points: vec![2, 5, 8, 6] },
            ],
        }
    }

    /// Returns a mesh with four qua8 cells
    ///
    /// ```text
    /// 2.0  14------16------13------20------18
    ///       |               |               |
    ///       |               |               |
    /// 1.5  17      [2]     15      [3]     19
    ///       |               |               |
    ///       |               |               |
    /// 1.0   3-------6-------2------12-------9
    ///       |               |               |
    ///       |               |               |
    /// 0.5   7      [0]      5      [1]     11
    ///       |               |               |
    ///       |               |               |
    /// 0.0   0-------4-------1------10-------8
    ///
    ///      0.0     0.5     1.0     1.5     2.0
    ///
    /// xmin = 0.0, xmax = 2.0
    /// ymin = 0.0, ymax = 2.0
    /// ```
    #[rustfmt::skip]
    pub fn block_2d_four_qua8() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id:  0, coords: vec![0.0, 0.0] },
                Point { id:  1, coords: vec![1.0, 0.0] },
                Point { id:  2, coords: vec![1.0, 1.0] },
                Point { id:  3, coords: vec![0.0, 1.0] },
                Point { id:  4, coords: vec![0.5, 0.0] },
                Point { id:  5, coords: vec![1.0, 0.5] },
                Point { id:  6, coords: vec![0.5, 1.0] },
                Point { id:  7, coords: vec![0.0, 0.5] },
                Point { id:  8, coords: vec![2.0, 0.0] },
                Point { id:  9, coords: vec![2.0, 1.0] },
                Point { id: 10, coords: vec![1.5, 0.0] },
                Point { id: 11, coords: vec![2.0, 0.5] },
                Point { id: 12, coords: vec![1.5, 1.0] },
                Point { id: 13, coords: vec![1.0, 2.0] },
                Point { id: 14, coords: vec![0.0, 2.0] },
                Point { id: 15, coords: vec![1.0, 1.5] },
                Point { id: 16, coords: vec![0.5, 2.0] },
                Point { id: 17, coords: vec![0.0, 1.5] },
                Point { id: 18, coords: vec![2.0, 2.0] },
                Point { id: 19, coords: vec![2.0, 1.5] },
                Point { id: 20, coords: vec![1.5, 2.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua8, points: vec![0, 1,  2,  3,  4,  5,  6,  7] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Qua8, points: vec![1, 8,  9,  2, 10, 11, 12,  5] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Qua8, points: vec![3, 2, 13, 14,  6, 15, 16, 17] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Qua8, points: vec![2, 9, 18, 13, 12, 19, 20, 15] },
            ],
        }
    }

    /// Returns a mesh with four qua9 cells
    ///
    /// ```text
    /// cell ids as in block_2d_four_qua8
    ///
    /// 2.0  16------18------15------23------21
    ///       |               |               |
    ///       |               |               |
    /// 1.5  19      20      17      24      22
    ///       |               |               |
    ///       |               |               |
    /// 1.0   3-------6-------2------13------10
    ///       |               |               |
    ///       |               |               |
    /// 0.5   7       8       5      14      12
    ///       |               |               |
    ///       |               |               |
    /// 0.0   0-------4-------1------11-------9
    ///
    ///      0.0     0.5     1.0     1.5     2.0
    ///
    /// xmin = 0.0, xmax = 2.0
    /// ymin = 0.0, ymax = 2.0
    /// ```
    #[rustfmt::skip]
    pub fn block_2d_four_qua9() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id:  0, coords: vec![0.0, 0.0] },
                Point { id:  1, coords: vec![1.0, 0.0] },
                Point { id:  2, coords: vec![1.0, 1.0] },
                Point { id:  3, coords: vec![0.0, 1.0] },
                Point { id:  4, coords: vec![0.5, 0.0] },
                Point { id:  5, coords: vec![1.0, 0.5] },
                Point { id:  6, coords: vec![0.5, 1.0] },
                Point { id:  7, coords: vec![0.0, 0.5] },
                Point { id:  8, coords: vec![0.5, 0.5] },
                Point { id:  9, coords: vec![2.0, 0.0] },
                Point { id: 10, coords: vec![2.0, 1.0] },
                Point { id: 11, coords: vec![1.5, 0.0] },
                Point { id: 12, coords: vec![2.0, 0.5] },
                Point { id: 13, coords: vec![1.5, 1.0] },
                Point { id: 14, coords: vec![1.5, 0.5] },
                Point { id: 15, coords: vec![1.0, 2.0] },
                Point { id: 16, coords: vec![0.0, 2.0] },
                Point { id: 17, coords: vec![1.0, 1.5] },
                Point { id: 18, coords: vec![0.5, 2.0] },
                Point { id: 19, coords: vec![0.0, 1.5] },
                Point { id: 20, coords: vec![0.5, 1.5] },
                Point { id: 21, coords: vec![2.0, 2.0] },
                Point { id: 22, coords: vec![2.0, 1.5] },
                Point { id: 23, coords: vec![1.5, 2.0] },
                Point { id: 24, coords: vec![1.5, 1.5] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua9, points: vec![0,  1,  2,  3,  4,  5,  6,  7,  8] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Qua9, points: vec![1,  9, 10,  2, 11, 12, 13,  5, 14] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Qua9, points: vec![3,  2, 15, 16,  6, 17, 18, 19, 20] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Qua9, points: vec![2, 10, 21, 15, 13, 22, 23, 17, 24] },
            ],
        }
    }

    /// Returns a mesh with four qua12 cells
    ///
    /// ```text
    /// cell ids as in block_2d_four_qua8
    ///
    /// 3.0  21---26---23----20---32---30----28
    ///       |               |               |
    /// 2.5  24              25              31
    ///       |               |               |
    /// 2.0  27              22              29
    ///       |               |               |
    /// 1.5   3---10-----6----2---19---16----13
    ///       |               |               |
    /// 1.0   7               9              18
    ///       |               |               |
    /// 0.5  11               5              15
    ///       |               |               |
    /// 0.0   0----4-----8----1---14---17----12
    ///
    ///      0.0  0.5   1.0  1.5  2.0  2.5   3.0
    ///
    /// xmin = 0.0, xmax = 3.0
    /// ymin = 0.0, ymax = 3.0
    /// ```
    #[rustfmt::skip]
    pub fn block_2d_four_qua12() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id:  0, coords: vec![0.0, 0.0] },
                Point { id:  1, coords: vec![1.5, 0.0] },
                Point { id:  2, coords: vec![1.5, 1.5] },
                Point { id:  3, coords: vec![0.0, 1.5] },
                Point { id:  4, coords: vec![0.5, 0.0] },
                Point { id:  5, coords: vec![1.5, 0.5] },
                Point { id:  6, coords: vec![1.0, 1.5] },
                Point { id:  7, coords: vec![0.0, 1.0] },
                Point { id:  8, coords: vec![1.0, 0.0] },
                Point { id:  9, coords: vec![1.5, 1.0] },
                Point { id: 10, coords: vec![0.5, 1.5] },
                Point { id: 11, coords: vec![0.0, 0.5] },
                Point { id: 12, coords: vec![3.0, 0.0] },
                Point { id: 13, coords: vec![3.0, 1.5] },
                Point { id: 14, coords: vec![2.0, 0.0] },
                Point { id: 15, coords: vec![3.0, 0.5] },
                Point { id: 16, coords: vec![2.5, 1.5] },
                Point { id: 17, coords: vec![2.5, 0.0] },
                Point { id: 18, coords: vec![3.0, 1.0] },
                Point { id: 19, coords: vec![2.0, 1.5] },
                Point { id: 20, coords: vec![1.5, 3.0] },
                Point { id: 21, coords: vec![0.0, 3.0] },
                Point { id: 22, coords: vec![1.5, 2.0] },
                Point { id: 23, coords: vec![1.0, 3.0] },
                Point { id: 24, coords: vec![0.0, 2.5] },
                Point { id: 25, coords: vec![1.5, 2.5] },
                Point { id: 26, coords: vec![0.5, 3.0] },
                Point { id: 27, coords: vec![0.0, 2.0] },
                Point { id: 28, coords: vec![3.0, 3.0] },
                Point { id: 29, coords: vec![3.0, 2.0] },
                Point { id: 30, coords: vec![2.5, 3.0] },
                Point { id: 31, coords: vec![3.0, 2.5] },
                Point { id: 32, coords: vec![2.0, 3.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua12, points: vec![0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Qua12, points: vec![1, 12, 13,  2, 14, 15, 16,  9, 17, 18, 19,  5] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Qua12, points: vec![3,  2, 20, 21, 10, 22, 23, 24,  6, 25, 26, 27] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Qua12, points: vec![2, 13, 28, 20, 19, 29, 30, 25, 16, 31, 32, 22] },
            ],
        }
    }

    /// Returns a mesh with four qua16 cells
    ///
    /// ```text
    /// cell ids as in block_2d_four_qua8
    ///
    /// 3.0  29---34----31---28---44---42----40
    ///       |               |               |
    /// 2.5  32   39    38   33   48   47    43
    ///       |               |               |
    /// 2.0  35   36    37   30   45   46    41
    ///       |               |               |
    /// 1.5   3---10-----6----2---23---20----17
    ///       |               |               |
    /// 1.0   7   15    14    9   27   26    22
    ///       |               |               |
    /// 0.5  11   12    13    5   24   25    19
    ///       |               |               |
    /// 0.0   0----4-----8----1---18---21----16
    ///
    ///      0.0  0.5   1.0  1.5  2.0  2.5   3.0
    ///
    /// xmin = 0.0, xmax = 3.0
    /// ymin = 0.0, ymax = 3.0
    /// ```
    #[rustfmt::skip]
    pub fn block_2d_four_qua16() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id:  0, coords: vec![0.0, 0.0] },
                Point { id:  1, coords: vec![1.5, 0.0] },
                Point { id:  2, coords: vec![1.5, 1.5] },
                Point { id:  3, coords: vec![0.0, 1.5] },
                Point { id:  4, coords: vec![0.5, 0.0] },
                Point { id:  5, coords: vec![1.5, 0.5] },
                Point { id:  6, coords: vec![1.0, 1.5] },
                Point { id:  7, coords: vec![0.0, 1.0] },
                Point { id:  8, coords: vec![1.0, 0.0] },
                Point { id:  9, coords: vec![1.5, 1.0] },
                Point { id: 10, coords: vec![0.5, 1.5] },
                Point { id: 11, coords: vec![0.0, 0.5] },
                Point { id: 12, coords: vec![0.5, 0.5] },
                Point { id: 13, coords: vec![1.0, 0.5] },
                Point { id: 14, coords: vec![1.0, 1.0] },
                Point { id: 15, coords: vec![0.5, 1.0] },
                Point { id: 16, coords: vec![3.0, 0.0] },
                Point { id: 17, coords: vec![3.0, 1.5] },
                Point { id: 18, coords: vec![2.0, 0.0] },
                Point { id: 19, coords: vec![3.0, 0.5] },
                Point { id: 20, coords: vec![2.5, 1.5] },
                Point { id: 21, coords: vec![2.5, 0.0] },
                Point { id: 22, coords: vec![3.0, 1.0] },
                Point { id: 23, coords: vec![2.0, 1.5] },
                Point { id: 24, coords: vec![2.0, 0.5] },
                Point { id: 25, coords: vec![2.5, 0.5] },
                Point { id: 26, coords: vec![2.5, 1.0] },
                Point { id: 27, coords: vec![2.0, 1.0] },
                Point { id: 28, coords: vec![1.5, 3.0] },
                Point { id: 29, coords: vec![0.0, 3.0] },
                Point { id: 30, coords: vec![1.5, 2.0] },
                Point { id: 31, coords: vec![1.0, 3.0] },
                Point { id: 32, coords: vec![0.0, 2.5] },
                Point { id: 33, coords: vec![1.5, 2.5] },
                Point { id: 34, coords: vec![0.5, 3.0] },
                Point { id: 35, coords: vec![0.0, 2.0] },
                Point { id: 36, coords: vec![0.5, 2.0] },
                Point { id: 37, coords: vec![1.0, 2.0] },
                Point { id: 38, coords: vec![1.0, 2.5] },
                Point { id: 39, coords: vec![0.5, 2.5] },
                Point { id: 40, coords: vec![3.0, 3.0] },
                Point { id: 41, coords: vec![3.0, 2.0] },
                Point { id: 42, coords: vec![2.5, 3.0] },
                Point { id: 43, coords: vec![3.0, 2.5] },
                Point { id: 44, coords: vec![2.0, 3.0] },
                Point { id: 45, coords: vec![2.0, 2.0] },
                Point { id: 46, coords: vec![2.5, 2.0] },
                Point { id: 47, coords: vec![2.5, 2.5] },
                Point { id: 48, coords: vec![2.0, 2.5] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua16, points: vec![0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Qua16, points: vec![1, 16, 17,  2, 18, 19, 20,  9, 21, 22, 23,  5, 24, 25, 26, 27] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Qua16, points: vec![3,  2, 28, 29, 10, 30, 31, 32,  6, 33, 34, 35, 36, 37, 38, 39] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Qua16, points: vec![2, 17, 40, 28, 23, 41, 42, 33, 20, 43, 44, 30, 45, 46, 47, 48] },
            ],
        }
    }

    /// Returns a mesh with four qua17 cells
    ///
    /// ```text
    /// cell ids as in block_2d_four_qua8
    ///
    /// 4.0  30---38---32---37---29---48---43---47---41
    ///       |                   |                   |
    /// 3.5  39                  36                  46
    ///       |                   |                   |
    /// 3.0  33        34        31        44        42
    ///       |                   |                   |
    /// 2.5  40                  35                  45
    ///       |                   |                   |
    /// 2.0   3---14----6---13----2---28---21---27---18
    ///       |                   |                   |
    /// 1.5  15                  12                  26
    ///       |                   |                   |
    /// 1.0   7         8         5        22        20
    ///       |                   |                   |
    /// 0.5  16                  11                  25
    ///       |                   |                   |
    /// 0.0   0----9----4---10----1---23---19---24---17
    ///
    ///      0.0  0.5  1.0  1.5  2.0  2.5  3.0  3.5  4.0
    ///
    /// xmin = 0.0, xmax = 4.0
    /// ymin = 0.0, ymax = 4.0
    /// ```
    #[rustfmt::skip]
    pub fn block_2d_four_qua17() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id:  0, coords: vec![0.0, 0.0] },
                Point { id:  1, coords: vec![2.0, 0.0] },
                Point { id:  2, coords: vec![2.0, 2.0] },
                Point { id:  3, coords: vec![0.0, 2.0] },
                Point { id:  4, coords: vec![1.0, 0.0] },
                Point { id:  5, coords: vec![2.0, 1.0] },
                Point { id:  6, coords: vec![1.0, 2.0] },
                Point { id:  7, coords: vec![0.0, 1.0] },
                Point { id:  8, coords: vec![1.0, 1.0] },
                Point { id:  9, coords: vec![0.5, 0.0] },
                Point { id: 10, coords: vec![1.5, 0.0] },
                Point { id: 11, coords: vec![2.0, 0.5] },
                Point { id: 12, coords: vec![2.0, 1.5] },
                Point { id: 13, coords: vec![1.5, 2.0] },
                Point { id: 14, coords: vec![0.5, 2.0] },
                Point { id: 15, coords: vec![0.0, 1.5] },
                Point { id: 16, coords: vec![0.0, 0.5] },
                Point { id: 17, coords: vec![4.0, 0.0] },
                Point { id: 18, coords: vec![4.0, 2.0] },
                Point { id: 19, coords: vec![3.0, 0.0] },
                Point { id: 20, coords: vec![4.0, 1.0] },
                Point { id: 21, coords: vec![3.0, 2.0] },
                Point { id: 22, coords: vec![3.0, 1.0] },
                Point { id: 23, coords: vec![2.5, 0.0] },
                Point { id: 24, coords: vec![3.5, 0.0] },
                Point { id: 25, coords: vec![4.0, 0.5] },
                Point { id: 26, coords: vec![4.0, 1.5] },
                Point { id: 27, coords: vec![3.5, 2.0] },
                Point { id: 28, coords: vec![2.5, 2.0] },
                Point { id: 29, coords: vec![2.0, 4.0] },
                Point { id: 30, coords: vec![0.0, 4.0] },
                Point { id: 31, coords: vec![2.0, 3.0] },
                Point { id: 32, coords: vec![1.0, 4.0] },
                Point { id: 33, coords: vec![0.0, 3.0] },
                Point { id: 34, coords: vec![1.0, 3.0] },
                Point { id: 35, coords: vec![2.0, 2.5] },
                Point { id: 36, coords: vec![2.0, 3.5] },
                Point { id: 37, coords: vec![1.5, 4.0] },
                Point { id: 38, coords: vec![0.5, 4.0] },
                Point { id: 39, coords: vec![0.0, 3.5] },
                Point { id: 40, coords: vec![0.0, 2.5] },
                Point { id: 41, coords: vec![4.0, 4.0] },
                Point { id: 42, coords: vec![4.0, 3.0] },
                Point { id: 43, coords: vec![3.0, 4.0] },
                Point { id: 44, coords: vec![3.0, 3.0] },
                Point { id: 45, coords: vec![4.0, 2.5] },
                Point { id: 46, coords: vec![4.0, 3.5] },
                Point { id: 47, coords: vec![3.5, 4.0] },
                Point { id: 48, coords: vec![2.5, 4.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua17, points: vec![0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Qua17, points: vec![1, 17, 18,  2, 19, 20, 21,  5, 22, 23, 24, 25, 26, 27, 28, 12, 11] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Qua17, points: vec![3,  2, 29, 30,  6, 31, 32, 33, 34, 14, 13, 35, 36, 37, 38, 39, 40] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Qua17, points: vec![2, 18, 41, 29, 21, 42, 43, 31, 44, 28, 27, 45, 46, 47, 48, 36, 35] },
            ],
        }
    }

    /// Returns a mesh with eight hex8 cells
    ///
    /// ```text
    ///              18------------------21------------------25
    ///              /.                  /.                  /|
    ///             / .                 / .                 / |
    ///            /  .                /  .                /  |
    ///           /   .               /   .               /   |
    ///          /    .              /    .              /    |
    ///        19------------------20------------------24     |
    ///        /.     .            /.     .            /|     |
    ///       / .     .           / .     .           / |     |
    ///      /  .     .          /  .     .          /  |     |
    ///     /   .     .         /   .     .         /   |     |
    ///    /    .     .        /    .     .        /    |     |
    ///  22==================23==================26     |     |
    ///   |     .     .       |     .     .       |     |     |
    ///   |     .     .       |     .     .       |     |     |
    ///   |     .     4 - - - | - - . - - 7 - - - | - - | - -15
    ///   |     .    /.       |     .    /.       |     |    /|
    ///   |     .   / .       |     .   / .       |     |   / |
    ///   |     .  /  .       |     .  /  .       |     |  /  |
    ///   |     . /   .       |     . /   .       |     | /   |
    ///   |     ./    .       |     ./    .       |     |/    |
    ///   |     5 - - - - - - | - - 6 - - - - - - | - -14     |
    ///   |    /.     .       |    /.     .       |    /|     |
    ///   |   / .     .       |   / .     .       |   / |     |
    ///   |  /  .     .       |  /  .     .       |  /  |     |
    ///   | /   .     .       | /   .     .       | /   |     |
    ///   |/    .     .       |/    .     .       |/    |     |
    ///  10==================11==================17     |     |
    ///   |     .     .       |     .     .       |     |     |
    ///   |     .     .       |     .     .       |     |     |
    ///   |     .     0 - - - | - - . - - 3 - - - | - - | - -13
    ///   |     .    /        |     .    /        |     |    /
    ///   |     .   /         |     .   /         |     |   /
    ///   |     .  /          |     .  /          |     |  /
    ///   |     . /           |     . /           |     | /
    ///   |     ./            |     ./            |     |/
    ///   |     1 - - - - - - | - - 2 - - - - - - | - -12
    ///   |    /              |    /              |    /
    ///   |   /               |   /               |   /
    ///   |  /                |  /                |  /
    ///   | /                 | /                 | /
    ///   |/                  |/                  |/
    ///   8===================9==================16
    ///
    /// xmin = 0.0, xmax = 2.0
    /// ymin = 0.0, ymax = 2.0
    /// zmin = 0.0, zmax = 4.0
    /// ```
    #[rustfmt::skip]
    pub fn block_3d_eight_hex8() -> Mesh {
        Mesh {
            ndim: 3,
            points: vec![
                Point { id:  0, coords: vec![0.0, 0.0, 0.0] },
                Point { id:  1, coords: vec![1.0, 0.0, 0.0] },
                Point { id:  2, coords: vec![1.0, 1.0, 0.0] },
                Point { id:  3, coords: vec![0.0, 1.0, 0.0] },
                Point { id:  4, coords: vec![0.0, 0.0, 2.0] },
                Point { id:  5, coords: vec![1.0, 0.0, 2.0] },
                Point { id:  6, coords: vec![1.0, 1.0, 2.0] },
                Point { id:  7, coords: vec![0.0, 1.0, 2.0] },
                Point { id:  8, coords: vec![2.0, 0.0, 0.0] },
                Point { id:  9, coords: vec![2.0, 1.0, 0.0] },
                Point { id: 10, coords: vec![2.0, 0.0, 2.0] },
                Point { id: 11, coords: vec![2.0, 1.0, 2.0] },
                Point { id: 12, coords: vec![1.0, 2.0, 0.0] },
                Point { id: 13, coords: vec![0.0, 2.0, 0.0] },
                Point { id: 14, coords: vec![1.0, 2.0, 2.0] },
                Point { id: 15, coords: vec![0.0, 2.0, 2.0] },
                Point { id: 16, coords: vec![2.0, 2.0, 0.0] },
                Point { id: 17, coords: vec![2.0, 2.0, 2.0] },
                Point { id: 18, coords: vec![0.0, 0.0, 4.0] },
                Point { id: 19, coords: vec![1.0, 0.0, 4.0] },
                Point { id: 20, coords: vec![1.0, 1.0, 4.0] },
                Point { id: 21, coords: vec![0.0, 1.0, 4.0] },
                Point { id: 22, coords: vec![2.0, 0.0, 4.0] },
                Point { id: 23, coords: vec![2.0, 1.0, 4.0] },
                Point { id: 24, coords: vec![1.0, 2.0, 4.0] },
                Point { id: 25, coords: vec![0.0, 2.0, 4.0] },
                Point { id: 26, coords: vec![2.0, 2.0, 4.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Hex8, points: vec![0,  1,  2,  3,  4,  5,  6,  7] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Hex8, points: vec![1,  8,  9,  2,  5, 10, 11,  6] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Hex8, points: vec![3,  2, 12, 13,  7,  6, 14, 15] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Hex8, points: vec![2,  9, 16, 12,  6, 11, 17, 14] },
                Cell { id: 4, attribute_id: 1, kind: GeoKind::Hex8, points: vec![4,  5,  6,  7, 18, 19, 20, 21] },
                Cell { id: 5, attribute_id: 1, kind: GeoKind::Hex8, points: vec![5, 10, 11,  6, 19, 22, 23, 20] },
                Cell { id: 6, attribute_id: 1, kind: GeoKind::Hex8, points: vec![7,  6, 14, 15, 21, 20, 24, 25] },
                Cell { id: 7, attribute_id: 1, kind: GeoKind::Hex8, points: vec![6, 11, 17, 14, 20, 23, 26, 24] },
            ],
        }
    }

    /// Returns a mesh with eight hex20 cells
    ///
    /// ```text
    ///              51--------58--------54--------74--------71
    ///              /.                  /.                  /|
    ///             / .                 / .                 / |
    ///           55  .               57  .               73  |
    ///           /   .               /   .               /   |
    ///          /    .              /    .              /    |
    ///        52--------56--------53--------72--------70     |
    ///        /.     .            /.     .            /|     |
    ///       / .    59           / .    62           / |    76
    ///     65  .     .         67  .     .         79  |     |
    ///     /   .     .         /   .     .         /   |     |
    ///    /    .     .        /    .     .        /    |     |
    ///  63========66========64========78========77     |     |
    ///   |     .     .       |     .     .       |     |     |
    ///   |    60     .       |    61     .       |    75     |
    ///   |     .     4 - - - |15 - . - - 7 - - - |41 - | - -35
    ///   |     .    /.       |     .    /.       |     |    /|
    ///   |     .   / .       |     .   / .       |     |   / |
    ///   |     . 12  .       |     . 14  .       |     | 40  |
    ///   |     . /   .       |     . /   .       |     | /   |
    ///  68     ./    .      69     ./    .      80     |/    |
    ///   |     5 - - - -13 - | - - 6 - - - -39 - | - -34     |
    ///   |    /.     .       |    /.     .       |    /|     |
    ///   |   / .    16       |   / .    19       |   / |    43
    ///   | 27  .     .       | 29  .     .       | 49  |     |
    ///   | /   .     .       | /   .     .       | /   |     |
    ///   |/    .     .       |/    .     .       |/    |     |
    ///  22========28========23========48========45     |     |
    ///   |     .     .       |     .     .       |     |     |
    ///   |    17     .       |    18     .       |    42     |
    ///   |     .     0 - - - |11 - . - - 3 - - - |38 - | - -33
    ///   |     .    /        |     .    /        |     |    /
    ///   |     .   /         |     .   /         |     |   /
    ///   |     .  8          |     . 10          |     | 37
    ///   |     . /           |     . /           |     | /
    ///  30     ./           31     ./           50     |/
    ///   |     1 - - - - 9 - | - - 2 - - - -36 - | - -32
    ///   |    /              |    /              |    /
    ///   |   /               |   /               |   /
    ///   | 24                | 26                | 47
    ///   | /                 | /                 | /
    ///   |/                  |/                  |/
    ///  20========25========21========46========44
    ///
    /// xmin = 0.0, xmax = 2.0
    /// ymin = 0.0, ymax = 2.0
    /// zmin = 0.0, zmax = 4.0
    /// ```
    #[rustfmt::skip]
    pub fn block_3d_eight_hex20() -> Mesh {
        Mesh {
            ndim: 3,
            points: vec![
                Point { id:  0, coords: vec![0.0, 0.0, 0.0] },
                Point { id:  1, coords: vec![1.0, 0.0, 0.0] },
                Point { id:  2, coords: vec![1.0, 1.0, 0.0] },
                Point { id:  3, coords: vec![0.0, 1.0, 0.0] },
                Point { id:  4, coords: vec![0.0, 0.0, 2.0] },
                Point { id:  5, coords: vec![1.0, 0.0, 2.0] },
                Point { id:  6, coords: vec![1.0, 1.0, 2.0] },
                Point { id:  7, coords: vec![0.0, 1.0, 2.0] },
                Point { id:  8, coords: vec![0.5, 0.0, 0.0] },
                Point { id:  9, coords: vec![1.0, 0.5, 0.0] },
                Point { id: 10, coords: vec![0.5, 1.0, 0.0] },
                Point { id: 11, coords: vec![0.0, 0.5, 0.0] },
                Point { id: 12, coords: vec![0.5, 0.0, 2.0] },
                Point { id: 13, coords: vec![1.0, 0.5, 2.0] },
                Point { id: 14, coords: vec![0.5, 1.0, 2.0] },
                Point { id: 15, coords: vec![0.0, 0.5, 2.0] },
                Point { id: 16, coords: vec![0.0, 0.0, 1.0] },
                Point { id: 17, coords: vec![1.0, 0.0, 1.0] },
                Point { id: 18, coords: vec![1.0, 1.0, 1.0] },
                Point { id: 19, coords: vec![0.0, 1.0, 1.0] },
                Point { id: 20, coords: vec![2.0, 0.0, 0.0] },
                Point { id: 21, coords: vec![2.0, 1.0, 0.0] },
                Point { id: 22, coords: vec![2.0, 0.0, 2.0] },
                Point { id: 23, coords: vec![2.0, 1.0, 2.0] },
                Point { id: 24, coords: vec![1.5, 0.0, 0.0] },
                Point { id: 25, coords: vec![2.0, 0.5, 0.0] },
                Point { id: 26, coords: vec![1.5, 1.0, 0.0] },
                Point { id: 27, coords: vec![1.5, 0.0, 2.0] },
                Point { id: 28, coords: vec![2.0, 0.5, 2.0] },
                Point { id: 29, coords: vec![1.5, 1.0, 2.0] },
                Point { id: 30, coords: vec![2.0, 0.0, 1.0] },
                Point { id: 31, coords: vec![2.0, 1.0, 1.0] },
                Point { id: 32, coords: vec![1.0, 2.0, 0.0] },
                Point { id: 33, coords: vec![0.0, 2.0, 0.0] },
                Point { id: 34, coords: vec![1.0, 2.0, 2.0] },
                Point { id: 35, coords: vec![0.0, 2.0, 2.0] },
                Point { id: 36, coords: vec![1.0, 1.5, 0.0] },
                Point { id: 37, coords: vec![0.5, 2.0, 0.0] },
                Point { id: 38, coords: vec![0.0, 1.5, 0.0] },
                Point { id: 39, coords: vec![1.0, 1.5, 2.0] },
                Point { id: 40, coords: vec![0.5, 2.0, 2.0] },
                Point { id: 41, coords: vec![0.0, 1.5, 2.0] },
                Point { id: 42, coords: vec![1.0, 2.0, 1.0] },
                Point { id: 43, coords: vec![0.0, 2.0, 1.0] },
                Point { id: 44, coords: vec![2.0, 2.0, 0.0] },
                Point { id: 45, coords: vec![2.0, 2.0, 2.0] },
                Point { id: 46, coords: vec![2.0, 1.5, 0.0] },
                Point { id: 47, coords: vec![1.5, 2.0, 0.0] },
                Point { id: 48, coords: vec![2.0, 1.5, 2.0] },
                Point { id: 49, coords: vec![1.5, 2.0, 2.0] },
                Point { id: 50, coords: vec![2.0, 2.0, 1.0] },
                Point { id: 51, coords: vec![0.0, 0.0, 4.0] },
                Point { id: 52, coords: vec![1.0, 0.0, 4.0] },
                Point { id: 53, coords: vec![1.0, 1.0, 4.0] },
                Point { id: 54, coords: vec![0.0, 1.0, 4.0] },
                Point { id: 55, coords: vec![0.5, 0.0, 4.0] },
                Point { id: 56, coords: vec![1.0, 0.5, 4.0] },
                Point { id: 57, coords: vec![0.5, 1.0, 4.0] },
                Point { id: 58, coords: vec![0.0, 0.5, 4.0] },
                Point { id: 59, coords: vec![0.0, 0.0, 3.0] },
                Point { id: 60, coords: vec![1.0, 0.0, 3.0] },
                Point { id: 61, coords: vec![1.0, 1.0, 3.0] },
                Point { id: 62, coords: vec![0.0, 1.0, 3.0] },
                Point { id: 63, coords: vec![2.0, 0.0, 4.0] },
                Point { id: 64, coords: vec![2.0, 1.0, 4.0] },
                Point { id: 65, coords: vec![1.5, 0.0, 4.0] },
                Point { id: 66, coords: vec![2.0, 0.5, 4.0] },
                Point { id: 67, coords: vec![1.5, 1.0, 4.0] },
                Point { id: 68, coords: vec![2.0, 0.0, 3.0] },
                Point { id: 69, coords: vec![2.0, 1.0, 3.0] },
                Point { id: 70, coords: vec![1.0, 2.0, 4.0] },
                Point { id: 71, coords: vec![0.0, 2.0, 4.0] },
                Point { id: 72, coords: vec![1.0, 1.5, 4.0] },
                Point { id: 73, coords: vec![0.5, 2.0, 4.0] },
                Point { id: 74, coords: vec![0.0, 1.5, 4.0] },
                Point { id: 75, coords: vec![1.0, 2.0, 3.0] },
                Point { id: 76, coords: vec![0.0, 2.0, 3.0] },
                Point { id: 77, coords: vec![2.0, 2.0, 4.0] },
                Point { id: 78, coords: vec![2.0, 1.5, 4.0] },
                Point { id: 79, coords: vec![1.5, 2.0, 4.0] },
                Point { id: 80, coords: vec![2.0, 2.0, 3.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Hex20, points: vec![0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Hex20, points: vec![1, 20, 21,  2,  5, 22, 23,  6, 24, 25, 26,  9, 27, 28, 29, 13, 17, 30, 31, 18] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Hex20, points: vec![3,  2, 32, 33,  7,  6, 34, 35, 10, 36, 37, 38, 14, 39, 40, 41, 19, 18, 42, 43] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Hex20, points: vec![2, 21, 44, 32,  6, 23, 45, 34, 26, 46, 47, 36, 29, 48, 49, 39, 18, 31, 50, 42] },
                Cell { id: 4, attribute_id: 1, kind: GeoKind::Hex20, points: vec![4,  5,  6,  7, 51, 52, 53, 54, 12, 13, 14, 15, 55, 56, 57, 58, 59, 60, 61, 62] },
                Cell { id: 5, attribute_id: 1, kind: GeoKind::Hex20, points: vec![5, 22, 23,  6, 52, 63, 64, 53, 27, 28, 29, 13, 65, 66, 67, 56, 60, 68, 69, 61] },
                Cell { id: 6, attribute_id: 1, kind: GeoKind::Hex20, points: vec![7,  6, 34, 35, 54, 53, 70, 71, 14, 39, 40, 41, 57, 72, 73, 74, 62, 61, 75, 76] },
                Cell { id: 7, attribute_id: 1, kind: GeoKind::Hex20, points: vec![6, 23, 45, 34, 53, 64, 77, 70, 29, 48, 49, 39, 67, 78, 79, 72, 61, 69, 80, 75] },
            ],
        }
    }

    /// Returns a mesh of a quarter of a ring with eight Qua8 (radius=1 thickness=1)
    ///
    /// ```text
    /// 2.0   14---36--,__11
    ///        |          `,-..33
    /// 1.75  24   [7]   22     `-,
    ///        |         ,  [5]    ,8.
    /// 1.5   13--35--10/        20   `*
    ///        |       ,`*32    ,'      30
    /// 1.25  23 [6] 21     *.7     [3]   *
    ///        |     ,  [4]  , *.          5
    /// 1.0   12-34-9      19    29     18' *
    ///              `31. ,' [2]   *  _,     *
    ///                  6.       _.4'        *
    ///                   28  _.17   *   [1]  27
    ///                     3'  [0]  26        *
    ///                     25        *        *
    ///        +             0---15---1---16---2
    ///
    ///                     1.0 1.25  1.5 1.75  2.0
    /// ```
    #[rustfmt::skip]
    pub fn ring_eight_qua8_rad1_thick1() -> Mesh {
        Mesh {
            ndim: 2,
            points: vec![
                Point { id:  0, coords: vec![1.000000000000000e+00, 0.000000000000000e+00] },
                Point { id:  1, coords: vec![1.500000000000000e+00, 0.000000000000000e+00] },
                Point { id:  2, coords: vec![2.000000000000000e+00, 0.000000000000000e+00] },
                Point { id:  3, coords: vec![9.127002149692066e-01, 4.086298050744269e-01] },
                Point { id:  4, coords: vec![1.357995128834866e+00, 6.079951288348657e-01] },
                Point { id:  5, coords: vec![1.825400429938413e+00, 8.172596101488537e-01] },
                Point { id:  6, coords: vec![7.071067811865476e-01, 7.071067811865475e-01] },
                Point { id:  7, coords: vec![1.060660171779821e+00, 1.060660171779821e+00] },
                Point { id:  8, coords: vec![1.414213562373095e+00, 1.414213562373095e+00] },
                Point { id:  9, coords: vec![4.086298050744270e-01, 9.127002149692066e-01] },
                Point { id: 10, coords: vec![6.079951288348662e-01, 1.357995128834866e+00] },
                Point { id: 11, coords: vec![8.172596101488541e-01, 1.825400429938413e+00] },
                Point { id: 12, coords: vec![6.123233995736766e-17, 1.000000000000000e+00] },
                Point { id: 13, coords: vec![9.184850993605148e-17, 1.500000000000000e+00] },
                Point { id: 14, coords: vec![1.224646799147353e-16, 2.000000000000000e+00] },
                Point { id: 15, coords: vec![1.250000000000000e+00, 0.000000000000000e+00] },
                Point { id: 16, coords: vec![1.750000000000000e+00, 0.000000000000000e+00] },
                Point { id: 17, coords: vec![1.131662607362388e+00, 5.066626073623883e-01] },
                Point { id: 18, coords: vec![1.584327650307344e+00, 7.093276503073436e-01] },
                Point { id: 19, coords: vec![8.838834764831844e-01, 8.838834764831842e-01] },
                Point { id: 20, coords: vec![1.237436867076458e+00, 1.237436867076458e+00] },
                Point { id: 21, coords: vec![5.066626073623884e-01, 1.131662607362388e+00] },
                Point { id: 22, coords: vec![7.093276503073440e-01, 1.584327650307344e+00] },
                Point { id: 23, coords: vec![7.654042494670958e-17, 1.250000000000000e+00] },
                Point { id: 24, coords: vec![1.071565949253934e-16, 1.750000000000000e+00] },
                Point { id: 25, coords: vec![9.759662299218728e-01, 2.179218163747866e-01] },
                Point { id: 26, coords: vec![1.448413825153672e+00, 3.234138251536718e-01] },
                Point { id: 27, coords: vec![1.951932459843746e+00, 4.358436327495732e-01] },
                Point { id: 28, coords: vec![8.212291819630956e-01, 5.705984846564395e-01] },
                Point { id: 29, coords: vec![1.228743911043582e+00, 8.537439110435825e-01] },
                Point { id: 30, coords: vec![1.642458363926191e+00, 1.141196969312879e+00] },
                Point { id: 31, coords: vec![5.705984846564395e-01, 8.212291819630954e-01] },
                Point { id: 32, coords: vec![8.537439110435825e-01, 1.228743911043582e+00] },
                Point { id: 33, coords: vec![1.141196969312879e+00, 1.642458363926191e+00] },
                Point { id: 34, coords: vec![2.179218163747867e-01, 9.759662299218727e-01] },
                Point { id: 35, coords: vec![3.234138251536719e-01, 1.448413825153672e+00] },
                Point { id: 36, coords: vec![4.358436327495734e-01, 1.951932459843745e+00] },
            ],
            cells: vec![
                Cell { id: 0, attribute_id: 1, kind: GeoKind::Qua8, points: vec![ 0, 1, 4, 3,15,26,17,25] },
                Cell { id: 1, attribute_id: 1, kind: GeoKind::Qua8, points: vec![ 1, 2, 5, 4,16,27,18,26] },
                Cell { id: 2, attribute_id: 1, kind: GeoKind::Qua8, points: vec![ 3, 4, 7, 6,17,29,19,28] },
                Cell { id: 3, attribute_id: 1, kind: GeoKind::Qua8, points: vec![ 4, 5, 8, 7,18,30,20,29] },
                Cell { id: 4, attribute_id: 1, kind: GeoKind::Qua8, points: vec![ 6, 7,10, 9,19,32,21,31] },
                Cell { id: 5, attribute_id: 1, kind: GeoKind::Qua8, points: vec![ 7, 8,11,10,20,33,22,32] },
                Cell { id: 6, attribute_id: 1, kind: GeoKind::Qua8, points: vec![ 9,10,13,12,21,35,23,34] },
                Cell { id: 7, attribute_id: 1, kind: GeoKind::Qua8, points: vec![10,11,14,13,22,36,24,35] },
            ],
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Samples;
    use crate::mesh::check_all;
    use crate::StrError;

    #[test]
    fn samples_work() -> Result<(), StrError> {
        let mesh = Samples::two_quads_horizontal();
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 6);
        assert_eq!(mesh.cells.len(), 2);
        check_all(&mesh)?;

        let mesh = Samples::two_cubes_vertical();
        assert_eq!(mesh.ndim, 3);
        assert_eq!(mesh.points.len(), 12);
        assert_eq!(mesh.cells.len(), 2);
        check_all(&mesh)?;

        let mesh = Samples::four_cubes_wall();
        assert_eq!(mesh.ndim, 3);
        assert_eq!(mesh.points.len(), 18);
        assert_eq!(mesh.cells.len(), 4);
        check_all(&mesh)?;

        let mesh = Samples::mixed_shapes_2d();
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 6);
        assert_eq!(mesh.cells.len(), 3);
        check_all(&mesh)?;

        let mesh = Samples::mixed_shapes_3d();
        assert_eq!(mesh.ndim, 3);
        assert_eq!(mesh.points.len(), 13);
        assert_eq!(mesh.cells.len(), 5);
        check_all(&mesh)?;

        let mesh = Samples::block_2d_four_qua4();
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 9);
        assert_eq!(mesh.cells.len(), 4);
        check_all(&mesh)?;

        let mesh = Samples::block_2d_four_qua8();
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 21);
        assert_eq!(mesh.cells.len(), 4);
        check_all(&mesh)?;

        let mesh = Samples::block_2d_four_qua9();
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 25);
        assert_eq!(mesh.cells.len(), 4);
        check_all(&mesh)?;

        let mesh = Samples::block_2d_four_qua12();
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 33);
        assert_eq!(mesh.cells.len(), 4);
        check_all(&mesh)?;

        let mesh = Samples::block_2d_four_qua16();
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 49);
        assert_eq!(mesh.cells.len(), 4);
        check_all(&mesh)?;

        let mesh = Samples::block_2d_four_qua17();
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 49);
        assert_eq!(mesh.cells.len(), 4);
        check_all(&mesh)?;

        let mesh = Samples::block_3d_eight_hex8();
        assert_eq!(mesh.ndim, 3);
        assert_eq!(mesh.points.len(), 27);
        assert_eq!(mesh.cells.len(), 8);
        check_all(&mesh)?;

        let mesh = Samples::block_3d_eight_hex20();
        assert_eq!(mesh.ndim, 3);
        assert_eq!(mesh.points.len(), 81);
        assert_eq!(mesh.cells.len(), 8);
        check_all(&mesh)?;

        let mesh = Samples::ring_eight_qua8_rad1_thick1();
        assert_eq!(mesh.ndim, 2);
        assert_eq!(mesh.points.len(), 37);
        assert_eq!(mesh.cells.len(), 8);
        check_all(&mesh)?;
        Ok(())
    }
}
