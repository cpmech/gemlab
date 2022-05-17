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

    /// Returns a mesh with four quad4 cells
    ///
    /// ```text
    ///   7---------------6---------------8
    ///   |               |               |
    ///   |               |               |
    ///   |      [2]      |      [3]      |
    ///   |               |               |
    ///   |               |               |
    ///   3---------------2---------------5
    ///   |               |               |
    ///   |               |               |
    ///   |      [0]      |      [1]      |
    ///   |               |               |
    ///   |               |               |
    ///   0---------------1---------------4
    ///
    /// xmin = 0.0, xmax = 2.0
    /// ymin = 0.0, ymax = 2.0
    /// ```
    #[rustfmt::skip]
    pub fn block_2d_four_quad4() -> Mesh {
        Mesh {
            space_ndim: 2,
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
                Cell { id: 0, attribute_id: 1, geo_ndim: 2, points: vec![0, 1, 2, 3] },
                Cell { id: 1, attribute_id: 1, geo_ndim: 2, points: vec![1, 4, 5, 2] },
                Cell { id: 2, attribute_id: 1, geo_ndim: 2, points: vec![3, 2, 6, 7] },
                Cell { id: 3, attribute_id: 1, geo_ndim: 2, points: vec![2, 5, 8, 6] },
            ],
        }
    }

    /// Returns a mesh with four quad8 cells
    ///
    /// ```text
    ///  14------16------13------20------18
    ///   |               |               |
    ///   |               |               |
    ///  17      [2]     15      [3]     19
    ///   |               |               |
    ///   |               |               |
    ///   3-------6-------2------12-------9
    ///   |               |               |
    ///   |               |               |
    ///   7      [0]      5      [1]     11
    ///   |               |               |
    ///   |               |               |
    ///   0-------4-------1------10-------8
    ///
    /// xmin = 0.0, xmax = 2.0
    /// ymin = 0.0, ymax = 2.0
    /// ```
    #[rustfmt::skip]
    pub fn block_2d_four_quad8() -> Mesh {
        Mesh {
            space_ndim: 2,
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
                Cell { id: 0, attribute_id: 1, geo_ndim: 2, points: vec![0, 1,  2,  3,  4,  5,  6,  7] },
                Cell { id: 1, attribute_id: 1, geo_ndim: 2, points: vec![1, 8,  9,  2, 10, 11, 12,  5] },
                Cell { id: 2, attribute_id: 1, geo_ndim: 2, points: vec![3, 2, 13, 14,  6, 15, 16, 17] },
                Cell { id: 3, attribute_id: 1, geo_ndim: 2, points: vec![2, 9, 18, 13, 12, 19, 20, 15] },
            ],
        }
    }

    /// Returns a mesh with four quad9 cells
    ///
    /// ```text
    /// cell ids as in block_2d_four_quad8
    ///
    ///  16------18------15------23------21
    ///   |               |               |
    ///   |               |               |
    ///  19      20      17      24      22
    ///   |               |               |
    ///   |               |               |
    ///   3-------6-------2------13------10
    ///   |               |               |
    ///   |               |               |
    ///   7       8       5      14      12
    ///   |               |               |
    ///   |               |               |
    ///   0-------4-------1------11-------9
    ///
    /// xmin = 0.0, xmax = 2.0
    /// ymin = 0.0, ymax = 2.0
    /// ```
    #[rustfmt::skip]
    pub fn block_2d_four_quad9() -> Mesh {
        Mesh {
            space_ndim: 2,
            points: vec![
                Point { id:  0, coords: vec![] },
                Point { id:  1, coords: vec![] },
                Point { id:  2, coords: vec![] },
                Point { id:  3, coords: vec![] },
                Point { id:  4, coords: vec![] },
                Point { id:  5, coords: vec![] },
                Point { id:  6, coords: vec![] },
                Point { id:  7, coords: vec![] },
                Point { id:  8, coords: vec![] },
                Point { id:  9, coords: vec![] },
                Point { id: 10, coords: vec![] },
                Point { id: 11, coords: vec![] },
                Point { id: 12, coords: vec![] },
                Point { id: 13, coords: vec![] },
                Point { id: 14, coords: vec![] },
                Point { id: 15, coords: vec![] },
                Point { id: 16, coords: vec![] },
                Point { id: 17, coords: vec![] },
                Point { id: 18, coords: vec![] },
                Point { id: 19, coords: vec![] },
                Point { id: 20, coords: vec![] },
                Point { id: 21, coords: vec![] },
                Point { id: 22, coords: vec![] },
                Point { id: 23, coords: vec![] },
                Point { id: 24, coords: vec![] },
            ],
            cells: vec![
                Cell { id: 4, attribute_id: 3, geo_ndim: 1, points: vec![] },
            ],
        }
    }

    /// Returns a mesh with four quad12 cells
    ///
    /// ```text
    /// cell ids as in block_2d_four_quad8
    ///
    ///  21---26---23----20---32---30----28
    ///   |               |               |
    ///  24              25              31
    ///   |               |               |
    ///  27              22              29
    ///   |               |               |
    ///   3---10-----6----2---19---16----13
    ///   |               |               |
    ///   7               9              18
    ///   |               |               |
    ///  11               5              15
    ///   |               |               |
    ///   0----4-----8----1---14---17----12
    ///
    /// xmin = 0.0, xmax = 3.0
    /// ymin = 0.0, ymax = 3.0
    /// ```
    #[rustfmt::skip]
    pub fn block_2d_four_quad12() -> Mesh {
        Mesh {
            space_ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![] },
            ],
            cells: vec![
                Cell { id: 4, attribute_id: 3, geo_ndim: 1, points: vec![] },
            ],
        }
    }

    /// Returns a mesh with four quad16 cells
    ///
    /// ```text
    /// cell ids as in block_2d_four_quad8
    ///
    ///  29---34----31---28---44---42----40
    ///   |               |               |
    ///  32   39    38   33   48   47    43
    ///   |               |               |
    ///  35   36    37   30   45   46    41
    ///   |               |               |
    ///   3---10-----6----2---23---20----17
    ///   |               |               |
    ///   7   15    14    9   27   26    22
    ///   |               |               |
    ///  11   12    13    5   24   25    19
    ///   |               |               |
    ///   0----4-----8----1---18---21----16
    ///
    /// xmin = 0.0, xmax = 3.0
    /// ymin = 0.0, ymax = 3.0
    /// ```
    #[rustfmt::skip]
    pub fn block_2d_four_quad16() -> Mesh {
        Mesh {
            space_ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![] },
            ],
            cells: vec![
                Cell { id: 4, attribute_id: 3, geo_ndim: 1, points: vec![] },
            ],
        }
    }

    /// Returns a mesh with four quad17 cells
    ///
    /// ```text
    /// cell ids as in block_2d_four_quad8
    ///
    ///  30---38---35---32---29---47---45---43---41
    ///   |                   |                   |
    ///  33                  37                  46
    ///   |                   |                   |
    ///  36        40        34        48        44
    ///   |                   |                   |
    ///  39                  31                  42
    ///   |                   |                   |
    ///   3---14---10----6----2---27---24---21---18
    ///   |                   |                   |
    ///   7                  13                  26
    ///   |                   |                   |
    ///  11        16         9        28        23
    ///   |                   |                   |
    ///  15                   5                  20
    ///   |                   |                   |
    ///   0----4----8---12----1---19---22---25---17
    ///
    /// xmin = 0.0, xmax = 4.0
    /// ymin = 0.0, ymax = 4.0
    /// ```
    #[rustfmt::skip]
    pub fn block_2d_four_quad13() -> Mesh {
        Mesh {
            space_ndim: 2,
            points: vec![
                Point { id: 0, coords: vec![] },
            ],
            cells: vec![
                Cell { id: 4, attribute_id: 3, geo_ndim: 1, points: vec![] },
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
            space_ndim: 3,
            points: vec![
                Point { id: 0, coords: vec![] },
            ],
            cells: vec![
                Cell { id: 4, attribute_id: 3, geo_ndim: 1, points: vec![] },
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
            space_ndim: 3,
            points: vec![
                Point { id: 0, coords: vec![] },
            ],
            cells: vec![
                Cell { id: 4, attribute_id: 3, geo_ndim: 1, points: vec![] },
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
