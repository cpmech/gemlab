use gemlab::mesh::Block;
use gemlab::StrError;

fn main() -> Result<(), StrError> {
    #[rustfmt::skip]
    let mut block = Block::new(&[
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 2.0],
        [0.0, 2.0],
    ])?;
    let mesh = block.subdivide(4)?;

    // will produce the following mesh:
    //
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

    println!("{}", mesh);
    /*
    # header
    # space_ndim npoint ncell
    2 25 4

    # points
    # id x y {z}
    0 0 0
    1 1 0
    2 1 1
    3 0 1
    4 0.5 0
    5 1 0.5
    6 0.5 1
    7 0 0.5
    8 0.5 0.5
    9 2 0
    10 2 1
    11 1.5 0
    12 2 0.5
    13 1.5 1
    14 1.5 0.5
    15 1 2
    16 0 2
    17 1 1.5
    18 0.5 2
    19 0 1.5
    20 0.5 1.5
    21 2 2
    22 2 1.5
    23 1.5 2
    24 1.5 1.5

    # cells
    # id attribute_id geo_ndim nnode  points
    0 1 2 9  0 1 2 3 4 5 6 7 8
    1 1 2 9  1 9 10 2 11 12 13 5 14
    2 1 2 9  3 2 15 16 6 17 18 19 20
    3 1 2 9  2 10 21 15 13 22 23 17 24
    */

    Ok(())
}
