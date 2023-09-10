use gemlab::mesh::{draw_mesh, Block};
use gemlab::shapes::GeoKind;
use gemlab::StrError;

fn main() -> Result<(), StrError> {
    #[rustfmt::skip]
    let mut block = Block::new(&[
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 2.0],
        [0.0, 2.0],
    ])?;
    let mesh = block.subdivide(GeoKind::Qua8)?;

    // will produce the following mesh:
    //
    // 2.0  14------16------13------20------18
    //       |               |               |
    //       |               |               |
    // 1.5  17      [2]     15      [3]     19
    //       |               |               |
    //       |               |               |
    // 1.0   3-------6-------2------12-------9
    //       |               |               |
    //       |               |               |
    // 0.5   7      [0]      5      [1]     11
    //       |               |               |
    //       |               |               |
    // 0.0   0-------4-------1------10-------8
    //
    //      0.0     0.5     1.0     1.5     2.0

    println!("{}", mesh);

    let correct = "# header
# ndim npoint ncell
2 21 4

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
8 2 0
9 2 1
10 1.5 0
11 2 0.5
12 1.5 1
13 1 2
14 0 2
15 1 1.5
16 0.5 2
17 0 1.5
18 2 2
19 2 1.5
20 1.5 2

# cells
# id att kind  points
0 1 qua8  0 1 2 3 4 5 6 7
1 1 qua8  1 8 9 2 10 11 12 5
2 1 qua8  3 2 13 14 6 15 16 17
3 1 qua8  2 9 18 13 12 19 20 15
";

    assert_eq!(format!("{}", mesh), correct);
    draw_mesh(&mesh, true, false, false, "/tmp/gemlab/example_block_2d_qua8.svg")?;
    Ok(())
}
