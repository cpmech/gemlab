use gemlab::mesh::{Features, Samples};
use gemlab::StrError;
use std::collections::HashMap;

#[test]
fn test_normal_vectors_2d() -> Result<(), StrError> {
    //  3---------2---------5
    //  |         |         |
    //  |   [0]   |   [1]   |
    //  |         |         |
    //  0---------1---------4
    let mesh = Samples::two_qua4();
    let features = Features::new(&mesh, false);

    // the magnitude (l) of the normal vector should be equal to
    // 0.5 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
    let l = 0.5; // magnitude of normal vector

    // edge keys and correct normal vectors (solutions)
    let solutions = HashMap::from([
        ((0, 1), (l, [0.0, -1.0])), // bottom
        ((1, 4), (l, [0.0, -1.0])), // bottom
        ((4, 5), (l, [1.0, 0.0])),  // right
        ((2, 3), (l, [0.0, 1.0])),  // top
        ((2, 5), (l, [0.0, 1.0])),  // top
        ((0, 3), (l, [-1.0, 0.0])), // left
    ]);

    // check if the normal vectors at boundary are outward
    mesh.check_2d_edge_normals(&features.edges, &solutions, 1e-15)
        .expect("ok");
    Ok(())
}

#[test]
fn test_normal_vectors_2d_qua8() -> Result<(), StrError> {
    // 14------16------13------20------18
    //  |               |               |
    //  |               |               |
    // 17              15              19
    //  |               |               |
    //  |               |               |
    //  3-------6-------2------12-------9
    //  |               |               |
    //  |               |               |
    //  7               5              11
    //  |               |               |
    //  |               |               |
    //  0-------4-------1------10-------8
    let mesh = Samples::block_2d_four_qua8();
    let features = Features::new(&mesh, false);

    // the magnitude (l) of the normal vector should be equal to
    // 0.5 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
    let l = 0.5; // magnitude of normal vector

    // edge keys and correct normal vectors (solutions)
    let solutions = HashMap::from([
        ((0, 1), (l, [0.0, -1.0])),  // bottom
        ((1, 8), (l, [0.0, -1.0])),  // bottom
        ((8, 9), (l, [1.0, 0.0])),   // right
        ((9, 18), (l, [1.0, 0.0])),  // right
        ((13, 14), (l, [0.0, 1.0])), // top
        ((13, 18), (l, [0.0, 1.0])), // top
        ((0, 3), (l, [-1.0, 0.0])),  // left
        ((3, 14), (l, [-1.0, 0.0])), // left
    ]);

    // check if the normal vectors at boundary are outward
    mesh.check_2d_edge_normals(&features.edges, &solutions, 1e-15)
        .expect("ok");
    Ok(())
}

#[test]
fn test_normal_vectors_2d_qua9() -> Result<(), StrError> {
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
    let mesh = Samples::block_2d_four_qua9();
    let features = Features::new(&mesh, false);

    // the magnitude (l) of the normal vector should be equal to
    // 0.5 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
    let l = 0.5; // magnitude of normal vector

    // edge keys and correct normal vectors (solutions)
    let solutions = HashMap::from([
        ((0, 1), (l, [0.0, -1.0])),  // bottom
        ((1, 9), (l, [0.0, -1.0])),  // bottom
        ((9, 10), (l, [1.0, 0.0])),  // right
        ((10, 21), (l, [1.0, 0.0])), // right
        ((15, 16), (l, [0.0, 1.0])), // top
        ((15, 21), (l, [0.0, 1.0])), // top
        ((0, 3), (l, [-1.0, 0.0])),  // left
        ((3, 16), (l, [-1.0, 0.0])), // left
    ]);

    // check if the normal vectors at boundary are outward
    mesh.check_2d_edge_normals(&features.edges, &solutions, 1e-15)
        .expect("ok");
    Ok(())
}

#[test]
fn test_normal_vectors_2d_qua12() -> Result<(), StrError> {
    // 21---26---23----20---32---30----28
    //  |               |               |
    // 24              25              31
    //  |               |               |
    // 27              22              29
    //  |               |               |
    //  3---10-----6----2---19---16----13
    //  |               |               |
    //  7               9              18
    //  |               |               |
    // 11               5              15
    //  |               |               |
    //  0----4-----8----1---14---17----12
    let mesh = Samples::block_2d_four_qua12();
    let features = Features::new(&mesh, false);

    // the magnitude (l) of the normal vector should be equal to
    // 0.75 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
    let l = 0.75; // magnitude of normal vector

    // edge keys and correct normal vectors (solutions)
    let solutions = HashMap::from([
        ((0, 1), (l, [0.0, -1.0])),  // bottom
        ((1, 12), (l, [0.0, -1.0])), // bottom
        ((12, 13), (l, [1.0, 0.0])), // right
        ((13, 28), (l, [1.0, 0.0])), // right
        ((20, 21), (l, [0.0, 1.0])), // top
        ((20, 28), (l, [0.0, 1.0])), // top
        ((0, 3), (l, [-1.0, 0.0])),  // left
        ((3, 21), (l, [-1.0, 0.0])), // left
    ]);

    // check if the normal vectors at boundary are outward
    mesh.check_2d_edge_normals(&features.edges, &solutions, 1e-15)
        .expect("ok");
    Ok(())
}

#[test]
fn test_normal_vectors_2d_qua16() -> Result<(), StrError> {
    // 29---34----31---28---44---42----40
    //  |               |               |
    // 32   39    38   33   48   47    43
    //  |               |               |
    // 35   36    37   30   45   46    41
    //  |               |               |
    //  3---10-----6----2---23---20----17
    //  |               |               |
    //  7   15    14    9   27   26    22
    //  |               |               |
    // 11   12    13    5   24   25    19
    //  |               |               |
    //  0----4-----8----1---18---21----16
    let mesh = Samples::block_2d_four_qua16();
    let features = Features::new(&mesh, false);

    // the magnitude (l) of the normal vector should be equal to
    // 0.75 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
    let l = 0.75; // magnitude of normal vector

    // edge keys and correct normal vectors (solutions)
    let solutions = HashMap::from([
        ((0, 1), (l, [0.0, -1.0])),  // bottom
        ((1, 16), (l, [0.0, -1.0])), // bottom
        ((16, 17), (l, [1.0, 0.0])), // right
        ((17, 40), (l, [1.0, 0.0])), // right
        ((28, 29), (l, [0.0, 1.0])), // top
        ((28, 40), (l, [0.0, 1.0])), // top
        ((0, 3), (l, [-1.0, 0.0])),  // left
        ((3, 29), (l, [-1.0, 0.0])), // left
    ]);

    // check if the normal vectors at boundary are outward
    mesh.check_2d_edge_normals(&features.edges, &solutions, 1e-15)
        .expect("ok");
    Ok(())
}

#[test]
fn test_normal_vectors_2d_qua17() -> Result<(), StrError> {
    // 30---38---35---32---29---47---45---43---41
    //  |                   |                   |
    // 33                  37                  46
    //  |                   |                   |
    // 36        40        34        48        44
    //  |                   |                   |
    // 39                  31                  42
    //  |                   |                   |
    //  3---14---10----6----2---27---24---21---18
    //  |                   |                   |
    //  7                  13                  26
    //  |                   |                   |
    // 11        16         9        28        23
    //  |                   |                   |
    // 15                   5                  20
    //  |                   |                   |
    //  0----4----8---12----1---19---22---25---17
    let mesh = Samples::block_2d_four_qua17();
    let features = Features::new(&mesh, false);

    // the magnitude (l) of the normal vector should be equal to
    // 1.0 = edge_length / 2.0 where 2.0 corresponds to the edge_length in the reference system
    let l = 1.0; // magnitude of normal vector

    // edge keys and correct normal vectors (solutions)
    let solutions = HashMap::from([
        ((0, 1), (l, [0.0, -1.0])),  // bottom
        ((1, 17), (l, [0.0, -1.0])), // bottom
        ((17, 18), (l, [1.0, 0.0])), // right
        ((18, 41), (l, [1.0, 0.0])), // right
        ((29, 30), (l, [0.0, 1.0])), // top
        ((29, 41), (l, [0.0, 1.0])), // top
        ((0, 3), (l, [-1.0, 0.0])),  // left
        ((3, 30), (l, [-1.0, 0.0])), // left
    ]);

    // check if the normal vectors at boundary are outward
    mesh.check_2d_edge_normals(&features.edges, &solutions, 1e-15)
        .expect("ok");
    Ok(())
}

#[test]
fn test_normal_vectors_3d() -> Result<(), StrError> {
    //      8-------------11
    //     /.             /|
    //    / .            / |
    //   /  .           /  |
    //  /   .          /   |
    // 9-------------10    |
    // |    .         |    |
    // |    4---------|----7
    // |   /.         |   /|
    // |  / .         |  / |
    // | /  .         | /  |
    // |/   .         |/   |
    // 5--------------6    |
    // |    .         |    |
    // |    0---------|----3
    // |   /          |   /
    // |  /           |  /
    // | /            | /
    // |/             |/
    // 1--------------2
    let mesh = Samples::two_hex8();
    let features = Features::new(&mesh, false);

    // the magnitude (l) of the normal vector should be equal to
    // 0.25 = face_area / 4.0 where 4.0 corresponds to the face_area in the reference system
    let l = 0.25; // magnitude of normal vector

    // face keys and correct normal vectors (solutions)
    let solutions = HashMap::from([
        ((0, 3, 4, 7), (l, [-1.0, 0.0, 0.0])),  // behind
        ((4, 7, 8, 11), (l, [-1.0, 0.0, 0.0])), // behind
        ((1, 2, 5, 6), (l, [1.0, 0.0, 0.0])),   // front
        ((5, 6, 9, 10), (l, [1.0, 0.0, 0.0])),  // front
        ((0, 1, 4, 5), (l, [0.0, -1.0, 0.0])),  // left
        ((4, 5, 8, 9), (l, [0.0, -1.0, 0.0])),  // left
        ((2, 3, 6, 7), (l, [0.0, 1.0, 0.0])),   // right
        ((6, 7, 10, 11), (l, [0.0, 1.0, 0.0])), // right
        ((0, 1, 2, 3), (l, [0.0, 0.0, -1.0])),  // bottom
        ((8, 9, 10, 11), (l, [0.0, 0.0, 1.0])), // top
    ]);

    // check if the normal vectors at boundary are outward
    mesh.check_face_normals(&features.faces, &solutions, 1e-15).expect("ok");
    Ok(())
}

#[test]
fn test_normal_vectors_3d_hex20() -> Result<(), StrError> {
    //              51--------58--------54--------74--------71
    //              /.                  /.                  /|
    //             / .                 / .                 / |
    //           55  .               57  .               73  |
    //           /   .               /   .               /   |
    //          /    .              /    .              /    |
    //        52--------56--------53--------72--------70     |
    //        /.     .            /.     .            /|     |
    //       / .    59           / .    62           / |    76
    //     65  .     .         67  .     .         79  |     |
    //     /   .     .         /   .     .         /   |     |
    //    /    .     .        /    .     .        /    |     |
    //  63========66========64========78========77     |     |
    //   |     .     .       |     .     .       |     |     |
    //   |    60     .       |    61     .       |    75     |
    //   |     .     4 - - - |15 - . - - 7 - - - |41 - | - -35
    //   |     .    /.       |     .    /.       |     |    /|
    //   |     .   / .       |     .   / .       |     |   / |
    //   |     . 12  .       |     . 14  .       |     | 40  |
    //   |     . /   .       |     . /   .       |     | /   |
    //  68     ./    .      69     ./    .      80     |/    |
    //   |     5 - - - -13 - | - - 6 - - - -39 - | - -34     |
    //   |    /.     .       |    /.     .       |    /|     |
    //   |   / .    16       |   / .    19       |   / |    43
    //   | 27  .     .       | 29  .     .       | 49  |     |
    //   | /   .     .       | /   .     .       | /   |     |
    //   |/    .     .       |/    .     .       |/    |     |
    //  22========28========23========48========45     |     |
    //   |     .     .       |     .     .       |     |     |
    //   |    17     .       |    18     .       |    42     |
    //   |     .     0 - - - |11 - . - - 3 - - - |38 - | - -33
    //   |     .    /        |     .    /        |     |    /
    //   |     .   /         |     .   /         |     |   /
    //   |     .  8          |     . 10          |     | 37
    //   |     . /           |     . /           |     | /
    //  30     ./           31     ./           50     |/
    //   |     1 - - - - 9 - | - - 2 - - - -36 - | - -32
    //   |    /              |    /              |    /
    //   |   /               |   /               |   /
    //   | 24                | 26                | 47
    //   | /                 | /                 | /
    //   |/                  |/                  |/
    //  20========25========21========46========44
    let mesh = Samples::block_3d_eight_hex20();
    let features = Features::new(&mesh, false);

    // the magnitude (l) of the normal vector should be equal to
    // face_area / 4.0 where 4.0 corresponds to the face_area in the reference system

    // face keys and correct normal vectors (solutions)
    let solutions = HashMap::from([
        ((0, 1, 2, 3), (0.25, [0.0, 0.0, -1.0])),
        ((52, 53, 63, 64), (0.25, [0.0, 0.0, 1.0])),
        ((4, 7, 51, 54), (0.5, [-1.0, 0.0, 0.0])),
        ((34, 45, 70, 77), (0.5, [0.0, 1.0, 0.0])),
    ]);

    // check if the normal vectors at boundary are outward
    mesh.check_face_normals(&features.faces, &solutions, 1e-15).expect("ok");
    Ok(())
}
