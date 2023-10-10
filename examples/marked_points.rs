use gemlab::prelude::*;
use gemlab::StrError;

fn main() -> Result<(), StrError> {
    // ```text
    // -500        -400
    // 8------7------6._
    // |       [3](3)|  '-.5
    // |  [0]        |     '-._
    // 9  (1)       10  [1]    '4 -300
    // |             |  (2)  .-'
    // |       [2](3)|   _.3'
    // 0------1------2.-'
    // -100         -200
    // ```
    let h = 0.866; // ~ SQRT_3 / 2
    let m = h / 2.0;
    #[rustfmt::skip]
    let mesh = Mesh {
        ndim: 2,
        points: vec![
            Point { id:  0, marker: -100, coords: vec![0.0,   0.0 ] },
            Point { id:  1, marker:    0, coords: vec![0.5,   0.0 ] },
            Point { id:  2, marker: -200, coords: vec![1.0,   0.0 ] },
            Point { id:  3, marker:    0, coords: vec![1.0+m, 0.25] },
            Point { id:  4, marker: -300, coords: vec![1.0+h, 0.5 ] },
            Point { id:  5, marker:    0, coords: vec![1.0+m, 0.75] },
            Point { id:  6, marker: -400, coords: vec![1.0,   1.0 ] },
            Point { id:  7, marker:    0, coords: vec![0.5,   1.0 ] },
            Point { id:  8, marker: -500, coords: vec![0.0,   1.0 ] },
            Point { id:  9, marker:    0, coords: vec![0.0,   0.5 ] },
            Point { id: 10, marker:    0, coords: vec![1.0,   0.5 ] },
        ],
        cells: vec![
            Cell { id: 0, attribute: 1, kind: GeoKind::Qua8, points: vec![0, 2, 6, 8, 1, 10, 7, 9] },
            Cell { id: 1, attribute: 2, kind: GeoKind::Tri6, points: vec![2, 4, 6, 3, 5, 10] },
            Cell { id: 2, attribute: 3, kind: GeoKind::Lin2, points: vec![2, 10] },
            Cell { id: 3, attribute: 3, kind: GeoKind::Lin2, points: vec![10, 6] },
        ],
    };
    assert_eq!(mesh.search_marked_points(-200, |_| true)?, &[2]);
    assert_eq!(mesh.search_marked_points(-400, |_| true)?, &[6]);
    assert_eq!(mesh.search_marked_points(0, |x| x[1] > 0.49 && x[1] < 0.51)?, &[9, 10]);
    Ok(())
}
