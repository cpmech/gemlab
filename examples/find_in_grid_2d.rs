use gemlab::geometry::Circle;
use gemlab::util::{GridSearch, GsNdiv, GsTol, PI};
use gemlab::StrError;
use std::collections::HashSet;

struct Segment {
    a: [f64; 2],
    b: [f64; 2],
}

fn main() -> Result<(), StrError> {
    let min = &[0.0, 0.0];
    let max = &[10.0, 10.0];
    let mut grid = GridSearch::new(min, max, GsNdiv::Default, GsTol::Default)?;

    // number of points in each entity
    const NPOINT: usize = 12;

    // circles
    let circles = vec![
        Circle {
            center: [5.0, 2.0],
            radius: 1.0,
        },
        Circle {
            center: [7.5, 5.0],
            radius: 1.2,
        },
        Circle {
            center: [5.0, 7.5],
            radius: 1.4,
        },
        Circle {
            center: [2.5, 5.0],
            radius: 1.2,
        },
    ];

    // add points to grid
    let mut id = 0;
    let mut x = vec![0.0; 2];
    for c in &circles {
        for n in 0..NPOINT {
            let alpha = (n as f64) * 2.0 * PI / (NPOINT as f64);
            x[0] = c.center[0] + c.radius * f64::cos(alpha);
            x[1] = c.center[1] + c.radius * f64::sin(alpha);
            grid.insert(id, &x)?;
            id += 1;
        }
    }

    // segments
    let segments = vec![
        Segment {
            a: [0.0, 0.0],
            b: [10.0, 10.0],
        },
        Segment {
            a: [0.0, 10.0],
            b: [10.0, 0.0],
        },
    ];
    for s in &segments {
        for n in 0..NPOINT {
            x[0] = s.a[0] + (n as f64) * (s.b[0] - s.a[0]) / ((NPOINT - 1) as f64);
            x[1] = s.a[1] + (n as f64) * (s.b[1] - s.a[1]) / ((NPOINT - 1) as f64);
            grid.insert(id, &x)?;
            id += 1;
        }
    }

    // plot
    let mut plot = grid.plot()?;
    plot.set_equal_axes(true)
        .set_figure_size_points(600.0, 600.0)
        .save("/tmp/gemlab/find_in_grid_2d.svg")?;

    // find points on circles
    let mut start_id = 0;
    let mut end_id = NPOINT;
    for c in &circles {
        let res = grid.find_on_circle(&c.center, c.radius)?;
        check(&res, &(start_id..end_id).collect::<Vec<_>>());
        start_id += NPOINT;
        end_id += NPOINT;
    }

    // find points on segments
    for s in &segments {
        let res = grid.find_on_line(&s.a, &s.b)?;
        check(&res, &(start_id..end_id).collect::<Vec<_>>());
        start_id += NPOINT;
        end_id += NPOINT;
    }

    Ok(())
}

fn check<T>(found: &HashSet<T>, correct: &[T])
where
    T: Copy + Ord + std::fmt::Debug,
{
    let mut ids: Vec<T> = found.iter().copied().collect();
    ids.sort();
    assert_eq!(ids, correct);
}
