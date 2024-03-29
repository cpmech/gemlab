use gemlab::geometry::Circle;
use gemlab::util::{any_x, GridSearch};
use gemlab::StrError;
use plotpy::Plot;
use russell_lab::math::PI;

#[test]
fn test_grid_search_2d() -> Result<(), StrError> {
    let xmin = &[0.0, 0.0];
    let xmax = &[1.0, 2.0];
    let mut grid = GridSearch::new(xmin, xmax, None, None, Some(0.0))?;

    // id
    let mut id = 0;

    // circle at y=0.5
    let bot = Circle {
        center: [0.5, 0.5],
        radius: 0.3,
    };
    let npoint = 12;
    let mut x = vec![0.0; 2];
    for n in 0..npoint {
        let alpha = (n as f64) * 2.0 * PI / (npoint as f64);
        x[0] = bot.center[0] + bot.radius * f64::cos(alpha);
        x[1] = bot.center[1] + bot.radius * f64::sin(alpha);
        grid.insert(id, &x)?;
        id += 1;
    }

    // circle at y=1.5
    let top = Circle {
        center: [0.5, 1.5],
        radius: 0.2,
    };
    for n in 0..npoint {
        let alpha = (n as f64) * 2.0 * PI / (npoint as f64);
        x[0] = top.center[0] + top.radius * f64::cos(alpha);
        x[1] = top.center[1] + top.radius * f64::sin(alpha);
        grid.insert(id, &x)?;
        id += 1;
    }

    // line
    let a = &[0.0, 0.5];
    let b = &[1.0, 1.5];
    for n in 0..npoint {
        x[0] = a[0] + (n as f64) * (b[0] - a[0]) / ((npoint - 1) as f64);
        x[1] = a[1] + (n as f64) * (b[1] - a[1]) / ((npoint - 1) as f64);
        grid.insert(id, &x)?;
        id += 1;
    }

    // draw grid
    if false {
        let mut plot = Plot::new();
        grid.draw(&mut plot, false)?;
        plot.set_equal_axes(true)
            .set_xrange(-0.2, 1.2)
            .set_figure_size_points(400.0, 800.0)
            .save("/tmp/gemlab/test_grid_search_2d.svg")?;
    }

    // search points on bot circle
    let res = grid.search_on_circle(&bot.center, bot.radius, any_x)?;
    let mut points: Vec<_> = res.iter().copied().collect();
    points.sort();
    let correct = (0..npoint).collect::<Vec<_>>();
    assert_eq!(points, correct);

    // search points on top circle
    let res = grid.search_on_circle(&top.center, top.radius, any_x)?;
    let mut points: Vec<_> = res.iter().copied().collect();
    points.sort();
    let correct = (npoint..npoint * 2).collect::<Vec<_>>();
    assert_eq!(points, correct);

    // search points on line
    let res = grid.search_on_line(a, b, any_x)?;
    let mut points: Vec<_> = res.iter().copied().collect();
    points.sort();
    let correct = (npoint * 2..npoint * 3).collect::<Vec<_>>();
    assert_eq!(points, correct);
    Ok(())
}
