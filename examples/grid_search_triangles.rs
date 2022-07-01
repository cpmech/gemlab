use gemlab::geometry::is_point_inside_triangle;
use gemlab::util::GridSearchCell;
use gemlab::StrError;
use plotpy::{Canvas, Plot, PolyCode, Text};

fn brute_force_search(triangles: &[[[f64; 2]; 3]], x: &[f64]) -> Option<usize> {
    for i in 0..triangles.len() {
        let t = &triangles[i];
        if is_point_inside_triangle(&t[0], &t[1], &t[2], x) {
            return Some(i);
        }
    }
    None
}

fn main() -> Result<(), StrError> {
    // [num_triangle][nnode=3][ndim=2]
    #[rustfmt::skip]
    const TRIS: [[[f64; 2]; 3]; 8] = [
        [[0.0, 0.0],  [1.0, 0.0],  [0.5, 0.85]],
        [[1.0, 0.0],  [2.0, 0.0],  [1.5, 0.85]],
        [[0.5, 0.85], [1.0, 0.0],  [1.5, 0.85]],
        [[1.5, 0.85], [2.0, 0.0],  [2.5, 0.85]],
        [[0.5, 0.85], [1.5, 0.85], [1.0, 1.7]],
        [[1.5, 0.85], [2.5, 0.85], [2.0, 1.7]],
        [[1.0, 1.7],  [1.5, 0.85], [2.0, 1.7]],
        [[2.0, 1.7],  [2.5, 0.85], [3.0, 1.7]],
    ];

    // closure that returns the number of nodes of a cell `t`
    let get_nnode = |_t| Ok(3);

    // closure that returns the coordinates of point `m` of cell `t`
    let get_x = |t: usize, m: usize| Ok(&TRIS[t][m][..]);

    // allocate grid search tool
    let ndim = 2;
    let grid = GridSearchCell::new(ndim, TRIS.len(), get_nnode, get_x, None, None)?;

    // print information
    println!("{}", grid);

    // closure that tells whether the point is in the cell or not
    let is_in_cell = |t: usize, x: &[f64]| Ok(is_point_inside_triangle(&TRIS[t][0], &TRIS[t][1], &TRIS[t][2], x));

    // find triangle given coords
    let x = &[1.0, 0.5];
    let id = grid.find_cell(x, is_in_cell)?;
    println!("\nwith x = {:?}", x);
    println!("found triangle with id = {:?} | {:?}", id, brute_force_search(&TRIS, x));
    assert_eq!(id, Some(2));

    // find with another triangle
    let x = &[2.9, 1.6];
    let id = grid.find_cell(x, is_in_cell)?;
    println!("\nwith x = {:?}", x);
    println!("found triangle with id = {:?} | {:?}", id, brute_force_search(&TRIS, x));
    assert_eq!(id, Some(7));

    // maybe find with another triangle
    let x = &[3.0, 1.0];
    let id = grid.find_cell(x, is_in_cell)?;
    println!("\nwith x = {:?}", x);
    println!("found triangle with id = {:?} | {:?}", id, brute_force_search(&TRIS, x));
    assert_eq!(id, None);

    // draw triangles and grid
    let mut plot = Plot::new();
    grid.draw(&mut plot, false)?;
    draw_triangles(&mut plot, &TRIS);
    plot.set_equal_axes(true)
        .set_figure_size_points(600.0, 600.0)
        .grid_and_labels("x", "y")
        .save("/tmp/gemlab/example_grid_search_triangles.svg")?;
    Ok(())
}

fn draw_triangles(plot: &mut Plot, triangles: &[[[f64; 2]; 3]]) {
    let mut ids = Text::new();
    ids.set_color("#cd0000")
        .set_bbox(true)
        .set_bbox_style("circle,pad=0.3")
        .set_bbox_edgecolor("None")
        .set_bbox_facecolor("white")
        .set_bbox_alpha(0.8);
    let mut canvas = Canvas::new();
    let (mut xc, mut yc) = (0.0, 0.0);
    canvas.set_face_color("#fefddc").set_edge_color("#fcb827");
    for t in 0..triangles.len() {
        canvas.polycurve_begin();
        for m in 0..3 {
            let code = if m == 0 {
                xc = triangles[t][m][0];
                yc = triangles[t][m][1];
                PolyCode::MoveTo
            } else {
                xc += triangles[t][m][0];
                yc += triangles[t][m][1];
                PolyCode::LineTo
            };
            canvas.polycurve_add(&triangles[t][m][0], &triangles[t][m][1], code);
        }
        canvas.polycurve_end(true);
        xc /= 3.0;
        yc /= 3.0;
        ids.draw(xc, yc, format!("{}", t).as_str());
    }
    plot.add(&canvas).add(&ids);
}
