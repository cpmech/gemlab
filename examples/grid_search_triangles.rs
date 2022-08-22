use gemlab::geometry::{in_triangle, triangle_coords};
use gemlab::util::GridSearchCell;
use gemlab::StrError;
use plotpy::{Canvas, Plot, PolyCode, Text};
use std::fs::File;
use std::io::{BufRead, BufReader};

fn brute_force_search(triangles: &Vec<Vec<Vec<f64>>>, x: &[f64]) -> Option<usize> {
    let mut zeta = vec![0.0; 3]; // 3 triangle coords
    for i in 0..triangles.len() {
        let t = &triangles[i];
        triangle_coords(&mut zeta, &t[0], &t[1], &t[2], x);
        if in_triangle(&zeta) {
            return Some(i);
        }
    }
    None
}

fn main() -> Result<(), StrError> {
    // select mesh for testing
    let num_triangle = 200;
    let path = format!("data/triangles/example_grid_search_gen_triangles_{}.dat", num_triangle);
    let tris = read_data(path.as_str())?;

    // closure that returns the number of nodes of a cell `t`
    let get_nnode = |_t| Ok(3);

    // closure that returns the coordinates of point `m` of cell `t`
    let get_x = |t: usize, m: usize| Ok(&tris[t][m][..]);

    // allocate grid search tool
    let ndim = 2;
    let grid = GridSearchCell::new(ndim, tris.len(), get_nnode, get_x, None, None)?;

    // print information
    // println!("{}", grid);
    grid.print_stat();

    // closure that tells whether the point is in the cell or not
    let is_in_cell = |t: usize, x: &[f64]| {
        let mut zeta = vec![0.0; 3];
        triangle_coords(&mut zeta, &tris[t][0], &tris[t][1], &tris[t][2], x);
        Ok(in_triangle(&zeta))
    };

    // find triangle given coords
    let x = &[1.0, 0.5];
    let id = grid.find_cell(x, is_in_cell)?;
    println!("\nwith x = {:?}", x);
    println!("found triangle with id = {:?} | {:?}", id, brute_force_search(&tris, x));

    // find with another triangle
    let x = &[2.9, 1.6];
    let id = grid.find_cell(x, is_in_cell)?;
    println!("\nwith x = {:?}", x);
    println!("found triangle with id = {:?} | {:?}", id, brute_force_search(&tris, x));

    // maybe find with another triangle
    let x = &[3.0, 1.0];
    let id = grid.find_cell(x, is_in_cell)?;
    println!("\nwith x = {:?}", x);
    println!("found triangle with id = {:?} | {:?}", id, brute_force_search(&tris, x));

    // draw triangles and grid
    let mut plot = Plot::new();
    grid.draw(&mut plot, false)?;
    draw_triangles(&mut plot, &tris, false);
    let path = format!("/tmp/gemlab/example_grid_search_triangles_{}.svg", tris.len());
    plot.set_equal_axes(true)
        .set_figure_size_points(800.0, 800.0)
        .grid_and_labels("x", "y")
        .save(path.as_str())?;
    Ok(())
}

fn draw_triangles(plot: &mut Plot, triangles: &Vec<Vec<Vec<f64>>>, with_ids: bool) {
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
        if with_ids {
            ids.draw(xc, yc, format!("{}", t).as_str());
        }
    }
    plot.add(&canvas);
    if with_ids {
        plot.add(&ids);
    }
}

// Returns [num_triangle][nnode=3][ndim=2]
fn read_data(path: &str) -> Result<Vec<Vec<Vec<f64>>>, StrError> {
    let mut data = Vec::new();
    let input = File::open(path).map_err(|_| "cannot open file")?;
    let buffered = BufReader::new(input);
    let mut lines_iter = buffered.lines();
    let num_triangle: usize = match lines_iter.next() {
        Some(v) => {
            let line = v.unwrap();
            line.parse().map_err(|_| "cannot parse num_triangle")?
        }
        None => return Err("cannot read num_triangle"),
    };
    for _ in 0..num_triangle {
        let mut vertices = Vec::new();
        match lines_iter.next() {
            Some(v) => {
                let line = v.unwrap();
                let mut array = line.split_whitespace();
                for _ in 0..3 {
                    let mut x = vec![0.0; 2];
                    for i in 0..2 {
                        x[i] = match array.next() {
                            Some(w) => w.parse().map_err(|_| "cannot parse vertex coordinate")?,
                            None => return Err("cannot read vertex coordinate"),
                        };
                    }
                    vertices.push(x);
                }
            }
            None => return Err("cannot read num_triangle"),
        }
        data.push(vertices);
    }
    Ok(data)
}
