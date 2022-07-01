use gemlab::util::SQRT_3;
use gemlab::StrError;
use plotpy::{Canvas, Plot, PolyCode};

const APPROXIMATED: bool = true;

// Equilateral triangle
//
//           /
//        2   \        0-----------2
//       / \   \        \         /
//      / ↑ \   l        \       /
//     /  h  \   \        \     /
//    /   ↓   \   \        \   /
//   /         \   /        \ /
//  0-----------1            1
//
//  |--s--|--s--|
//
//  |-----l-----|
//
// area = l * h / 2.0;
fn generate_equilateral_triangle(l: f64, x0: f64, y0: f64, flipped: bool) -> Vec<[f64; 2]> {
    let s = l / 2.0;
    let h = if APPROXIMATED { l * 1.7 / 2.0 } else { l * SQRT_3 / 2.0 };
    if flipped {
        return vec![[x0, y0 + h], [x0 + s, y0], [x0 + l, y0 + h]];
    } else {
        return vec![[x0, y0], [x0 + l, y0], [x0 + s, y0 + h]];
    }
}

fn main() -> Result<(), StrError> {
    // constants
    const L: f64 = 1.0;
    const NX: usize = 2;
    const NY: usize = 2;
    const X0: f64 = 0.0;
    const Y0: f64 = 0.0;
    let h = if APPROXIMATED { L * 1.7 / 2.0 } else { L * SQRT_3 / 2.0 };

    // generate triangles
    let mut y0 = Y0;
    let mut triangles: Vec<Vec<[f64; 2]>> = Vec::new();
    for j in 0..NY {
        let mut x0 = X0 + (j as f64) * L / 2.0;
        for _ in 0..NX {
            let tri = generate_equilateral_triangle(L, x0, y0, false);
            triangles.push(tri);
            x0 += L;
        }
        x0 = X0 + (j as f64) * L / 2.0 + L / 2.0;
        for _ in 0..NX {
            let tri = generate_equilateral_triangle(L, x0, y0, true);
            triangles.push(tri);
            x0 += L;
        }
        y0 += h;
    }

    // draw
    let mut plot = Plot::new();
    draw_triangles(&mut plot, &triangles);
    plot.set_equal_axes(true)
        .set_figure_size_points(600.0, 600.0)
        .grid_and_labels("x", "y")
        .save("/tmp/gemlab/example_grid_search_gen_triangles.svg")?;

    // print code to be used in doc example
    println!("const TRIS: [[[f64; 2]; 3]; {}] = [", triangles.len());
    for t in &triangles {
        println!(
            "[[{:?},{:?}], [{:?},{:?}], [{:?},{:?}]],",
            t[0][0], t[0][1], t[1][0], t[1][1], t[2][0], t[2][1]
        );
    }
    println!("];");
    Ok(())
}

fn draw_triangles(plot: &mut Plot, triangles: &Vec<Vec<[f64; 2]>>) {
    let mut xmin = vec![f64::MAX; 2];
    let mut xmax = vec![f64::MIN; 2];
    let mut canvas = Canvas::new();
    canvas.set_face_color("#fefddc").set_edge_color("#fcb827");
    for t in 0..triangles.len() {
        canvas.polycurve_begin();
        for m in 0..3 {
            for i in 0..2 {
                xmin[i] = f64::min(xmin[i], triangles[t][m][i]);
                xmax[i] = f64::max(xmax[i], triangles[t][m][i]);
            }
            let code = if m == 0 { PolyCode::MoveTo } else { PolyCode::LineTo };
            canvas.polycurve_add(&triangles[t][m][0], &triangles[t][m][1], code);
        }
        canvas.polycurve_end(true);
    }
    plot.add(&canvas);
    plot.set_range(xmin[0], xmax[0], xmin[1], xmax[1]);
}
