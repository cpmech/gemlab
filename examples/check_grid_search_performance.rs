use gemlab::prelude::*;
use gemlab::util::GridSearch;
use russell_lab::*;

const THICKNESS: f64 = 1.0; // out-of-plane thickness
const R1: f64 = 3.0; // inner radius
const R2: f64 = 6.0; // outer radius

fn run(grid_search_ndiv: usize, nr: usize, na: usize) -> Result<(), StrError> {
    println!("\ngrid.ndiv       = {}", grid_search_ndiv);

    let mut sw = Stopwatch::new("");

    let den = 6.0;

    sw.reset();
    let delta_x = (R2 - R1) / (nr as f64);
    let global_max_volume = Some(delta_x * delta_x * delta_x / den);
    println!("max vol         = {:?}", global_max_volume);
    let mesh = Unstructured::quarter_ring_3d(R1, R2, THICKNESS, nr, na, GeoKind::Tet10, global_max_volume)?;
    println!("mesh generation : {}", format_nanoseconds(sw.stop()));

    sw.reset();
    mesh.check_all()?;
    println!("check_all       : {}", format_nanoseconds(sw.stop()));

    // // too slow
    // sw.reset();
    // mesh.check_overlapping_points(0.001)?;
    // println!("check_overlapping_points: {}", format_nanoseconds(sw.stop()));

    let tol = 0.001;

    let (xmin, xmax) = mesh.get_limits();
    sw.reset();
    let mut grid = GridSearch::new(&xmin, &xmax, Some(grid_search_ndiv), Some(tol), None)?;
    println!("grid.new        : {}", format_nanoseconds(sw.stop()));

    sw.reset();
    for point in &mesh.points {
        grid.insert(point.id, &point.coords)?;
    }
    println!("grid.insert     : {}", format_nanoseconds(sw.stop()));

    sw.reset();
    for point in &mesh.points {
        match grid.search(&point.coords)? {
            Some(id) => {
                if id != point.id {
                    println!("found overlapping points: {} => {}", id, point.id);
                    return Err("found overlapping points");
                }
            }
            None => (),
        }
    }
    println!("grid.search     : {}", format_nanoseconds(sw.stop()));
    Ok(())
}

fn main() -> Result<(), StrError> {
    for (nr, na) in [(1, 2), (4, 8), (50, 100)] {
        println!(
            "\n......................... nr = {} na = {} .........................",
            nr, na
        );
        run(5, nr, na)?;
        run(20, nr, na)?;
        run(40, nr, na)?;
        run(100, nr, na)?;
    }
    Ok(())
}
