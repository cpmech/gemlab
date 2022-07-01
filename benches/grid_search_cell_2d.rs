use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, SamplingMode, Throughput};
use gemlab::geometry::is_point_inside_triangle;
use gemlab::util::GridSearchCell;
use russell_lab::{generate2d, StrError};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Duration;

fn brute_force_search(triangles: &Vec<Vec<Vec<f64>>>, x: &[f64]) -> Option<usize> {
    for i in 0..triangles.len() {
        let t = &triangles[i];
        if is_point_inside_triangle(&t[0], &t[1], &t[2], x) {
            return Some(i);
        }
    }
    None
}

const NP: usize = 100;

fn my_benchmark(crit: &mut Criterion) {
    // let sizes = &[200];
    let sizes = &[200, 800, 3200, 7200, 20000];
    let mut meshes = HashMap::new();
    let mut grids = HashMap::new();
    let mut xxs = HashMap::new();
    let mut yys = HashMap::new();
    for size in sizes {
        let path = format!("data/triangles/example_grid_search_gen_triangles_{}.dat", size);
        let (tris, xmin, xmax) = read_data(path.as_str()).unwrap();
        let get_x = |t: usize, m: usize| Ok(&tris[t][m][..]);
        let grid = GridSearchCell::new(2, tris.len(), |_| Ok(3), get_x, None, None).unwrap();
        let (xx, yy) = generate2d(xmin[0], xmax[0], xmin[1], xmax[1], NP, NP);
        meshes.insert(size, tris);
        grids.insert(size, grid);
        xxs.insert(size, xx);
        yys.insert(size, yy);
    }
    let mut group = crit.benchmark_group("my_benchmark");
    for size in sizes {
        group.sampling_mode(SamplingMode::Flat);
        group.throughput(Throughput::Elements(*size as u64));
        group.bench_with_input(BenchmarkId::new("Grid", size), size, |b, &size| {
            let tris = meshes.get(&size).unwrap();
            let grid = grids.get(&size).unwrap();
            // println!("{}", grids.get(&size).unwrap());
            let in_cell = |t: usize, x: &[f64]| Ok(is_point_inside_triangle(&tris[t][0], &tris[t][1], &tris[t][2], x));
            let xx = xxs.get(&size).unwrap();
            let yy = yys.get(&size).unwrap();
            b.iter(|| {
                for i in 0..NP {
                    for j in 0..NP {
                        let x = &[xx[i][j], yy[i][j]];
                        grid.find_cell(x, in_cell).unwrap();
                    }
                }
            });
        });
        group.bench_with_input(BenchmarkId::new("Brute", size), size, |b, &size| {
            let tris = meshes.get(&size).unwrap();
            let xx = xxs.get(&size).unwrap();
            let yy = yys.get(&size).unwrap();
            b.iter(|| {
                for i in 0..NP {
                    for j in 0..NP {
                        let x = &[xx[i][j], yy[i][j]];
                        brute_force_search(&tris, x);
                    }
                }
            });
        });
    }
    group.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(100));
    targets= my_benchmark);
criterion_main!(benches);

// Returns [num_triangle][nnode=3][ndim=2]
fn read_data(path: &str) -> Result<(Vec<Vec<Vec<f64>>>, Vec<f64>, Vec<f64>), StrError> {
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
    let mut xmin = vec![f64::MAX; 2];
    let mut xmax = vec![f64::MIN; 2];
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
                        xmin[i] = f64::min(xmin[i], x[i]);
                        xmax[i] = f64::max(xmax[i], x[i]);
                    }
                    vertices.push(x);
                }
            }
            None => return Err("cannot read num_triangle"),
        }
        data.push(vertices);
    }
    Ok((data, xmin, xmax))
}
