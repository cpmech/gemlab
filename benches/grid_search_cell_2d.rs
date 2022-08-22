use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, SamplingMode, Throughput};
use gemlab::geometry::{in_triangle, triangle_coords};
use gemlab::util::GridSearchCell;
use russell_lab::{generate2d, Matrix, StrError};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Duration;

const NP: usize = 100; // search 100*100 points

struct Analysis {
    tris: Vec<Vec<Vec<f64>>>,
    grid: GridSearchCell,
    xx: Matrix,
    yy: Matrix,
}

fn load_analyses(sizes: &[usize]) -> HashMap<usize, Analysis> {
    let mut analyses = HashMap::new();
    for size in sizes {
        let path = format!("data/triangles/example_grid_search_gen_triangles_{}.dat", size);
        let (tris, xmin, xmax) = read_data(path.as_str()).unwrap();
        let get_x = |t: usize, m: usize| Ok(&tris[t][m][..]);
        let grid = GridSearchCell::new(2, tris.len(), |_| Ok(3), get_x, None, None).unwrap();
        let (xx, yy) = generate2d(xmin[0], xmax[0], xmin[1], xmax[1], NP, NP);
        analyses.insert(*size, Analysis { tris, grid, xx, yy });
    }
    analyses
}

fn bench_new_grid(crit: &mut Criterion) {
    // let sizes = &[200, 800];
    let sizes = &[200, 800, 3200, 7200, 20000];
    let mut all_tris = HashMap::new();
    for size in sizes {
        let path = format!("data/triangles/example_grid_search_gen_triangles_{}.dat", size);
        let (tris, _, _) = read_data(path.as_str()).unwrap();
        all_tris.insert(size, tris);
    }
    let mut group = crit.benchmark_group("bench_new_grid");
    for size in sizes {
        group.throughput(Throughput::Elements(*size as u64));
        group.bench_with_input(BenchmarkId::new("NewGrid", size), size, |b, &size| {
            let tris = all_tris.get(&size).unwrap();
            let get_x = |t: usize, m: usize| Ok(&tris[t][m][..]);
            b.iter(|| {
                GridSearchCell::new(2, tris.len(), |_| Ok(3), get_x, None, None).unwrap();
            });
        });
    }
    group.finish();
}

fn bench_grid_vs_brute(crit: &mut Criterion) {
    // let sizes = &[200, 800];
    let sizes = &[200, 800, 3200, 7200, 20000];
    let analyses = load_analyses(sizes);
    let mut group = crit.benchmark_group("bench_grid_vs_brute");
    for size in sizes {
        group.sampling_mode(SamplingMode::Flat);
        group.throughput(Throughput::Elements(*size as u64));
        group.bench_with_input(BenchmarkId::new("Grid", size), size, |b, &size| {
            let ana = analyses.get(&size).unwrap();
            let is_in_cell = |t: usize, x: &[f64]| {
                let mut zeta = vec![0.0; 3];
                triangle_coords(&mut zeta, &ana.tris[t][0], &ana.tris[t][1], &ana.tris[t][2], x);
                Ok(in_triangle(&zeta))
            };
            b.iter(|| {
                for i in 0..NP {
                    for j in 0..NP {
                        let x = &[ana.xx[i][j], ana.yy[i][j]];
                        ana.grid.find_cell(x, is_in_cell).unwrap();
                    }
                }
            });
        });
        group.bench_with_input(BenchmarkId::new("Brute", size), size, |b, &size| {
            let ana = analyses.get(&size).unwrap();
            b.iter(|| {
                for i in 0..NP {
                    for j in 0..NP {
                        let x = &[ana.xx[i][j], ana.yy[i][j]];
                        brute_force_search(&ana.tris, x);
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
    targets = bench_new_grid, bench_grid_vs_brute);
criterion_main!(benches);

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
