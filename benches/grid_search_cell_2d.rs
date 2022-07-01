use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use gemlab::geometry::is_point_inside_triangle;
use gemlab::util::GridSearchCell;
use russell_lab::generate2d;

fn brute_force_search(triangles: &[[[f64; 2]; 3]], x: &[f64]) -> Option<usize> {
    for i in 0..triangles.len() {
        let t = &triangles[i];
        if is_point_inside_triangle(&t[0], &t[1], &t[2], x) {
            return Some(i);
        }
    }
    None
}

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

fn my_benchmark(crit: &mut Criterion) {
    let get_x = |t: usize, m: usize| Ok(&TRIS[t][m][..]);
    let grid = GridSearchCell::new(2, TRIS.len(), |_| Ok(3), get_x, None, None).unwrap();
    let in_cell = |t: usize, x: &[f64]| Ok(is_point_inside_triangle(&TRIS[t][0], &TRIS[t][1], &TRIS[t][2], x));
    let sizes = &[2]; // min must be 2 because of generate2d
    let mut group = crit.benchmark_group("my_benchmark");
    for size in sizes {
        group.throughput(Throughput::Elements(*size as u64));
        group.bench_with_input(BenchmarkId::new("Grid", size), size, |b, &size| {
            let (xmin, xmax) = grid.limits();
            let (xx, yy) = generate2d(xmin[0], xmax[0], xmin[1], xmax[1], size, size);
            b.iter(|| {
                for i in 0..size {
                    for j in 0..size {
                        let x = &[xx[i][j], yy[i][j]];
                        grid.find_cell(x, in_cell).unwrap();
                    }
                }
            });
        });
        group.bench_with_input(BenchmarkId::new("Brute", size), size, |b, &size| {
            let (xmin, xmax) = grid.limits();
            let (xx, yy) = generate2d(xmin[0], xmax[0], xmin[1], xmax[1], size, size);
            b.iter(|| {
                for i in 0..size {
                    for j in 0..size {
                        let x = &[xx[i][j], yy[i][j]];
                        brute_force_search(&TRIS, x);
                    }
                }
            });
        });
    }
    group.finish();
}

criterion_group!(benches, my_benchmark);
criterion_main!(benches);
