# Benchmarks

We run the benchmarks using the following systems:

**System-1**

`cat /proc/cpuinfo`

```text
processor	: 11
vendor_id	: GenuineIntel
cpu family	: 6
model		: 158
model name	: Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz
```

`lsb_release -a`

```text
Distributor ID:	Ubuntu
Description:	Ubuntu 20.04.4 LTS
Release:	20.04
Codename:	focal
```

`grep MemTotal /proc/meminfo`

```
MemTotal:       49048068 kB
```

## GridSearchCell

### Grid allocation times

First, we benchmark the grid initialization which loops over all triangles, finds the bounding boxes, and finds the largest triangle.

```text
bench_new_grid/NewGrid/200                                                                             
                        time:   [93.762 us 93.817 us 93.881 us]
                        thrpt:  [2.1304 Melem/s 2.1318 Melem/s 2.1331 Melem/s]
bench_new_grid/NewGrid/800                                                                             
                        time:   [370.13 us 370.33 us 370.56 us]
                        thrpt:  [2.1589 Melem/s 2.1603 Melem/s 2.1614 Melem/s]
bench_new_grid/NewGrid/3200                                                                            
                        time:   [1.4503 ms 1.4509 ms 1.4517 ms]
                        thrpt:  [2.2044 Melem/s 2.2055 Melem/s 2.2064 Melem/s]
bench_new_grid/NewGrid/7200                                                                            
                        time:   [3.2617 ms 3.2630 ms 3.2645 ms]
                        thrpt:  [2.2056 Melem/s 2.2066 Melem/s 2.2075 Melem/s]
bench_new_grid/NewGrid/20000                                                                            
                        time:   [9.2292 ms 9.2349 ms 9.2412 ms]
                        thrpt:  [2.1642 Melem/s 2.1657 Melem/s 2.1670 Melem/s]
```

This chart below shows the mean measured time as the number of triangles increases.

![bench_new_grid_1](https://github.com/cpmech/gemlab/raw/main/benches/figures/bench_new_grid_1.svg)

### Grid search versus brute-force search times

Next, we perform the search of 100 x 100 points over the 2D domain of a mesh of triangles. We use the GridSearchCell tool and a naive (brute-force) implementation which simply loops over all triangles for each point.

```text
bench_grid_vs_brute/Grid/200                                                                             
                        time:   [613.52 us 613.81 us 614.12 us]
                        thrpt:  [325.67 Kelem/s 325.84 Kelem/s 325.99 Kelem/s]
bench_grid_vs_brute/Brute/200                                                                             
                        time:   [10.973 ms 10.978 ms 10.983 ms]
                        thrpt:  [18.210 Kelem/s 18.219 Kelem/s 18.226 Kelem/s]
bench_grid_vs_brute/Grid/800                                                                             
                        time:   [722.73 us 723.04 us 723.39 us]
                        thrpt:  [1.1059 Melem/s 1.1064 Melem/s 1.1069 Melem/s]
bench_grid_vs_brute/Brute/800                                                                             
                        time:   [39.245 ms 39.262 ms 39.281 ms]
                        thrpt:  [20.366 Kelem/s 20.376 Kelem/s 20.385 Kelem/s]
bench_grid_vs_brute/Grid/3200                                                                             
                        time:   [831.14 us 831.46 us 831.83 us]
                        thrpt:  [3.8469 Melem/s 3.8486 Melem/s 3.8501 Melem/s]
bench_grid_vs_brute/Brute/3200                                                                            
                        time:   [150.57 ms 150.64 ms 150.72 ms]
                        thrpt:  [21.232 Kelem/s 21.242 Kelem/s 21.252 Kelem/s]
bench_grid_vs_brute/Grid/7200                                                                             
                        time:   [917.37 us 917.79 us 918.28 us]
                        thrpt:  [7.8408 Melem/s 7.8449 Melem/s 7.8485 Melem/s]
bench_grid_vs_brute/Brute/7200                                                                            
                        time:   [331.90 ms 332.05 ms 332.21 ms]
                        thrpt:  [21.673 Kelem/s 21.683 Kelem/s 21.693 Kelem/s]
bench_grid_vs_brute/Grid/20000                                                                            
                        time:   [1.0091 ms 1.0096 ms 1.0102 ms]
                        thrpt:  [19.798 Melem/s 19.810 Melem/s 19.820 Melem/s]
bench_grid_vs_brute/Brute/20000                                                                            
                        time:   [893.07 ms 893.54 ms 894.05 ms]
                        thrpt:  [22.370 Kelem/s 22.383 Kelem/s 22.395 Kelem/s]
```

This chart below shows the mean measured time as the number of triangles increases.

![bench_grid_vs_brute_1](https://github.com/cpmech/gemlab/raw/main/benches/figures/bench_grid_vs_brute_1.svg)

### Discussion

For 20000 triangles, the grid search average time is 9.2349 ms (initialization) plus 1.0096 ms (search) and the brute force time is 893.54 ms. Therefore, including the initialization time, `GridSearchCell` is roughly 87 times faster than the brute force algorithm.
