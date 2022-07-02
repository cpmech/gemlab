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
                        time:   [94.311 us 94.405 us 94.514 us]
                        thrpt:  [2.1161 Melem/s 2.1185 Melem/s 2.1206 Melem/s]
bench_new_grid/NewGrid/800                                                                             
                        time:   [373.66 us 373.86 us 374.08 us]
                        thrpt:  [2.1386 Melem/s 2.1398 Melem/s 2.1410 Melem/s]
bench_new_grid/NewGrid/3200                                                                            
                        time:   [1.4656 ms 1.4663 ms 1.4670 ms]
                        thrpt:  [2.1813 Melem/s 2.1823 Melem/s 2.1834 Melem/s]
bench_new_grid/NewGrid/7200                                                                            
                        time:   [3.3000 ms 3.3018 ms 3.3037 ms]
                        thrpt:  [2.1794 Melem/s 2.1806 Melem/s 2.1818 Melem/s]
bench_new_grid/NewGrid/20000                                                                            
                        time:   [9.3440 ms 9.3493 ms 9.3550 ms]
                        thrpt:  [2.1379 Melem/s 2.1392 Melem/s 2.1404 Melem/s]
```

This chart below shows the mean measured time as the number of triangles increases.

![bench_new_grid_1](https://github.com/cpmech/gemlab/raw/main/benches/figures/bench_new_grid_1.svg)

### Grid search versus brute-force search times

Next, we perform the search of 100 x 100 points over the 2D domain of a mesh of triangles. We use the GridSearchCell tool and a naive (brute-force) implementation which simply loops over all triangles for each point.


```text
bench_grid_vs_brute/Grid/200                                                                             
                        time:   [776.10 us 776.49 us 776.91 us]
                        thrpt:  [257.43 Kelem/s 257.57 Kelem/s 257.70 Kelem/s]
bench_grid_vs_brute/Brute/200                                                                             
                        time:   [10.524 ms 10.529 ms 10.534 ms]
                        thrpt:  [18.987 Kelem/s 18.995 Kelem/s 19.003 Kelem/s]
bench_grid_vs_brute/Grid/800                                                                             
                        time:   [877.44 us 877.84 us 878.27 us]
                        thrpt:  [910.88 Kelem/s 911.33 Kelem/s 911.74 Kelem/s]
bench_grid_vs_brute/Brute/800                                                                             
                        time:   [37.116 ms 37.136 ms 37.159 ms]
                        thrpt:  [21.529 Kelem/s 21.542 Kelem/s 21.554 Kelem/s]
bench_grid_vs_brute/Grid/3200                                                                             
                        time:   [983.19 us 984.27 us 985.75 us]
                        thrpt:  [3.2463 Melem/s 3.2511 Melem/s 3.2547 Melem/s]
bench_grid_vs_brute/Brute/3200                                                                            
                        time:   [144.54 ms 144.65 ms 144.76 ms]
                        thrpt:  [22.105 Kelem/s 22.123 Kelem/s 22.139 Kelem/s]
bench_grid_vs_brute/Grid/7200                                                                            
                        time:   [1.1606 ms 1.1661 ms 1.1715 ms]
                        thrpt:  [6.1457 Melem/s 6.1743 Melem/s 6.2035 Melem/s]
bench_grid_vs_brute/Brute/7200                                                                            
                        time:   [331.91 ms 333.21 ms 334.50 ms]
                        thrpt:  [21.525 Kelem/s 21.608 Kelem/s 21.692 Kelem/s]
bench_grid_vs_brute/Grid/20000                                                                            
                        time:   [1.3010 ms 1.3124 ms 1.3241 ms]
                        thrpt:  [15.105 Melem/s 15.239 Melem/s 15.373 Melem/s]
bench_grid_vs_brute/Brute/20000                                                                            
                        time:   [879.77 ms 886.01 ms 892.72 ms]
                        thrpt:  [22.403 Kelem/s 22.573 Kelem/s 22.733 Kelem/s]
```

This chart below shows the mean measured time as the number of triangles increases.

![bench_grid_vs_brute_1](https://github.com/cpmech/gemlab/raw/main/benches/figures/bench_grid_vs_brute_1.svg)

### Discussion

For 20000 triangles, the grid search average time is 9.3493 ms (initialization) plus 1.3124 ms (search) and the brute force time is 886.01 ms. Therefore, including the initialization time, `GridSearchCell` is roughly 80 times faster than the brute force algorithm.
