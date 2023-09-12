# Gemlab examples

* [block_2d_qua8](#block_2d_qua8) -- Shows how to subdivide a rectangle (Block) into Qua8 cells
* [find_in_grid_2d](#find_in_grid_2d) -- Uses GridSearch to assist in finding points in the x-y space
* [integ_mom_inertia_disk](#integ_mom_inertia_disk) -- Approximates the second moment of inertia of a disk (circle)
* [integ_mom_inertia_ring](#integ_mom_inertia_ring) -- Approximates the second moment of inertia of a ring
* [integration_tet4](#integration_tet4) -- Computes the stiffness matrix of a Tet4 element
* [integration_tri3](#integration_tri3) -- Performs numerical integrations with a Tri3 element
* [isoparametric_qua4](#isoparametric_qua4) -- Illustrates the isoparametric formula with a Qua4
* [mesh_2d_tri3_qua4](#mesh_2d_tri3_qua4) -- Converts a string into a Mesh of Tri3 and Qua4
* [grid_search_triangles](#grid_search_triangles) -- Fast-search triangles in a mesh using the GridSearchCell tool

## block_2d_qua8

Shows how to subdivide a rectangle (Block) into Qua8 cells.

[Source code](https://github.com/cpmech/gemlab/blob/main/examples/block_2d_qua8.rs)

The output is:

```text
# header
# ndim npoint ncell
2 21 4

# points
# id marker x y {z}
0 0 0.0 0.0
1 0 1.0 0.0
2 0 1.0 1.0
3 0 0.0 1.0
4 0 0.5 0.0
5 0 1.0 0.5
6 0 0.5 1.0
7 0 0.0 0.5
8 0 2.0 0.0
9 0 2.0 1.0
10 0 1.5 0.0
11 0 2.0 0.5
12 0 1.5 1.0
13 0 1.0 2.0
14 0 0.0 2.0
15 0 1.0 1.5
16 0 0.5 2.0
17 0 0.0 1.5
18 0 2.0 2.0
19 0 2.0 1.5
20 0 1.5 2.0

# cells
# id attribute kind points
0 1 qua8 0 1 2 3 4 5 6 7
1 1 qua8 1 8 9 2 10 11 12 5
2 1 qua8 3 2 13 14 6 15 16 17
3 1 qua8 2 9 18 13 12 19 20 15
```

![mesh](https://github.com/cpmech/gemlab/raw/main/data/figures/example_block_2d_qua8.svg)

## find_in_grid_2d

[Source code](https://github.com/cpmech/gemlab/blob/main/examples/find_in_grid_2d.rs)

Uses GridSearch to assist in finding points in the x-y space.

![find_in_grid_2d](https://github.com/cpmech/gemlab/raw/main/data/figures/example_find_in_grid_2d.svg)

## integ_mom_inertia_disk

Approximates the second moment of inertia of a disk (circle).

[Source code](https://github.com/cpmech/gemlab/blob/main/examples/integ_mom_inertia_disk.rs)

First mesh:

![mesh 1](https://github.com/cpmech/gemlab/raw/main/data/figures/example_mom_inertia_disk_1.svg)

Second mesh:

![mesh 2](https://github.com/cpmech/gemlab/raw/main/data/figures/example_mom_inertia_disk_2.svg)

The output is:

```text
mesh 1: second_mom_inertia = 63.61725744867274 (63.61725123519331): err = 6.21e-6
mesh 2: second_mom_inertia = 63.61725006590863 (63.61725123519331): err = 1.17e-6
```

## integ_mom_inertia_ring

[Source code](https://github.com/cpmech/gemlab/blob/main/examples/integ_mom_inertia_ring.rs)

Approximates the second moment of inertia of a ring.

![mesh](https://github.com/cpmech/gemlab/raw/main/data/figures/example_mom_inertia_ring.svg)

The output is:

```text
second_mom_inertia = 62.83185310895357 (62.83185307179586): err = 3.72e-8
```

## integration_tet4

Computes the stiffness matrix of a Tet4 element.

[Source code](https://github.com/cpmech/gemlab/blob/main/examples/integration_tet4.rs)

The output is:

```text
┌                                                             ┐
│  149  108   24   -1    6   12  -54  -48    0  -94  -66  -36 │
│  108  344   54  -24  104   42  -24 -216  -12  -60 -232  -84 │
│   24   54  113    0   30   35    0  -24  -54  -24  -60  -94 │
│   -1  -24    0   29  -18  -12  -18   24    0  -10   18   12 │
│    6  104   30  -18   44   18   12  -72  -12    0  -76  -36 │
│   12   42   35  -12   18   29    0  -24  -18    0  -36  -46 │
│  -54  -24    0  -18   12    0   36    0    0   36   12    0 │
│  -48 -216  -24   24  -72  -24    0  144    0   24  144   48 │
│    0  -12  -54    0  -12  -18    0    0   36    0   24   36 │
│  -94  -60  -24  -10    0    0   36   24    0   68   36   24 │
│  -66 -232  -60   18  -76  -36   12  144   24   36  164   72 │
│  -36  -84  -94   12  -36  -46    0   48   36   24   72  104 │
└                                                             ┘
```

## integration_tri3

Performs numerical integrations with a Tri3 element.

[Source code](https://github.com/cpmech/gemlab/blob/main/examples/integration_tri3.rs)

```text
a =
┌                   ┐
│ 36.00000000000001 │
│                36 │
│                36 │
└                   ┘
b =
┌                    ┐
│ 24.000000000000007 │
│ 24.000000000000007 │
│                 24 │
│                 24 │
│                 24 │
│                 24 │
└                    ┘
c =
┌    ┐
│ -2 │
│ -4 │
│  6 │
└    ┘
d =
┌     ┐
│ -15 │
│ -10 │
│  12 │
│   4 │
│   3 │
│   6 │
└     ┘
```

## isoparametric_qua4

Illustrates the isoparametric formula with a Qua4

[Source code](https://github.com/cpmech/gemlab/blob/main/examples/isoparametric_qua4.rs)

```text
xm = 4, ym = 4.5
x_interpolated =
┌     ┐
│   4 │
│ 4.5 │
└     ┘
```

## mesh_2d_tri3_qua4

Converts a string into a Mesh of Tri3 and Qua4

[Source code](https://github.com/cpmech/gemlab/blob/main/examples/mesh_2d_tri3_qua4.rs)

Draws the following mesh:

![mesh](https://github.com/cpmech/gemlab/raw/main/data/figures/example_mesh_2d_tri3_qua4.svg)

## grid_search_triangles

Fast-search triangles in a mesh using the GridSearchCell tool.

[Source code](https://github.com/cpmech/gemlab/blob/main/examples/grid_search_triangles.rs)

We run the code with 200, 800, 3200, 7200, and 20000 triangles. The next sections show the results.
The majority of the containers are occupied by 8 triangles.

### 200 triangles

```text
Summary
=======
ncell = 200
xmin = [-0.01, -0.01]
xmax = [15.9932, 8.9918]
side_length = 1.0002
num of non-empty containers = 109
max container num items = 12

Histogram of container num items
================================
[ 0, 1) |  0 
[ 1, 2) | 11 ***********
[ 2, 3) |  0 
[ 3, 4) |  9 *********
[ 4, 5) |  8 ********
[ 5, 6) |  7 *******
[ 6, 7) |  2 **
[ 7, 8) |  7 *******
[ 8, 9) | 56 ************************************************************
[ 9,10) |  0 
[10,11) |  0 
[11,12) |  2 **
[12,13) |  7 *******
   sum = 109
```

![example_grid_search_triangles_200](https://github.com/cpmech/gemlab/raw/main/data/triangles/example_grid_search_triangles_200.png)

### 800 triangles

```text
Summary
=======
ncell = 800
xmin = [-0.01, -0.01]
xmax = [30.996199999999998, 17.993599999999997]
side_length = 1.0002
num of non-empty containers = 397
max container num items = 12

Histogram of container num items
================================
[ 0, 1) |   0 
[ 1, 2) |  19 ****
[ 2, 3) |   0 
[ 3, 4) |  19 ****
[ 4, 5) |  18 ****
[ 5, 6) |  15 ***
[ 6, 7) |   2 
[ 7, 8) |  15 ***
[ 8, 9) | 270 ************************************************************
[ 9,10) |   2 
[10,11) |   0 
[11,12) |   2 
[12,13) |  35 *******
    sum = 397
```

![example_grid_search_triangles_800](https://github.com/cpmech/gemlab/raw/main/data/triangles/example_grid_search_triangles_800.png)

### 3200 triangles

```text
Summary
=======
ncell = 3200
xmin = [-0.01, -0.01]
xmax = [61.0022, 34.997]
side_length = 1.0002
num of non-empty containers = 1463
max container num items = 15

Histogram of container num items
================================
[ 0, 1) |    0 
[ 1, 2) |   35 *
[ 2, 3) |    1 
[ 3, 4) |   34 *
[ 4, 5) |   37 **
[ 5, 6) |   29 *
[ 6, 7) |    6 
[ 7, 8) |   29 *
[ 8, 9) | 1083 ************************************************************
[ 9,10) |    5 
[10,11) |   12 
[11,12) |    6 
[12,13) |  184 **********
[13,14) |    0 
[14,15) |    0 
[15,16) |    2 
    sum = 1463
```

![example_grid_search_triangles_3200](https://github.com/cpmech/gemlab/raw/main/data/triangles/example_grid_search_triangles_3200.png)

### 7200 triangles

```text
Summary
=======
ncell = 7200
xmin = [-0.01, -0.01]
xmax = [91.00819999999999, 52.0004]
side_length = 1.0002
num of non-empty containers = 3140
max container num items = 15

Histogram of container num items
================================
[ 0, 1) |    0 
[ 1, 2) |   50 *
[ 2, 3) |    0 
[ 3, 4) |   52 *
[ 4, 5) |    0 
[ 5, 6) |   42 *
[ 6, 7) |    8 
[ 7, 8) |   44 *
[ 8, 9) | 2429 ************************************************************
[ 9,10) |    8 
[10,11) |   43 *
[11,12) |    8 
[12,13) |  448 ***********
[13,14) |    0 
[14,15) |    0 
[15,16) |    8 
    sum = 3140
```

![example_grid_search_triangles_7200](https://github.com/cpmech/gemlab/raw/main/data/triangles/example_grid_search_triangles_7200.png)

### 20000 triangles

```text
Summary
=======
ncell = 20000
xmin = [-0.01, -0.01]
xmax = [150.02, 86.0072]
side_length = 1.0002
num of non-empty containers = 8634
max container num items = 15

Histogram of container num items
================================
[ 0, 1) |    0 
[ 1, 2) |   84 
[ 2, 3) |    0 
[ 3, 4) |   86 
[ 4, 5) |    0 
[ 5, 6) |   70 
[ 6, 7) |   14 
[ 7, 8) |   71 
[ 8, 9) | 6852 ************************************************************
[ 9,10) |   15 
[10,11) |   70 
[11,12) |   14 
[12,13) | 1344 ***********
[13,14) |    0 
[14,15) |    0 
[15,16) |   14 
    sum = 8634
```

![example_grid_search_triangles_20000](https://github.com/cpmech/gemlab/raw/main/data/triangles/example_grid_search_triangles_20000.png)
