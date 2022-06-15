# Gemlab examples

* [block_2d_qua8](#block_2d_qua8) -- Shows how to subdivide a rectangle (Block) into Qua8 cells
* [find_in_grid_2d](#find_in_grid_2d) -- Uses GridSearch to assist in finding points in the x-y space
* [integ_mom_inertia_disk](#integ_mom_inertia_disk) -- Approximates the second moment of inertia of a disk (circle)
* [integ_mom_inertia_ring](#integ_mom_inertia_ring) -- Approximates the second moment of inertia of a ring
* [integration_tet4](#integration_tet4) -- Computes the stiffness matrix of a Tet4 element
* [integration_tri3](#integration_tri3) -- Performs numerical integrations with a Tri3 element
* [isoparametric_qua4](#isoparametric_qua4) -- Illustrates the isoparametric formula with a Qua4
* [mesh_2d_tri3_qua4](#mesh_2d_tri3_qua4) -- Converts a string into a Mesh of Tri3 and Qua4

## block_2d_qua8

Shows how to subdivide a rectangle (Block) into Qua8 cells.

[Source code](https://github.com/cpmech/gemlab/blob/main/examples/integ_mom_inertia_disk.rs)

The output is:

```text
# header
# ndim npoint ncell
2 21 4

# points
# id x y {z}
0 0 0
1 1 0
2 1 1
3 0 1
4 0.5 0
5 1 0.5
6 0.5 1
7 0 0.5
8 2 0
9 2 1
10 1.5 0
11 2 0.5
12 1.5 1
13 1 2
14 0 2
15 1 1.5
16 0.5 2
17 0 1.5
18 2 2
19 2 1.5
20 1.5 2

# cells
# id att kind  points
0 1 qua8  0 1 2 3 4 5 6 7
1 1 qua8  1 8 9 2 10 11 12 5
2 1 qua8  3 2 13 14 6 15 16 17
3 1 qua8  2 9 18 13 12 19 20 15
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
