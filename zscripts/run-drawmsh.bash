#!/bin/bash

cargo run --bin drawmsh -- data/meshes/rectangle_tris_quads.msh \
    -p -v -c -a -n -e -w 800 -h 800 \
    --m-normal-vector 0.07 \
    --m-normal-vector-marker 0.05
