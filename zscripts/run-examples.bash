#!/bin/bash

for example in examples/*.rs; do
    filename="$(basename "$example")"
    filekey="${filename%%.*}"
    if ! [ "$filekey" = "check_grid_search_performance" ]; then
        cargo run --example $filekey
    fi
done
