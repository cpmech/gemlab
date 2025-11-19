#!/bin/bash

cargo clean
cargo build --bin drawmsh --release
cargo build --bin hex2msh --release
cargo build --bin msh2tri --release
cargo build --bin qua2msh --release

sudo cp ~/rust_modules/release/drawmsh /usr/local/bin/
sudo cp ~/rust_modules/release/hex2msh /usr/local/bin/
sudo cp ~/rust_modules/release/msh2tri /usr/local/bin/
sudo cp ~/rust_modules/release/qua2msh /usr/local/bin/

echo "Installed tools to /usr/local/bin/"
