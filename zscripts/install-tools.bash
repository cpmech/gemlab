#!/bin/bash

cargo clean
cargo build --bin hex2msh --release
cargo build --bin qua2msh --release

sudo cp ~/rust_modules/release/hex2msh /usr/local/bin/
sudo cp ~/rust_modules/release/qua2msh /usr/local/bin/

echo "Installed hex2msh and qua2msh to /usr/local/bin/"
