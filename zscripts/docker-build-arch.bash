#!/bin/bash

set -euo pipefail

DOCKERFILE="zdocker/Dockerfile.Arch"
NAME="cpmech/gemlab_arch"

# build Docker image
docker build \
    -f "${DOCKERFILE}" \
    -t "${NAME}" \
    .

echo
echo "... SUCCESS: image ${NAME} created ..."
echo
