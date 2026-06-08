#!/bin/bash

set -euo pipefail

ALTERNATIVE="${1:-0}"

# image name and Dockerfile
if [ "${ALTERNATIVE}" = "1" ]; then
    DOCKERFILE="zdocker/Dockerfile.Arch.Mkl.Local"
    NAME="cpmech/gemlab_arch_mkl_local"
else
    DOCKERFILE="zdocker/Dockerfile.Arch"
    NAME="cpmech/gemlab_arch"
fi

# build Docker image
docker build \
    -f "${DOCKERFILE}" \
    -t "${NAME}" \
    .

echo
echo "... SUCCESS: image ${NAME} created ..."
echo
