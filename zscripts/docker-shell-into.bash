#!/bin/bash

set -euo pipefail

NAME="cpmech/gemlab_arch"

docker run --rm -it "${NAME}:latest" /bin/bash
