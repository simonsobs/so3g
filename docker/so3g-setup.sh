#!/usr/bin/env bash

np_incl=$(python3 -c 'import numpy as np; print(np.get_include())')

mkdir -p build
cd build
cmake \
    -DCMAKE_VERBOSE_MAKEFILE=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DPython_EXECUTABLE=$(which python3) \
    -DPython_NumPy_INCLUDE_DIRS="${np_incl}" \
    ..
make -j 2
make install
