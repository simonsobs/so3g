#!/usr/bin/env bash

mkdir -p build
cd build
export PKG_CONFIG_PATH=$(numpy-config --pkgconfigdir)
cmake \
    -DCMAKE_VERBOSE_MAKEFILE=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DPython_EXECUTABLE=$(which python3) \
    ..
make -j 2
make install
