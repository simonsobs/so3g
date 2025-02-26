#!/usr/bin/env bash

mkdir -p build
cd build
cmake \
    -DCMAKE_VERBOSE_MAKEFILE=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DPython_EXECUTABLE=$(which python3) \
    -DCMAKE_MODULE_PATH=../cmake
    ..
make -j 2
make install
