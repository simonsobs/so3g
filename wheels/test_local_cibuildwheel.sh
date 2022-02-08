#!/bin/bash

# Before running this from the so3g git checkout directory,
# you should pip install cibuildwheel

export CIBW_BUILD="cp38-manylinux_x86_64"
export CIBW_MANYLINUX_X86_64_IMAGE="manylinux2014"
export CIBW_BUILD_VERBOSITY=3
export CIBW_ENVIRONMENT_LINUX="CC=gcc CXX=g++ CFLAGS='-O3 -fPIC' CXXFLAGS='-O3 -fPIC -std=c++11' BOOST_ROOT=/usr/local"
export CIBW_ENVIRONMENT_MACOS="CC=clang CXX=clang++ CFLAGS='-O3 -fPIC' CXXFLAGS='-O3 -fPIC -std=c++11 -stdlib=libc++' BOOST_ROOT=/usr/local"
export CIBW_BEFORE_BUILD_LINUX="./wheels/install_deps_linux.sh"
export CIBW_REPAIR_WHEEL_COMMAND_LINUX="./wheels/repair_wheel_linux.sh {dest_dir} {wheel}"
export CIBW_REPAIR_WHEEL_COMMAND_MACOS="./wheels/repair_wheel_macos.sh {dest_dir} {wheel} {delocate_archs}"
export CIBW_BEFORE_TEST=
export CIBW_TEST_COMMAND="python -c 'import so3g.smurf.reader; from spt3g import core'"


# Get the current date for logging
now=$(date "+%Y-%m-%d_%H:%M:%S")

# Run it
cibuildwheel --platform linux --archs x86_64 --output-dir wheelhouse . 2>&1 | tee log_${now}
