#!/bin/bash
#
# This script is designed to run within a container managed by cibuildwheel.
# This will use a recent version of OS X.
#

set -e

# Location of this script
pushd $(dirname $0) >/dev/null 2>&1
scriptdir=$(pwd)
popd >/dev/null 2>&1
echo "Wheel script directory = ${scriptdir}"

# Build options.

# FIXME:  would be nice to switch to clang once spt3g / cereal
# runtime registration works.
use_gcc=yes
# use_gcc=no
gcc_version=14

if [ "x${use_gcc}" = "xyes" ]; then
    CC=gcc-${gcc_version}
    CXX=g++-${gcc_version}
    FC=gfortran-${gcc_version}
    CFLAGS="-O3 -fPIC"
    CXXFLAGS="-O3 -fPIC -std=c++14"
    FCFLAGS="-O3 -fPIC -pthread"
    OMPFLAGS="-fopenmp"
else
    export MACOSX_DEPLOYMENT_TARGET=$(python3 -c "import sysconfig as s; print(s.get_config_vars()['MACOSX_DEPLOYMENT_TARGET'])")
    echo "Using MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET}"
    CC=clang
    CXX=clang++
    #FC=""
    CFLAGS="-O3 -fPIC"
    CXXFLAGS="-O3 -fPIC -std=c++14 -stdlib=libc++"
    #FCFLAGS=""
    #OMPFLAGS=""
fi

MAKEJ=2

PREFIX=/usr/local

# Workaround permissions on macos-14 github runner
# https://github.com/actions/runner-images/issues/9272
sudo chown -R runner:admin /usr/local

# Install library dependencies with homebrew
brew install sqlite3 netcdf

# Optionally install gcc.
if [ "x${use_gcc}" = "xyes" ]; then
    brew install gcc@${gcc_version}
fi

# Update pip
python3 -m pip install --upgrade pip

# Install a couple of base packages that are always required
python3 -m pip install cmake wheel setuptools

pyver=$(python3 --version 2>&1 | awk '{print $2}' | sed -e "s#\(.*\)\.\(.*\)\..*#\1.\2#")

# Install build requirements.
CC="${CC}" CFLAGS="${CFLAGS}" python3 -m pip install -v -r "${scriptdir}/../requirements.txt" --prefer-binary

# Install Openblas

openblas_version=0.3.29
openblas_dir=OpenBLAS-${openblas_version}
openblas_pkg=${openblas_dir}.tar.gz

if [ ! -e ${openblas_pkg} ]; then
    echo "Fetching OpenBLAS..."
    curl -SL https://github.com/xianyi/OpenBLAS/archive/v${openblas_version}.tar.gz -o ${openblas_pkg}
fi

echo "Building OpenBLAS..."

omp="USE_OPENMP=1"
if [ "x${OMPFLAGS}" = "x" ]; then
    omp="USE_OPENMP=0"
fi

rm -rf ${openblas_dir}
tar xzf ${openblas_pkg} \
    && pushd ${openblas_dir} >/dev/null 2>&1 \
    && make ${omp} NO_STATIC=1 \
    MAKE_NB_JOBS=${MAKEJ} \
    CC="${CC}" FC="${FC}" DYNAMIC_ARCH=1 TARGET=GENERIC \
    COMMON_OPT="${CFLAGS}" FCOMMON_OPT="${FCFLAGS}" \
    EXTRALIB="${OMPFLAGS}" libs netlib shared \
    && make NO_STATIC=1 DYNAMIC_ARCH=1 TARGET=GENERIC PREFIX="${PREFIX}" install \
    && popd >/dev/null 2>&1

# Build GSL

gsl_version=2.8
gsl_dir=gsl-${gsl_version}
gsl_pkg=${gsl_dir}.tar.gz

if [ ! -e ${gsl_pkg} ]; then
    echo "Fetching GSL..."
    curl -SL https://ftp.gnu.org/gnu/gsl/gsl-${gsl_version}.tar.gz -o ${gsl_pkg}
fi

echo "Building GSL..."

rm -rf ${gsl_dir}
tar xzf ${gsl_pkg} \
    && pushd ${gsl_dir} >/dev/null 2>&1 \
    && mkdir -p build \
    && pushd build >/dev/null 2>&1 \
    && CC="${CC}" CFLAGS="-O3 -fPIC" ../configure --prefix="${PREFIX}" \
    && make -j ${MAKEJ} \
    && make install \
    && popd >/dev/null 2>&1 \
    && popd >/dev/null 2>&1

# Build Eigen

eigen_version=3.4.0
eigen_dir=eigen-${eigen_version}
eigen_pkg=${eigen_dir}.tar.gz

echo "Fetching Eigen..."

if [ ! -e ${eigen_pkg} ]; then
    curl -SL "https://gitlab.com/libeigen/eigen/-/archive/${eigen_version}/eigen-${eigen_version}.tar.bz2" -o "${eigen_pkg}"
fi

echo "Building Eigen..."

rm -rf ${eigen_dir}
tar xjf ${eigen_pkg} \
    && pushd ${eigen_dir} >/dev/null 2>&1 \
    && mkdir -p build \
    && pushd build >/dev/null 2>&1 \
    && cmake \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    .. \
    && make install \
    && popd >/dev/null 2>&1 \
    && popd >/dev/null 2>&1

# Build GLOG

glog_version=0.7.1
glog_dir=glog-${glog_version}
glog_pkg=${glog_dir}.tar.gz

echo "Fetching GLOG..."

if [ ! -e ${glog_pkg} ]; then
    curl -SL "https://github.com/google/glog/archive/refs/tags/v${glog_version}.tar.gz" -o "${glog_pkg}"
fi

echo "Building GLOG..."

rm -rf ${glog_dir}
tar xzf ${glog_pkg} \
    && pushd ${glog_dir} >/dev/null 2>&1 \
    && mkdir -p build \
    && pushd build >/dev/null 2>&1 \
    && cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="${CC}" \
    -DCMAKE_C_FLAGS="${CFLAGS}" \
    -DCMAKE_CXX_COMPILER="${CXX}" \
    -DCMAKE_CXX_FLAGS="${CXXFLAGS}" \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DWITH_GFLAGS:BOOL=OFF \
    -DWITH_GTEST:BOOL=OFF \
    -DBUILD_SHARED_LIBS:BOOL=OFF \
    .. \
    && make -j ${MAKEJ} install \
    && popd >/dev/null 2>&1 \
    && popd >/dev/null 2>&1

# Build Ceres

ceres_version=2.2.0
ceres_dir=ceres-solver-${ceres_version}
ceres_pkg=${ceres_dir}.tar.gz

echo "Fetching ceres-solver..."

if [ ! -e ${ceres_pkg} ]; then
    curl -SL "http://ceres-solver.org/ceres-solver-${ceres_version}.tar.gz" -o "${ceres_pkg}"
fi

echo "Building ceres-solver..."

rm -rf ${ceres_dir}
tar xzf ${ceres_pkg} \
    && pushd ${ceres_dir} >/dev/null 2>&1 \
    && mkdir -p build_dir \
    && pushd build_dir >/dev/null 2>&1 \
    && cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="${CC}" \
    -DCMAKE_C_FLAGS="${CFLAGS}" \
    -DCMAKE_CXX_COMPILER="${CXX}" \
    -DCMAKE_CXX_FLAGS="${CXXFLAGS}" \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DBUILD_EXAMPLES=OFF \
    -DBUILD_BENCHMARKS=OFF \
    -DBUILD_SHARED_LIBS=OFF \
    -DBUILD_TESTING=OFF \
    -DGFLAGS=OFF \
    -DSUITESPARSE=OFF \
    -DBLAS_LIBRARIES='/usr/local/lib/libopenblas.dylib' \
    .. \
    && make -j ${MAKEJ} install \
    && popd >/dev/null 2>&1 \
    && popd >/dev/null 2>&1

# Astropy caching...

echo "Attempting to trigger astropy IERS download..."

python3 -c '
from astropy.utils.iers import IERS_Auto
columns = ["year", "month", "day", "MJD", "PM_x", "PM_y", "UT1_UTC"]
iers_table = IERS_Auto.open()[columns].as_array()
'

echo "Done with IERS download."
