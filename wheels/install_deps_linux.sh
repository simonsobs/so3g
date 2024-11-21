#!/bin/bash
#
# This script is designed to run within a container managed by cibuildwheel.
#

set -e

# Location of this script
pushd $(dirname $0) >/dev/null 2>&1
scriptdir=$(pwd)
popd >/dev/null 2>&1
echo "Wheel script directory = ${scriptdir}"

# Install library dependencies
yum update -y
yum install -y flac-devel bzip2-devel zlib-devel sqlite-devel netcdf-devel

# Build options

CC=gcc
CXX=g++
FC=gfortran

CFLAGS="-O3 -fPIC -pthread"
CXXFLAGS="-O3 -fPIC -pthread -std=c++17"
FCFLAGS="-O3 -fPIC -pthread"

MAKEJ=2

PREFIX=/usr/local

# Update pip
pip install --upgrade pip

# Install a couple of base packages that are always required
pip install -v cmake wheel setuptools

# In order to maximize ABI compatibility with numpy, build with the newest numpy
# version containing the oldest ABI version compatible with the python we are using.
pyver=$(python3 --version 2>&1 | awk '{print $2}' | sed -e "s#\(.*\)\.\(.*\)\..*#\1.\2#")
if [ ${pyver} == "3.7" ]; then
    numpy_ver="1.20"
fi
if [ ${pyver} == "3.8" ]; then
    numpy_ver="1.20"
fi
if [ ${pyver} == "3.9" ]; then
    numpy_ver="1.20"
fi
if [ ${pyver} == "3.10" ]; then
    numpy_ver="1.22"
fi
if [ ${pyver} == "3.11" ]; then
    numpy_ver="1.24"
fi

# Install build requirements.
CC="${CC}" CFLAGS="${CFLAGS}" pip install -v -r "${scriptdir}/build_requirements.txt" "numpy<${numpy_ver}"

# Install Openblas

openblas_version=0.3.28
openblas_dir=OpenBLAS-${openblas_version}
openblas_pkg=${openblas_dir}.tar.gz

if [ ! -e ${openblas_pkg} ]; then
    echo "Fetching OpenBLAS..."
    curl -SL https://github.com/xianyi/OpenBLAS/archive/v${openblas_version}.tar.gz -o ${openblas_pkg}
fi

echo "Building OpenBLAS..."

rm -rf ${openblas_dir}
tar xzf ${openblas_pkg} \
    && pushd ${openblas_dir} >/dev/null 2>&1 \
    && make USE_OPENMP=1 NO_SHARED=1 \
    MAKE_NB_JOBS=${MAKEJ} \
    CC="${CC}" FC="${FC}" DYNAMIC_ARCH=1 TARGET=GENERIC \
    COMMON_OPT="${CFLAGS}" FCOMMON_OPT="${FCFLAGS}" \
    LDFLAGS="-fopenmp -lm" libs netlib \
    && make NO_SHARED=1 DYNAMIC_ARCH=1 TARGET=GENERIC PREFIX="${PREFIX}" install \
    && popd >/dev/null 2>&1

# Install boost

boost_version=1_86_0
boost_dir=boost_${boost_version}
boost_pkg=${boost_dir}.tar.bz2

echo "Fetching boost..."

if [ ! -e ${boost_pkg} ]; then
    curl -SL "https://archives.boost.io/release/1.86.0/source/${boost_pkg}" -o "${boost_pkg}"
fi

echo "Building boost..."

pyincl=$(for d in $(python3-config --includes | sed -e 's/-I//g'); do echo "include=${d}"; done | xargs)

rm -rf ${boost_dir}
tar xjf ${boost_pkg} \
    && pushd ${boost_dir} \
    && echo "using gcc : : ${CXX} ;" > tools/build/user-config.jam \
    && echo "option jobs : ${MAKEJ} ;" >> tools/build/user-config.jam \
    && BOOST_BUILD_USER_CONFIG=tools/build/user-config.jam \
    ./bootstrap.sh \
    --with-python=python3 \
    --prefix=${PREFIX} \
    && ./b2 --layout=tagged --user-config=./tools/build/user-config.jam \
    ${pyincl} cxxflags="${CXXFLAGS}" variant=release threading=multi link=shared runtime-link=shared install \
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
    && popd >/dev/null 2>&1

# Astropy caching...

echo "Attempting to trigger astropy IERS download..."

python3 -c '
from astropy.utils.iers import IERS_Auto
columns = ["year", "month", "day", "MJD", "PM_x", "PM_y", "UT1_UTC"]
iers_table = IERS_Auto.open()[columns].as_array()
'

echo "Done with IERS download."
