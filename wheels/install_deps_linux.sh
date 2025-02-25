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
yum install -y bzip2-devel zlib-devel sqlite-devel

# Build options

CC=gcc
CXX=g++
FC=gfortran

CFLAGS="-O3 -fPIC -pthread"
CXXFLAGS="-O3 -fPIC -pthread -std=c++17"
FCFLAGS="-O3 -fPIC -pthread"

MAKEJ=2

PREFIX=/usr/local

# Link lib64 directory to lib
echo "Symlink /usr/local/lib64 to /usr/local/lib"
rm -rf /usr/local/lib64
ln -s /usr/local/lib /usr/local/lib64

# Update pip
python3 -m pip install --upgrade pip

# Install a couple of base packages that are always required
python3 -m pip install -v cmake wheel setuptools

pyver=$(python3 --version 2>&1 | awk '{print $2}' | sed -e "s#\(.*\)\.\(.*\)\..*#\1.\2#")

# Install build requirements.
CC="${CC}" CFLAGS="${CFLAGS}" python3 -m pip install -v -r "${scriptdir}/../requirements.txt"

# Install Openblas

openblas_version=0.3.29
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

boost_version=1_87_0
boost_dir=boost_${boost_version}
boost_pkg=${boost_dir}.tar.bz2

echo "Fetching boost..."

if [ ! -e ${boost_pkg} ]; then
    curl -SL "https://archives.boost.io/release/1.87.0/source/${boost_pkg}" -o "${boost_pkg}"
fi

echo "Building boost..."

pyincl=$(for d in $(python3-config --includes | sed -e 's/-I//g'); do echo "include=${d}"; done | xargs)

rm -rf ${boost_dir}
tar xjf ${boost_pkg} \
    && pushd ${boost_dir} \
    && echo "using gcc : : ${CXX} ;" > tools/build/user-config.jam \
    && echo "option jobs : ${MAKEJ} ;" >> tools/build/user-config.jam \
    && BOOST_BUILD_PATH=tools/build \
    ./bootstrap.sh \
    --with-python=python3 \
    --prefix=${PREFIX} \
    --with-libraries="iostreams,python,regex" \
    && ./b2 --layout=tagged --user-config=./tools/build/user-config.jam \
    ${pyincl} cxxflags="${CXXFLAGS}" variant=release threading=multi link=shared runtime-link=shared install \
    && popd >/dev/null 2>&1

# Install libFLAC

flac_version=1.5.0
flac_dir=flac-${flac_version}
flac_pkg=${flac_dir}.tar.gz

echo "Fetching libFLAC..."

if [ ! -e ${flac_pkg} ]; then
    curl -SL "https://github.com/xiph/flac/archive/refs/tags/${flac_version}.tar.gz" -o "${flac_pkg}"
fi

echo "Building libFLAC..."

rm -rf ${flac_dir}
tar xzf ${flac_pkg} \
    && pushd ${flac_dir} >/dev/null 2>&1 \
    && mkdir -p build \
    && pushd build >/dev/null 2>&1 \
    && cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="${CC}" \
    -DCMAKE_C_FLAGS="${CFLAGS}" \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -DBUILD_DOCS=OFF \
    -DWITH_OGG=OFF \
    -DBUILD_CXXLIBS=OFF \
    -DBUILD_PROGRAMS=OFF \
    -DBUILD_UTILS=OFF \
    -DBUILD_TESTING=OFF \
    -DBUILD_EXAMPLES=OFF \
    -DBUILD_SHARED_LIBS=ON \
    -DINSTALL_MANPAGES=OFF \
    -DENABLE_MULTITHREADING=ON \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    .. \
    && make -j ${MAKEJ} install \
    && popd >/dev/null 2>&1 \
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
