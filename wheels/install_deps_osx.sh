#!/bin/bash
#
# This script is designed to run within a container managed by cibuildwheel.
# This will use a recent version of OS X.
#

set -e

# Location of this script
pushd $(dirname $0) >/dev/null 2>&1
topdir=$(pwd)
popd >/dev/null 2>&1

# Build options.  If we use clang, then use accelerate framework.  Otherwise
# build and use OpenBLAS.

use_gcc=yes

CC=gcc-11
CXX=g++-11
FC=gfortran-11
#CC=clang
#CXX=clang++
#FC=

CFLAGS="-O3 -fPIC"
FCFLAGS="-O3 -fPIC"
# Use the second when building with clang
CXXFLAGS="-O3 -fPIC -std=c++11"
#CXXFLAGS="-O3 -fPIC -std=c++11 -stdlib=libc++"

MAKEJ=2

PREFIX=/usr/local

# Install library dependencies with homebrew
brew install flac
brew install bzip2
brew install netcdf
brew install sqlite3

# Optionally install gcc
if [ "x${use_gcc}" = "xyes" ]; then
    brew install gcc
fi

# Update pip
pip install --upgrade pip

# Install a couple of base packages that are always required
pip install pyaml numpy cmake

# Install requirements
pip install -r requirements.txt
pip install pytest

# Optionally Install Openblas

if [ "x${use_gcc}" = "xyes" ]; then
    openblas_version=0.3.13
    openblas_dir=OpenBLAS-${openblas_version}
    openblas_pkg=${openblas_dir}.tar.gz

    echo "Fetching OpenBLAS..."

    if [ ! -e ${openblas_pkg} ]; then
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
        LDFLAGS="-fopenmp -lm" \
        && make NO_SHARED=1 DYNAMIC_ARCH=1 TARGET=GENERIC PREFIX="${PREFIX}" install \
        && popd >/dev/null 2>&1
fi

# Install boost

boost_version=1_76_0
boost_dir=boost_${boost_version}
boost_pkg=${boost_dir}.tar.bz2

echo "Fetching boost..."

if [ ! -e ${boost_pkg} ]; then
    curl -SL "https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/${boost_pkg}" -o "${boost_pkg}"
fi

echo "Building boost..."

pyincl=$(for d in $(python3-config --includes | sed -e 's/-I//g'); do echo "include=${d}"; done | xargs)

use_line="using darwin : : ${CXX} ;"
extra_link="linkflags=\"-stdlib=libc++\""
if [ "x${use_gcc}" = "xyes" ]; then
    use_line="using gcc : : ${CXX} ;"
    extra_link=""
fi

rm -rf ${boost_dir}
tar xjf ${boost_pkg} \
    && pushd ${boost_dir} \
    && echo ${use_line} > tools/build/user-config.jam \
    && echo "option jobs : ${MAKEJ} ;" >> tools/build/user-config.jam \
    && BOOST_BUILD_USER_CONFIG=tools/build/user-config.jam \
    ./bootstrap.sh \
    --with-python=python3 \
    --prefix=${PREFIX} \
    && ./b2 --layout=tagged --user-config=./tools/build/user-config.jam \
    ${pyincl} -sNO_LZMA=1 -sNO_ZSTD=1 \
    cxxflags="${CXXFLAGS}" ${extra_link} \
    variant=release threading=multi link=shared runtime-link=shared install \
    && popd >/dev/null 2>&1

# Install qpoint

qpoint_version=828126de9f195f88bfaf1996527f633382457461
qpoint_dir="qpoint"

echo "Fetching qpoint..."

if [ ! -d ${qpoint_dir} ]; then
    git clone https://github.com/arahlin/qpoint.git ${qpoint_dir}
fi

echo "Building qpoint..."

pushd ${qpoint_dir} \
    && git checkout master \
    && if [ "x$(git branch -l | grep so3g)" != x ]; then \
        git branch -D so3g; fi \
    && git fetch \
    && git checkout -b so3g ${qpoint_version} \
    && CC="${CC}" CFLAGS="${CFLAGS}" \
    python3 setup.py install --prefix "${PREFIX}" \
    && popd > /dev/null
