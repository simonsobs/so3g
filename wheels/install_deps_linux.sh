#!/bin/bash
#
# This script is designed to run within a container managed by cibuildwheel.
#

set -e

# Location of this script
pushd $(dirname $0) >/dev/null 2>&1
topdir=$(pwd)
popd >/dev/null 2>&1

# Install library dependencies
yum update -y
yum install -y flac-devel bzip2-devel zlib-devel sqlite-devel netcdf-devel

#export PATH=/opt/python/cp38-cp38/bin:${PATH}

# Update pip
pip install --upgrade pip

# Install a couple of base packages that are always required
pip install pyaml numpy cmake

# Install requirements
pip install -r requirements.txt

# Build options

CC=gcc
CXX=g++
FC=gfortran

CFLAGS="-O3 -fPIC -pthread"
CXXFLAGS="-O3 -fPIC -pthread -std=c++11"
FCFLAGS="-O3 -fPIC -pthread"

MAKEJ=2

PREFIX=/usr/local

# Install Openblas

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
