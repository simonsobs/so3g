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
yum install -y flac-devel bzip2-devel zlib-devel sqlite3-devel netcdf-devel

#export PATH=/opt/python/cp38-cp38/bin:${PATH}

# Install requirements
pip install -r requirements.txt

# Build options

CC=gcc
CXX=g++

CFLAGS="-O3 -fPIC -pthread"
CXXFLAGS="-O3 -fPIC -pthread -std=c++11"

MAKEJ=2

PREFIX=/usr

# Install boost

boost_version=1_72_0
boost_dir=boost_${boost_version}
boost_pkg=${boost_dir}.tar.bz2

echo "Fetching boost..."

if [ ! -e ${boost_pkg} ]; then
    curl -SL "https://dl.bintray.com/boostorg/release/1.72.0/source/${boost_pkg}" -o "${boost_pkg}"
fi

echo "Building boost..."

rm -rf ${boost_dir}
tar xjf ${boost_pkg} \
    && pushd ${boost_dir} \
    && echo "option jobs : ${MAKEJ} ;" >> tools/build/user-config.jam \
    && BOOST_BUILD_USER_CONFIG=tools/build/user-config.jam \
    ./bootstrap.sh \
    --with-toolset=gcc \
    --with-python=python3 \
    --prefix=${PREFIX} \
    && ./b2 --layout=tagged --user-config=./tools/build/user-config.jam \
    $(python3-config --includes | sed -e 's/-I//g' -e 's/\([^[:space:]]\+\)/ include=\1/g') \
    variant=release threading=multi link=shared runtime-link=shared install \
    && popd >/dev/null 2>&1
