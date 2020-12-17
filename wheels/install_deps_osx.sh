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

# Install library dependencies with homebrew
brew install flac
brew install bzip2
brew install netcdf
brew install sqlite3

# Install requirements
pip install -r requirements.txt

# Build options

CC=clang
CXX=clang++

CFLAGS="-O3 -fPIC"
CXXFLAGS="-O3 -fPIC -std=c++11 -stdlib=libc++"

MAKEJ=2

PREFIX=/usr/local

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
    --with-toolset=clang \
    --with-python=python3 \
    --prefix=${PREFIX} \
    && ./b2 --layout=tagged --user-config=./tools/build/user-config.jam \
    $(python3-config --includes | sed -e 's/-I//g' -e 's/\([^[:space:]]\+\)/ include=\1/g') \
    toolset=clang cxxflags="-stdlib=libc++" linkflags="-stdlib=libc++" \
    variant=release threading=multi link=shared runtime-link=shared install \
    && popd >/dev/null 2>&1
