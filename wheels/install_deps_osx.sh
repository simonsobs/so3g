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

# Update pip
pip install --upgrade pip

# Install a couple of base packages that are always required
pip install pyaml numpy cmake

# Install requirements
pip install -r requirements.txt

# Build options

CC=clang
CXX=clang++

CFLAGS="-O3 -fPIC"
CXXFLAGS="-O3 -fPIC -std=c++11"
#CXXFLAGS="-O3 -fPIC -std=c++11 -stdlib=libc++"

MAKEJ=2

PREFIX=/usr/local

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
    && echo "using darwin : : ${CXX} ;" > tools/build/user-config.jam \
    && echo "option jobs : ${MAKEJ} ;" >> tools/build/user-config.jam \
    && BOOST_BUILD_USER_CONFIG=tools/build/user-config.jam \
    ./bootstrap.sh \
    --with-python=python3 \
    --prefix=${PREFIX} \
    && ./b2 --layout=tagged --user-config=./tools/build/user-config.jam \
    ${pyincl} -sNO_LZMA=1 -sNO_ZSTD=1 \
    cxxflags="${CXXFLAGS} -stdlib=libc++" linkflags="-stdlib=libc++" \
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
