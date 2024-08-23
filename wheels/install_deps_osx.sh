#!/bin/bash
#
# This script is designed to run within a container managed by cibuildwheel.
# This will use a recent version of OS X.
#

set -e

# Note:  we are not cross-compiling here.  This is for selecting the
# openblas tarball to fetch.
arch=$1

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

if [ "x${use_gcc}" = "xyes" ]; then
    CC=gcc-12
    CXX=g++-12
    CFLAGS="-O3 -fPIC"
    CXXFLAGS="-O3 -fPIC -std=c++14"
else
    CC=clang
    CXX=clang++
    CFLAGS="-O3 -fPIC"
    CXXFLAGS="-O3 -fPIC -std=c++11 -stdlib=libc++"
fi

MAKEJ=2

PREFIX=/usr/local

# Install library dependencies with homebrew
brew install netcdf
brew install sqlite3
brew install flac

# Optionally install gcc
if [ "x${use_gcc}" = "xyes" ]; then
    brew install gcc@12
fi

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

# Install openblas from the multilib package- the same one numpy uses.

if [ "${arch}" = "macosx_arm64" ]; then
    openblas_pkg="openblas-v0.3.27-macosx_11_0_arm64-gf_5272328.tar.gz"
else
    openblas_pkg="openblas-v0.3.27-macosx_10_9_x86_64-gf_c469a42.tar.gz"
fi
openblas_url="https://anaconda.org/multibuild-wheels-staging/openblas-libs/v0.3.27/download/${openblas_pkg}"

if [ ! -e ${openblas_pkg} ]; then
    echo "Fetching OpenBLAS..."
    curl -SL ${openblas_url} -o ${openblas_pkg}
fi

echo "Extracting OpenBLAS"
tar -x -z -v -C "${PREFIX}" --strip-components 2 -f ${openblas_pkg}

# Install the gfortran (and libgcc) that was used for openblas compilation

curl -L https://github.com/MacPython/gfortran-install/raw/master/archives/gfortran-4.9.0-Mavericks.dmg -o gfortran.dmg
GFORTRAN_SHA256=$(shasum -a 256 gfortran.dmg)
KNOWN_SHA256="d2d5ca5ba8332d63bbe23a07201c4a0a5d7e09ee56f0298a96775f928c3c4b30  gfortran.dmg"
if [ "$GFORTRAN_SHA256" != "$KNOWN_SHA256" ]; then
    echo sha256 mismatch
    exit 1
fi

hdiutil attach -mountpoint /Volumes/gfortran gfortran.dmg
sudo installer -pkg /Volumes/gfortran/gfortran.pkg -target /
otool -L /usr/local/gfortran/lib/libgfortran.3.dylib

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

# Astropy caching...

echo "Attempting to trigger astropy IERS download..."

python3 -c '
from astropy.utils.iers import IERS_Auto
columns = ["year", "month", "day", "MJD", "PM_x", "PM_y", "UT1_UTC"]
iers_table = IERS_Auto.open()[columns].as_array()
'

echo "Done with IERS download."
