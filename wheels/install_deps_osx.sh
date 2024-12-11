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
    CC=gcc-14
    CXX=g++-14
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
    brew install gcc@14
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
CC="${CC}" CFLAGS="${CFLAGS}" pip install -v -r "${scriptdir}/build_requirements.txt" "numpy<${numpy_ver}" "scipy_openblas32"

# We use the scipy openblas wheel to get the openblas to use.

# First ensure that pkg-config is set to search somewhere
if [ -z "${PKG_CONFIG_PATH}" ]; then
    export PKG_CONFIG_PATH="/usr/local/lib/pkgconfig"
fi

python3 -c "import scipy_openblas32; print(scipy_openblas32.get_pkg_config())" > ${PKG_CONFIG_PATH}/scipy-openblas.pc

# To help delocate find the libraries, we copy them into /usr/local
python3 <<EOF
import os, scipy_openblas32, shutil
srcdir = os.path.dirname(scipy_openblas32.__file__)
incdir = os.path.join(srcdir, "include")
libdir = os.path.join(srcdir, "lib")
shutil.copytree(libdir, "/usr/local/lib", dirs_exist_ok=True)
shutil.copytree(incdir, "/usr/local/include", dirs_exist_ok=True)
EOF

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
    && CC="${CC}" CFLAGS="-O3 -fPIC" ./configure --prefix="${PREFIX}" \
    && make -j ${MAKEJ} install \
    && popd >/dev/null 2>&1

# Astropy caching...

echo "Attempting to trigger astropy IERS download..."

python3 -c '
from astropy.utils.iers import IERS_Auto
columns = ["year", "month", "day", "MJD", "PM_x", "PM_y", "UT1_UTC"]
iers_table = IERS_Auto.open()[columns].as_array()
'

echo "Done with IERS download."
