#!/bin/bash

# This script tests the local use of setup.py on MacOS.  It assumes
# that you have homebrew installed and have the "brew" command in
# your PATH.
#
# This script takes one argument, which is the path to the virtualenv
# to use for installation.  If that path does not exist, it will be
# created.
#

set -e

# Location of the source tree
pushd $(dirname $0) >/dev/null 2>&1
topdir=$(dirname $(pwd))
popd >/dev/null 2>&1

brew_com=$(which brew)
if [ "x${brew_com}" = "x" ]; then
    echo "Homebrew must be installed and the brew command available"
    exit 1
fi

venv_path=$1
if [ "x${venv_path}" = "x" ]; then
    echo "Usage:  $0 <virtualenv path>"
    echo "  If the path to the virtualenv exists, it will be activated."
    echo "  Otherwise it will be created."
    exit 1
fi

# Deactivate any current venv
if [ "x$(type -t deactivate)" != "x" ]; then
    deactivate
fi

# Path to homebrew
brew_root=$(dirname $(dirname ${brew_com}))
echo "Using homebrew installation in ${brew_root}"

# Export compiler information
use_gcc=yes
export CC=gcc-14
export CXX=g++-14
export FC=gfortran-14
# export CC=clang
# export CXX=clang++
# export FC=
export CFLAGS="-O3 -fPIC"
export FCFLAGS="-O3 -fPIC"
# Use the second when building with clang
CXXFLAGS="-O3 -fPIC -std=c++14"
#CXXFLAGS="-O3 -fPIC -std=c++14 -stdlib=libc++"

# Install most dependencies with homebrew, including python-3.9
eval ${brew_com} install flac
eval ${brew_com} install bzip2
eval ${brew_com} install netcdf
eval ${brew_com} install sqlite3
eval ${brew_com} install openblas
eval ${brew_com} install gsl
eval ${brew_com} install boost-python3

if [ "x${use_gcc}" = "xyes" ]; then
    eval ${brew_com} install gcc
    # eval ${brew_com} install openblas
fi

# Add this homebrew python to our path
export PATH="${brew_root}/opt/python@3.9/bin:$PATH"
brew_py_ver=$(python3 --version | awk '{print $2}')
echo "Using python $(which python3), version ${brew_py_ver}"

if [ -d ${venv_path} ]; then
    echo "Virtualenv \"${venv_path}\" already exists, activating"
    source "${venv_path}/bin/activate"
    venv_py_ver=$(python3 --version | awk '{print $2}')
    if [ "${venv_py_ver}" != "${brew_py_ver}" ]; then
        echo "Virtualenv python version ${venv_py_ver} does not match"
        echo "homebrew python version ${brew_py_ver}.  Remove this"
        echo "virtualenv or use a different path."
        exit 1
    fi
else
    echo "Creating virtualenv \"${venv_path}\""
    python3 -m venv "${venv_path}"
    source "${venv_path}/bin/activate"
fi

# Update pip
python3 -m pip install --upgrade pip setuptools wheel

# Install a couple of base packages that are always required
python3 -m pip install pyaml numpy cmake

# Install requirements
python3 -m pip install -r requirements.txt
python3 -m pip install pytest
python3 -m pip install delocate

# Package location
export BOOST_ROOT="${brew_root}"
export LD_LIBRARY_PATH="${brew_root}/lib"
export DYLD_LIBRARY_PATH="${brew_root}/lib"
export CPATH="${brew_root}/include:${brew_root}/opt/openblas/include"

# Tell setup.py to look in the homebrew prefix for libraries when
# building spt3g and so3g.
export SPT3G_BUILD_CMAKE_INCLUDE_PATH="${brew_root}/include"
export SPT3G_BUILD_CMAKE_LIBRARY_PATH="${brew_root}/lib"
export SO3G_BUILD_CMAKE_INCLUDE_PATH="${brew_root}/include"
export SO3G_BUILD_CMAKE_LIBRARY_PATH="${brew_root}/lib"

export SO3G_BUILD_BLA_VENDOR="OpenBLAS"
export SO3G_BUILD_BLAS_LIBRARIES="${brew_root}/opt/openblas/lib/libopenblas.dylib"

# Now build a wheel
python3 setup.py clean
python3 -m pip wheel ${topdir} --wheel-dir=build/temp_wheels --no-deps -vvv

# The file
input_wheel=$(ls ${topdir}/build/temp_wheels/*.whl)
wheel_file=$(basename ${input_wheel})

# Repair it
spt3g_libdir=$(ls -d ${topdir}/build/temp.macosx*/spt3g_install/lib)
export DYLD_LIBRARY_PATH=${spt3g_libdir}:${DYLD_LIBRARY_PATH}

delocate-listdeps ${input_wheel} \
&& delocate-wheel -w ${topdir} ${input_wheel}

# Install it
python3 -m pip install ${topdir}/${wheel_file}
