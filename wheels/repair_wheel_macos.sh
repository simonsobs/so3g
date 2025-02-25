#!/bin/bash
#
# This script runs the "usual" command to repair wheels, but adds the
# build directory to the library search path so that the spt3g / so3g
# libraries can be found
#

set -e

dest_dir=$1
wheel=$2
delocate_archs=$3

# Location of this script
pushd $(dirname $0) >/dev/null 2>&1
scriptdir=$(pwd)
popd >/dev/null 2>&1

spt3g_install=$(ls -d ${scriptdir}/../build/temp.*/spt3g_install)
dylib_path="/usr/local/lib":"${spt3g_install}/lib":"${spt3g_install}/lib64":${DYLD_LIBRARY_PATH}

# DYLD_LIBRARY_PATH must be set on the command invocation line.  See:
# https://github.com/pypa/cibuildwheel/issues/816
#
DYLD_LIBRARY_PATH=${dylib_path} delocate-listdeps --all ${wheel} && \
DYLD_LIBRARY_PATH=${dylib_path} delocate-wheel -v --require-archs ${delocate_archs} -w ${dest_dir} ${wheel}

