#!/bin/bash
#
# This script runs the "usual" command to repair wheels, but adds the
# build directory to the library search path so that the spt3g / so3g
# libraries can be found
#

set -e

dest_dir=$1
wheel=$2

# Location of this script
pushd $(dirname $0) >/dev/null 2>&1
scriptdir=$(pwd)
popd >/dev/null 2>&1

# On Linux, we need to add this to LD_LIBRARY_PATH
spt3g_install=$(ls -d ${scriptdir}/../build/temp.*/spt3g_install)
export LD_LIBRARY_PATH="/usr/local/lib":"${spt3g_install}/lib":"${spt3g_install}/lib64":${LD_LIBRARY_PATH}

auditwheel repair -w ${dest_dir} ${wheel}
