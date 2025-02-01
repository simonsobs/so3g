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

spt3g_install=$(ls -d /project/build/temp.*/spt3g_install)

export DYLD_LIBRARY_PATH="/usr/local/lib":"${spt3g_install}/lib":"${spt3g_install}/lib64":${DYLD_LIBRARY_PATH}

delocate-listdeps ${wheel} \
&& delocate-wheel --require-archs ${delocate_archs} -w ${dest_dir} ${wheel}
