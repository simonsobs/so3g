#!/bin/bash
#
# This script creates a source distribution.
#

# Location of this script
pushd $(dirname $0) >/dev/null 2>&1
topdir=$(pwd)
popd >/dev/null 2>&1

mkdir -p dist
rm -f dist/*

python setup.py sdist
