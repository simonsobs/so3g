====
so3g
====

.. image:: https://img.shields.io/github/workflow/status/simonsobs/so3g/Build%20Official%20Docker%20Images/master
    :target: https://github.com/simonsobs/so3g/actions?query=workflow%3A%22Build+Official+Docker+Images%22
    :alt: GitHub Workflow Status (branch)

.. image:: https://readthedocs.org/projects/so3g/badge/?version=latest
    :target: https://so3g.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://coveralls.io/repos/github/simonsobs/so3g/badge.svg?branch=master
    :target: https://coveralls.io/github/simonsobs/so3g?branch=master

.. image:: https://img.shields.io/badge/dockerhub-latest-blue
    :target: https://hub.docker.com/r/simonsobs/so3g/tags

Glue functions and new classes for SO work in the spt3g paradigm.

Environment Setup
=================

Before installing SO software, make sure you know what python
environment you will be using:

- Are you working on Linux or MacOS?

- Are you using a "python3" executable provided by your OS or one
  provided by another source (Anaconda, conda-forge / miniforge,
  homebrew, macports, etc)?

- Are you going to be actively developing so3g or just installing
  and using it?

After you have determined the answers to these questions, you can
set up your working environment.

Using Conda
-----------

If you already have a conda installation (a.k.a. conda "base" or "root"
environment) that is recent, then you can use that to create an
environment.  First, verify some info about your installation::

  which python3
  python3 --version
  which conda

Your python version should be at least 3.7.0.  Does the location of python3
match the location of the conda command (are they in the same bin
directory)?  If so, then you are ready.  If you do not have conda installed
but would like to use it, you might consider installing the "miniforge"
root environment (https://github.com/conda-forge/miniforge) which is
configured to get packages from the conda-forge channel.

The next step is to create a dedicated conda environment for your SO work.
This allows us to create and delete these environments as needed without
messing up the root environment::

  conda create -n simons # <- Only do this once
  conda activate simons

Now install as many dependencies as possible from conda packages.  These
are listed in a text file in the top of this git repo::

  conda install --file conda_deps.txt

Using a Virtualenv
------------------

If you are using a python3 provided by your OS or somewhere else, you
can work inside a "virtualenv".  This is like a sandbox where you can
install packages for working on one project and you can always just
wipe the directory and start over if something gets messed up.  We
will create a virtualenv in our home directory::

  python3 -m venv ${HOME}/simons # <- Just do this once
  source ${HOME}/simons/bin/activate

Now install some basic packages and then all of our requirements::

  python3 -m pip install --upgrade pip setuptools wheel
  python3 -m pip install -r requirements.txt

Other Python Packages
----------------------

If you will be using the pointing code in so3g, install pixell and qpoint
with pip (regardless of whether you are using a conda env or virtualenv)::

  pip install pixell
  pip install git+https://github.com/arahlin/qpoint


Installing Pre-Built Wheels
===========================

If you are just using so3g (not modifying or developing it), then you can
install the latest release of the package with::

  pip install so3g

This command should be used regardless of whether you are working in a
conda env or a virtualenv.

Now any time you activate your virtualenv / conda env, you can use so3g.


Building From Source
====================

If you will be developing so3g or want more control over the build, then
you should build from source.  You will need to install boost, FLAC, and
a BLAS/LAPACK library.

Prerequisites on Linux
----------------------

The easiest approach in this case is to use your OS package manager.  For
example::

  apt install \
  libboost-all-dev \
  libopenblas-openmp-dev \
  libflac-dev

Next, download and install spt3g_software.  Check the major / minor version
of your python (e.g. 3.7, 3.8 or 3.9).  We use that information to install
spt3g into our virtualenv or conda environment.  Below we assume that our
environment is in our home directory in a folder called "simons" and that
we are using python3.9::

  cd spt3g_software
  mkdir -p build
  cd build
  cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="gcc" \
    -DCMAKE_CXX_COMPILER="g++" \
    -DCMAKE_C_FLAGS="-O3 -g -fPIC" \
    -DCMAKE_CXX_FLAGS="-O3 -g -fPIC -std=c++11" \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -DPython_EXECUTABLE:FILEPATH=$(which python3) \
    -DPYTHON_MODULE_DIR="${HOME}/simons/lib/python3.9/site-packages" \
    -DCMAKE_INSTALL_PREFIX="${HOME}/simons"
  make -j 2 install

Prerequisites on MacOS
----------------------

The so3g / spt3g_software does not seem to run on MacOS when built with the
clang++ compiler (unit tests fail with cereal error).  Instead, we will use
homebrew to install our dependencies and the latest gcc compiler tools.

  brew install \
  flac \
  bzip2 \
  netcdf \
  sqlite3 \
  boost-python3 \
  gcc

Next, download and install spt3g_software.  Check the major / minor version
of your python (e.g. 3.7, 3.8 or 3.9).  We use that information to install
spt3g into our virtualenv or conda environment.  Below we assume that our
environment is in our home directory in a folder called "simons" and that
we are using python3.9.  We further assume that the homebrew gcc version
is called "gcc-11".  Also, this assumes that homebrew is installing things
to /usr/local::

  cd spt3g_software
  mkdir -p build
  cd build
  cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="gcc-11" \
    -DCMAKE_CXX_COMPILER="g++-11" \
    -DCMAKE_C_FLAGS="-O3 -g -fPIC" \
    -DCMAKE_CXX_FLAGS="-O3 -g -fPIC -std=c++11" \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -DBOOST_ROOT="/usr/local" \
    -DPython_EXECUTABLE:FILEPATH=$(which python3) \
    -DPYTHON_MODULE_DIR="${HOME}/simons/lib/python3.9/site-packages" \
    -DCMAKE_INSTALL_PREFIX="${HOME}/simons"
  make -j 2 install

Compilation and Installation
----------------------------

To compile and install the so3g package to our virtualenv / conda env, run::

  cd so3g
  mkdir -p build
  cd build
  cmake \
    -DCMAKE_PREFIX_PATH=/path/to/spt3g_software/build \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -DPython_EXECUTABLE:FILEPATH=$(which python3) \
    -DCMAKE_INSTALL_PREFIX="${HOME}/simons" \
    ..
  make -j 2 install

The definition of `CMAKE_PREFIX_PATH` must point to the build
directory for `spt3g`, because cmake output there will be used to
generate best compilation and/or linking instructions for Boost and
other dependencies of spt3g/so3g.


Local configuration through local.cmake
---------------------------------------

Optional, site-specific parameters may be set in the file local.cmake.
Lines declaring set(VARIABLE, value) should have the same effect as
passing -DVARIABLE=value to the cmake invocation.

To change the destination directory for the installation, add a line
like this one::

  set(PYTHON_INSTALL_DEST $ENV{HOME}/.local/lib/python3.7/site-packages/)

To point cmake to the spt3g build directory, add a line like this
one::

  set(CMAKE_PREFIX_PATH $ENV{HOME}/code/spt3g_software/build)


Testing
=======

The unit tests are not installed with the so3g package, so in order to run
them you must have a git checkout of so3g (even if you installed so3g from
a pre-built wheel).

After installing the so3g package, you can run the unit tests by passing the
path to the test directory to the pytest command.  You should do this from
**outside** the so3g git checkout, since otherwise the "so3g" in the source
tree may be imported rather than the installed package::

  pytest /path/to/so3g/test

You can run specific tests by calling them directly::

  python3 -m unittest /path/to/so3g/test/test_indexed
