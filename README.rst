====
so3g
====

.. image:: https://img.shields.io/github/actions/workflow/status/simonsobs/so3g/official-docker-images.yml?branch=master
    :target: https://github.com/simonsobs/so3g/actions?query=workflow%3A%22Build+Official+Docker+Images%22
    :alt: GitHub Workflow Status (branch)

.. image:: https://readthedocs.org/projects/so3g/badge/?version=latest
    :target: https://so3g.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://coveralls.io/repos/github/simonsobs/so3g/badge.svg?branch=master
    :target: https://coveralls.io/github/simonsobs/so3g?branch=master

.. image:: https://img.shields.io/badge/dockerhub-latest-blue
    :target: https://hub.docker.com/r/simonsobs/so3g/tags

.. image:: https://img.shields.io/pypi/v/so3g
   :target: https://pypi.org/project/so3g/
   :alt: PyPI Package

Glue functions and new classes for SO work in the spt3g paradigm.

Installation from Binary Packages
===================================

If you are just "using" `so3g` and not actively modifying the source, simply install the binary wheels from PyPI::

    pip install so3g

Building from Source
======================

When developing the `so3g` code, you will need to build from source.  There are two methods documented here:  (1) using a conda environment to provide python and all compiled dependencies and (2) using a virtualenv for python and OS packages for compiled dependencies.  In both cases, the compiled dependencies include:

- A C++ compiler supporting the c++17 standard

- BLAS / LAPACK

- Boost (at least version 1.87 for numpy-2 compatibility)

- GSL

- libFLAC

Building with Conda Tools
----------------------------

This method is the most reliable, since we will be using a self-consistent set of dependencies and the same compilers that were used to build those.  First, ensure that you have a conda base environment that uses the conda-forge channels.  The easiest way to get this is to use the "mini-forge" installer (https://github.com/conda-forge/miniforge).

Once you have the conda "base" environment installed, create a new environment for Simons Observatory work.  We force the python version to 3.12, since the default (3.13) is still missing some of our dependencies::

    conda create -n simons python==3.12 # <- Only do this once
    conda activate simons

Now install all of our dependencies (except for spt3g)::

    conda install --file conda_dev_requirements.txt

Next, choose how to install spt3g.

Bundled SPT3G
~~~~~~~~~~~~~~~~~

If you are just testing a quick change, you can use `pip` to install so3g.  This will download a copy of spt3g and bundle it into the the installed package.  The downside is that **every time** you run pip, it will re-build all of spt3g and so3g under the hood with cmake::

    pip install -vv .

Separate SPT3G
~~~~~~~~~~~~~~~~~

If you are going to be developing so3g and repeatedly building it, you probably want to install spt3g once.  See the `instructions from that package <https://github.com/CMB-S4/spt3g_software>`_ to download and install.  When building, you can install into your conda environment like this::

    cd spt3g_software
    mkdir -p build
    cd build
    cmake \
        -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} \
        -DCMAKE_C_COMPILER=${CC} \
        -DCMAKE_CXX_COMPILER=${CXX} \
        -DPython_ROOT_DIR=${CONDA_PREFIX} \
        ..
    make -j 4 install
    # Copy the python package into place
    cp -r ./spt3g ${CONDA_PREFIX}/lib/python3.12/site-packages/

When building `so3g` against a stand-alone version of `spt3g`, you need to use cmake directly::

    cd so3g
    mkdir -p build
    cd build
    cmake \
        -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} \
        -DCMAKE_C_COMPILER=${CC} \
        -DCMAKE_CXX_COMPILER=${CXX} \
        -DPython_ROOT_DIR=${CONDA_PREFIX} \
        -DBLAS_LIBRARIES='-L${CONDA_PREFIX}/lib -lopenblas -fopenmp' \
        ..
    make -j 4 install


Building with OS Packages
----------------------------

Another option is to use a virtualenv for python packages and use the compilers and
libraries from your OS to provide so3g dependencies. Install dependencies, for example::

    apt install \
        libboost-all-dev \
        libopenblas-openmp-dev \
        libflac-dev \
        libgsl-dev \
        libnetcdf-dev

Then activate your virtualenv. Next you should install to someplace in your library
search path. Note that the commands below will not work unless you change the install
prefix to a user-writable directory (or make install with sudo). You should decide where
you want to install and make sure that the location is in your PATH and
LD_LIBRARY_PATH::

    cd spt3g_software
    mkdir -p build
    cd build
    cmake \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        ..
    make -j 4 install
    #

And similarly for so3g::

    cd so3g
    mkdir -p build
    cd build
    cmake \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBLAS_LIBRARIES='-lopenblas -fopenmp' \
        ..
    make -j 4 install


Testing
=======

The unit tests are not installed with the so3g package, so in order to run
them you must have a git checkout of so3g (even if you installed so3g from
a pre-built wheel).

After installing the so3g package, you can run the unit tests by passing the
path to the test directory to the pytest command::

  pytest /path/to/so3g/test

You can run specific tests by calling them directly::

  python3 -m unittest /path/to/so3g/test/test_indexed
