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

If you are just "using" `so3g` and not actively modifying the source, simply install
the binary wheels from PyPI::

    pip install so3g

Building from Source
======================

When developing the `so3g` code, you will need to build from source.  There are two
methods documented here:  (1) using a conda environment to provide python and all
compiled dependencies and (2) using a virtualenv for python and OS packages for
compiled dependencies.  In both cases, the compiled dependencies include:

- A C++ compiler supporting the c++17 standard

- BLAS / LAPACK

- Pybind11

- GSL

- Ceres Solver / Eigen 3

- CMake + scikit_build_core

Building with Conda Tools
----------------------------

This method is the most reliable, since we will be using a self-consistent set
of dependencies and the same compilers that were used to build those. First,
ensure that you have a conda base environment that uses the conda-forge
channels. The easiest way to get this is to use the "mini-forge" installer
(https://github.com/conda-forge/miniforge).

Creating the Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~

When building within a conda environment, there is some one-time setup to do in order
to have all dependencies ready.  Once you have the conda "base" environment installed,
create a new environment for Simons Observatory work.  We force the python version to
3.13, since the default (3.14) is still missing some of our dependencies::

    conda create -n simons python==3.13
    conda activate simons

Now install all of our dependencies except for spt3g, which is not yet on conda-forge::

    conda install --file conda_dev_requirements.txt

Some of the above dependencies (compilers) will not be available until re-activating
the conda environment::

    conda deactivate
    conda activate simons

Next, install spt3g with pip::

    pip install spt3g

Installing SO3G
~~~~~~~~~~~~~~~~~

The so3g package now uses scikit_build_core, which runs cmake "under the hood".
**You should no longer run cmake directly**.  If your dependencies are in place and
your conda environment is activated, you can install so3g with::

    pip install -v .

If you are actively hacking on so3g, then you can install the package in "editable"
mode.  This will install symlinks that point back to your build directory.  If you edit
the source files in this mode, cmake will be triggered to rebuild on the next import of
so3g.  To install in editable mode run::

    pip install --no-build-isolation -v -e .


Building with OS Packages
----------------------------

Another option is to use a virtualenv for python packages and use the compilers and
libraries from your OS to provide so3g dependencies. Install dependencies, for example::

    apt install \
        libopenblas-openmp-dev \
        libgsl-dev \
        libceres-dev \
        libeigen3-dev

**NOTE:  Ubuntu 22.04 (for example) has a version of Ceres that is too old.**  Then
create and activate a virtualenv.  For example::

    python3 -m venv ~/env_simons
    source ~/env_simons/bin/activate

Installing SO3G
~~~~~~~~~~~~~~~~~

The so3g package now uses scikit_build_core, which runs cmake "under the hood".
**You should no longer run cmake directly**.  With your virtualenv activated you can
install so3g with::

    pip install -v .

If you are actively hacking on so3g, then you can install the package in "editable"
mode.  This will install symlinks that point back to your build directory.  If you edit
the source files in this mode, cmake will be triggered to rebuild on the next import of
so3g.  To install in editable mode run::

    pip install --no-build-isolation -v -e .

Customizing the Build
-------------------------

Build options can be changed by editing pyproject.toml, or by overriding those same
options on the command line.  For example, when debugging you might want to use Debug
mode and turn off compiler optimizations::

    pip install \
        --no-build-isolation \
        -Ccmake.build-type=Debug \
        -Ccmake.args="-DCMAKE_CXX_FLAGS='-O0 -g'" \
        -v -e .

Testing
=======

The unit tests are not installed with the so3g package, so in order to run
them you must have a git checkout of so3g (even if you installed so3g from
a pre-built wheel).

After installing the so3g package, you can run the unit tests by passing the
path to the test directory to the pytest command::

  pytest /path/to/so3g/test

You can run specific tests by calling them directly::

  pytest /path/to/so3g/test/test_indexed.py
