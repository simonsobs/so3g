====
so3g
====

.. image:: https://travis-ci.com/simonsobs/so3g.svg?branch=master
    :target: https://travis-ci.com/simonsobs/so3g

.. image:: https://coveralls.io/repos/github/simonsobs/so3g/badge.svg?branch=master
    :target: https://coveralls.io/github/simonsobs/so3g?branch=master

.. image:: https://img.shields.io/badge/dockerhub-latest-blue
    :target: https://hub.docker.com/r/simonsobs/so3g/tags

Glue functions and new classes for SO work in the spt3g paradigm.
This might turn into a permanent part of the SO software stack... but
for now let's treat it as an experimentation and development area for
things that require a high level of integration with spt3g_software.

Requirements
============

- spt3g_software and its dependencies.

  - N.B. When compiling ``spt3g_software``, you may have to explicitly
    tell it to use python3 via ``cmake`` in order for it to run
    properly with `so3g`. The invocation is:
    ``cmake .. -DPYTHON_EXECUTABLE=`which python3```

- (Optional) ``environment-modules`` - Load ``spt3g`` environment
  automatically. For details see the README in `modules/`_

.. _modules/: ./modules


Compilation and Installation
============================

The installation system is a bit different from spt3g in that there is
a separate "install" target defined.  The "install" step is required;
the install destination is customizable through the local.cmake file;
see below.

To compile the library run::

  mkdir build
  cd build
  cmake ..
  make
  make install

The build process will try to find boost, python, and spt3g.  For
spt3g, environment variable SPT3G_SOFTWARE_PATH and
SPT3G_SOFTWARE_BUILD_PATH must both be defined.  The first variable
should point to the root of spt3g repository, for finding header
files.  The second variable should point to the cmake build directory,
for finding libraries.


Local configuration through local.cmake
---------------------------------------

Optional, site-specific parameters may be set in the file local.cmake.
Lines declaring set(VARIABLE, value) should have the same effect as
passing -DVARIABLE=value to the cmake invocation.

To change the destination directory for the installation, add a line
like this one::

  set(PYTHON_INSTALL_DEST $ENV{HOME}/.local/lib/python3.7/site-packages/)

If you need to hard-code the boost python package name, add a line
like this one::

  set(Boost_PYTHON_TYPE python-py35)

Testing
=======
We use the built-in python unittest for testing. To run the tests install so3g
and run::

    cd test/
    python3 -m unittest

You can run specific tests by calling them directly::

    python3 -m unittest test_indexed
