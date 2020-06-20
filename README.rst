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
    tell it to prefer Python3, rather than Python2, when you run
    ``cmake``.  It may be enough to do
    ``cmake .. -DPYTHON_EXECUTABLE=`which python3```.

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
  cmake .. -DCMAKE_PREFIX_PATH=${SPT3G_SOFTWARE_BUILD_PATH}
  make
  make install

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

We use the built-in python unittest for testing. To run the tests install so3g
and run::

    cd test/
    python3 -m unittest

You can run specific tests by calling them directly::

    python3 -m unittest test_indexed
