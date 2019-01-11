====
so3g
====

Glue functions and new classes for SO work in the spt3g paradigm.
This might turn into a permanent part of the SO software stack... but
for now let's treat it as an experimentation and development area for
things that require a high level of integration with spt3g_software.

Requirements
============

- spt3g_software and its dependencies.
- (Optional) ``environment-modules`` - Load ``spt3g`` environment automatically. For details see the README in `modules/`_

.. _modules/: ./modules


Compilation and Installation
============================

The installation system is a bit different from spt3g in that there is
a separate "install" target defined.  The "install" step is
optional, and the install destination is customizable through the
local.cmake file; see below.

To compile the library run::

  mkdir build
  cd build
  cmake ..
  make

The build process will try to find boost, python, and spt3g.  For
spt3g, environment variable SPT3G_SOFTWARE_PATH should give the path
to the spt3g_software source tree (or, possibly, anywhere where the
relevants headers and libraries could plausibly be found after digging).


Local configuration through local.cmake
=======================================

Optional, site-specific parameters may be set in the file local.cmake.

To change the destination directory for the installation, add a line
like this one::

  set(PYTHON_INSTALL_DEST $ENV{HOME}/.local/lib/python3.7/site-packages/)

