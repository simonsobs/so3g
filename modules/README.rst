Environment Modules
===================

The spt3g_ installation instructions describe how to build the software, which
creates a script that initializes a set of environment variables allowing you
to import spt3g.

To avoid needing to run this script every time you want to use spt3g we
recommend using environment modules to setup the environment and loading the
modules in your ``.bashrc`` file.

Installation
------------

In addition to needing the ``spt3g`` dependencies you also need to install the
``environment-modules`` package, as well as the ``tcl`` package::

    $ sudo apt-get update
    $ sudo apt-get install -y tcl environment-modules

Setup
-----

To setup, copy this modules directory somewhere you would like to install it,
for the example we will use your home directory. You should then edit the first
uncommented line in ``modules/spt3g_shared`` to set the ``g3root`` to the
location you have cloned the ``spt3g_software`` to. For a user called "vagrant"
with ``spt3g_softare`` in their home directory it will look like this::

    set g3root /home/vagrant/spt3g_software

Once you have updated this, add the following lines to your ``.bashrc`` file::

    # load spt3g using environment modules
    module use --append /home/vagrant/modules
    module load spt3g_shared

Replace ``/vagrant/modules`` with the location you have copied this directory
to. You can then test by sourcing your ``.bashrc`` and trying to load spt3g::
  
    $ source ~/.bashrc
    $ python3 -c "import spt3g.core"

Shared Installation
-------------------

This can be used to make a shared installation of ``spt3g``, just have other
users add the same lines to their ``.bashrc`` file and make sure ``spt3g`` is
installed somewhere they have permissions to read it.

.. _spt3g: https://github.com/CMB-S4/spt3g_software
