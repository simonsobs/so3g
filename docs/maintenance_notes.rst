=================
Maintenance Notes
=================

Triggering readthedocs updates
==============================

readthedocs, understandably, does not make it easy to compile the
whole C++ library in order to extract docstrings.  The solution we
have developed is to use a compiled version of so3g to make mock
versions of the C++ classes, and push that mocked library to the
``readthedocs`` branch of so3g on github.

I know this is a messy tangle.  If you know how to improve it, go
ahead.

Here is how to update readthedocs:

- You have to be working in a build directory that is within the
  cloned repo, on master branch, with the up-to-date so3g installed.

  - The repo should be clean, because commands like ``git checkout``
    and ``git merge`` will be run on this tree, and those aren't
    guaranteed to work if there are local modifications.
  - The installed version of so3g must be available (in the Python
    search path) because the code will ```import so3g``` in order to
    find and extract the docstrings.

- From that build directory, run "make prep-readthedocs".  The
  following will (should) happen:

  - You will be asked to confirm merge of master into the local
    readthedocs branch.
  - Code from docs/extract_docstrings.py will be run, which will scan
    the docstrings of the installed so3g library and create a file
    python/_libso3g_docstring_shells.py that contains dummy versions
    of the C++ classes with docstrings populated.
  - The dummy class will be committed to the local branch of
    readthedocs.

- If that went ok, run ``git push readthedocs``.  github will notify
  readthedocs of the new action on our readthedocs branch, and
  readthedocs will try to build the docs.  At this point you are done
  and can return to master, ``git checkout master``.

- *Before pushing the readthedocs branch*, if you want, you can try to
  build the docs locally.  Make sure you set READTHEDOCS=True so that
  the sphinx conf.py knows to load the dummy class definitions::

    cd docs
    READTHEDOCS=True make html

  Note that building this locally will leave some annoying garbage
  lying around.  This can cause future builds/imports to fail if you
  don't clean it up.  Specifically the sphinx code will create
  importable ``so3g`` (a symlink to ``python``) and ``spt3g`` (a pure
  python tree containing a couple of spt3g classes used by so3g python
  code) at root level of the repo.  Those should be removed!::

    rm so3g
    rm -r spt3g

Python Wheels
=============

Binary wheels are built using the included setup.py, which downloads
a copy of spt3g_software and uses cmake to build both that package and
so3g.  Although you can build wheels locally, usually these are built
using github workflows.  The "wheels" directory contains scripts which
are used by the cibuildwheel package running in the github workflows.

Testing Wheels on Github CI
---------------------------

The wheels.yml workflow builds the wheels for a variety of python versions
and platforms.  It normally runs once per day and takes about an hour.
After building each wheel, it installs the wheel into a virtualenv and
runs the tests.  If you are debugging the builds, open a PR with your
changes.  Then comment out the cron entry at the top of this workflow and
un-comment the other entry which enables building for pull requests against
master.  Push commits to your PR to get things working, then swap these
entries back.

If you make changes to wheels.yml, make the same changes to deploy.yml,
which builds the wheels and uploads to PyPI whenever a git tag is made.

Testing Wheels Locally
----------------------

On linux, if you have docker installed, you can test the wheel building
locally, which saves a lot of time and commit noise.  Create a virtualenv
for testing and pip install the "cibuildwheel" package.  Then run
``./wheels/test_local_cibuildwheel.sh`` from the top of the source tree.
This will run cibuildwheel, which builds the wheels inside the same
manylinux docker container with the same scripts as used during the CI
builds.

Deploying New Wheels
--------------------

Just make a tag, and ensure that the version string of the package conforms
to the python package versioning scheme (e.g. 1.2.3, 2.3.4rc1, etc).  This
should happen as long as the git tag starts with a "v".  The deploy.yml
workflow will use an API token stored in the github repo "secret" parameters
to upload the resulting wheels and source dist with twine.
