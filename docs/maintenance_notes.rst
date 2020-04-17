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

