#
# Based on spt3g.core.load_bindings.
#
import platform, sys, imp, os
import glob

if platform.system().startswith('freebsd') or platform.system().startswith('FreeBSD'):
    # C++ modules are extremely fragile when loaded with RTLD_LOCAL,
    # which is what Python uses on FreeBSD by default, and maybe other
    # systems. Convince it to use RTLD_GLOBAL.

    # See thread by Abrahams et al:
    # http://mail.python.org/pipermail/python-dev/2002-May/024074.html
    sys.setdlopenflags(0x102)

def load_pybindings(paths, name=None, lib_suffix=None):
    """
    Load all non-private items from the libraries in the list "paths".
    Provide the full path to each library, but without extension.  The
    .so or .dylib will be appended depending on the system
    architecture.  The namespace into which the items are imported
    will be determined from the first path, unless name= is explicitly
    provided.
    """
    if lib_suffix is None:
        if platform.system().startswith('Darwin'):
            # OSX compatibility requires .dylib suffix
            lib_suffix = ".dylib"
        else:
            lib_suffix = ".so"
    for path in paths:
        if name is None:
            name = os.path.split(path)[1]
        # The shared library may have other identifiers if built with setuptools
        fullpath = glob.glob("{}*{}".format(path, lib_suffix))[0]
        # Save copy of current module def
        mod = sys.modules[name]
        m = imp.load_dynamic(name, fullpath)
        sys.modules[name] = mod # Don't override Python mod with C++

        for (k,v) in m.__dict__.items():
            if not k.startswith("_"):
                mod.__dict__[k] = v
