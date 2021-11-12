# This setup.py file simply builds so3g using the underlying cmake build
# system.  This is only preferred in certain cases where the automation is
# easier from a setup.py (e.g. readthedocs, pip, etc).

import os
import sys
import re
import subprocess as sp
import glob
import shutil
from pathlib import Path
import datetime

import numpy

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.command.clean import clean

try:
    from numpy.distutils.ccompiler import CCompiler_compile
    import distutils.ccompiler

    distutils.ccompiler.CCompiler.compile = CCompiler_compile
except ImportError:
    print("Numpy distutils not found, parallel compile not available")

# Absolute path to the directory with this file
topdir = Path(__file__).resolve().parent

# The version of spt3g we will be installing
upstream_spt3g_version = "1341ea5fc1584f2fce2454a9032ad0abe03bfe89"

# The name of the spt3g source and package dirs
spt3g_pkg_dir = os.path.join(topdir, "python", "spt3g_internal")
spt3g_src_dir = os.path.join(topdir, "spt3g_software")


def get_version():
    # Call the same python function used by cmake to get the version
    ver = None
    try:
        sys.path.insert(0, os.path.abspath(topdir))
        from version_h import get_versions as ver_function

        ver_info = ver_function()
        ver = ver_info["version"]
        sys.path.pop(0)
    except:
        raise RuntimeError("Cannot call get_versions() from version_h.py!")
    return ver


def get_spt3g():
    if not os.path.isdir(spt3g_pkg_dir):
        # First, download a tarball with our desired version.
        ziproot = "spt3g_software-{}".format(upstream_spt3g_version)
        zipfile = "{}.zip".format(ziproot)
        if not os.path.isfile(zipfile):
            sp.check_call(
                [
                    "curl",
                    "-SL",
                    "https://github.com/CMB-S4/spt3g_software/archive/{}.zip".format(
                        upstream_spt3g_version
                    ),
                    "-o",
                    zipfile,
                ]
            )
        sp.check_call(["unzip", "-x", zipfile])
        os.rename(ziproot, spt3g_src_dir)

        # Apply a patch with any changes
        start_dir = os.getcwd()
        os.chdir(spt3g_src_dir)
        sp.check_call(
            ["patch", "-p1", "-i", os.path.join("..", "wheels", "spt3g.patch")]
        )
        os.chdir(start_dir)

        # The code organization and build of spt3g does not follow python norms.
        # Instead, we create a package directory containing all the pure python code.
        # Then we separately build the compiled extensions.
        os.makedirs(spt3g_pkg_dir)
        with open(os.path.join(spt3g_pkg_dir, "__init__.py"), "w") as f:
            f.write("# Bundled spt3g package.  This software was obtained from:\n")
            f.write("#\n")
            f.write("#   https://github.com/CMB-S4/spt3g_software\n")
            f.write("#\n")
            f.write(
                "# At {}, and uses git version '{}'\n".format(
                    datetime.datetime.now(), upstream_spt3g_version
                )
            )
            f.write("#\n")
        for subpkg in ["calibration", "core", "dfmux", "gcp", "maps"]:
            subsrc = os.path.join(spt3g_src_dir, subpkg, "python")
            subdest = os.path.join(spt3g_pkg_dir, subpkg)
            shutil.copytree(subsrc, subdest)


# The spt3g directory needs to be in place before defining extensions below.
get_spt3g()


class RealClean(clean):
    """Really clean up.

    Delete all temporary build directories when running `python setup.py clean`.
    """

    def run(self):
        super().run()
        clean_files = [
            "./build",
            "./dist",
            "./__pycache__",
            "./*.egg-info",
            spt3g_pkg_dir,
            spt3g_src_dir,
            "./include/_version.h",
        ]
        for cf in clean_files:
            # Make paths absolute and relative to this path
            apaths = glob.glob(os.path.abspath(cf))
            for path in apaths:
                if os.path.isdir(path):
                    shutil.rmtree(path)
                elif os.path.isfile(path):
                    os.remove(path)
        return


def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile

    devnull = None
    oldstderr = None
    try:
        with tempfile.NamedTemporaryFile("w", suffix=".cpp") as f:
            f.write("int main (int argc, char **argv) { return 0; }")
            try:
                devnull = open("/dev/null", "w")
                oldstderr = os.dup(sys.stderr.fileno())
                os.dup2(devnull.fileno(), sys.stderr.fileno())
                compiler.compile([f.name], extra_postargs=[flagname])
            except:
                return False
            return True
    finally:
        if oldstderr is not None:
            os.dup2(oldstderr, sys.stderr.fileno())
        if devnull is not None:
            devnull.close()


def cpp11_flag(compiler):
    """Return the -std=c++11 compiler flag."""
    if has_flag(compiler, "-std=c++11"):
        return "-std=c++11"
    else:
        raise RuntimeError("Unsupported compiler -- at least C++11 support is needed!")


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""

    c_opts = {
        "msvc": ["/EHsc"],
        "unix": ["-DBOOST_PYTHON_MAX_ARITY=20"],
    }

    if sys.platform.lower() == "darwin":
        c_opts["unix"] += ["-stdlib=libc++", "-mmacosx-version-min=10.7"]

    def run(self):
        """Set up version header before building any extensions."""
        # get_spt3g()
        sp.check_call(
            [
                sys.executable,
                os.path.join(topdir, "version_h.py"),
                "SO3G_VERSION_STRING",
                os.path.join(topdir, "include", "_version.h"),
            ]
        )
        super().run()

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        pyver = "{}.{}".format(sys.version_info[0], sys.version_info[1])
        linkopts = [
            "-lboost_system",
            "-lboost_iostreams",
            "-lboost_filesystem",
            "-lboost_python38",
            "-lboost_regex",
            "-lFLAC",
            "-lnetcdf",
            # "-lpython{}".format(pyver),
        ]

        if ct == "unix":
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp11_flag(self.compiler))
            if has_flag(self.compiler, "-fvisibility=hidden"):
                opts.append("-fvisibility=hidden")
            if has_flag(self.compiler, "-fopenmp"):
                opts.append("-fopenmp")
                linkopts.append("-fopenmp")
            if sys.platform.lower() == "darwin":
                linkopts.append("-stdlib=libc++")
                linkopts.append("-framework Accelerate")
            else:
                linkopts.append("-lopenblas")
        elif ct == "msvc":
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())

        # Build our extensions
        for ext in self.extensions:
            ext.extra_compile_args.extend(opts)
            ext.extra_link_args.extend(linkopts)
        build_ext.build_extensions(self)


# Define all our extensions, including the spt3g ones.

# The various spt3g subdirectories have cross dependencies on included headers.
# Here we just make a list of all of them and use those for all extensions.
spt3g_includes = [numpy.get_include()]
for g3sub in ["core", "dfmux", "gcp", "maps", "calibration"]:
    spt3g_includes.append(os.path.join(spt3g_src_dir, g3sub, "include"))
    spt3g_includes.append(os.path.join(spt3g_src_dir, g3sub, "include", g3sub))

spt3g_sources = dict()
for g3sub in ["core", "dfmux", "gcp", "maps", "calibration"]:
    spt3g_sources[g3sub] = glob.glob(os.path.join(spt3g_src_dir, g3sub, "src", "*.c*"))
    if g3sub == "core" and sys.platform.lower() != "darwin":
        # Special exception for an Apple-specific source file...
        spt3g_sources[g3sub].remove(
            os.path.join(spt3g_src_dir, "core", "src", "ApplePthreadBarrier.cxx")
        )

ext_modules = list()
for g3sub in ["core", "dfmux", "gcp", "maps", "calibration"]:
    ext_modules.append(
        Extension(
            "so3g.spt3g_internal.{}.lib{}".format(g3sub, g3sub),
            spt3g_sources[g3sub],
            include_dirs=[os.path.join(spt3g_src_dir, g3sub, "src")] + spt3g_includes,
            language="c++",
        )
    )

# For libso3g, we include the spt3g/core objects directly, rather than trying to link
# to a compiled extension.
libso3g_sources = glob.glob(os.path.join("src", "*.cxx"))
libso3g_sources.extend(spt3g_sources["core"])

ext_modules.append(
    Extension(
        "so3g.libso3g",
        libso3g_sources,
        include_dirs=["include"] + spt3g_includes,
        language="c++",
    ),
)

# Install the python scripts from spt3g
scripts = glob.glob(os.path.join(spt3g_src_dir, "*", "bin", "*"))


def readme():
    with open("README.rst") as f:
        return f.read()


conf = dict()
conf["name"] = "so3g"
conf["description"] = "Tools for Simons Observatory work with spt3g_software"
conf["long_description"] = readme()
conf["long_description_content_type"] = "text/rst"
conf["author"] = "Simons Observatory Collaboration"
conf["author_email"] = "so_software@simonsobservatory.org"
conf["license"] = "MIT"
conf["url"] = "https://github.com/simonsobs/so3g"
conf["version"] = get_version()
conf["python_requires"] = ">=3.6.0"
conf["setup_requires"] = (["wheel"],)
conf["install_requires"] = [
    "numpy",
    "astropy",
    "matplotlib",
    "ephem",
    "pytz",
    "pyaml",
    "sqlalchemy",
    "pysqlite3",
]

# Since the so3g python package is in a directory called "python", we can't use the
# normal find_packages() function to recursively set these up.  Instead we specify them
# manually.

conf["packages"] = ["so3g", "so3g.spt3g_internal"]
conf["package_dir"] = {
    "so3g": "python",
    "so3g.spt3g_internal": os.path.join("python", "spt3g_internal"),
}
for sub in ["core", "dfmux", "gcp", "maps", "calibration"]:
    psub = "so3g.spt3g_internal.{}".format(sub)
    pdir = os.path.join("python", "spt3g_internal", sub)
    conf["packages"].append(psub)
    conf["package_dir"][psub] = pdir
for sub in ["hk", "proj", "smurf"]:
    psub = "so3g.{}".format(sub)
    pdir = os.path.join("python", sub)
    conf["packages"].append(psub)
    conf["package_dir"][psub] = pdir

# conf["packages"] = find_packages("so3g")
# conf["package_dir"] = {"": "so3g"}
conf["ext_modules"] = ext_modules
conf["scripts"] = scripts
conf["cmdclass"] = {"build_ext": BuildExt, "clean": RealClean}
conf["zip_safe"] = False
conf["classifiers"] = [
    "Development Status :: 5 - Production/Stable",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Operating System :: POSIX",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering :: Astronomy",
]

setup(**conf)
