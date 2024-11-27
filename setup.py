# This setup.py file simply builds so3g using the underlying cmake build
# system.  This is only preferred in certain cases where the automation is
# easier from a setup.py (e.g. readthedocs, pip, etc).

import os
import sys
import sysconfig
import re
import subprocess as sp
import glob
import shutil
from pathlib import Path

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.command.clean import clean

import numpy as np

# Absolute path to the directory with this file
topdir = Path(__file__).resolve().parent

# The version of spt3g we will be installing.  Get this from the
# Dockerfile for consistency.
def get_spt3g_version():
    dockerfile = os.path.join(topdir, "Dockerfile")
    ver = None
    linepat = re.compile(r".*simonsobs/spt3g:(.*)\s*")
    verpat = re.compile(r".*-g(.*)")
    with open(dockerfile, "r") as f:
        for line in f:
            mat = linepat.match(line)
            if mat is not None:
                fullver = mat.group(1)
                vermat = verpat.match(fullver)
                if vermat is None:
                    # This must be an actual tag
                    ver = fullver
                else:
                    # Extract the short hash
                    ver = vermat.group(1)
    return ver

upstream_spt3g_version = get_spt3g_version()
print(f"Using upstream spt3g_software version {upstream_spt3g_version}")

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


# Define some helper functions to do the actual fetch / build / install


def get_spt3g():
    # if os.path.isdir(spt3g_pkg_dir):
    #     return
    # We use git to get the repo, since spt3g uses git to get its version
    # information.
    print("DEBUG get_spt3g start", flush=True)
    if not os.path.isdir(spt3g_src_dir):
        print("DEBUG get_spt3g downloading", flush=True)
        sp.check_call(
            [
                "git",
                "clone",
                "https://github.com/CMB-S4/spt3g_software",
                spt3g_src_dir,
            ]
        )
        sp.check_call(
            [
                "git",
                "-C",
                spt3g_src_dir,
                "checkout",
                "-b",
                upstream_spt3g_version,
                upstream_spt3g_version,
            ]
        )
        # Apply patches with any changes
        patches = glob.glob(f"{topdir}/wheels/spt3g*.patch")
        for patch_file in patches:
            if os.path.isfile(patch_file):
                start_dir = os.getcwd()
                os.chdir(spt3g_src_dir)
                sp.check_call(["patch", "-p1", "-i", patch_file])
                os.chdir(start_dir)
    else:
        print("DEBUG get_spt3g found existing source tree, skipping", flush=True)
    print("DEBUG get_spt3g end", flush=True)


def extract_cmake_env(varprefix):
    cmake_opts = list()
    cpat = re.compile(f"{varprefix}_(.*)")
    for k, v in os.environ.items():
        mat = cpat.match(k)
        if mat is not None:
            cmake_opts.append(f"-D{mat.group(1)}={v}")
    return cmake_opts


def build_common(src_dir, build_dir, install_dir, cmake_extra, debug, pkg, version):
    print(f"DEBUG build_common start ({src_dir}, {build_dir}, {install_dir}, {cmake_extra}, {debug}, {pkg}, {version})", flush=True)
    cmake_args = list()
    cfg = "Debug" if debug else "Release"
    cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
    cmake_args += ["-DCMAKE_VERBOSE_MAKEFILE=ON"]
    cmake_args += [f"-DCMAKE_INSTALL_PREFIX={install_dir}"]
    cmake_args.extend(extract_cmake_env(f"{pkg}_BUILD"))
    cmake_args.extend(cmake_extra)

    build_args = ["--config", cfg]

    # Make a copy of the environment so that we can modify it
    env = os.environ.copy()

    ccomp = env.get("CC", None)
    cxxcomp = env.get("CXX", None)
    cflags = env.get("CFLAGS", None)
    cxxflags = env.get("CXXFLAGS", "")
    cxxflags = f"{cxxflags} -DVERSION_INFO='{version}'"
    # if sys.platform.lower() == "darwin":
    #     cmake_args += ["-DCMAKE_SHARED_LINKER_FLAGS='-undefined dynamic_lookup'"]

    # Add numpy includes
    numpy_inc = np.get_include()
    cxxflags += f" -I{numpy_inc}"

    env["CXXFLAGS"] = cxxflags

    if ccomp is not None:
        cmake_args += [f"-DCMAKE_C_COMPILER={ccomp}"]
    if cxxcomp is not None:
        cmake_args += [f"-DCMAKE_CXX_COMPILER={cxxcomp}"]
    if cflags is not None:
        cmake_args += [f"-DCMAKE_C_FLAGS={cflags}"]
    cmake_args += [f"-DCMAKE_CXX_FLAGS={cxxflags}"]

    if not os.path.exists(build_dir):
        os.makedirs(build_dir)

    # CMakeLists.txt is in the source dir
    cmake_list_dir = os.path.abspath(src_dir)
    print("-" * 10, f"Running {pkg} CMake", "-" * 40)
    print(f"cmake {cmake_list_dir} {' '.join(cmake_args)}")
    sp.check_call(["cmake", cmake_list_dir] + cmake_args, cwd=build_dir, env=env)

    make_j = 2
    if "CPU_COUNT" in os.environ:
        make_j = int(os.environ["CPU_COUNT"])

    print("-" * 10, f"Building {pkg}", "-" * 40)
    cmake_cmd = ["cmake", "--build", "."] + build_args + ["--", f"-j{make_j}"]
    sp.check_call(cmake_cmd, cwd=build_dir)
    cmake_cmd = ["cmake", "--install", "."] + build_args
    sp.check_call(cmake_cmd, cwd=build_dir)
    print("DEBUG build_common end {pkg}", flush=True)


def build_spt3g(src_dir, build_dir, install_dir, cmake_extra, debug):
    print("DEBUG build_spt3g start", flush=True)
    # Build spt3g with cmake, using any customizations passed through
    # environment variables named SPT3G_BUILD_*.  For example, the value
    # of "SPT3G_BUILD_BLAH" is passed to cmake as "-DBLAH=<value>".
    build_common(
        src_dir, build_dir, install_dir, cmake_extra, debug, "SPT3G", upstream_spt3g_version
    )
    print("DEBUG build_spt3g end", flush=True)


def build_so3g(src_dir, build_dir, install_dir, cmake_extra, debug):
    # Build so3g with cmake, using any customizations passed through
    # environment variables named SO3G_BUILD_*.  For example, the value
    # of "SO3G_BUILD_BLAH" is passed to cmake as "-DBLAH=<value>".
    print("DEBUG build_so3g start", flush=True)
    build_common(src_dir, build_dir, install_dir, cmake_extra, debug, "SO3G", get_version())
    print("DEBUG build_so3g end", flush=True)


# The spt3g directory needs to be in place before we start.
get_spt3g()


class RealClean(clean):
    """Really clean up.

    Delete all temporary build directories when running `python setup.py clean`.
    """

    def run(self):
        print("DEBUG RealClean.run()", flush=True)
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
            if not os.path.exists(cf):
                continue
            # Make paths absolute and relative to this path
            apaths = glob.glob(os.path.abspath(cf))
            for path in apaths:
                if os.path.isdir(path):
                    shutil.rmtree(path)
                elif os.path.isfile(path):
                    os.remove(path)
        return


class CMakeExtension(Extension):
    """
    This overrides the built-in extension class and essentially does nothing,
    since all extensions are compiled in one go by the custom build_ext class.
    """

    def __init__(self, name, sources=[]):
        print("DEBUG CMakeExtension ctor {name}, {sources}", flush=True)
        super().__init__(name=name, sources=sources)


class CMakeBuild(build_ext):
    """
    Builds the full package using CMake.
    """

    def initialize_options(self):
        super().initialize_options()
        self.cmake_build_done = False

    def run(self):
        """
        Perform build_cmake before doing the 'normal' stuff
        """
        print("DEBUG CMakeBuild.run()", flush=True)
        for extension in self.extensions:
            print("DEBUG CMakeBuild.run() extension {extension.name}", flush=True)
            if extension.name == "so3g._libso3g":
                # We just trigger this on one of the extensions.  build_cmake()
                # will actually build everything.
                self.build_cmake()
        # We DO NOT want to run the base class method.  It will try to copy files
        # to places that we don't want.
        # super().run()

    def build_cmake(self):
        print("DEBUG build_cmake start", flush=True)
        if self.cmake_build_done:
            return
        try:
            out = sp.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        # Path to build/temp.<arch>
        temp_build = Path(self.build_temp).resolve()

        # CMake build directory for so3g
        temp_so3g = os.path.join(temp_build, "so3g")

        # CMake build directory for spt3g
        temp_spt3g = os.path.join(temp_build, "spt3g")

        # Use CMake to install to the distutils build location
        install_so3g = os.path.dirname(
            Path(self.get_ext_fullpath("so3g._libso3g")).resolve().parents[0]
        )

        # Use CMake to install spt3g python code into a subdirectory of so3g, but
        # install the headers and other files to a separate location.
        install_spt3g_fake = os.path.join(temp_build, "spt3g_install")
        install_spt3g_py = install_so3g

        # # By default, the spt3g build system attempts to link to libpython, which
        # # should never be done when building wheels.  This link resolution should
        # # only be done at runtime on the target system after installation.  We
        # # have patched spt3g to not look for the python "Development" target, so
        # # here we specify the associated CMake variables directly.
        # py_exe = sys.executable
        # py_maj = sys.version_info[0]
        # py_min = sys.version_info[1]
        # # The includes vary slightly between builds and versions, so we call out
        # # to the python-config script for this.
        # out = sp.check_output(
        #     ["python3-config", "--includes"],
        #     universal_newlines=True,
        # )
        # raw_incl = out.split()[0]
        # py_incl = re.sub("-I", "", raw_incl)
        # dlist3g = [
        #     f"-DPython_EXECUTABLE={py_exe}",
        #     f"-DPython_INCLUDE_DIRS={py_incl}",
        #     f"-DPython_LIBRARIES=''",
        #     f"-DPython_RUNTIME_LIBRARY_DIRS=''",
        #     f"-DPython_LIBRARY_DIRS=''",
        #     f"-DPython_VERSION_MAJOR={py_maj}",
        #     f"-DPython_VERSION_MINOR={py_min}",
        #     "-DBoost_ARCHITECTURE=-x64",
        #     f"-DBoost_PYTHON_TYPE=python{py_maj}{py_min}",
        #     "-DBoost_DEBUG=ON",
        #     f"-DPYTHON_MODULE_DIR={install_spt3g_py}",
        # ]
        dlist3g = list()
        if "BOOST_ROOT" in os.environ:
            dlist3g.append(f"-DBOOST_ROOT={os.environ['BOOST_ROOT']}")
        if "FLAC_ROOT" in os.environ:
            # The spt3g package uses a custom FindFLAC.cmake, while so3g uses
            # the built-in one.  Override the spt3g detection.
            flcroot = os.environ["FLAC_ROOT"]
            flcext = "so"
            if sys.platform.lower() == "darwin":
                flcext = "dylib"
            dlist3g.extend(
                [
                    f"-DFLAC_LIBRARIES={flcroot}/lib/libFLAC.{flcext}",
                    f"-DFLAC_INCLUDE_DIR={flcroot}/include",
                    f"-DFLAC_FOUND=1",
                ]
            )

        print("DEBUG build_cmake calling spt3g", flush=True)

        build_spt3g(
            spt3g_src_dir,
            temp_spt3g,
            install_spt3g_fake,
            dlist3g,
            self.debug,
        )

        # Move spt3g python directory into place.  Remove any stale copy of the
        # directory.
        sh_ext = os.path.splitext(sysconfig.get_config_var("EXT_SUFFIX"))[1]

        install_spt3g_internal = os.path.join(install_so3g, "so3g", "spt3g_internal")
        if os.path.isdir(install_spt3g_internal):
            print(f"DEBUG rmtree {install_spt3g_internal}", flush=True)
            shutil.rmtree(install_spt3g_internal)
        print(f"DEBUG mv/rename {os.path.join(install_spt3g_py, 'spt3g')}, {install_spt3g_internal}", flush=True)
        os.rename(os.path.join(install_spt3g_py, "spt3g"), install_spt3g_internal)

        print("DEBUG build_cmake calling so3g", flush=True)

        build_so3g(
            topdir,
            temp_so3g,
            install_so3g,
            [
                f"-DPYTHON_INSTALL_DEST={install_so3g}",
                f"-DCMAKE_PREFIX_PATH={install_spt3g_fake}",
            ],
            self.debug,
        )
        self.cmake_build_done = True
        print("DEBUG build_cmake end", flush=True)


ext_modules = [
    CMakeExtension("so3g._libso3g"),
    # CMakeExtension("so3g.spt3g_internal._libcore"),
    # CMakeExtension("so3g.spt3g_internal._libdfmux"),
    # CMakeExtension("so3g.spt3g_internal._libcalibration"),
    # CMakeExtension("so3g.spt3g_internal._libgcp"),
    # CMakeExtension("so3g.spt3g_internal._libmaps"),
]


# Install the python scripts from spt3g
scripts = glob.glob(os.path.join(spt3g_src_dir, "*", "bin", "*"), root_dir=topdir)
scripts = [x.removeprefix(f"{topdir}/") for x in scripts]
print(f"DEBUG scripts = {scripts}", flush=True)

def readme():
    with open("README.rst") as f:
        return f.read()


conf = dict()
conf["name"] = "so3g"
conf["version"] = get_version()

# Since the so3g python package is in a directory called "python", we can't use the
# normal find_packages() function to recursively set these up.  Instead we specify them
# manually.

conf["packages"] = ["so3g",]
conf["package_dir"] = {
    "so3g": "python",
}

for sub in ["hk", "proj", "smurf"]:
    psub = f"so3g.{sub}"
    pdir = os.path.join("python", sub)
    conf["packages"].append(psub)
    conf["package_dir"][psub] = pdir

conf["ext_modules"] = ext_modules
conf["scripts"] = scripts
conf["cmdclass"] = {"build_ext": CMakeBuild, "clean": RealClean}
conf["zip_safe"] = False

print(f"DEBUG calling setup({conf})", flush=True)
setup(**conf)
