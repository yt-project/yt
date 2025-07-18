import glob
import os
import sys
from collections import defaultdict
from distutils.ccompiler import get_default_compiler
from importlib import resources as importlib_resources

from setuptools import Distribution, setup

# ensure enclosing directory is in PYTHON_PATH to allow importing from setupext.py
if (script_dir := os.path.dirname(__file__)) not in sys.path:
    sys.path.insert(0, script_dir)

from setupext import (
    NUMPY_MACROS,
    check_CPP14_flags,
    check_for_openmp,
    check_for_pyembree,
    create_build_ext,
    get_python_include_dirs,
    install_ccompiler,
)

install_ccompiler()

if os.path.exists("MANIFEST"):
    os.remove("MANIFEST")

with open("README.md") as file:
    long_description = file.read()

CPP14_CONFIG = defaultdict(
    lambda: check_CPP14_flags(["-std=c++14", "-std=c++1y", "-std=gnu++0x"]),
    {"msvc": ["/std:c++14"]},
)
CPP11_CONFIG = defaultdict(lambda: ["-std=c++11"], {"msvc": ["/std:c++11"]})

_COMPILER = get_default_compiler()

omp_args, _ = check_for_openmp()

if os.name == "nt":
    std_libs = []
else:
    std_libs = ["m"]

CPP14_FLAG = CPP14_CONFIG[_COMPILER]
CPP11_FLAG = CPP11_CONFIG[_COMPILER]

FIXED_INTERP = "fixed_interpolator"

cythonize_aliases = {
    "LIB_DIR": "yt/utilities/lib/",
    "LIB_DIR_GEOM": ["yt/utilities/lib/", "yt/geometry/"],
    "LIB_DIR_GEOM_ARTIO": [
        "yt/utilities/lib/",
        "yt/geometry/",
        "yt/frontends/artio/artio_headers/",
    ],
    "STD_LIBS": std_libs,
    "EWAH_LIBS": std_libs
    + [os.path.abspath(importlib_resources.files("ewah_bool_utils"))],
    "OMP_ARGS": omp_args,
    "FIXED_INTERP": FIXED_INTERP,
    "ARTIO_SOURCE": sorted(glob.glob("yt/frontends/artio/artio_headers/*.c")),
    "CPP14_FLAG": CPP14_FLAG,
    "CPP11_FLAG": CPP11_FLAG,
}

lib_exts = [
    "yt/geometry/*.pyx",
    "yt/utilities/cython_fortran_utils.pyx",
    "yt/frontends/ramses/io_utils.pyx",
    "yt/frontends/gamer/cfields.pyx",
    "yt/utilities/lib/cykdtree/kdtree.pyx",
    "yt/utilities/lib/cykdtree/utils.pyx",
    "yt/frontends/artio/_artio_caller.pyx",
    "yt/utilities/lib/*.pyx",
]

embree_libs, embree_aliases = check_for_pyembree(std_libs)
cythonize_aliases.update(embree_aliases)
lib_exts += embree_libs

# This overrides using lib_exts, so it has to happen after lib_exts is fully defined
build_ext, sdist = create_build_ext(lib_exts, cythonize_aliases)


# Force setuptools to consider that there are ext modules, even if empty.
# See https://github.com/yt-project/yt/issues/2922 and
# https://stackoverflow.com/a/62668026/2601223 for the fix.
class BinaryDistribution(Distribution):
    """Distribution which always forces a binary package with platform name."""

    def has_ext_modules(self):
        return True


if __name__ == "__main__":
    # Avoid a race condition on fixed_interpolator.o during parallel builds by
    # building it only once and storing it in a static library.
    # See https://github.com/yt-project/yt/issues/4278 and
    # https://github.com/pypa/setuptools/issues/3119#issuecomment-2076922303
    # for the inspiration for this fix.

    # build_clib doesn't add the Python include dirs (for Python.h) by default,
    # as opposed to build_ext, so we need to add them manually.
    clib_include_dirs = get_python_include_dirs()

    # fixed_interpolator.cpp uses Numpy types
    import numpy

    clib_include_dirs.append(numpy.get_include())

    fixed_interp_lib = (
        FIXED_INTERP,
        {
            "sources": ["yt/utilities/lib/fixed_interpolator.cpp"],
            "include_dirs": clib_include_dirs,
            "define_macros": NUMPY_MACROS,
        },
    )

    setup(
        cmdclass={"sdist": sdist, "build_ext": build_ext},
        distclass=BinaryDistribution,
        libraries=[fixed_interp_lib],
        ext_modules=[],  # !!! We override this inside build_ext above
    )
