import glob
import os
from collections import defaultdict
from distutils.ccompiler import get_default_compiler

from setuptools import Distribution, setup

from setupext import (
    check_CPP14_flags,
    check_for_openmp,
    check_for_pyembree,
    create_build_ext,
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
CPP03_CONFIG = defaultdict(lambda: ["-std=c++03"], {"msvc": ["/std:c++03"]})

_COMPILER = get_default_compiler()

omp_args, _ = check_for_openmp()

if os.name == "nt":
    std_libs = []
else:
    std_libs = ["m"]

CPP14_FLAG = CPP14_CONFIG[_COMPILER]
CPP03_FLAG = CPP03_CONFIG[_COMPILER]

cythonize_aliases = {
    "LIB_DIR": "yt/utilities/lib/",
    "LIB_DIR_EWAH": ["yt/utilities/lib/", "yt/utilities/lib/ewahboolarray/"],
    "LIB_DIR_GEOM": ["yt/utilities/lib/", "yt/geometry/"],
    "LIB_DIR_GEOM_ARTIO": [
        "yt/utilities/lib/",
        "yt/geometry/",
        "yt/frontends/artio/artio_headers/",
    ],
    "STD_LIBS": std_libs,
    "OMP_ARGS": omp_args,
    "FIXED_INTERP": "yt/utilities/lib/fixed_interpolator.cpp",
    "ARTIO_SOURCE": glob.glob("yt/frontends/artio/artio_headers/*.c"),
    "CPP14_FLAG": CPP14_FLAG,
    "CPP03_FLAG": CPP03_FLAG,
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
    setup(
        cmdclass={"sdist": sdist, "build_ext": build_ext},
        distclass=BinaryDistribution,
        ext_modules=[],  # !!! We override this inside build_ext above
    )
