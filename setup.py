import glob
import os
import sys
from distutils.ccompiler import get_default_compiler
from distutils.version import LooseVersion

import pkg_resources
from setuptools import Distribution, find_packages, setup

from setupext import (
    check_for_openmp,
    check_for_pyembree,
    create_build_ext,
    install_ccompiler,
)

install_ccompiler()

try:
    distribute_ver = LooseVersion(pkg_resources.get_distribution("distribute").version)
    if distribute_ver < LooseVersion("0.7.3"):
        print("Distribute is a legacy package obsoleted by setuptools.")
        print("We strongly recommend that you just uninstall it.")
        print("If for some reason you cannot do it, you'll need to upgrade it")
        print("to latest version before proceeding:")
        print("    pip install -U distribute")
        sys.exit(1)
except pkg_resources.DistributionNotFound:
    pass  # yay!

VERSION = "4.0.dev0"

if os.path.exists("MANIFEST"):
    os.remove("MANIFEST")

with open("README.md") as file:
    long_description = file.read()

if check_for_openmp():
    omp_args = ["-fopenmp"]
else:
    omp_args = []

if os.name == "nt":
    std_libs = []
else:
    std_libs = ["m"]

if get_default_compiler() == "msvc":
    CPP14_FLAG = ["/std:c++14"]
else:
    CPP14_FLAG = ["--std=c++14"]

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
}

lib_exts = [
    "yt/geometry/*.pyx",
    "yt/utilities/cython_fortran_utils.pyx",
    "yt/frontends/ramses/io_utils.pyx",
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
        name="yt",
        version=VERSION,
        description="An analysis and visualization toolkit for volumetric data",
        long_description=long_description,
        long_description_content_type="text/markdown",
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: POSIX :: AIX",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: C",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Topic :: Scientific/Engineering :: Astronomy",
            "Topic :: Scientific/Engineering :: Physics",
            "Topic :: Scientific/Engineering :: Visualization",
            "Framework :: Matplotlib",
        ],
        keywords="astronomy astrophysics visualization " + "amr adaptivemeshrefinement",
        entry_points={
            "console_scripts": ["yt = yt.utilities.command_line:run_main",],
            "nose.plugins.0.10": [
                "answer-testing = yt.utilities.answer_testing.framework:AnswerTesting"
            ],
        },
        packages=find_packages(),
        include_package_data=True,
        install_requires=[
            "matplotlib>=1.5.3",
            "setuptools>=19.6",
            "sympy>=1.2",
            "numpy>=1.10.4",
            "IPython>=1.0",
            "unyt>=2.7.2",
        ],
        extras_require={"hub": ["girder_client"], "mapserver": ["bottle"]},
        cmdclass={"sdist": sdist, "build_ext": build_ext},
        author="The yt project",
        author_email="yt-dev@python.org",
        url="https://github.com/yt-project/yt",
        project_urls={
            "Homepage": "https://yt-project.org/",
            "Documentation": "https://yt-project.org/doc/",
            "Source": "https://github.com/yt-project/yt/",
            "Tracker": "https://github.com/yt-project/yt/issues",
        },
        license="BSD 3-Clause",
        zip_safe=False,
        scripts=["scripts/iyt"],
        distclass=BinaryDistribution,
        ext_modules=[],  # !!! We override this inside build_ext above
        python_requires=">=3.6",
    )
