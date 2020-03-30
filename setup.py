import os
import platform
from concurrent.futures import ThreadPoolExecutor as Pool
import glob
import sys
from sys import platform as _platform
from setuptools import setup, find_packages
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist
from setupext import \
    check_for_openmp, \
    check_for_pyembree, \
    read_embree_location, \
    in_conda_env
from distutils.version import LooseVersion
from distutils.ccompiler import CCompiler
import pkg_resources


def _get_cpu_count():
    if platform.system() != "Windows":
        return os.cpu_count()
    return 0


def _compile(
    self, sources, output_dir=None, macros=None, include_dirs=None,
    debug=0, extra_preargs=None, extra_postargs=None, depends=None,
):
    """Function to monkey-patch distutils.ccompiler.CCompiler"""
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(
        output_dir, macros, include_dirs, sources, depends, extra_postargs
    )
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)

    for obj in objects:
        try:
            src, ext = build[obj]
        except KeyError:
            continue
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

    # Return *all* object filenames, not just the ones we just built.
    return objects

CCompiler.compile = _compile

if sys.version_info < (3, 5):
    print("yt currently supports versions newer than Python 3.5")
    print("certain features may fail unexpectedly and silently with older "
          "versions.")
    sys.exit(1)

try:
    distribute_ver = \
        LooseVersion(pkg_resources.get_distribution("distribute").version)
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

if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

with open('README.md') as file:
    long_description = file.read()

if check_for_openmp() is True:
    omp_args = ['-fopenmp']
else:
    omp_args = None

if os.name == "nt":
    std_libs = []
else:
    std_libs = ["m"]

from Cython.Build import cythonize
cythonize_aliases = {'LIB_DIR': 'yt/utilities/lib/',
                     'LIB_DIR_EWAH': ['yt/utilities/lib/', 'yt/utilities/lib/ewahboolarray/'],
                     'LIB_DIR_GEOM': ['yt/utilities/lib/', 'yt/geometry/'],
                     'LIB_DIR_GEOM_ARTIO': ['yt/utilities/lib/', 'yt/geometry/',
                                            'yt/frontends/artio/artio_headers/'],
                     'STD_LIBS': std_libs,
                     'OMP_ARGS': omp_args,
                     'FIXED_INTERP': 'yt/utilities/lib/fixed_interpolator.c',
                     'ARTIO_SOURCE': glob.glob("yt/frontends/artio/artio_headers/*.c"),
}

cython_extensions = []

lib_exts = ["yt/geometry/*.pyx",
            "yt/utilities/cython_fortran_utils.pyx",
            "yt/frontends/ramses/io_utils.pyx",
            "yt/utilities/lib/cykdtree/kdtree.pyx",
            "yt/utilities/lib/cykdtree/utils.pyx",
            "yt/frontends/artio/_artio_caller.pyx",
            "yt/utilities/lib/*.pyx",
]

# EMBREE
if check_for_pyembree() is not None:

    embree_prefix = os.path.abspath(read_embree_location())
    embree_inc_dir = os.path.join(embree_prefix, 'include')
    embree_lib_dir = os.path.join(embree_prefix, 'lib')
    if in_conda_env():
        conda_basedir = os.path.dirname(os.path.dirname(sys.executable))
        embree_inc_dir.append(os.path.join(conda_basedir, 'include'))
        embree_lib_dir.append(os.path.join(conda_basedir, 'lib'))

    if _platform == "darwin":
        embree_lib_name = "embree.2"
    else:
        embree_lib_name = "embree"

    cythonize_aliases['EMBREE_INC_DIR'] = ['yt/utilities/lib/', embree_inc_dir]
    cythonize_aliases['EMBREE_LIB_DIR'] = [embree_lib_dir]
    cythonize_aliases['EMBREE_LIBS'] = std_libs + [embree_lib_name]
    lib_exts += ['yt/utilities/lib/embree_mesh/*.pyx']

cython_extensions += cythonize(lib_exts, aliases = cythonize_aliases,
                               language_level = 2)


if os.environ.get("GPERFTOOLS", "no").upper() != "NO":
    gpd = os.environ["GPERFTOOLS"]
    idir = os.path.join(gpd, "include")
    ldir = os.path.join(gpd, "lib")
    print(("INCLUDE AND LIB DIRS", idir, ldir))
    cython_extensions.append(
        Extension("yt.utilities.lib.perftools_wrap",
                  ["yt/utilities/lib/perftools_wrap.pyx"],
                  libraries=["profiler"],
                  library_dirs=[ldir],
                  include_dirs=[idir]))

class build_ext(_build_ext):
    # subclass setuptools extension builder to avoid importing cython and numpy
    # at top level in setup.py. See http://stackoverflow.com/a/21621689/1382869
    def finalize_options(self):
        try:
            import cython
            import numpy
        except ImportError:
            raise ImportError(
"""Could not import cython or numpy. Building yt from source requires
cython and numpy to be installed. Please install these packages using
the appropriate package manager for your python environment.""")
        if LooseVersion(cython.__version__) < LooseVersion('0.26.1'):
            raise RuntimeError(
"""Building yt from source requires Cython 0.26.1 or newer but
Cython %s is installed. Please update Cython using the appropriate
package manager for your python environment.""" %
                cython.__version__)
        if LooseVersion(numpy.__version__) < LooseVersion('1.13.3'):
            raise RuntimeError(
"""Building yt from source requires NumPy 1.13.3 or newer but
NumPy %s is installed. Please update NumPy using the appropriate
package manager for your python environment.""" %
                numpy.__version__)
        from Cython.Build import cythonize
        self.distribution.ext_modules[:] = cythonize(
            self.distribution.ext_modules,
            compiler_directives={'language_level': 2},
            nthreads=_get_cpu_count(),
        )
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process
        # see http://stackoverflow.com/a/21621493/1382869
        if isinstance(__builtins__, dict):
            # sometimes this is a dict so we need to check for that
            # https://docs.python.org/3/library/builtins.html
            __builtins__["__NUMPY_SETUP__"] = False
        else:
            __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

    def build_extensions(self):
        self.check_extensions_list(self.extensions)

        ncpus = _get_cpu_count()
        if ncpus > 0:
            with Pool(ncpus) as pool:
                pool.map(self.build_extension, self.extensions)
        else:
            super().build_extensions()


class sdist(_sdist):
    # subclass setuptools source distribution builder to ensure cython
    # generated C files are included in source distribution.
    # See http://stackoverflow.com/a/18418524/1382869
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        from Cython.Build import cythonize
        cythonize(
            cython_extensions,
            compiler_directives={'language_level': 2},
            nthreads=_get_cpu_count(),
        )
        _sdist.run(self)

setup(
    name="yt",
    version=VERSION,
    description="An analysis and visualization toolkit for volumetric data",
    long_description = long_description,
    long_description_content_type='text/markdown',
    classifiers=["Development Status :: 5 - Production/Stable",
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
                 "Topic :: Scientific/Engineering :: Astronomy",
                 "Topic :: Scientific/Engineering :: Physics",
                 "Topic :: Scientific/Engineering :: Visualization"],
    keywords='astronomy astrophysics visualization ' +
    'amr adaptivemeshrefinement',
    entry_points={'console_scripts': [
        'yt = yt.utilities.command_line:run_main',
    ],
        'nose.plugins.0.10': [
            'answer-testing = yt.utilities.answer_testing.framework:AnswerTesting'
    ]
    },
    packages=find_packages(),
    include_package_data = True,
    install_requires=[
        'matplotlib>=1.5.3',
        'setuptools>=19.6',
        'sympy>=1.2',
        'numpy>=1.10.4',
        'IPython>=1.0',
        'unyt>=2.2.2',
    ],
    extras_require = {
        'hub':  ["girder_client"],
        'mapserver': ["bottle"]
    },
    cmdclass={'sdist': sdist, 'build_ext': build_ext},
    author="The yt project",
    author_email="yt-dev@python.org",
    url="https://github.com/yt-project/yt",
    project_urls={
        'Homepage': 'https://yt-project.org/',
        'Documentation': 'https://yt-project.org/doc/',
        'Source': 'https://github.com/yt-project/yt/',
        'Tracker': 'https://github.com/yt-project/yt/issues'
    },
    license="BSD 3-Clause",
    zip_safe=False,
    scripts=["scripts/iyt"],
    ext_modules=cython_extensions,
    python_requires='>=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*,!=3.4.*'
)
