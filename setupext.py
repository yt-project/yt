import contextlib
import glob
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor as Pool
from distutils import log
from distutils.ccompiler import CCompiler, new_compiler
from distutils.errors import CompileError, LinkError
from distutils.sysconfig import customize_compiler
from distutils.version import LooseVersion
from subprocess import PIPE, Popen
from sys import platform as _platform

from pkg_resources import resource_filename
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist

CCODE = """
#include <omp.h>
#include <stdio.h>
int main() {
  omp_set_num_threads(2);
  #pragma omp parallel
  printf("nthreads=%d\\n", omp_get_num_threads());
  return 0;
}
"""


@contextlib.contextmanager
def stdchannel_redirected(stdchannel, dest_filename):
    """
    A context manager to temporarily redirect stdout or stderr

    e.g.:

    with stdchannel_redirected(sys.stderr, os.devnull):
        if compiler.has_function('clock_gettime', libraries=['rt']):
            libraries.append('rt')

    Code adapted from https://stackoverflow.com/a/17752455/1382869
    """

    try:
        oldstdchannel = os.dup(stdchannel.fileno())
        dest_file = open(dest_filename, "w")
        os.dup2(dest_file.fileno(), stdchannel.fileno())

        yield
    finally:
        if oldstdchannel is not None:
            os.dup2(oldstdchannel, stdchannel.fileno())
        if dest_file is not None:
            dest_file.close()


def check_for_openmp():
    """Returns True if local setup supports OpenMP, False otherwise

    Code adapted from astropy_helpers, originally written by Tom
    Robitaille and Curtis McCully.
    """

    # Create a temporary directory
    ccompiler = new_compiler()
    customize_compiler(ccompiler)

    tmp_dir = tempfile.mkdtemp()
    start_dir = os.path.abspath(".")

    if os.name == "nt":
        # TODO: make this work with mingw
        # AFAICS there's no easy way to get the compiler distutils
        # will be using until compilation actually happens
        compile_flag = "-openmp"
        link_flag = ""
    else:
        compile_flag = "-fopenmp"
        link_flag = "-fopenmp"

    try:
        os.chdir(tmp_dir)

        with open("test_openmp.c", "w") as f:
            f.write(CCODE)

        os.mkdir("objects")

        # Compile, link, and run test program
        with stdchannel_redirected(sys.stderr, os.devnull):
            ccompiler.compile(
                ["test_openmp.c"], output_dir="objects", extra_postargs=[compile_flag]
            )
            ccompiler.link_executable(
                glob.glob(os.path.join("objects", "*")),
                "test_openmp",
                extra_postargs=[link_flag],
            )
            output = (
                subprocess.check_output("./test_openmp")
                .decode(sys.stdout.encoding or "utf-8")
                .splitlines()
            )

        if "nthreads=" in output[0]:
            nthreads = int(output[0].strip().split("=")[1])
            if len(output) == nthreads:
                using_openmp = True
            else:
                log.warn(
                    "Unexpected number of lines from output of test "
                    "OpenMP program (output was %s)",
                    output,
                )
                using_openmp = False
        else:
            log.warn(
                "Unexpected output from test OpenMP program (output was %s)", output
            )
            using_openmp = False

    except (CompileError, LinkError):
        using_openmp = False
    finally:
        os.chdir(start_dir)

    if using_openmp:
        log.warn("Using OpenMP to compile parallel extensions")
    else:
        log.warn(
            "Unable to compile OpenMP test program so Cython\n"
            "extensions will be compiled without parallel support"
        )

    return using_openmp


def check_for_pyembree(std_libs):
    embree_libs = []
    embree_aliases = {}
    try:
        _ = resource_filename("pyembree", "rtcore.pxd")
    except ImportError:
        return embree_libs, embree_aliases

    embree_prefix = os.path.abspath(read_embree_location())
    embree_inc_dir = os.path.join(embree_prefix, "include")
    embree_lib_dir = os.path.join(embree_prefix, "lib")

    if _platform == "darwin":
        embree_lib_name = "embree.2"
    else:
        embree_lib_name = "embree"

    embree_aliases["EMBREE_INC_DIR"] = ["yt/utilities/lib/", embree_inc_dir]
    embree_aliases["EMBREE_LIB_DIR"] = [embree_lib_dir]
    embree_aliases["EMBREE_LIBS"] = std_libs + [embree_lib_name]
    embree_libs += ["yt/utilities/lib/embree_mesh/*.pyx"]

    if in_conda_env():
        conda_basedir = os.path.dirname(os.path.dirname(sys.executable))
        embree_aliases["EMBREE_INC_DIR"].append(os.path.join(conda_basedir, "include"))
        embree_aliases["EMBREE_LIB_DIR"].append(os.path.join(conda_basedir, "lib"))

    return embree_libs, embree_aliases


def in_conda_env():
    return any(s in sys.version for s in ("Anaconda", "Continuum", "conda-forge"))


def read_embree_location():
    """

    Attempts to locate the embree installation. First, we check for an
    EMBREE_DIR environment variable. If one is not defined, we look for
    an embree.cfg file in the root yt source directory. Finally, if that
    is not present, we default to /usr/local. If embree is installed in a
    non-standard location and none of the above are set, the compile will
    not succeed. This only gets called if check_for_pyembree() returns
    something other than None.

    """

    rd = os.environ.get("EMBREE_DIR")
    if rd is None:
        try:
            rd = open("embree.cfg").read().strip()
        except IOError:
            rd = "/usr/local"

    fail_msg = (
        "I attempted to find Embree headers in %s. \n"
        "If this is not correct, please set your correct embree location \n"
        "using EMBREE_DIR environment variable or your embree.cfg file. \n"
        "Please see http://yt-project.org/docs/dev/visualizing/unstructured_mesh_rendering.html "
        "for more information. \n" % rd
    )

    # Create a temporary directory
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()

    try:
        os.chdir(tmpdir)

        # Get compiler invocation
        compiler = os.getenv("CXX", "c++")
        compiler = compiler.split(" ")

        # Attempt to compile a test script.
        filename = r"test.cpp"
        file = open(filename, "wt", 1)
        file.write('#include "embree2/rtcore.h"\n' "int main() {\n" "return 0;\n" "}")
        file.flush()
        p = Popen(
            compiler + ["-I%s/include/" % rd, filename],
            stdin=PIPE,
            stdout=PIPE,
            stderr=PIPE,
        )
        output, err = p.communicate()
        exit_code = p.returncode

        if exit_code != 0:
            log.warn(
                "Pyembree is installed, but I could not compile Embree " "test code."
            )
            log.warn("The error message was: ")
            log.warn(err)
            log.warn(fail_msg)

        # Clean up
        file.close()

    except OSError:
        log.warn(
            "read_embree_location() could not find your C compiler. "
            "Attempted to use '%s'.",
            compiler,
        )
        return False

    finally:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)

    return rd


def get_cpu_count():
    if platform.system() == "Windows":
        return 0

    cpu_count = os.cpu_count()
    try:
        user_max_cores = int(os.getenv("MAX_BUILD_CORES", cpu_count))
    except ValueError as e:
        raise ValueError(
            "MAX_BUILD_CORES must be set to an integer. "
            + "See above for original error."
        ) from e
    max_cores = min(cpu_count, user_max_cores)
    return max_cores


def install_ccompiler():
    def _compile(
        self,
        sources,
        output_dir=None,
        macros=None,
        include_dirs=None,
        debug=0,
        extra_preargs=None,
        extra_postargs=None,
        depends=None,
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


def create_build_ext(lib_exts, cythonize_aliases):
    class build_ext(_build_ext):
        # subclass setuptools extension builder to avoid importing cython and numpy
        # at top level in setup.py. See http://stackoverflow.com/a/21621689/1382869
        def finalize_options(self):
            try:
                import cython
                import numpy
            except ImportError as e:
                raise ImportError(
                    """Could not import cython or numpy. Building yt from source requires
    cython and numpy to be installed. Please install these packages using
    the appropriate package manager for your python environment."""
                ) from e
            if LooseVersion(cython.__version__) < LooseVersion("0.26.1"):
                raise RuntimeError(
                    """Building yt from source requires Cython 0.26.1 or newer but
    Cython %s is installed. Please update Cython using the appropriate
    package manager for your python environment."""
                    % cython.__version__
                )
            if LooseVersion(numpy.__version__) < LooseVersion("1.13.3"):
                raise RuntimeError(
                    """Building yt from source requires NumPy 1.13.3 or newer but
    NumPy %s is installed. Please update NumPy using the appropriate
    package manager for your python environment."""
                    % numpy.__version__
                )
            from Cython.Build import cythonize

            # Override the list of extension modules
            self.distribution.ext_modules[:] = cythonize(
                lib_exts,
                aliases=cythonize_aliases,
                compiler_directives={"language_level": 2},
                nthreads=get_cpu_count(),
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

            ncpus = get_cpu_count()
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
                lib_exts,
                aliases=cythonize_aliases,
                compiler_directives={"language_level": 2},
                nthreads=get_cpu_count(),
            )
            _sdist.run(self)

    return build_ext, sdist
