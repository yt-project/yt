import contextlib
import glob
import logging
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from textwrap import dedent
from concurrent.futures import ThreadPoolExecutor
from distutils.ccompiler import CCompiler, new_compiler
from distutils.sysconfig import customize_compiler
from subprocess import PIPE, Popen
from sys import platform as _platform

from pkg_resources import resource_filename
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist
from setuptools.errors import CompileError, LinkError

log = logging.getLogger("setupext")

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
    """Returns OpenMP compiler and linker flags if local setup supports
    OpenMP or [], [] otherwise

    Code adapted from astropy_helpers, originally written by Tom
    Robitaille and Curtis McCully.
    """

    # Create a temporary directory
    ccompiler = new_compiler()
    customize_compiler(ccompiler)

    tmp_dir = tempfile.mkdtemp()
    start_dir = os.path.abspath(".")

    CCODE = dedent("""\
        #include <omp.h>
        #include <stdio.h>
        int main() {
            omp_set_num_threads(2);
            #pragma omp parallel
            printf("nthreads=%d\\n", omp_get_num_threads());
            return 0;
        }"""
    )

    # TODO: test more known compilers:
    # MinGW, AppleClang with libomp, MSVC, ICC, XL, PGI, ...
    if os.name == "nt":
        # TODO: make this work with mingw
        # AFAICS there's no easy way to get the compiler distutils
        # will be using until compilation actually happens
        compile_flags = ["-openmp"]
        link_flags = [""]
    else:
        compile_flags = ["-fopenmp"]
        link_flags = ["-fopenmp"]

    try:
        os.chdir(tmp_dir)

        with open("test_openmp.c", "w") as f:
            f.write(CCODE)

        os.mkdir("objects")

        # Compile, link, and run test program
        with stdchannel_redirected(sys.stderr, os.devnull):
            ccompiler.compile(
                ["test_openmp.c"], output_dir="objects", extra_postargs=compile_flags
            )
            ccompiler.link_executable(
                glob.glob(os.path.join("objects", "*")),
                "test_openmp",
                extra_postargs=link_flags,
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
                log.warning(
                    "Unexpected number of lines from output of test "
                    "OpenMP program (output was %s)",
                    output,
                )
                using_openmp = False
        else:
            log.warning(
                "Unexpected output from test OpenMP program (output was %s)", output
            )
            using_openmp = False

    except (CompileError, LinkError):
        using_openmp = False
    finally:
        os.chdir(start_dir)

    if using_openmp:
        log.warning("Using OpenMP to compile parallel extensions")
    else:
        log.warning(
            "Unable to compile OpenMP test program so Cython\n"
            "extensions will be compiled without parallel support"
        )

    if using_openmp:
        return compile_flags, link_flags
    else:
        return [], []


def check_CPP14_flag(compile_flags):
    # Create a temporary directory
    ccompiler = new_compiler()
    customize_compiler(ccompiler)

    tmp_dir = tempfile.mkdtemp()
    start_dir = os.path.abspath(".")

    # Note: This code requires C++14 functionalities (also required to compile yt)
    # It compiles on gcc 4.7.4 (together with the entirety of yt) with the flag "-std=gnu++0x".
    # It does not compile on gcc 4.6.4 (neither does yt).
    CPPCODE = dedent("""\
        #include <vector>

        struct node {
            std::vector<int> vic;
            bool visited = false;
        };

        int main() {
            return 0;
        }"""
    )

    os.chdir(tmp_dir)
    try:
        with open("test_cpp14.cpp", "w") as f:
            f.write(CPPCODE)

        os.mkdir("objects")

        # Compile, link, and run test program
        with stdchannel_redirected(sys.stderr, os.devnull):
            ccompiler.compile(
                ["test_cpp14.cpp"], output_dir="objects", extra_postargs=compile_flags
            )
        return True
    except CompileError:
        return False
    finally:
        os.chdir(start_dir)


def check_CPP14_flags(possible_compile_flags):
    for flags in possible_compile_flags:
        if check_CPP14_flag([flags]):
            return flags

    log.warning(
        "Your compiler seems to be too old to support C++14. "
        "yt may not be able to compile. Please use a newer version."
    )
    return []


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
        CCODE = dedent("""\
            #include "embree2/rtcore.h
            int main() {
                return 0;
            }"""
        )
        file.write(CCODE)
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
            log.warning(
                "Pyembree is installed, but I could not compile Embree test code."
            )
            log.warning("The error message was: ")
            log.warning(err)
            log.warning(fail_msg)

        # Clean up
        file.close()

    except OSError:
        log.warning(
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
        # NOTE: this is likely not necessary anymore since
        # pyproject.toml was introduced in the project

        def finalize_options(self):
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
                with ThreadPoolExecutor(ncpus) as executor:
                    results = {
                        executor.submit(self.build_extension, extension): extension
                        for extension in self.extensions
                    }
                for result in results:
                    result.result()
            else:
                super().build_extensions()

        def build_extension(self, extension):
            try:
                super().build_extension(extension)
            except CompileError as exc:
                print(f"While building '{extension.name}' following error was raised:\n {exc}")
                raise

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
