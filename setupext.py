import os
from pkg_resources import resource_filename
import shutil
from subprocess import Popen, PIPE
import sys
import tempfile


def check_for_openmp():
    """Returns True if local setup supports OpenMP, False otherwise"""

    # See https://bugs.python.org/issue25150
    if sys.version_info[:3] == (3, 5, 0) or os.name == 'nt':
        return False

    # Create a temporary directory
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    exit_code = 1

    try:
        os.chdir(tmpdir)

        # Get compiler invocation
        compiler = os.getenv('CC', 'cc')
        compiler = compiler.split(' ')

        # Attempt to compile a test script.
        # See http://openmp.org/wp/openmp-compilers/
        filename = r'test.c'
        file = open(filename, 'wt', 1)
        file.write(
            "#include <omp.h>\n"
            "#include <stdio.h>\n"
            "int main() {\n"
            "#pragma omp parallel\n"
            "printf(\"Hello from thread %d, nthreads %d\\n\", omp_get_thread_num(), omp_get_num_threads());\n"
            "}"
        )
        file.flush()
        p = Popen(compiler + ['-fopenmp', filename],
                  stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate()
        exit_code = p.returncode
        
        if exit_code != 0:
            print("Compilation of OpenMP test code failed with the error: ")
            print(err.decode('utf8'))
            print("Disabling OpenMP support. ")

        # Clean up
        file.close()
    except OSError:
        print("check_for_openmp() could not find your C compiler. "
              "Attempted to use '%s'. " % compiler)
        return False
    finally:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)

    return exit_code == 0


def check_for_pyembree():
    try:
        fn = resource_filename("pyembree", "rtcore.pxd")
    except ImportError:
        return None
    return os.path.dirname(fn)

def in_conda_env():
    return any(s in sys.version for s in ("Anaconda", "Continuum"))

def read_embree_location():
    '''

    Attempts to locate the embree installation. First, we check for an
    EMBREE_DIR environment variable. If one is not defined, we look for
    an embree.cfg file in the root yt source directory. Finally, if that
    is not present, we default to /usr/local. If embree is installed in a
    non-standard location and none of the above are set, the compile will
    not succeed. This only gets called if check_for_pyembree() returns
    something other than None.

    '''

    rd = os.environ.get('EMBREE_DIR')
    if rd is None:
        try:
            rd = open("embree.cfg").read().strip()
        except IOError:
            rd = '/usr/local'

    fail_msg = ("I attempted to find Embree headers in %s. \n"
               "If this is not correct, please set your correct embree location \n"
               "using EMBREE_DIR environment variable or your embree.cfg file. \n"
               "Please see http://yt-project.org/docs/dev/visualizing/unstructured_mesh_rendering.html "
                "for more information. \n" % rd)

    # Create a temporary directory
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()

    try:
        os.chdir(tmpdir)

        # Get compiler invocation
        compiler = os.getenv('CXX', 'c++')
        compiler = compiler.split(' ')

        # Attempt to compile a test script.
        filename = r'test.cpp'
        file = open(filename, 'wt', 1)
        file.write(
            '#include "embree2/rtcore.h"\n'
            'int main() {\n'
            'return 0;\n'
            '}'
        )
        file.flush()
        p = Popen(compiler + ['-I%s/include/' % rd, filename], 
                  stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate()
        exit_code = p.returncode

        if exit_code != 0:
            print("Pyembree is installed, but I could not compile Embree test code.")
            print("The error message was: ")
            print(err)
            print(fail_msg)

        # Clean up
        file.close()

    except OSError:
        print("read_embree_location() could not find your C compiler. "
              "Attempted to use '%s'. " % compiler)
        return False

    finally:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)

    return rd


def get_mercurial_changeset_id(target_dir):
    '''
    Returns changeset and branch using hglib
    '''
    try:
        import hglib
    except ImportError:
        return None
    try:
        with hglib.open(target_dir) as repo:
            changeset = repo.identify(
                id=True, branch=True).strip().decode('utf8')
    except hglib.error.ServerError:
        return None
    return changeset
