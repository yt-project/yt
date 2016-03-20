import os
from pkg_resources import resource_filename
import shutil
import subprocess
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
        with open(os.devnull, 'w') as fnull:
            exit_code = subprocess.call(compiler + ['-fopenmp', filename],
                                        stdout=fnull, stderr=fnull)

        # Clean up
        file.close()
    except OSError:
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

    fail_msg = "Pyembree is installed, but I could not compile Embree test code. \n" + \
               "Attempted to find Embree headers in %s. \n" % rd + \
               "If this is not correct, please set your correct embree location \n" + \
               "using EMBREE_DIR environment variable or your embree.cfg file. \n" + \
               "Please see http://yt-project.org/docs/dev/visualizing/unstructured_mesh_rendering.html " + \
               "for more information."

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
        with open(os.devnull, 'w') as fnull:
            exit_code = subprocess.call(compiler + ['-I%s/include/' % rd, filename],
                             stdout=fnull, stderr=fnull)

        # Clean up
        file.close()

    except OSError:
        print fail_msg

    finally:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)

    if exit_code != 0:
        print fail_msg

    return rd


def get_mercurial_changeset_id(target_dir):
    '''
    Returns changeset and branch using hglib
    '''
    try:
        import hglib
    except ImportError:
        return None
    with hglib.open(target_dir) as repo:
        changeset = repo.identify(id=True, branch=True).strip().decode('utf8')
    return changeset
