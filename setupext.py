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
    if rd is not None:
        return rd
    print("EMBREE_DIR not set. Attempting to read embree.cfg")
    try:
        rd = open("embree.cfg").read().strip()
        return rd
    except IOError:
        print("Reading Embree location from embree.cfg failed.")
        print("If compilation fails, please place the base directory")
        print("of your Embree install in embree.cfg and restart.")
        return '/usr/local'


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
