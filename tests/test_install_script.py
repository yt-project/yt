import contextlib
import glob
import os
import shutil
import subprocess
import sys
import tempfile

# dependencies that are always installed
REQUIRED_DEPS = [
    'mercurial',
    'jupyter',
    'numpy',
    'matplotlib',
    'h5py',
    'cython',
    'nose',
    'sympy',
    'setuptools',
    'hglib'
]

# dependencies that aren't installed by default
OPTIONAL_DEPS = [
    'unstructured',
    'pyx',
    'rockstar',
    'scipy',
    'astropy',
]

# dependencies that are only installable when yt is built from source
YT_SOURCE_ONLY_DEPS = [
    'unstructured',
    'rockstar'
]

# dependencies that are only installable when yt is built from source under conda
YT_SOURCE_CONDA_ONLY_DEPS = [
    'unstructured'
]

DEPENDENCY_IMPORT_TESTS = {
    'unstructured': "from yt.utilities.lib import mesh_traversal",
    'rockstar': ("from yt.analysis_modules.halo_finding.rockstar "
                 "import rockstar_interface")
}


def call_unix_command(command):
    print ('Running "%s" in %s' % (command, os.getcwd()))
    output = ''
    try:
        output = subprocess.check_output(command, stderr=subprocess.STDOUT,
                                         shell=True)
    except subprocess.CalledProcessError as er:
        raise RuntimeError(
            'Command \'%s\' failed with return code \'%s\' and error:\n\n%s' %
            (command, er.returncode, er.output))
    finally:
        if len(output.splitlines()) > 25:
            print ('truncated output:')
            print ('\n'.join((output.splitlines())[-25:]))
        else:
            print (output)


@contextlib.contextmanager
def working_directory(path):
    """A context manager which changes the working directory to the given
    path, and then changes it back to its previous value on exit.

    """
    prev_cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)
        shutil.rmtree(path)


def run_install_script(install_script_path, inst_py3,
                       conda=False, binary_yt=False):
    msg = 'Testing installation with conda={}, inst_py3={}, and binary_yt={}'
    print (msg.format(conda, inst_py3, binary_yt))
    shutil.copy(install_script_path, os.curdir)
    with open('install_script.sh', 'r') as source:
        with open('install_script_edited.sh', 'w') as target:
            data = source.read()
            for dep in OPTIONAL_DEPS:
                if binary_yt is True and dep in YT_SOURCE_ONLY_DEPS:
                    continue
                if conda is False and dep in YT_SOURCE_CONDA_ONLY_DEPS:
                    continue
                dname = 'INST_%s' % dep.upper()
                data = data.replace(dname + '=0', dname + '=1')
            if inst_py3 is True:
                data = data.replace('INST_PY3=0', 'INST_PY3=1')
            if conda is False:
                data = data.replace('INST_CONDA=1', 'INST_CONDA=0')
            if binary_yt is False:
                data = data.replace('INST_YT_SOURCE=0', 'INST_YT_SOURCE=1')
            target.write(data)
    shutil.copyfile('install_script_edited.sh', 'install_script.sh')
    call_unix_command('bash install_script.sh --yes')


def verify_yt_installation(binary_yt, conda):
    yt_dir = glob.glob('yt-*')
    ndirs = len(yt_dir)
    if ndirs != 1:
        raise RuntimeError(
            'A yt installation was not properly cleaned up, exiting')
    yt_dir = yt_dir[0]
    python_path = os.sep.join([yt_dir, 'bin', 'python'])
    for dep in OPTIONAL_DEPS + REQUIRED_DEPS:
        if binary_yt is True and dep in YT_SOURCE_ONLY_DEPS:
            continue
        if conda is False and dep in YT_SOURCE_CONDA_ONLY_DEPS:
            continue
        elif dep == 'mercurial':
            hg_path = os.sep.join([yt_dir, 'bin', 'hg'])
            call_unix_command('{} --version'.format(hg_path))
        elif dep in DEPENDENCY_IMPORT_TESTS:
            cmd = "{} -c '{}'"
            if dep == 'rockstar':
                cmd = 'LD_LIBRARY_PATH={} '.format(
                    os.sep.join([os.curdir, yt_dir, 'lib'])) + cmd
                if sys.platform == 'darwin':
                    cmd = 'DY' + cmd
            call_unix_command(cmd.format(
                python_path, DEPENDENCY_IMPORT_TESTS[dep]))
        else:
            call_unix_command("{} -c 'import {}'".format(python_path, dep))
    return yt_dir


if __name__ == '__main__':
    install_script_path = os.path.abspath(os.path.sep.join(
        [os.getcwd(), os.pardir, 'doc', 'install_script.sh']))
    for inst_py3 in [False, True]:
        tmpdir = tempfile.mkdtemp()
        with working_directory(tmpdir):
            run_install_script(
                install_script_path, inst_py3, conda=True, binary_yt=True)
            conda_binary_path = verify_yt_installation(
                binary_yt=True, conda=True)
            shutil.rmtree(conda_binary_path)

            run_install_script(
                install_script_path, inst_py3, conda=True, binary_yt=False)
            conda_source_path = verify_yt_installation(
                binary_yt=False, conda=True)
            shutil.rmtree(conda_source_path)

            run_install_script(
                install_script_path, inst_py3, conda=False, binary_yt=False)
            source_path = verify_yt_installation(
                binary_yt=False, conda=False)
            shutil.rmtree(source_path)
