import contextlib
import glob
import os
import platform
import shutil
import subprocess
import sys
import tempfile

# dependencies that are always installed
REQUIRED_DEPS = [
    "git",
    "jupyter",
    "numpy",
    "matplotlib",
    "h5py",
    "cython",
    "nose",
    "sympy",
    "setuptools",
]

# dependencies that aren't installed by default
OPTIONAL_DEPS = [
    "embree",
    "pyx",
    "scipy",
    "astropy",
    "cartopy",
    "pooch",
]

# dependencies that are only installable when yt is built from source
YT_SOURCE_ONLY_DEPS = [
    "embree",
]

DEPENDENCY_IMPORT_TESTS = {
    "embree": "from yt.utilities.lib.embree_mesh import mesh_traversal",
}


def call_unix_command(command):
    print(f'Running "{command}" in {os.getcwd()}')
    output = ""
    try:
        output = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError as er:
        raise RuntimeError(
            "Command '%s' failed with return code '%s' and error:\n\n%s"
            % (command, er.returncode, er.output.decode("utf-8"))
        )
    finally:
        if len(output.splitlines()) > 25:
            print("truncated output:")
            print("\n".join((output.decode().splitlines())[-25:]))
        else:
            print(output)


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


def run_install_script(install_script_path, inst_py3, binary_yt=False):
    msg = "Testing installation with inst_py3={} and binary_yt={}"
    print(msg.format(inst_py3, binary_yt))
    shutil.copy(install_script_path, os.curdir)
    with open("install_script.sh") as source:
        with open("install_script_edited.sh", "w") as target:
            data = source.read()
            for dep in OPTIONAL_DEPS:
                if binary_yt and dep in YT_SOURCE_ONLY_DEPS:
                    continue
                if dep == "rockstar":
                    # compiling rockstar is broken on newer MacOS releases
                    if platform.mac_ver()[0].startswith(("10.12", "10.13")):
                        continue
                dname = f"INST_{dep.upper()}"
                data = data.replace(dname + "=0", dname + "=1")
            if inst_py3:
                data = data.replace("INST_PY3=0", "INST_PY3=1")
            if not binary_yt:
                data = data.replace("INST_YT_SOURCE=0", "INST_YT_SOURCE=1")
            target.write(data)
    shutil.copyfile("install_script_edited.sh", "install_script.sh")
    call_unix_command("bash install_script.sh --yes")


def verify_yt_installation(binary_yt):
    yt_dir = glob.glob("yt-*")
    ndirs = len(yt_dir)
    if ndirs != 1:
        raise RuntimeError("A yt installation was not properly cleaned up, exiting")
    yt_dir = yt_dir[0]
    python_path = os.sep.join([yt_dir, "bin", "python"])
    for dep in OPTIONAL_DEPS + REQUIRED_DEPS:
        if binary_yt and dep in YT_SOURCE_ONLY_DEPS:
            continue
        if dep == "git":
            git_path = os.sep.join([yt_dir, "bin", "git"])
            call_unix_command(f"{git_path} --version")
        if dep in DEPENDENCY_IMPORT_TESTS:
            cmd = "{} -c '{}'"
            if dep == "rockstar":
                # compiling rockstar is broken on newer MacOS releases
                if platform.mac_ver()[0].startswith(("10.12", "10.13")):
                    continue
                cmd = (
                    f"LD_LIBRARY_PATH={os.sep.join([os.curdir, yt_dir, 'lib'])} " + cmd
                )
                if sys.platform == "darwin":
                    cmd = "DY" + cmd
            call_unix_command(cmd.format(python_path, DEPENDENCY_IMPORT_TESTS[dep]))
        else:
            call_unix_command(f"{python_path} -c 'import {dep}'")
    return yt_dir


if __name__ == "__main__":
    install_script_path = os.path.abspath(
        os.path.sep.join([os.getcwd(), os.pardir, "doc", "install_script.sh"])
    )
    for inst_py3 in [False, True]:
        tmpdir = tempfile.mkdtemp()
        with working_directory(tmpdir):
            run_install_script(install_script_path, inst_py3, binary_yt=True)
            conda_binary_path = verify_yt_installation(binary_yt=True)
            shutil.rmtree(conda_binary_path)

            run_install_script(install_script_path, inst_py3, binary_yt=False)
            conda_source_path = verify_yt_installation(binary_yt=False)
            shutil.rmtree(conda_source_path)
