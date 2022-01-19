"""Module for cookbook testing


This test should be run from main yt directory.

Example:

      $ sed -e '/where/d' -i nose.cfg setup.cfg
      $ nosetests doc/source/cookbook/tests/test_cookbook.py -P -v
"""
import glob
import os
import subprocess
import sys


def run_with_capture(*args, **kwargs):
    sp = subprocess.Popen(
        *args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs
    )
    out, err = sp.communicate()
    if out:
        sys.stdout.write(out.decode("UTF-8"))
    if err:
        sys.stderr.write(err.decode("UTF-8"))

    if sp.returncode != 0:
        retstderr = " ".join(args[0])
        retstderr += "\n\nTHIS IS THE REAL CAUSE OF THE FAILURE:\n"
        retstderr += err.decode("UTF-8") + "\n"
        raise subprocess.CalledProcessError(sp.returncode, retstderr)

    return sp.returncode


PARALLEL_TEST = {"rockstar_nest.py": "3"}
BLACKLIST = [
    "opengl_ipython.py",
    "opengl_vr.py",
    "matplotlib-animation.py",
    "rockstar_nest.py",
]


def test_recipe():
    """Dummy test grabbing all cookbook's recipes"""
    for fname in glob.glob("doc/source/cookbook/*.py"):
        recipe = os.path.basename(fname)
        if recipe in BLACKLIST:
            continue
        check_recipe.description = f"Testing recipe: {recipe}"
        if recipe in PARALLEL_TEST:
            yield check_recipe, [
                "mpiexec",
                "-n",
                PARALLEL_TEST[recipe],
                "python",
                fname,
            ]
        else:
            yield check_recipe, ["python", fname]


def check_recipe(cmd):
    """Run single recipe"""
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    if out:
        sys.stdout.write(out.decode("utf8"))
    if err:
        sys.stderr.write(err.decode("utf8"))

    if proc.returncode != 0:
        retstderr = " ".join(cmd)
        retstderr += "\n\nTHIS IS THE REAL CAUSE OF THE FAILURE:\n"
        retstderr += err.decode("UTF-8") + "\n"
        raise subprocess.CalledProcessError(proc.returncode, retstderr)
