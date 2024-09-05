"""Module for cookbook testing


This test should be run from main yt directory.

Example:

      $ sed -e '/where/d' -i nose.cfg setup.cfg
      $ nosetests doc/source/cookbook/tests/test_cookbook.py -P -v
"""

import subprocess
import sys

from pathlib import Path


BLACKLIST = [
    "matplotlib-animation.py",
]


def test_recipe():
    """Dummy test grabbing all cookbook's recipes"""
    COOKBOOK_DIR = Path("doc", "source", "cookbook")
    for fname in sorted(COOKBOOK_DIR.glob("*.py")):
        if fname.name in BLACKLIST:
            continue
        check_recipe.description = f"Testing recipe: {fname.name}"
        yield check_recipe, ["python", str(fname)]


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
