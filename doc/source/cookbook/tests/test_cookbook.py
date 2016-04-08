# -*- coding: utf-8 -*-
"""Module for cookbook testing


This test should be run from main yt directory.

Example:

      $ sed -e '/where/d' -i nose.cfg setup.cfg
      $ nosetests doc/source/cookbook/tests/test_cookbook.py -P -v
"""
import glob
import os
import sys
import subprocess


PARALLEL_TEST = {"rockstar_nest.py": "3"}
BLACKLIST = ["opengl_ipython.py", "opengl_vr.py"]

if sys.version_info >= (3,0,0):
    BLACKLIST.append("rockstar_nest.py")

def test_recipe():
    '''Dummy test grabbing all cookbook's recipes'''
    for fname in glob.glob("doc/source/cookbook/*.py"):
        recipe = os.path.basename(fname)
        if recipe in BLACKLIST:
            continue
        check_recipe.description = "Testing recipe: %s" % recipe
        if recipe in PARALLEL_TEST:
            yield check_recipe, \
                ["mpiexec", "-n", PARALLEL_TEST[recipe], "python", fname]
        else:
            yield check_recipe, ["python", fname]


def check_recipe(cmd):
    '''Run single recipe'''
    try:
        subprocess.check_call(cmd)
        result = True
    except subprocess.CalledProcessError as e:
        print(("Stdout output:\n", e.output))
        result = False
    assert result
