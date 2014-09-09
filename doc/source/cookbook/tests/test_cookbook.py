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

sys.path.append(os.path.join(os.getcwd(), "doc/source/cookbook"))


def test_recipe():
    '''Dummy test grabbing all cookbook's recipes'''
    for fname in glob.glob("doc/source/cookbook/*.py"):
        module_name = os.path.splitext(os.path.basename(fname))[0]
        yield check_recipe, module_name


def check_recipe(module_name):
    '''Run single recipe'''
    __import__(module_name)
    assert True
