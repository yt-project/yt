# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 18:24:09 2015

@author: stuart
"""

from yt.utilities.hierarchy_inspection import find_lowest_subclass

from nose.tools import assert_raises

class level1(object):
    pass

class level1a(object):
    pass

class level2(level1):
    pass

class level3(level2):
    pass

class level4(level3):
    pass




def test_two_classes():
    result = find_lowest_subclass([level1, level2])
    assert result is level2

def test_four_deep():
    result = find_lowest_subclass([level1, level2, level3, level4])
    assert result is level4

def test_four_deep_outoforder():
    result = find_lowest_subclass([level2, level3, level1, level4])
    assert result is level4

def test_diverging_tree():
    with assert_raises(TypeError):
        find_lowest_subclass([level1, level2, level3, level1a])
