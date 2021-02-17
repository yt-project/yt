"""
Created on Wed Feb 18 18:24:09 2015

@author: stuart
"""

from ..hierarchy_inspection import find_lowest_subclasses


class level1:
    pass


class level1a:
    pass


class level2(level1):
    pass


class level3(level2):
    pass


class level4(level3):
    pass


def test_single():
    result = find_lowest_subclasses([level2])
    assert len(result) == 1
    assert result[0] is level2


def test_two_classes():
    result = find_lowest_subclasses([level1, level2])
    assert len(result) == 1
    assert result[0] is level2


def test_four_deep():
    result = find_lowest_subclasses([level1, level2, level3, level4])
    assert len(result) == 1
    assert result[0] is level4


def test_four_deep_outoforder():
    result = find_lowest_subclasses([level2, level3, level1, level4])
    assert len(result) == 1
    assert result[0] is level4


def test_diverging_tree():
    result = find_lowest_subclasses([level1, level2, level3, level1a])
    assert len(result) == 2
    assert level1a in result and level3 in result
