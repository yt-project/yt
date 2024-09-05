from numpy.testing import assert_equal

from yt.utilities.periodic_table import _elements, periodic_table


def test_element_accuracy():
    for num, w, name, sym in _elements:
        e0 = periodic_table[num]
        e1 = periodic_table[name]
        e2 = periodic_table[sym]
        # If num == -1, then we are in one of the things like Deuterium or El
        # that are not elements by themselves.
        if num == -1:
            e0 = e1
        assert_equal(id(e0), id(e1))
        assert_equal(id(e0), id(e2))
        assert_equal(e0.num, num)
        assert_equal(e0.weight, w)
        assert_equal(e0.name, name)
        assert_equal(e0.symbol, sym)
