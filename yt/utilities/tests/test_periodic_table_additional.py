import numpy as np
from numpy.testing import assert_equal

from yt.utilities.periodic_table import periodic_table


def test_periodic_table_numeric_lookup_and_repr():
    oxygen = periodic_table[np.int64(8)]

    assert_equal(oxygen.num, 8)
    assert_equal(oxygen.weight, 15.9994000000)
    assert_equal(oxygen.name, "Oxygen")
    assert_equal(oxygen.symbol, "O")
    assert_equal(repr(oxygen), "Element: O (Oxygen)")


def test_periodic_table_three_letter_symbol_lookup():
    ununoctium = periodic_table["Uuo"]

    assert_equal(ununoctium.num, 118)
    assert_equal(ununoctium.weight, 294.0000000000)
    assert_equal(ununoctium.name, "Ununoctium")
    assert_equal(ununoctium.symbol, "Uuo")


def test_periodic_table_invalid_key_type_raises_keyerror():
    invalid_key = ("O",)

    try:
        periodic_table[invalid_key]
    except KeyError as exc:
        assert_equal(exc.args[0], invalid_key)
    else:
        raise AssertionError("Expected KeyError for tuple lookup in periodic_table")
