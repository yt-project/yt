"""
Test all (reasonable) operation between floats (native or numpy), quantities,
numpy ndarrays, and yt arrays.

2013, Casey W. Stark <caseywstark@gmail.com>.

"""

from nose import with_setup
from nose.plugins.skip import SkipTest
import numpy as np

from yt.data_objects.yt_array import YTArray
from yt.utilities.quantities import Quantity
from yt.utilities.units import Unit, InvalidUnitOperationError

from yt.config import ytcfg

ignore_setting = "ignore_invalid_unit_operation_error"


def create_operands():
    # Make units less objects first
    s = [1.1, 2.2]
    a = [np.array([6.0, 4.3]), np.array([0.4, 7.9])]

    # Now make version with units.
    u = [Unit("pc"), Unit("Gyr")]
    q = [Quantity(s[0], u[0]), Quantity(s[1], u[1])]
    y = [YTArray(a[0], u[0]), Quantity(a[1], u[1])]

    operand_pair_list = [ (s[0], q[1]), (q[0], s[1]),
                          (a[0], q[1]), (q[0], a[1]),
                          (s[0], y[1]), (y[0], s[1]),
                          (a[0], y[1]), (y[0], a[1]),
                          (q[0], y[1]), (y[0], q[1]),
                          (q[0], q[1]), (y[0], y[1]) ]


def empty_teardown():
    pass

#@with_setup(create_operands, empty_teardown)
def test_add_sub_float_quantity():
    raise SkipTest

    f1 = 1.1
    f2 = 2.2
    f3 = 1.0
    a1 = np.array([4.6, 3.3])
    a2 = np.array([9.1, 2.7])
    a3 = np.array([2.1, 9.2])
    u1 = Unit("pc")
    u2 = Unit("yr")
    u3 = Unit("m")
    y1 = YTArray(a1, u1)
    y2 = YTArray(a2, u2)
    y3 = YTArray(a3, u3)

    # first the version with the rc ignore off.
    ytcfg[ignore_setting] = False

    # test operand pairs that should fail with the ignore off
    operand_pairs = [(f1, y2), (y1, f2),
                     (n1, y2), (y2, n1),
                     (y1, y2)]
    for o1, o2 in operand_pairs:
        try:
            res = o1 + o2
        except InvalidUnitOperationError:
            pass
        else:
            assert False

        try:
            res = o1 - o2
        except InvalidUnitOperationError:
            pass
        else:
            assert False

    # test pairs that should be good.
    res = f1 + y3
    assert isinstance(res, YTArray)
    assert res == a1 + a3 * u3.cgs_value / u1.cgs_value
    assert res.units == u1

    # now with rc ignore on
    ytcfg[ignore_setting] = True

    res = f1 + y2
    assert isistance(res, Quantity)
    assert res[0] == f1 + f2
    assert res.units == u2
