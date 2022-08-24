import numpy as np

from yt.frontends.moab.api import MoabHex8Dataset
from yt.testing import (
    assert_almost_equal,
    assert_equal,
    requires_file,
    units_override_check,
)
from yt.utilities.answer_testing.framework import (
    FieldValuesTest,
    data_dir_load,
    requires_ds,
)

_fields = (("moab", "flux"),)

c5 = "c5/c5.h5m"


@requires_ds(c5)
def test_cantor_5():
    np.random.seed(0x4D3D3D3)
    ds = data_dir_load(c5)
    assert_equal(str(ds), "c5")
    dso = [
        None,
        ("sphere", ("c", (0.1, "unitary"))),
        ("sphere", ("c", (0.2, "unitary"))),
    ]
    dd = ds.all_data()
    assert_almost_equal(ds.index.get_smallest_dx(), 0.00411522633744843, 10)
    assert_equal(dd[("gas", "x")].shape[0], 63 * 63 * 63)
    assert_almost_equal(
        dd[("index", "cell_volume")].in_units("code_length**3").sum(dtype="float64").d,
        1.0,
        10,
    )
    for offset_1 in [1e-9, 1e-4, 0.1]:
        for offset_2 in [1e-9, 1e-4, 0.1]:
            DLE = ds.domain_left_edge
            DRE = ds.domain_right_edge
            ray = ds.ray(DLE + offset_1 * DLE.uq, DRE - offset_2 * DRE.uq)
            assert_almost_equal(ray["dts"].sum(dtype="float64"), 1.0, 8)
    for p1 in np.random.random((5, 3)):
        for p2 in np.random.random((5, 3)):
            ray = ds.ray(p1, p2)
            assert_almost_equal(ray["dts"].sum(dtype="float64"), 1.0, 8)
    for field in _fields:
        for dobj_name in dso:
            yield FieldValuesTest(c5, field, dobj_name)


@requires_file(c5)
def test_MoabHex8Dataset():
    assert isinstance(data_dir_load(c5), MoabHex8Dataset)


@requires_file(c5)
def test_units_override():
    units_override_check(c5)
