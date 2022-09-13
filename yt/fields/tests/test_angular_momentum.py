import numpy as np

from yt.testing import assert_allclose_units, fake_amr_ds


def test_AM_value():
    ds = fake_amr_ds(
        fields=("Density", "velocity_x", "velocity_y", "velocity_z"),
        units=("g/cm**3", "cm/s", "cm/s", "cm/s"),
        length_unit=0.5,
    )

    sp = ds.sphere([0.5] * 3, (0.1, "code_length"))

    x0 = sp.center
    v0 = ds.arr([1, 2, 3], "km/s")

    sp.set_field_parameter("bulk_velocity", v0)

    X = (ds.arr([sp[("index", k)] for k in "xyz"]) - x0[:, None]).T
    V = (ds.arr([sp[("gas", f"velocity_{k}")] for k in "xyz"]) - v0[:, None]).T

    sAM_manual = ds.arr(np.cross(X, V), X.units * V.units)
    sAM = ds.arr([sp[("gas", f"specific_angular_momentum_{k}")] for k in "xyz"]).T

    assert_allclose_units(sAM_manual, sAM)
