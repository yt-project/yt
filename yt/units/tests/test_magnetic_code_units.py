import numpy as np

from yt.frontends.stream.api import load_uniform_grid
from yt.testing import assert_allclose


def test_magnetic_code_units():

    sqrt4pi = np.sqrt(4.0 * np.pi)
    ddims = (16,) * 3
    data = {"density": (np.random.uniform(size=ddims), "g/cm**3")}

    ds1 = load_uniform_grid(
        data, ddims, magnetic_unit=(sqrt4pi, "gauss"), unit_system="cgs"
    )

    assert_allclose(ds1.magnetic_unit.value, sqrt4pi)
    assert str(ds1.magnetic_unit.units) == "G"

    mucu = ds1.magnetic_unit.to("code_magnetic")
    assert_allclose(mucu.value, 1.0)
    assert str(mucu.units) == "code_magnetic"

    ds2 = load_uniform_grid(data, ddims, magnetic_unit=(1.0, "T"), unit_system="cgs")

    assert_allclose(ds2.magnetic_unit.value, 10000.0)
    assert str(ds2.magnetic_unit.units) == "G"

    mucu = ds2.magnetic_unit.to("code_magnetic")
    assert_allclose(mucu.value, 1.0)
    assert str(mucu.units) == "code_magnetic"

    ds3 = load_uniform_grid(data, ddims, magnetic_unit=(1.0, "T"), unit_system="mks")

    assert_allclose(ds3.magnetic_unit.value, 1.0)
    assert str(ds3.magnetic_unit.units) == "T"

    mucu = ds3.magnetic_unit.to("code_magnetic")
    assert_allclose(mucu.value, 1.0)
    assert str(mucu.units) == "code_magnetic"

    ds4 = load_uniform_grid(
        data, ddims, magnetic_unit=(1.0, "gauss"), unit_system="mks"
    )

    assert_allclose(ds4.magnetic_unit.value, 1.0e-4)
    assert str(ds4.magnetic_unit.units) == "T"

    mucu = ds4.magnetic_unit.to("code_magnetic")
    assert_allclose(mucu.value, 1.0)
    assert str(mucu.units) == "code_magnetic"
