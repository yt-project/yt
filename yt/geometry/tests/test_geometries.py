from importlib.util import find_spec

import pytest

import yt
from yt.testing import fake_amr_ds


@pytest.mark.parametrize(
    "geometry",
    [
        "cartesian",
        "polar",
        "cylindrical",
        "spherical",
        "geographic",
        "internal_geographic",
        "spectral_cube",
    ],
)
def test_testable_geometries(geometry):
    # check that initializing a simple fake dataset works in any geometry
    ds = fake_amr_ds(fields=[("gas", "density")], units=["g/cm**3"], geometry=geometry)
    # make sure basic plotting works
    for axis in range(3):
        if "geographic" in geometry and axis == 2 and find_spec("cartopy") is None:
            pytest.skip(
                reason=(
                    "cannot test this case with vanilla yt (requires cartopy) "
                    "see https://github.com/yt-project/yt/issues/4182"
                )
            )
        yt.SlicePlot(ds, axis, ("gas", "density"), buff_size=(8, 8))
