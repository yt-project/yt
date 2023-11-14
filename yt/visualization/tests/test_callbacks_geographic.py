import pytest

from yt import SlicePlot
from yt.testing import assert_fname, fake_amr_ds, requires_module


@requires_module("cartopy")
@pytest.mark.parametrize("geometry", ["geographic", "internal_geographic"])
def test_quiver_callback_geographic(geometry, tmp_path):
    flds = ("density", "velocity_ew", "velocity_ns")
    units = ("g/cm**3", "m/s", "m/s")

    ds = fake_amr_ds(fields=flds, units=units, geometry=geometry)
    tmpfi = str(tmp_path / "geo_quiver.png")

    for ax in ds.coordinates.axis_order:
        slc = SlicePlot(ds, ax, "density", buff_size=(50, 50))
        slc.annotate_quiver(("stream", "velocity_ew"), ("stream", "velocity_ns"))
        slc.save(tmpfi)
        assert_fname(tmpfi)
