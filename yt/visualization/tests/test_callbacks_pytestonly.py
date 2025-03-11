import pytest

import yt
from yt.testing import fake_amr_ds


@pytest.mark.parametrize("plot_args", [{"cmap": "jet"}, {"colors": "red"}])
def test_contour_callback_kwargs(plot_args):
    ds = fake_amr_ds()
    slc = yt.SlicePlot(ds, "x", ("stream", "Density"))
    slc.annotate_contour(("stream", "Density"), plot_args=plot_args)

    pargs = slc._callbacks[0].plot_args
    assert all(pargs[ky] == val for ky, val in plot_args.items())
    _ = slc.render()
