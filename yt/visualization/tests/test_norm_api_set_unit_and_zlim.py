import numpy.testing as npt

from yt.testing import fake_random_ds
from yt.visualization.api import SlicePlot


def test_sliceplot_set_unit_and_zlim_order():
    ds = fake_random_ds(16)
    field = ("gas", "density")

    p0 = SlicePlot(ds, "z", field)
    p0.set_unit(field, "kg/m**3")
    p0.set_zlim(field, zmin=0)

    # reversing order of operations
    p1 = SlicePlot(ds, "z", field)
    p1.set_zlim(field, zmin=0)
    p1.set_unit(field, "kg/m**3")

    p0._setup_plots()
    p1._setup_plots()

    im0 = p0.plots[field].image._A
    im1 = p1.plots[field].image._A

    npt.assert_allclose(im0, im1)
