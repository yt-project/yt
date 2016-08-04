from yt.testing import \
    fake_amr_ds, _geom_transforms
from yt.utilities.answer_testing.framework import \
    AxialPixelizationTest

def test_axial_pixelization():
    for geom in sorted(_geom_transforms):
        ds = fake_amr_ds(geometry=geom)
        yield AxialPixelizationTest(ds)
