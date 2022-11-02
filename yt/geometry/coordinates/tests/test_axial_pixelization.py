from yt.testing import _geom_transforms, fake_amr_ds
from yt.utilities.answer_testing.framework import AxialPixelizationTest


def test_axial_pixelization():
    for geom in sorted(_geom_transforms):
        if geom == "spectral_cube":
            # skip this case as it was added much later and we don't want to keep
            # adding yield-based tests during the nose->pytest migration
            continue
        ds = fake_amr_ds(geometry=geom)
        yield AxialPixelizationTest(ds)
