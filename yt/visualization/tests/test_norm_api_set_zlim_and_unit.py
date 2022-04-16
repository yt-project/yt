from nose.plugins.attrib import attr

from yt.testing import ANSWER_TEST_TAG, fake_random_ds
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.api import SlicePlot

ds = fake_random_ds(16)


@attr(ANSWER_TEST_TAG)
def test_sliceplot_set_zlim_and_unit():
    def create_image(filename_prefix):
        field = ("gas", "density")
        p = SlicePlot(ds, "z", field)
        p.set_zlim(field, zmin=0)
        p.set_unit(field, "kg/m**3")
        p.save(filename_prefix)

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_sliceplot_set_zlim_and_unit"
    test.answer_name = "sliceplot_set_zlim_and_unit"
    yield test
