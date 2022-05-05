from nose.plugins.attrib import attr

from yt.testing import ANSWER_TEST_TAG, fake_random_ds
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.api import SlicePlot


@attr(ANSWER_TEST_TAG)
def test_sliceplot_set_background_color():
    # see https://github.com/yt-project/yt/issues/3854
    ds = fake_random_ds(16)

    def create_image(filename_prefix):
        field = ("gas", "density")
        p = SlicePlot(ds, "z", field, width=1.5)
        p.set_background_color(field, color="C0")
        p.save(f"{filename_prefix}_log")
        p.set_log(("gas", "density"), False)
        p.save(f"{filename_prefix}_lin")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_sliceplot_set_background_color"
    test.answer_name = "sliceplot_set_background_color"
    yield test
