from nose.plugins.attrib import attr

from yt.testing import ANSWER_TEST_TAG, fake_random_ds
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.api import LinePlot


@attr(ANSWER_TEST_TAG)
def test_lineplot_set_axis_properties():
    ds = fake_random_ds(16)

    def create_image(filename_prefix):
        p = LinePlot(
            ds,
            ("gas", "density"),
            start_point=[0, 0, 0],
            end_point=[1, 1, 1],
            npoints=32,
        )
        p.set_x_unit("cm")
        p.save(f"{filename_prefix}_xunit")

        p.set_unit(("gas", "density"), "kg/cm**3")
        p.save(f"{filename_prefix}_xunit_zunit")

        p.set_log(("gas", "density"), False)
        p.save(f"{filename_prefix}_xunit_zunit_lin")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_lineplot_set_axis_properties"
    test.answer_name = "lineplot_set_axis_properties"
    yield test
