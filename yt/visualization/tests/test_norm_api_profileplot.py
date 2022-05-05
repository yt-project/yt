from nose.plugins.attrib import attr

from yt.testing import ANSWER_TEST_TAG, fake_random_ds
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.api import ProfilePlot


@attr(ANSWER_TEST_TAG)
def test_profileplot_set_axis_properties():
    ds = fake_random_ds(16)

    def create_image(filename_prefix):
        disk = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], (10, "m"), (1, "m"))
        p = ProfilePlot(disk, ("gas", "density"), [("gas", "velocity_x")])
        p.save(f"{filename_prefix}_defaults")

        p.set_unit(("gas", "density"), "kg/cm**3")
        p.save(f"{filename_prefix}_xunit")

        p.set_log(("gas", "density"), False)
        p.save(f"{filename_prefix}_xunit_xlin")

        p.set_unit(("gas", "velocity_x"), "mile/hour")
        p.save(f"{filename_prefix}_xunit_xlin_yunit")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_profileplot_set_axis_properties"
    test.answer_name = "profileplot_set_axis_properties"
    yield test
