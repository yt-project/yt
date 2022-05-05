from nose.plugins.attrib import attr

from yt.testing import ANSWER_TEST_TAG, add_noise_fields, fake_random_ds
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.api import PhasePlot


@attr(ANSWER_TEST_TAG)
def test_phaseplot_set_colorbar_properties_implicit():
    ds = fake_random_ds(16)
    add_noise_fields(ds)

    def create_image(filename_prefix):
        my_sphere = ds.sphere("c", 1)
        p = PhasePlot(
            my_sphere,
            ("gas", "noise1"),
            ("gas", "noise3"),
            [("gas", "density")],
            weight_field=None,
        )
        # using implicit units
        p.set_zlim(("gas", "density"), zmax=10)
        p.save(f"{filename_prefix}_set_zlim_implicit")

        # changing units should affect the colorbar and not the image
        p.set_unit(("gas", "density"), "kg/AU**3")
        p.save(f"{filename_prefix}_set_zlim_set_unit_implicit")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_phaseplot_set_colorbar_properties_implicit"
    test.answer_name = "phaseplot_set_colorbar_properties_implicit"
    yield test
