from nose.plugins.attrib import attr

from yt.testing import ANSWER_TEST_TAG, fake_particle_ds
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.api import ParticleProjectionPlot


@attr(ANSWER_TEST_TAG)
def test_particleprojectionplot_set_colorbar_properties():
    ds = fake_particle_ds(npart=100)

    def create_image(filename_prefix):
        field = ("all", "particle_mass")
        p = ParticleProjectionPlot(ds, 2, field)
        p.set_buff_size(10)

        p.set_unit(field, "Msun")
        p.save(f"{filename_prefix}_set_unit")

        p.set_zlim(field, zmax=1e-35)
        p.save(f"{filename_prefix}_set_unit_zlim")

        p.set_log(field, False)
        p.save(f"{filename_prefix}_set_unit_zlim_log")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_particleprojectionplot_set_colorbar_properties"
    test.answer_name = "particleprojectionplot_set_colorbar_properties"
    yield test
