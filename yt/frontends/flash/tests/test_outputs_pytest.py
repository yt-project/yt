# test_output.py requires pytest, and yield-based tests are only supported by nose
# this tests in this file need to be translated to pytest before they can reintegrate test_outputs.py

import pytest

from yt.utilities.answer_testing.framework import data_dir_load, requires_ds
from yt.visualization.api import ParticleProjectionPlot

fid_1to3_b1 = "fiducial_1to3_b1/fiducial_1to3_b1_hdf5_part_0080"


class TestFid1to3B1:
    @requires_ds(fid_1to3_b1, big_data=True)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load(fid_1to3_b1)
        cls.dd = cls.ds.all_data()

    def test_nparticles(self):
        expected = 6684119
        actual = sum(
            self.dd[ptype, "particle_position"].shape[0]
            for ptype in self.ds.particle_types_raw
        )
        assert actual == expected

    @pytest.mark.parametrize(
        "fields",
        [
            (("all", "particle_mass"), None),
            (("all", "particle_ones"), None),
            (("all", "particle_velocity_x"), ("all", "particle_mass")),
            (("all", "particle_velocity_y"), ("all", "particle_mass")),
            (("all", "particle_velocity_z"), ("all", "particle_mass")),
        ],
    )
    @pytest.mark.parametrize("axis", range(3))
    @pytest.mark.mpl_image_compare
    def test_plot(self, fields, axis):
        field, weight_field = fields
        assert field[0] in self.ds.particle_types

        pp = ParticleProjectionPlot(self.ds, axis, field, weight_field=weight_field)
        pp._setup_plots()
        return pp[field].figure
