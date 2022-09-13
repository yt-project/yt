# test_output.py requires pytest, and yield-based tests are only supported by nose
# this tests in this file need to be translated to pytest before they can reintegrate test_outputs.py

import pytest

from yt.utilities.answer_testing.framework import data_dir_load, requires_ds
from yt.visualization.api import ParticleProjectionPlot

pkdgrav = "halo1e11_run1.00400/halo1e11_run1.00400"
gasoline_dmonly = "agora_1e11.00400/agora_1e11.00400"
tipsy_gal = "TipsyGalaxy/galaxy.00300"

tested_fields = [
    (("all", "particle_mass"), None),
    (("all", "particle_ones"), None),
    (("all", "particle_velocity_x"), ("all", "particle_mass")),
    (("all", "particle_velocity_y"), ("all", "particle_mass")),
    (("all", "particle_velocity_z"), ("all", "particle_mass")),
]


class TestPKDGrav:
    @requires_ds(pkdgrav, big_data=True, file_check=True)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load(
            pkdgrav,
            kwargs=dict(
                field_dtypes={"Coordinates": "d"},
                cosmology_parameters=dict(
                    current_redshift=0.0,
                    omega_lambda=0.728,
                    omega_matter=0.272,
                    hubble_constant=0.702,
                ),
                unit_base={"length": (60.0, "Mpccm/h")},
            ),
        )
        cls.dd = cls.ds.all_data()

    def test_nparticles(self):
        expected = 26847360
        actual = sum(
            self.dd[ptype, "particle_position"].shape[0]
            for ptype in self.ds.particle_types_raw
        )
        assert actual == expected

    @pytest.mark.parametrize(
        "fields",
        tested_fields,
    )
    @pytest.mark.parametrize("axis", range(3))
    @pytest.mark.mpl_image_compare
    def test_plot(self, fields, axis):
        field, weight_field = fields
        assert field[0] in self.ds.particle_types

        pp = ParticleProjectionPlot(self.ds, axis, field, weight_field=weight_field)
        pp._setup_plots()
        return pp[field].figure


class TestGasolineDMOnly:
    @requires_ds(gasoline_dmonly, big_data=True, file_check=True)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load(
            gasoline_dmonly,
            kwargs=dict(
                cosmology_parameters=dict(
                    current_redshift=0.0,
                    omega_lambda=0.728,
                    omega_matter=0.272,
                    hubble_constant=0.702,
                ),
                unit_base={"length": (60.0, "Mpccm/h")},
            ),
        )
        cls.dd = cls.ds.all_data()

    def test_nparticles(self):
        expected = 10550576
        actual = sum(
            self.dd[ptype, "particle_position"].shape[0]
            for ptype in self.ds.particle_types_raw
        )
        assert actual == expected

    @pytest.mark.parametrize("fields", tested_fields)
    @pytest.mark.parametrize("axis", range(3))
    @pytest.mark.mpl_image_compare
    def test_plot(self, fields, axis):
        field, weight_field = fields
        assert field[0] in self.ds.particle_types

        pp = ParticleProjectionPlot(self.ds, axis, field, weight_field=weight_field)
        pp._setup_plots()
        return pp[field].figure


class TestTipsyGal:
    @requires_ds(tipsy_gal)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load(
            tipsy_gal,
            kwargs={"bounding_box": [[-2000, 2000], [-2000, 2000], [-2000, 2000]]},
        )
        cls.dd = cls.ds.all_data()

    def test_nparticles(self):
        expected = 315372
        actual = sum(
            self.dd[ptype, "particle_position"].shape[0]
            for ptype in self.ds.particle_types_raw
        )
        assert actual == expected

    @pytest.mark.parametrize(
        "fields",
        [
            (("Stars", "Metals"), None),
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
