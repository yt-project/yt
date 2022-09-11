# image tests using pytest-mpl
import sys
from typing import Dict

import numpy as np
import numpy.testing as npt
import pytest
from packaging.version import Version

from yt.data_objects.profiles import create_profile
from yt.loaders import load_uniform_grid
from yt.testing import add_noise_fields, fake_amr_ds, fake_particle_ds, fake_random_ds
from yt.visualization.api import (
    LinePlot,
    ParticleProjectionPlot,
    PhasePlot,
    ProfilePlot,
    SlicePlot,
)

if sys.version_info >= (3, 8):
    from importlib.metadata import version
else:
    from importlib_metadata import version

MPL_VERSION = Version(version("matplotlib"))


def test_sliceplot_set_unit_and_zlim_order():
    ds = fake_random_ds(16)
    field = ("gas", "density")

    p0 = SlicePlot(ds, "z", field)
    p0.set_unit(field, "kg/m**3")
    p0.set_zlim(field, zmin=0)

    # reversing order of operations
    p1 = SlicePlot(ds, "z", field)
    p1.set_zlim(field, zmin=0)
    p1.set_unit(field, "kg/m**3")

    p0._setup_plots()
    p1._setup_plots()

    im0 = p0.plots[field].image._A
    im1 = p1.plots[field].image._A

    npt.assert_allclose(im0, im1)


@pytest.mark.mpl_image_compare(filename="inf_and_finite_values_set_zlim.png")
def test_inf_and_finite_values_set_zlim():
    # see https://github.com/yt-project/yt/issues/3901
    shape = (32, 16, 1)
    a = np.ones(16)
    b = np.ones((32, 16))
    c = np.reshape(a * b, shape)

    # injecting an inf
    c[0, 0, 0] = np.inf

    field = ("gas", "density")
    data = {field: c}

    ds = load_uniform_grid(
        data,
        shape,
        bbox=np.array([[0, 1], [0, 1], [0, 1]]),
    )

    p = SlicePlot(ds, "z", field)

    # setting zlim manually
    p.set_zlim(field, -10, 10)

    p._setup_plots()
    return p.plots[field].figure


@pytest.mark.skipif(
    MPL_VERSION < Version("3.2"),
    reason=f"TwoSlopeNorm requires MPL 3.2, we have {MPL_VERSION}",
)
@pytest.mark.mpl_image_compare(filename="sliceplot_custom_norm.png")
def test_sliceplot_custom_norm():
    from matplotlib.colors import TwoSlopeNorm

    ds = fake_random_ds(16)
    field = ("gas", "density")
    p = SlicePlot(ds, "z", field)
    p.set_norm(field, norm=TwoSlopeNorm(vcenter=0, vmin=-0.5, vmax=1))

    # this shouldn't be necessary, we should have public api to do this
    p._setup_plots()
    return p.plots[field].figure


@pytest.mark.mpl_image_compare(filename="lineplot_set_axis_properties.png")
def test_lineplot_set_axis_properties():
    ds = fake_random_ds(16)
    field = ("gas", "density")
    p = LinePlot(
        ds,
        field,
        start_point=[0, 0, 0],
        end_point=[1, 1, 1],
        npoints=32,
    )
    p.set_x_unit("cm")
    p.set_unit(field, "kg/cm**3")
    p.set_log(field, False)

    p._setup_plots()
    return p.plots[field].figure


@pytest.mark.mpl_image_compare(filename="profileplot_set_axis_properties.png")
def test_profileplot_set_axis_properties():
    ds = fake_random_ds(16)

    disk = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], (10, "m"), (1, "m"))
    p = ProfilePlot(disk, ("gas", "density"), [("gas", "velocity_x")])
    p.set_unit(("gas", "density"), "kg/cm**3")
    p.set_log(("gas", "density"), False)
    p.set_unit(("gas", "velocity_x"), "mile/hour")

    p._setup_plots()
    return p.plots["gas", "velocity_x"].figure


@pytest.mark.mpl_image_compare(
    filename="particleprojectionplot_set_colorbar_properties.png"
)
def test_particleprojectionplot_set_colorbar_properties():
    ds = fake_particle_ds(npart=100)

    field = ("all", "particle_mass")
    p = ParticleProjectionPlot(ds, 2, field)
    p.set_buff_size(10)

    p.set_unit(field, "Msun")
    p.set_zlim(field, zmax=1e-35)
    p.set_log(field, False)

    p._setup_plots()
    return p.plots[field].figure


class TestProfilePlot:
    @classmethod
    def setup_class(cls):
        fields = ("density", "temperature", "velocity_x", "velocity_y", "velocity_z")
        units = ("g/cm**3", "K", "cm/s", "cm/s", "cm/s")
        cls.ds = fake_random_ds(16, fields=fields, units=units)
        regions = [cls.ds.region([0.5] * 3, [0.4] * 3, [0.6] * 3), cls.ds.all_data()]
        pr_fields = [
            [("gas", "density"), ("gas", "temperature")],
            [("gas", "density"), ("gas", "velocity_x")],
            [("gas", "temperature"), ("gas", "mass")],
            [("gas", "density"), ("index", "radius")],
            [("gas", "velocity_magnitude"), ("gas", "mass")],
        ]
        cls.profiles: Dict[str, ProfilePlot] = {}
        for i_reg, reg in enumerate(regions):
            id_prefix = f"region{i_reg}"
            for x_field, y_field in pr_fields:
                id_suffix = f"{x_field[0]}_{x_field[1]}_{y_field[0]}_{y_field[1]}"
                base_id = f"{id_prefix}_{id_suffix}"
                cls.profiles[base_id] = ProfilePlot(reg, x_field, y_field)
                cls.profiles[f"{base_id}_fractional_accumulation"] = ProfilePlot(
                    reg, x_field, y_field, fractional=True, accumulation=True
                )

                p1d = create_profile(reg, x_field, y_field)
                cls.profiles[f"{base_id}_from_profiles"] = ProfilePlot.from_profiles(
                    p1d
                )

        p1 = create_profile(
            cls.ds.all_data(), ("gas", "density"), ("gas", "temperature")
        )
        p2 = create_profile(
            cls.ds.all_data(), ("gas", "density"), ("gas", "velocity_x")
        )
        cls.profiles["from_multiple_profiles"] = ProfilePlot.from_profiles(
            [p1, p2], labels=["temperature", "velocity"]
        )

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_density_gas_temperature.png"
    )
    def test_profileplot_region0_gas_density_gas_temperature(self):
        plots = list(
            self.profiles["region0_gas_density_gas_temperature"].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_density_gas_temperature_fractional_accumulation.png"
    )
    def test_profileplot_region0_gas_density_gas_temperature_fractional_accumulation(
        self,
    ):
        plots = list(
            self.profiles[
                "region0_gas_density_gas_temperature_fractional_accumulation"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_density_gas_temperature_from_profiles.png"
    )
    def test_profileplot_region0_gas_density_gas_temperature_from_profiles(self):
        plots = list(
            self.profiles[
                "region0_gas_density_gas_temperature_from_profiles"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_density_gas_velocity_x.png"
    )
    def test_profileplot_region0_gas_density_gas_velocity_x(self):
        plots = list(self.profiles["region0_gas_density_gas_velocity_x"].plots.values())
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_density_gas_velocity_x_fractional_accumulation.png"
    )
    def test_profileplot_region0_gas_density_gas_velocity_x_fractional_accumulation(
        self,
    ):
        plots = list(
            self.profiles[
                "region0_gas_density_gas_velocity_x_fractional_accumulation"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_density_gas_velocity_x_from_profiles.png"
    )
    def test_profileplot_region0_gas_density_gas_velocity_x_from_profiles(self):
        plots = list(
            self.profiles[
                "region0_gas_density_gas_velocity_x_from_profiles"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_temperature_gas_mass.png"
    )
    def test_profileplot_region0_gas_temperature_gas_mass(self):
        plots = list(self.profiles["region0_gas_temperature_gas_mass"].plots.values())
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_temperature_gas_mass_fractional_accumulation.png"
    )
    def test_profileplot_region0_gas_temperature_gas_mass_fractional_accumulation(self):
        plots = list(
            self.profiles[
                "region0_gas_temperature_gas_mass_fractional_accumulation"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_temperature_gas_mass_from_profiles.png"
    )
    def test_profileplot_region0_gas_temperature_gas_mass_from_profiles(self):
        plots = list(
            self.profiles[
                "region0_gas_temperature_gas_mass_from_profiles"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_density_index_radius.png"
    )
    def test_profileplot_region0_gas_density_index_radius(self):
        plots = list(self.profiles["region0_gas_density_index_radius"].plots.values())
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_density_index_radius_fractional_accumulation.png"
    )
    def test_profileplot_region0_gas_density_index_radius_fractional_accumulation(self):
        plots = list(
            self.profiles[
                "region0_gas_density_index_radius_fractional_accumulation"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_density_index_radius_from_profiles.png"
    )
    def test_profileplot_region0_gas_density_index_radius_from_profiles(self):
        plots = list(
            self.profiles[
                "region0_gas_density_index_radius_from_profiles"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_velocity_magnitude_gas_mass.png"
    )
    def test_profileplot_region0_gas_velocity_magnitude_gas_mass(self):
        plots = list(
            self.profiles["region0_gas_velocity_magnitude_gas_mass"].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_velocity_magnitude_gas_mass_fractional_accumulation.png"
    )
    def test_profileplot_region0_gas_velocity_magnitude_gas_mass_fractional_accumulation(
        self,
    ):
        plots = list(
            self.profiles[
                "region0_gas_velocity_magnitude_gas_mass_fractional_accumulation"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region0_gas_velocity_magnitude_gas_mass_from_profiles.png"
    )
    def test_profileplot_region0_gas_velocity_magnitude_gas_mass_from_profiles(self):
        plots = list(
            self.profiles[
                "region0_gas_velocity_magnitude_gas_mass_from_profiles"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_density_gas_temperature.png"
    )
    def test_profileplot_region1_gas_density_gas_temperature(self):
        plots = list(
            self.profiles["region1_gas_density_gas_temperature"].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_density_gas_temperature_fractional_accumulation.png"
    )
    def test_profileplot_region1_gas_density_gas_temperature_fractional_accumulation(
        self,
    ):
        plots = list(
            self.profiles[
                "region1_gas_density_gas_temperature_fractional_accumulation"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_density_gas_temperature_from_profiles.png"
    )
    def test_profileplot_region1_gas_density_gas_temperature_from_profiles(self):
        plots = list(
            self.profiles[
                "region1_gas_density_gas_temperature_from_profiles"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_density_gas_velocity_x.png"
    )
    def test_profileplot_region1_gas_density_gas_velocity_x(self):
        plots = list(self.profiles["region1_gas_density_gas_velocity_x"].plots.values())
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_density_gas_velocity_x_fractional_accumulation.png"
    )
    def test_profileplot_region1_gas_density_gas_velocity_x_fractional_accumulation(
        self,
    ):
        plots = list(
            self.profiles[
                "region1_gas_density_gas_velocity_x_fractional_accumulation"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_density_gas_velocity_x_from_profiles.png"
    )
    def test_profileplot_region1_gas_density_gas_velocity_x_from_profiles(self):
        plots = list(
            self.profiles[
                "region1_gas_density_gas_velocity_x_from_profiles"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_temperature_gas_mass.png"
    )
    def test_profileplot_region1_gas_temperature_gas_mass(self):
        plots = list(self.profiles["region1_gas_temperature_gas_mass"].plots.values())
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_temperature_gas_mass_fractional_accumulation.png"
    )
    def test_profileplot_region1_gas_temperature_gas_mass_fractional_accumulation(self):
        plots = list(
            self.profiles[
                "region1_gas_temperature_gas_mass_fractional_accumulation"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_temperature_gas_mass_from_profiles.png"
    )
    def test_profileplot_region1_gas_temperature_gas_mass_from_profiles(self):
        plots = list(
            self.profiles[
                "region1_gas_temperature_gas_mass_from_profiles"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_density_index_radius.png"
    )
    def test_profileplot_region1_gas_density_index_radius(self):
        plots = list(self.profiles["region1_gas_density_index_radius"].plots.values())
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_density_index_radius_fractional_accumulation.png"
    )
    def test_profileplot_region1_gas_density_index_radius_fractional_accumulation(self):
        plots = list(
            self.profiles[
                "region1_gas_density_index_radius_fractional_accumulation"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_density_index_radius_from_profiles.png"
    )
    def test_profileplot_region1_gas_density_index_radius_from_profiles(self):
        plots = list(
            self.profiles[
                "region1_gas_density_index_radius_from_profiles"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_velocity_magnitude_gas_mass.png"
    )
    def test_profileplot_region1_gas_velocity_magnitude_gas_mass(self):
        plots = list(
            self.profiles["region1_gas_velocity_magnitude_gas_mass"].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_velocity_magnitude_gas_mass_fractional_accumulation.png"
    )
    def test_profileplot_region1_gas_velocity_magnitude_gas_mass_fractional_accumulation(
        self,
    ):
        plots = list(
            self.profiles[
                "region1_gas_velocity_magnitude_gas_mass_fractional_accumulation"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(
        filename="profile_plot_region1_gas_velocity_magnitude_gas_mass_from_profiles.png"
    )
    def test_profileplot_region1_gas_velocity_magnitude_gas_mass_from_profiles(self):
        plots = list(
            self.profiles[
                "region1_gas_velocity_magnitude_gas_mass_from_profiles"
            ].plots.values()
        )
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare(filename="profile_plot_from_multiple_profiles_0.png")
    def test_profileplot_from_multiple_profiles_0(self):
        plots = list(self.profiles["from_multiple_profiles"].plots.values())
        assert len(plots) == 2
        return plots[0].figure

    @pytest.mark.mpl_image_compare(filename="profile_plot_from_multiple_profiles_1.png")
    def test_profileplot_from_multiple_profiles_1(self):
        plots = list(self.profiles["from_multiple_profiles"].plots.values())
        assert len(plots) == 2
        return plots[1].figure


class TestPhasePlot:
    @classmethod
    def setup_class(cls):
        ds = fake_random_ds(16)
        add_noise_fields(ds)
        cls.ds = (
            ds  # create a persistent reference to avoid premature garbage collection
        )
        cls.data = ds.sphere("c", 1)

    @pytest.mark.mpl_image_compare(
        filename="phaseplot_set_zlim_with_implicit_units.png"
    )
    def test_phaseplot_set_zlim_with_implicit_units(self):
        p = PhasePlot(
            self.data,
            ("gas", "noise1"),
            ("gas", "noise3"),
            [("gas", "density")],
            weight_field=None,
        )
        field = ("gas", "density")
        p.set_zlim(field, zmax=10)
        p._setup_plots()
        return p.plots[field].figure

    @pytest.mark.mpl_image_compare(
        filename="phaseplot_set_zlim_with_explicit_units.png"
    )
    def test_phaseplot_set_zlim_with_explicit_units(self):
        p = PhasePlot(
            self.data,
            ("gas", "noise1"),
            ("gas", "noise3"),
            [("gas", "density")],
            weight_field=None,
        )
        field = ("gas", "density")
        # using explicit units, we expect the colorbar units to stay unchanged
        p.set_zlim(field, zmin=(1e36, "kg/AU**3"))
        p._setup_plots()
        return p.plots[field].figure

    @pytest.mark.mpl_image_compare(filename="phaseplot_set_unit.png")
    def test_phaseplot_set_unit(self):
        p = PhasePlot(
            self.data,
            ("gas", "noise1"),
            ("gas", "noise3"),
            [("gas", "density")],
            weight_field=None,
        )
        field = ("gas", "density")
        p.set_unit(field, "kg/AU**3")
        p._setup_plots()
        return p.plots[field].figure


class TestSetBackgroundColor:
    # see https://github.com/yt-project/yt/issues/3854

    @classmethod
    def setup_class(cls):
        cls.ds = fake_random_ds(16)

    @pytest.mark.mpl_image_compare(filename="sliceplot_set_background_color_lin.png")
    def test_sliceplot_set_background_color_lin(self):
        field = ("gas", "density")
        p = SlicePlot(self.ds, "z", field, width=1.5)
        p.set_background_color(field, color="C0")
        p.set_log(field, False)

        p._setup_plots()
        return p.plots[field].figure

    @pytest.mark.mpl_image_compare(filename="sliceplot_set_background_color_log.png")
    def test_sliceplot_set_background_color_log(self):
        field = ("gas", "density")
        p = SlicePlot(self.ds, "z", field, width=1.5)
        p.set_background_color(field, color="C0")

        p._setup_plots()
        return p.plots[field].figure


class TestCylindricalZSlicePlot:
    @classmethod
    def setup_class(cls):
        cls.ds = fake_amr_ds(geometry="cylindrical")
        add_noise_fields(cls.ds)
        fields = ["noise%d" % i for i in range(4)]
        cls.plot = SlicePlot(cls.ds, "z", fields)

    @pytest.mark.mpl_image_compare(filename="test_cylindrical_z_noise0_log.png")
    def test_noise0_log(self):
        return self.plot.plots["noise0"].figure

    @pytest.mark.mpl_image_compare(filename="test_cylindrical_z_noise1_log.png")
    def test_noise1_log(self):
        return self.plot.plots["noise1"].figure

    @pytest.mark.mpl_image_compare(filename="test_cylindrical_z_noise2_log.png")
    def test_noise2_log(self):
        return self.plot.plots["noise2"].figure

    @pytest.mark.mpl_image_compare(filename="test_cylindrical_z_noise3_log.png")
    def test_noise3_log(self):
        return self.plot.plots["noise3"].figure

    @pytest.mark.mpl_image_compare(filename="test_cylindrical_z_noise0_lin.png")
    def test_noise0_lin(self):
        self.plot.set_log("noise0", False)
        return self.plot.plots["noise0"].figure

    @pytest.mark.mpl_image_compare(filename="test_cylindrical_z_noise1_lin.png")
    def test_noise1_lin(self):
        self.plot.set_log("noise1", False)
        return self.plot.plots["noise1"].figure

    @pytest.mark.mpl_image_compare(filename="test_cylindrical_z_noise2_lin.png")
    def test_noise2_lin(self):
        self.plot.set_log("noise2", False)
        return self.plot.plots["noise2"].figure

    @pytest.mark.mpl_image_compare(filename="test_cylindrical_z_noise3_lin.png")
    def test_noise3_lin(self):
        self.plot.set_log("noise3", False)
        return self.plot.plots["noise3"].figure


class TestSphericalPhiSlicePlot:
    @classmethod
    def setup_class(cls):
        cls.ds = fake_amr_ds(geometry="spherical")
        add_noise_fields(cls.ds)
        fields = ["noise%d" % i for i in range(4)]
        cls.plot = SlicePlot(cls.ds, "phi", fields)

    @pytest.mark.mpl_image_compare(filename="test_sph_phi_noise0_log.png")
    def test_noise0_log(self):
        return self.plot.plots["noise0"].figure

    @pytest.mark.mpl_image_compare(filename="test_sph_phi_noise1_log.png")
    def test_noise1_log(self):
        return self.plot.plots["noise1"].figure

    @pytest.mark.mpl_image_compare(filename="test_sph_phi_noise2_log.png")
    def test_noise2_log(self):
        return self.plot.plots["noise2"].figure

    @pytest.mark.mpl_image_compare(filename="test_sph_phi_noise3_log.png")
    def test_noise3_log(self):
        return self.plot.plots["noise3"].figure


class TestSphericalThetaSlicePlot:
    @classmethod
    def setup_class(cls):
        cls.ds = fake_amr_ds(geometry="spherical")
        add_noise_fields(cls.ds)
        fields = ["noise%d" % i for i in range(4)]
        cls.plot = SlicePlot(cls.ds, "theta", fields)

    @pytest.mark.mpl_image_compare(filename="test_sph_theta_noise0_log.png")
    def test_noise0_log(self):
        return self.plot.plots["noise0"].figure

    @pytest.mark.mpl_image_compare(filename="test_sph_theta_noise1_log.png")
    def test_noise1_log(self):
        return self.plot.plots["noise1"].figure

    @pytest.mark.mpl_image_compare(filename="test_sph_theta_noise2_log.png")
    def test_noise2_log(self):
        return self.plot.plots["noise2"].figure

    @pytest.mark.mpl_image_compare(filename="test_sph_theta_noise3_log.png")
    def test_noise3_log(self):
        return self.plot.plots["noise3"].figure
