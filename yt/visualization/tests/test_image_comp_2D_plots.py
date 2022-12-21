# image tests using pytest-mpl
import re
from importlib.metadata import version
from itertools import chain
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

    p0.render()
    p1.render()

    im0 = p0.plots[field].image._A
    im1 = p1.plots[field].image._A

    npt.assert_allclose(im0, im1)


@pytest.mark.mpl_image_compare
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

    p.render()
    return p.plots[field].figure


@pytest.mark.skipif(
    MPL_VERSION < Version("3.5"),
    reason=f"Correct behaviour requires MPL 3.5 we have {MPL_VERSION}",
)
@pytest.mark.mpl_image_compare
def test_sliceplot_custom_norm():
    from matplotlib.colors import TwoSlopeNorm

    ds = fake_random_ds(16)
    field = ("gas", "density")
    p = SlicePlot(ds, "z", field)
    p.set_norm(field, norm=TwoSlopeNorm(vcenter=0, vmin=-0.5, vmax=1))

    p.render()
    return p.plots[field].figure


@pytest.mark.skipif(
    MPL_VERSION >= Version("3.5"),
    reason=f"Testing a warning that should only happen with MPL < 3.5, we have {MPL_VERSION}",
)
@pytest.mark.mpl_image_compare
def test_sliceplot_custom_norm_warning():
    from matplotlib.colors import TwoSlopeNorm

    ds = fake_random_ds(16)
    field = ("gas", "density")
    p = SlicePlot(ds, "z", field)
    p.set_norm(field, norm=TwoSlopeNorm(vcenter=0, vmin=-0.5, vmax=1))

    with pytest.warns(
        UserWarning,
        match=re.escape(
            "Failed to format colorbar ticks. "
            "This is expected when using the set_norm method "
            "with some matplotlib classes (e.g. TwoSlopeNorm) "
            "with matplotlib versions older than 3.5\n"
            "Please try upgrading matplotlib to a more recent version. "
            "If the problem persists, please file a report to "
            "https://github.com/yt-project/yt/issues/new"
        ),
    ):
        p.render()
    return p.plots[field].figure


@pytest.mark.mpl_image_compare
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

    p.render()
    return p.plots[field].figure


@pytest.mark.mpl_image_compare
def test_profileplot_set_axis_properties():
    ds = fake_random_ds(16)

    disk = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], (10, "m"), (1, "m"))
    p = ProfilePlot(disk, ("gas", "density"), [("gas", "velocity_x")])
    p.set_unit(("gas", "density"), "kg/cm**3")
    p.set_log(("gas", "density"), False)
    p.set_unit(("gas", "velocity_x"), "mile/hour")

    p.render()
    return p.plots["gas", "velocity_x"].figure


@pytest.mark.mpl_image_compare
def test_particleprojectionplot_set_colorbar_properties():
    ds = fake_particle_ds(npart=100)

    field = ("all", "particle_mass")
    p = ParticleProjectionPlot(ds, 2, field)
    p.set_buff_size(10)

    p.set_unit(field, "Msun")
    p.set_zlim(field, zmax=1e-35)
    p.set_log(field, False)

    p.render()
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
            id_prefix = str(i_reg)
            for x_field, y_field in pr_fields:
                id_suffix = "_".join([*x_field, *y_field])
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

    @pytest.mark.parametrize(
        "suffix",
        [None, "from_profiles", "fractional_accumulation"],
    )
    @pytest.mark.parametrize("region", ["0", "1"])
    @pytest.mark.parametrize(
        "xax, yax",
        [
            (("gas", "density"), ("gas", "temperature")),
            (("gas", "density"), ("gas", "velocity_x")),
            (("gas", "temperature"), ("gas", "mass")),
            (("gas", "density"), ("index", "radius")),
            (("gas", "velocity_magnitude"), ("gas", "mass")),
        ],
    )
    @pytest.mark.mpl_image_compare
    def test_profileplot_simple(self, region, xax, yax, suffix):
        key = "_".join(chain(region, xax, yax))
        if suffix is not None:
            key += f"_{suffix}"
        plots = list(self.profiles[key].plots.values())
        assert len(plots) == 1
        return plots[0].figure

    @pytest.mark.mpl_image_compare
    def test_profileplot_from_multiple_profiles_0(self):
        plots = list(self.profiles["from_multiple_profiles"].plots.values())
        assert len(plots) == 2
        return plots[0].figure

    @pytest.mark.mpl_image_compare
    def test_profileplot_from_multiple_profiles_1(self):
        plots = list(self.profiles["from_multiple_profiles"].plots.values())
        assert len(plots) == 2
        return plots[1].figure


class TestPhasePlot:
    @classmethod
    def setup_class(cls):
        fields = ("density", "temperature", "velocity_x", "velocity_y", "velocity_z")
        units = ("g/cm**3", "K", "cm/s", "cm/s", "cm/s")
        cls.ds = fake_random_ds(16, fields=fields, units=units)
        regions = [cls.ds.region([0.5] * 3, [0.4] * 3, [0.6] * 3), cls.ds.all_data()]
        pr_fields = [
            [("gas", "density"), ("gas", "temperature"), ("gas", "mass")],
            [("gas", "density"), ("gas", "velocity_x"), ("gas", "mass")],
            [
                ("index", "radius"),
                ("gas", "temperature"),
                ("gas", "velocity_magnitude"),
            ],
        ]
        cls.profiles: Dict[str, PhasePlot] = {}
        for i_reg, reg in enumerate(regions):
            id_prefix = str(i_reg)
            for x_field, y_field, z_field in pr_fields:
                id_suffix = "_".join([*x_field, *y_field, *z_field])
                base_id = f"{id_prefix}_{id_suffix}"
                cls.profiles[base_id] = PhasePlot(
                    reg, x_field, y_field, z_field, x_bins=16, y_bins=16
                )
                cls.profiles[f"{base_id}_fractional_accumulation"] = PhasePlot(
                    reg,
                    x_field,
                    y_field,
                    z_field,
                    fractional=True,
                    accumulation=True,
                    x_bins=16,
                    y_bins=16,
                )

                p2d = create_profile(reg, [x_field, y_field], z_field, n_bins=[16, 16])
                cls.profiles[f"{base_id}_from_profiles"] = PhasePlot.from_profile(p2d)

    @pytest.mark.parametrize(
        "suffix",
        [None, "from_profiles", "fractional_accumulation"],
    )
    @pytest.mark.parametrize("region", ["0", "1"])
    @pytest.mark.parametrize(
        "xax, yax, zax",
        [
            (("gas", "density"), ("gas", "temperature"), ("gas", "mass")),
            (("gas", "density"), ("gas", "velocity_x"), ("gas", "mass")),
            (
                ("index", "radius"),
                ("gas", "temperature"),
                ("gas", "velocity_magnitude"),
            ),
        ],
    )
    @pytest.mark.mpl_image_compare
    def test_phaseplot(self, region, xax, yax, zax, suffix):
        key = "_".join(chain(region, xax, yax, zax))
        if suffix is not None:
            key += f"_{suffix}"
        plots = list(self.profiles[key].plots.values())
        assert len(plots) == 1
        return plots[0].figure


class TestPhasePlotSetZlim:
    @classmethod
    def setup_class(cls):
        cls.ds = fake_random_ds(16)
        add_noise_fields(cls.ds)
        cls.data = cls.ds.sphere("c", 1)

    @pytest.mark.mpl_image_compare
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
        p.render()
        return p.plots[field].figure

    @pytest.mark.mpl_image_compare
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
        p.render()
        return p.plots[field].figure


class TestSetBackgroundColor:
    # see https://github.com/yt-project/yt/issues/3854

    @classmethod
    def setup_class(cls):
        cls.ds = fake_random_ds(16)

    @pytest.mark.mpl_image_compare
    def test_sliceplot_set_background_color_linear(self):
        field = ("gas", "density")
        p = SlicePlot(self.ds, "z", field, width=1.5)
        p.set_background_color(field, color="C0")
        p.set_log(field, False)

        p.render()
        return p.plots[field].figure

    @pytest.mark.mpl_image_compare
    def test_sliceplot_set_background_color_log(self):
        field = ("gas", "density")
        p = SlicePlot(self.ds, "z", field, width=1.5)
        p.set_background_color(field, color="C0")

        p.render()
        return p.plots[field].figure


class TestCylindricalZSlicePlot:
    @classmethod
    def setup_class(cls):
        cls.ds = fake_amr_ds(geometry="cylindrical")
        add_noise_fields(cls.ds)
        fields = ["noise%d" % i for i in range(4)]
        cls.plot = SlicePlot(cls.ds, "z", fields)

    @pytest.mark.parametrize("field", ["noise0", "noise1", "noise2", "noise3"])
    @pytest.mark.mpl_image_compare
    def test_cylindrical_z_log(self, field):
        return self.plot.plots[field].figure

    @pytest.mark.parametrize("field", ["noise0", "noise1", "noise2", "noise3"])
    @pytest.mark.mpl_image_compare
    def test_cylindrical_z_linear(self, field):
        self.plot.set_log("noise0", False)
        return self.plot.plots[field].figure


class TestSphericalPhiSlicePlot:
    @classmethod
    def setup_class(cls):
        cls.ds = fake_amr_ds(geometry="spherical")
        add_noise_fields(cls.ds)
        fields = ["noise%d" % i for i in range(4)]
        cls.plot = SlicePlot(cls.ds, "phi", fields)

    @pytest.mark.parametrize("field", ["noise0", "noise1", "noise2", "noise3"])
    @pytest.mark.mpl_image_compare
    def test_spherical_phi_log(self, field):
        return self.plot.plots[field].figure


class TestSphericalThetaSlicePlot:
    @classmethod
    def setup_class(cls):
        cls.ds = fake_amr_ds(geometry="spherical")
        add_noise_fields(cls.ds)
        fields = ["noise%d" % i for i in range(4)]
        cls.plot = SlicePlot(cls.ds, "theta", fields)

    @pytest.mark.parametrize("field", ["noise0", "noise1", "noise2", "noise3"])
    @pytest.mark.mpl_image_compare
    def test_spherical_theta_log(self, field):
        return self.plot.plots[field].figure
