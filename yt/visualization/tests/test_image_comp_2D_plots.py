# image tests using pytest-mpl
from itertools import chain

import matplotlib as mpl
import numpy as np
import numpy.testing as npt
import pytest
from matplotlib.colors import SymLogNorm

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


def test_annotation_parse_h():
    ds = fake_random_ds(16)

    # Make sure `h` (reduced Hubble constant) is not equal to 1
    ds.unit_registry.modify("h", 0.7)

    rad = ds.quan(0.1, "cm/h")
    center = ds.arr([0.5] * 3, "code_length")

    # Twice the same slice plot
    p1 = SlicePlot(ds, "x", "density")
    p2 = SlicePlot(ds, "x", "density")

    # But the *same* center is given in different units
    p1.annotate_sphere(center.to("cm"), rad, circle_args={"color": "black"})
    p2.annotate_sphere(center.to("cm/h"), rad, circle_args={"color": "black"})

    # Render annotations, and extract matplotlib image
    # as an RGB array
    p1.render()
    p1.plots["gas", "density"].figure.canvas.draw()
    img1 = p1.plots["gas", "density"].figure.canvas.renderer.buffer_rgba()

    p2.render()
    p2.plots["gas", "density"].figure.canvas.draw()
    img2 = p2.plots["gas", "density"].figure.canvas.renderer.buffer_rgba()

    # This should be the same image
    npt.assert_allclose(img1, img2)


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


@pytest.mark.mpl_image_compare
def test_sliceplot_custom_norm():
    from matplotlib.colors import TwoSlopeNorm

    ds = fake_random_ds(16)
    field = ("gas", "density")
    p = SlicePlot(ds, "z", field)
    p.set_norm(field, norm=TwoSlopeNorm(vcenter=0, vmin=-0.5, vmax=1))

    p.render()
    return p.plots[field].figure


@pytest.mark.mpl_image_compare
def test_sliceplot_custom_norm_symlog_int_base():
    ds = fake_random_ds(16)
    add_noise_fields(ds)
    field = "noise3"
    p = SlicePlot(ds, "z", field)

    # using integer base !=10 and >2 to exercise special case
    # for colorbar minor ticks
    p.set_norm(field, norm=SymLogNorm(linthresh=0.1, base=5))

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


class TestMultipanelPlot:
    @classmethod
    def setup_class(cls):
        cls.fields = [
            ("gas", "density"),
            ("gas", "velocity_x"),
            ("gas", "velocity_y"),
            ("gas", "velocity_magnitude"),
        ]
        cls.ds = fake_random_ds(16)

    @pytest.mark.skipif(
        mpl.__version_info__ < (3, 7),
        reason="colorbar cannot currently be set horizontal in multi-panel plot with matplotlib older than 3.7.0",
    )
    @pytest.mark.parametrize("cbar_location", ["top", "bottom", "left", "right"])
    @pytest.mark.mpl_image_compare
    def test_multipanelplot_colorbar_orientation_simple(self, cbar_location):
        p = SlicePlot(self.ds, "z", self.fields)
        return p.export_to_mpl_figure((2, 2), cbar_location=cbar_location)

    @pytest.mark.parametrize("cbar_location", ["top", "bottom"])
    def test_multipanelplot_colorbar_orientation_warning(self, cbar_location):
        p = SlicePlot(self.ds, "z", self.fields)
        if mpl.__version_info__ < (3, 7):
            with pytest.warns(
                UserWarning,
                match="Cannot properly set the orientation of colorbar.",
            ):
                p.export_to_mpl_figure((2, 2), cbar_location=cbar_location)
        else:
            p.export_to_mpl_figure((2, 2), cbar_location=cbar_location)


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
        cls.profiles: dict[str, ProfilePlot] = {}
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
        cls.profiles: dict[str, PhasePlot] = {}
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

        def some_nans_field(field, data):
            ret = data["gas", "density"]
            ret[::2] *= np.nan
            return ret

        cls.ds.add_field(
            name=("gas", "polluted_field"),
            function=some_nans_field,
            sampling_type="local",
        )

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

    @pytest.mark.mpl_image_compare
    def test_sliceplot_set_background_color_and_bad_value(self):
        # see https://github.com/yt-project/yt/issues/4639
        field = ("gas", "polluted_field")
        p = SlicePlot(self.ds, "z", field, width=1.5)
        p.set_background_color(field, color="black")

        # copy the default colormap
        cmap = mpl.colormaps["cmyt.arbre"]
        cmap.set_bad("red")
        p.set_cmap(field, cmap)

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

    @pytest.mark.parametrize(
        "theta_min, theta_max",
        [
            pytest.param(0, 2 * np.pi, id="full_azimuthal_domain"),
            pytest.param(3 / 4 * np.pi, 5 / 4 * np.pi, id="restricted_sector"),
        ],
    )
    @pytest.mark.mpl_image_compare
    def test_exclude_pixels_with_partial_bbox_intersection(self, theta_min, theta_max):
        rmin = 1.0
        rmax = 2.0
        ds = fake_amr_ds(
            geometry="cylindrical",
            domain_left_edge=[rmin, 0, theta_min],
            domain_right_edge=[rmax, 1, theta_max],
        )
        add_noise_fields(ds)
        plot = SlicePlot(ds, "z", ("gas", "noise0"))
        for radius in [rmin - 0.01, rmax]:
            plot.annotate_sphere(
                center=[0, 0, 0],
                radius=radius,
                circle_args={"color": "red", "alpha": 0.4, "linewidth": 3},
            )
        plot.annotate_title("all pixels beyond (or on) red lines should be white")
        plot.render()
        return plot.plots["gas", "noise0"].figure


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
